#include "deep_app.h"

//#define DEMO

#include <assert.h>
#include <limits.h>
#include <stdlib.h>

#include "audio.h"
#include "gamepad.h"
#include "gfx-pixslice.h"
#include "lcd.h"
#include "math-util.h"
#include "pam8019.h"
#include "text.h"

#ifdef DEMO
    #include "demo.h"
#else
    #include "synth.h"
#endif

#define GFX_PACK565(r, g, b) \
    (((r) & 0xF8) << (11 - 3) | (((g) & 0xFC) << (5 - 2)) | (((b) & 0xF8) >> 3))

static const gfx_rgb565 bg_color = GFX_PACK565(0, 0, 0); // black
static const gfx_rgb565 fg_color = GFX_PACK565(0, 0, 0xFF); // blue

#define NO_NOTE (-1)

static const int pad_to_note[16] = {
    [GAMEPAD_BLEFT_BIT]   =  0,  // left   => C
    [GAMEPAD_BUP_BIT]     =  1,  // up     => C#
    [GAMEPAD_BDOWN_BIT]   =  2,  // down   => D
    [GAMEPAD_BRIGHT_BIT]  =  3,  // right  => D#
    [GAMEPAD_BSELECT_BIT] =  4,  // select => E
    [GAMEPAD_BSTART_BIT]  =  5,  // start  => F
    [GAMEPAD_BY_BIT]      =  6,  // Y      => F#
    [GAMEPAD_BB_BIT]      =  7,  // B      => G
    [GAMEPAD_BX_BIT]      =  8,  // X      => G#
    [GAMEPAD_BA_BIT]      =  9,  // A      => A
    [GAMEPAD_BVOLP_BIT]   = 10,  // Vol+   => Bflat
    [GAMEPAD_BVOLM_BIT]   = 11,  // Vol-   => B
    [GAMEPAD_BL_BIT]      = NO_NOTE,
    [GAMEPAD_BR_BIT]      = NO_NOTE,
    [GAMEPAD_BRES1_BIT]   = NO_NOTE,
    [GAMEPAD_BRES2_BIT]   = NO_NOTE,
};

#define X_CENTER (LCD_WIDTH / 2)
#define SPARK_RAD 5
#define BLUE_SPREAD_L 1
#define BLUE_SPREAD_R 1
#define BLUE_SPREAD_U 1
#define BLUE_SPREAD_D 1

#define SEG0_X0  X_CENTER
#define SEG0_Y   (LCD_HEIGHT / 4)
#define SEG0_X1  (LCD_WIDTH * 6 / 7)
#define SEG0_LEN (SEG0_X1 - SEG0_X0)

#define SEG1_XL  (LCD_WIDTH - SEG0_X1)
#define SEG1_XR  SEG0_X1
#define SEG1_Y0  (LCD_HEIGHT * 1 / 4)
#define SEG1_Y1  (LCD_HEIGHT * 3 / 4)
#define SEG1_LEN (SEG1_Y1 - SEG1_Y0)

#define SEG2_XR0 (LCD_WIDTH * 6 / 7)
#define SEG2_XL0 (LCD_WIDTH - SEG2_XR0)
#define SEG2_Y   (LCD_HEIGHT * 3 / 4)
#define SEG2_LEN (SEG2_XR0 - X_CENTER + 1)

static gfx_ipoint spark_points[2];
static gfx_rgb565 bright_pix;
static int anim_frame;
static int anim_seg;
static int seg0_xmin;
static int seg0_xmax;
static int seg1_ymax;
static int seg2_xlmax;
static int seg2_xrmin;
static size_t spark_count;
static gfx_rgb565 txt_color;

#ifdef DEMO
    DEFINE_AUDIO_BUFFER(abuf, DEMO_BEST_FRAME_COUNT, ACC_STEREO, ASD_12BIT);
    static demo_state d_state;
#else
    DEFINE_AUDIO_BUFFER(abuf, SYNTH_BEST_FRAME_COUNT, ACC_STEREO, ASD_12BIT);
#endif

static void audio_callback(void *buf, size_t frame_count)
{
#ifdef DEMO
    assert(frame_count == DEMO_BEST_FRAME_COUNT);
    demo_frame *frames = buf;
    demo_gen_samples(&d_state, frames, frame_count);
#else
    assert(frame_count == SYNTH_BEST_FRAME_COUNT);
    synth_frame *frames = buf;
    synth_gen_samples(frames, frame_count);
#endif
}

void deep_app_init(void)
{
#ifdef DEMO
    demo_init_state(&d_state);
    audio_init(DEMO_SAMPLE_RATE, ACC_STEREO, ASD_12BIT, abuf, sizeof abuf);
#else
    synth_init_state();
    audio_init(SYNTH_SAMPLE_RATE, ACC_STEREO, ASD_12BIT, abuf, sizeof abuf);
#endif
    audio_register_callback(audio_callback);
    pam8019_set_mode(PM_NORMAL);
    audio_start();
    anim_frame = 0;
    anim_seg = 0;
}

void deep_app_end(void)
{
    audio_stop();
    pam8019_set_mode(PM_SHUTDOWN);
}

void deep_animate(void)
{
#ifndef DEMO
    // Read buttons.
    if (anim_seg > 2) {
        uint16_t curr_gamepad = gamepad_get();
        static uint16_t prev_gamepad;

        for (int pad = 0; pad < 16; pad++) {
            uint16_t mask = 1 << pad;
            int note = pad_to_note[pad];
            if (note == NO_NOTE)
                continue;
            if (curr_gamepad & mask && !(prev_gamepad & mask)) {
                synth_note_on(note);
            } else if (!(curr_gamepad & mask) && prev_gamepad & mask) {
                synth_note_off(note);
            }
        }
        prev_gamepad = curr_gamepad;
    }
#endif

    // video animation
    int seg_len;
    switch (anim_seg) {

    case 0:
        // segment 0: sparks from top-center to top corners
        {
            seg_len = SEG0_LEN;
            int x = anim_frame;
            seg0_xmin = X_CENTER - x;
            seg0_xmax = X_CENTER + x;
            seg1_ymax = SEG1_Y0;
            seg2_xlmax = SEG2_XL0;
            seg2_xrmin = SEG2_XR0;
            txt_color = 0;
            spark_count = 2;
            spark_points[0].x = seg0_xmin;
            spark_points[0].y = SEG0_Y;
            spark_points[1].x = seg0_xmax;
            spark_points[1].y = SEG0_Y;
        }
        break;

    case 1:
        // segment1: sparks run down sides.
        {
            seg_len = SEG1_LEN;
            int y = anim_frame;
            seg1_ymax = SEG1_Y0 + y;
            spark_count = 2;
            spark_points[0].x = SEG1_XL;
            spark_points[0].y = seg1_ymax;
            spark_points[1].x = SEG1_XR;
            spark_points[1].y = seg1_ymax;
        }
        break;

    case 2:
        // segment 2: sparks from bottom corners to bottom center
        {
            seg_len = SEG2_LEN;
            int x = anim_frame;
            seg2_xlmax = SEG2_XL0 + x;
            seg2_xrmin = SEG2_XR0 - x;
            spark_count = 2;
            spark_points[0].x = seg2_xlmax;
            spark_points[0].y = SEG2_Y;
            spark_points[1].x = seg2_xrmin;
            spark_points[1].y = SEG2_Y;
        }
        break;

    case 3:
        // segment 3: fade in text
        {
            seg_len = 16 + 20;
            int b = anim_frame;
            b -= 20;
            if (b < 0)
                b = 0;
            else {
                b *= b;
                if (b > 255)
                    b = 255;
            }
            txt_color = GFX_PACK565(b, b, b);
        }
        break;

    default:
        seg_len = INT_MAX;
        spark_count = 0;
        break;
    }
    if (++anim_frame >= seg_len) {
        anim_seg++;
        anim_frame = 0;
    }

    int ill = (rand() + 127) & 255;
    bright_pix = GFX_PACK565(ill, ill, ill);
}

static void deep_render_slice(gfx_pixslice *slice)
{
    // top blue line
    for (int y = SEG0_Y - BLUE_SPREAD_U; y < SEG0_Y + BLUE_SPREAD_D; y++) {
        gfx_rgb565 *line = gfx_pixel_address(slice, 0, y);
        if (line) {
            for (int x = seg0_xmin; x < seg0_xmax; x++)
                line[x] = fg_color;
        }
    }

    // side blue lines
    for (int y = SEG1_Y0; y < seg1_ymax; y++) {
        gfx_rgb565 *line = gfx_pixel_address(slice, 0, y);
        if (line) {
            for (int x = SEG1_XL - BLUE_SPREAD_L; x < SEG1_XL + BLUE_SPREAD_R; x++)
                line[x] = fg_color;
            for (int x = SEG1_XR - BLUE_SPREAD_L; x < SEG1_XR + BLUE_SPREAD_R; x++)
                line[x] = fg_color;
        }
    }

    // bottom blue lines
    for (int y = SEG2_Y - BLUE_SPREAD_U; y < SEG2_Y + BLUE_SPREAD_D; y++) {
        gfx_rgb565 *line = gfx_pixel_address(slice, 0, y);
        if (line) {
            for (int x = SEG2_XL0; x < seg2_xlmax; x++)
                line[x] = fg_color;
            for (int x = seg2_xrmin; x < SEG2_XR0; x++)
                line[x] = fg_color;
        }
    }

    // text
    if (txt_color) {
        const char msg[] = "Press and hold any key";
        const size_t nc = sizeof msg - 1;
        int x = (LCD_WIDTH - nc * 8) / 2 / 8;
        int y = LCD_HEIGHT * 8/15 / 16;
        text_draw_str16(slice, msg, x, y, txt_color);
    }

    // sparks
    int size = rand() % (SPARK_RAD - 1) + 1;
    if (anim_seg < 3) {
        for (int i = 0; i < spark_count; i++) {
            int cx = spark_points[i].x;
            int cy = spark_points[i].y;
            for (int y = -size; y <= +size; y++) {
                gfx_rgb565 *line = gfx_pixel_address(slice, 0, cy + y);
                if (line) {
                    int ay = abs(y);
                    for (int x = -SPARK_RAD + ay; x < SPARK_RAD - ay; x++) {
                        line[cx + x] = bright_pix;
                    }
                }
            }
        }
    }
}

void deep_render(void)
{
    if (lcd_bg_color() != bg_color) {
        lcd_set_bg_color(bg_color, true);
    }

    for (size_t y = 0, h; y < LCD_HEIGHT; y += h) {
        h = MIN(LCD_MAX_SLICE_ROWS, LCD_HEIGHT - y);
        gfx_pixslice *slice = lcd_alloc_pixslice(0, y, LCD_WIDTH, h);
        deep_render_slice(slice);
        lcd_send_pixslice(slice);
    }
}
