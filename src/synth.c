#include "synth.h"

//#define HURRY_UP 8

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "intr.h"

// XXX I like the part where all the configuration parameters
// are in a big struct.  Makes it easy to find and track them.

// XXX todo: randomize noise generators' starting phase.

#define VOICE_COUNT      3
#define NOTE_COUNT       12 
#define MIDI_NOTE_COUNT 128
#define NOTES_per_OCTAVE 12

#define Fs           ((float)SYNTH_SAMPLE_RATE)
#define Fk           500
#define MAX_ENV_SEGS 3
#define NO_NOTE      (-1)

// Voice Parameters
#define V_RAND_FREQ_MIN    200
#define V_RAND_FREQ_MAX    400
#define V_RAND_NOISE_FREQ  0.5
#define V_RAND_NOISE_MAX   180.0
#define V_DEST_OCTAVES     5
#define V_DEST_NOISE_FREQ  0.1
#define V_DEST_NOISE_MAX   10.0

#define V_SWEEP_FOCUS_DUR_MIN 3.0
#define V_SWEEP_FOCUS_DUR_MAX 4.0
#define V_SWEEP_FOCUS_EXP_MIN 3.0
#define V_SWEEP_FOCUS_EXP_MAX 5.0

// #define V_SWEEP_SLIDE_DUR_MIN 0.8
// #define V_SWEEP_SLIDE_DUR_MAX 1.2
// #define V_SWEEP_SLIDE_EXP_MIN 1.0
// #define V_SWEEP_SLIDE_EXP_MAX 2.0

#define V_SWEEP_BLUR_DUR_MIN  3.0
#define V_SWEEP_BLUR_DUR_MAX  3.5
#define V_SWEEP_BLUR_EXP_MIN  3.0
#define V_SWEEP_BLUR_EXP_MAX  5.0

#define V_LOWPASS_H        6.0  // number of harmonics
#define V_LOWPASS_Q        1.67
#define V_PAN_LEFT         (-0.5)
#define V_PAN_RIGHT        (+0.5)

// Global Parameters
#define G_LOWPASS_FREQ_MIN 2000.0
#define G_LOWPASS_FREQ_MAX 20000.0
#define G_LOWPASS_DUR1     8.0
#define G_LOWPASS_CURVE1   2.0
#define G_LOWPASS_Q        2.0
#define G_AMP_RISE_DUR     4.0
#define G_AMP_CURVE1       2.0

static float midi_to_freq[MIDI_NOTE_COUNT];

typedef float frame[2];

typedef struct noise_cfg {
    float             freq;
    float             amp;
} noise_cfg;

typedef struct env_cfg {
    size_t            seg_count;
    float             levels[MAX_ENV_SEGS + 1];
    float             durations[MAX_ENV_SEGS];
    float             curves[MAX_ENV_SEGS];
} env_cfg;

typedef struct sweeper_cfg {
    float             focus_dur;
    float             focus_exp;
    float             blur_dur;
    float             blur_exp;
} sweeper_cfg;

typedef struct filter_cfg {
    float             q;
} filter_cfg;

typedef struct voice_cfg {
    float             rand_freq;
    noise_cfg         rand_noise;
    noise_cfg         dest_noise;
    sweeper_cfg       sweep;
    filter_cfg        lowpass;
    float             amp;
    float             pan;
} voice_cfg;

typedef struct global_cfg {
    voice_cfg         voices[VOICE_COUNT];
    env_cfg           lowpass_env;
    filter_cfg        lowpass;
    env_cfg           amp_env;
} global_cfg;

//  --  //  --  //  --  //  --  //  --  //  --  //  --  //  --  //  --  //  --

typedef struct note_state {
    bool              is_on;
} note_state;

typedef struct noise_state {
    float             level;
    float             next_midpoint;
    float             next_level;
    float             slope;
    float             curve;
    int               count;
} noise_state;

typedef struct env_state {
    int               frame_number;
    int               end_fn;
    ssize_t           seg_index;
    float             curve;
    float             value;
    float             a2;
    float             b1;
    float             grow;
} env_state;

typedef enum direction {
    D_FOCUSING,
    D_BLURRING,
} direction;

typedef struct sweeper_state {
    enum direction    direction;
    float             value;
    float             delta;
    float             grow;
} sweeper_state;

typedef struct osc_state {
    float             phase;
    float             inc;
} osc_state;

typedef struct filter_state {
    int               samples;
    float             f;
    float             q;
    float             s0;
    float             s1;
} filter_state;


typedef struct voice_state {
    int               dest_note;
    float             dest_freq;
    int               prev_note;
    float             prev_freq;
    noise_state       rand_noise;
    noise_state       dest_noise;
    sweeper_state     sweep;
    osc_state         osc;
    filter_state      lowpass;
    float             left_gain;
    float             right_gain;
} voice_state;

typedef struct global_state {
    const global_cfg *cfg;
    note_state        notes[NOTE_COUNT];
    voice_state       voices[VOICE_COUNT];
    env_state         lowpass_env;
    filter_state      lowpass[2]; // stereo
    env_state         amp_env;
} global_state;

//  --  //  --  //  --  //  --  //  --  //  --  //  --  //  --  //  --  //  --

static global_cfg gcfg;
static global_state gstate;
static size_t saved_frame_count;
static unsigned active_note_count;

// MIDI note is a float because the original Deep Note was tune
//  halfway between D and E-Flat.
static float midicps(float midi_note)
{
    // Ref. to A4 (note 69) == 440 Hz.
    return 440.0 * powf(2.0, (midi_note - 69) / 12.0);
}

static void init_osc_tables(void)
{
    for (size_t i = 0; i < MIDI_NOTE_COUNT; i++)
        midi_to_freq[i] = midicps(i);
}

static float rrand(float a, float b)
{
    return a + drand48() * (b - a);
}

static int compare_voices(const void *v0, const void *v1)
{
    float f0 = ((const voice_cfg *)v0)->rand_freq;
    float f1 = ((const voice_cfg *)v1)->rand_freq;
    // N.B., sorts in descending order.
    if (f0 < f1)
        return +1;
    else if (f0 > f1)
        return -1;
    return 0;
}

static void init_cfg(global_cfg *g)
{
    memset(g, 0, sizeof *g);

    // Choose random frequencies, then sort in descending order.
    for (size_t i = 0; i < VOICE_COUNT; i++)
        g->voices[i].rand_freq = rrand(V_RAND_FREQ_MIN, V_RAND_FREQ_MAX);
    qsort(&g->voices, VOICE_COUNT, sizeof g->voices[0], compare_voices);

    for (size_t i = 0; i < VOICE_COUNT; i++) {
        voice_cfg *v = &g->voices[i];

        v->rand_noise.freq = V_RAND_NOISE_FREQ;
        v->rand_noise.amp = (VOICE_COUNT - i) * V_RAND_NOISE_MAX / VOICE_COUNT;

        v->dest_noise.freq = V_DEST_NOISE_FREQ;
        v->dest_noise.amp = (i + 1) * V_DEST_NOISE_MAX / VOICE_COUNT;

        v->sweep.focus_dur =
            rrand(V_SWEEP_FOCUS_DUR_MIN, V_SWEEP_FOCUS_DUR_MAX);
        v->sweep.focus_exp =
            rrand(V_SWEEP_FOCUS_EXP_MIN, V_SWEEP_FOCUS_EXP_MAX);
        v->sweep.blur_dur = rrand(V_SWEEP_BLUR_DUR_MIN, V_SWEEP_BLUR_DUR_MAX);
        v->sweep.blur_exp = rrand(V_SWEEP_BLUR_EXP_MIN, V_SWEEP_BLUR_EXP_MAX);

        v->lowpass.q = V_LOWPASS_Q;

        // 1/2, 2/3, 3/4, ... from lowest to highest
        v->amp = 1.0 - 1.0 / (i + 2) / VOICE_COUNT;

        v->pan = rrand(V_PAN_LEFT, V_PAN_RIGHT);

#ifdef HURRY_UP
        v->sweep.focus_dur /= HURRY_UP;
        v->sweep.blur_dur  /= HURRY_UP;
#endif /* HURRY_UP */
    }

    g->lowpass.q = V_LOWPASS_Q;

    g->lowpass_env.seg_count = 1;
    g->lowpass_env.levels[0] = 0;
    g->lowpass_env.levels[1] = 1;
    g->lowpass_env.durations[0] = G_LOWPASS_DUR1;
    g->lowpass_env.curves[0] = G_LOWPASS_CURVE1;

    g->lowpass.q = G_LOWPASS_Q;

    g->amp_env.seg_count = 1;
    g->amp_env.levels[0] = 0.0;
    g->amp_env.levels[1] = 1.0;
    g->amp_env.durations[0] = G_AMP_RISE_DUR;
    g->amp_env.curves[0] = G_AMP_CURVE1;

#ifdef HURRY_UP
    g->lowpass_env.durations[0] /= HURRY_UP;
    g->lowpass_env.durations[1] /= HURRY_UP;
    g->amp_env.durations[0] /= HURRY_UP;
    g->amp_env.durations[1] /= HURRY_UP;
    g->amp_env.durations[2] /= HURRY_UP;
#endif
}

//  --  //  --  //  --  //  --  //  --  //  --  //  --  //  --  //  --  //  --

static void init_noise_state(noise_state *nstate, const noise_cfg *ncfg)
{
    float amp = ncfg->amp;
    float next_level = rrand(-amp, +amp);

    nstate->level = 0;
    nstate->next_level = next_level;
    nstate->next_midpoint = next_level / 2;
    nstate->slope = 0;
    nstate->curve = 0;
    nstate->count = 0;
}

static void init_env_state(env_state *estate, const env_cfg *ecfg)
{
    estate->frame_number = 0;
    estate->end_fn = 0;
    estate->seg_index = -1;
    estate->curve = 0.0;
    estate->value = ecfg->levels[0];
}

static void init_sweeper_state(sweeper_state *sstate, const sweeper_cfg *scfg)
{
    sstate->direction = D_BLURRING;
    sstate->value = 0;
    sstate->delta = 0;
    sstate->grow = 1;
}

static void init_filter_state(filter_state *fstate, const filter_cfg *fcfg)
{
    fstate->q = 1 / fcfg->q;
    fstate->s0 = 0;
    fstate->s1 = 0;
}

static void init_state(global_state *state, const global_cfg *cfg)
{
    memset(state, 0, sizeof *state);
    state->cfg = cfg;
    for (size_t i = 0; i < VOICE_COUNT; i++) {
        voice_state *vstate = &state->voices[i];
        const voice_cfg *vcfg = &cfg->voices[i];
        init_noise_state(&vstate->rand_noise, &vcfg->rand_noise);
        init_noise_state(&vstate->dest_noise, &vcfg->dest_noise);
        init_sweeper_state(&vstate->sweep, &vcfg->sweep);
        init_filter_state(&vstate->lowpass, &vcfg->lowpass);
        float pos = (vcfg->pan + 1) * M_PI_4;
        vstate->left_gain = vcfg->amp * cosf(pos) / VOICE_COUNT;
        vstate->right_gain = vcfg->amp * sinf(pos) / VOICE_COUNT;
    }
    init_env_state(&state->lowpass_env, &cfg->lowpass_env);
    init_filter_state(&state->lowpass[0], &cfg->lowpass);
    init_filter_state(&state->lowpass[1], &cfg->lowpass);
    init_env_state(&state->amp_env, &cfg->amp_env);
}

static void update_noise_state(noise_state *nstate, const noise_cfg *ncfg)
{
    int count = nstate->count;
    float slope = nstate->slope;
    if (count <= 0) {
        float amp = ncfg->amp;
        float curr_level = nstate->next_level;
        float next_level = rrand(-amp, +amp);
        float level = nstate->next_midpoint;
        float midpoint = (next_level + curr_level) * 0.5;
        nstate->next_midpoint = midpoint;
        nstate->next_level = next_level;
        count = Fk / ncfg->freq;
        if (count < 2)
            count = 2;
        float fcount = (float)count;
        float curve_numer = 2 * (midpoint - level - fcount * slope);
        float curve_denom = fcount * fcount + fcount;
        nstate->curve = curve_numer / curve_denom;
        nstate->level = level;
    }
    nstate->count = count - 1;
    nstate->slope = slope += nstate->curve;
    nstate->level += slope;
}

static void update_env_state(env_state *estate, const env_cfg *ecfg)
{
    int fn = estate->frame_number;
    int seg = estate->seg_index;
    if (fn >= estate->end_fn) {
        seg++;
        if (seg >= (ssize_t)ecfg->seg_count) {
            // Stay at final level forever.
            estate->frame_number = 0;
            estate->end_fn = INT_MAX;
            estate->seg_index = seg = ecfg->seg_count;
            estate->curve = 0;
            estate->value = ecfg->levels[seg];
            estate->grow = 0;
        } else {
            estate->seg_index = seg;
            estate->frame_number = fn -= estate->end_fn;
            estate->end_fn = ecfg->durations[seg] * Fk;
            float level0 = ecfg->levels[seg];
            float level1 = ecfg->levels[seg + 1];
            float seg_dur = ecfg->durations[seg];
            float curve = ecfg->curves[seg];
            if (fabsf(curve) < 0.001) {
                estate->curve = 0.0;
                estate->value = level0;
                estate->grow = (level1 - level0) / (Fk * seg_dur);
            } else {
                double a1 = (level1 - level0) / (1.0 - expf(curve));
                estate->curve = curve;
                estate->a2 = level0 + a1;
                estate->b1 = a1;
                estate->grow = expf(estate->curve / (Fk * seg_dur));
            }
        }
    }
    if (estate->curve == 0.0) {
        estate->value += estate->grow;
    } else {
        double a2 = estate->a2;
        double b1 = estate->b1;
        double grow = estate->grow;
        b1 *= grow;
        estate->value = a2 - b1;
        estate->b1 = b1;
    }
    estate->frame_number++;
}

static void update_sweeper_state(sweeper_state     *sstate,
                                 sweeper_cfg const *scfg)
{

    switch (sstate->direction) {

        case D_FOCUSING:
            if (sstate->value >= 1) {
                sstate->value = 1;
                return;
            }
            break;

        case D_BLURRING:
            if (sstate->value <= 0) {
                sstate->value = 0;
                return;
            }
            break;
    }

    float value = sstate->value;
    float delta = sstate->delta;
    value += delta;
    delta *= sstate->grow;
    sstate->value = value;
    sstate->delta = delta;
}

static float sweeper_interpolate(sweeper_state const *sstate,
                                 float                blur_freq,
                                 float                dest_freq)
{
    float sweep = sstate->value;
    return sweep * dest_freq + (1 - sweep) * blur_freq;
}

static void sweeper_focus(sweeper_state *sstate, sweeper_cfg const *scfg)
{
    float dur = scfg->focus_dur;
    float exp = scfg->focus_exp;
    float value = sstate->value;
    float kdur = Fk * dur;
    float c = expf(exp) - 1;
    if (fabsf(c) < 0.005) {
        sstate->delta = 1.0 / kdur;
        sstate->grow = 1.0;
    } else {
        sstate->delta = exp * expf(exp * value) / c / kdur;
        sstate->grow = expf(exp / kdur);
    }
    sstate->direction = D_FOCUSING;
}

static void sweeper_blur(sweeper_state *sstate, sweeper_cfg const *scfg)
{
    float dur = scfg->blur_dur;
    float exp = scfg->blur_exp;
    float value = sstate->value;
    float kdur = Fk * dur;
    float c = expf(exp) - 1;
    if (fabsf(c) < 0.005) {
        sstate->delta = -1.0 / kdur;
        sstate->grow = 1.0;
    } else {
        sstate->delta = -exp * expf(exp * (1.0 - value)) / c / kdur;
        sstate->grow = expf(exp / kdur);
    }
    sstate->direction = D_BLURRING;
}

static void update_filter_state(filter_state     *fstate,
                                filter_cfg const *fcfg,
                                float             freq)
{
    // There is no theoretical basis for this calculation of the
    // oversampling rate, but it seems to work for the values of Q
    // we're using.  See Tim Stilson's PhD thesis pp 97-106 for more
    // info.
    //
    // https://ccrma.stanford.edu/~stilti/papers/TimStilsonPhDThesis2006.pdf

    int samples = (int)(freq * ((2 * M_PI / Fs) + (3 / Fs))) + 1;
    fstate->samples = samples;
    fstate->f = 2 * sinf(M_PI * freq / samples / Fs);
}

static void update_state(global_state *state)
{
    const global_cfg *cfg = state->cfg;
    for (size_t i = 0; i < VOICE_COUNT; i++) {
        voice_state *vstate = &state->voices[i];
        const voice_cfg *vcfg = &cfg->voices[i];
        update_noise_state(&vstate->rand_noise, &vcfg->rand_noise);
        update_noise_state(&vstate->dest_noise, &vcfg->dest_noise);
        update_sweeper_state(&vstate->sweep, &vcfg->sweep);
        float blur_freq = vcfg->rand_freq + vstate->rand_noise.level;
        // float prev_freq = vstate->prev_freq + vstate->dest_noise.level;
        float dest_freq = vstate->dest_freq + vstate->dest_noise.level;
        float freq = sweeper_interpolate(&vstate->sweep, blur_freq, dest_freq);
        float filt_freq = freq * V_LOWPASS_H;
        update_filter_state(&vstate->lowpass, &vcfg->lowpass, filt_freq);
        vstate->osc.inc = freq / Fs * 2;
    }
    update_env_state(&state->lowpass_env, &cfg->lowpass_env);
    update_env_state(&state->amp_env, &cfg->amp_env);
    float cutoff = G_LOWPASS_FREQ_MIN +
        (G_LOWPASS_FREQ_MAX - G_LOWPASS_FREQ_MIN) * state->lowpass_env.value;
    update_filter_state(&state->lowpass[0], &cfg->lowpass, cutoff);
    update_filter_state(&state->lowpass[1], &cfg->lowpass, cutoff);
}

static void osc_gen_chunk(osc_state *ostate,
                          float     *frames_out,
                          size_t     frame_count)
{
    float phase = ostate->phase;
    const float inc = ostate->inc;

    for (size_t i = 0; i < frame_count; i++) {
        frames_out[i] = phase;
        phase += inc;
        if (phase >= 1.0)
            phase -= 2.0;
    }
    ostate->phase = phase;
}

static void filter_chunk(float        *buf,
                         size_t        frame_count,
                         size_t        stride,
                         filter_state *fstate)
{
    const int samples = fstate->samples;
    const float f = fstate->f;
    const float q = fstate->q;
    float s0 = fstate->s0;
    float s1 = fstate->s1;

    for (size_t i = 0; i < stride * frame_count; i += stride) {
        for (int j = 0; j < samples; j++) {
            s0 += f * (buf[i] - q * s0 - s1);
            s1 += f * s0;
        }
        buf[i] = s1;
    }

    fstate->s0 = s0;
    fstate->s1 = s1;
}

static void filter_mono_chunk(float        *buf,
                              size_t        frame_count,
                              filter_state *fstate)
{
    filter_chunk(buf, frame_count, 1, fstate);
}

static void filter_stereo_chunk(frame       *buf,
                                size_t       frame_count,
                                filter_state fstate[2])
{
    filter_chunk(&buf[0][0], frame_count, 2, &fstate[0]);
    filter_chunk(&buf[0][1], frame_count, 2, &fstate[1]);
}

static void mix_chunk(frame       *frames_out,
                      size_t      frame_count,
                      float       *input,
                      voice_state *vstate)
{
    float left_gain = vstate->left_gain;
    float right_gain = vstate->right_gain;

    for (size_t i = 0; i < frame_count; i++) {
        float sample = input[i];
        frames_out[i][0] += left_gain * sample;
        frames_out[i][1] += right_gain * sample;
    }
}

static void amp_stereo_chunk(frame *buf, size_t frame_count, float amp)
{
    for (size_t i = 0; i < frame_count; i++) {
        buf[i][0] *= amp;
        buf[i][1] *= amp;
    }
}

static void gen_chunk(global_state *gstate,
                      frame        *frames_out,
                      size_t        frame_count)
{
    memset(frames_out, 0, frame_count * sizeof *frames_out);
    update_state(gstate);
    for (size_t i = 0; i < VOICE_COUNT; i++) {
        float voice_buf[frame_count];
        voice_state *vstate = &gstate->voices[i];
        osc_gen_chunk(&vstate->osc, voice_buf, frame_count);
        filter_mono_chunk(voice_buf, frame_count, &vstate->lowpass);
        mix_chunk(frames_out, frame_count, voice_buf, vstate);
    }
    filter_stereo_chunk(frames_out, frame_count, gstate->lowpass);
    float amp = (2 + gstate->lowpass_env.value) * gstate->amp_env.value / 2;
    amp_stereo_chunk(frames_out, frame_count, amp);
}

static void convert_frames(const frame *in, synth_frame *out, size_t count)
{
    const float SCALE = 2047.0f;
    const float OFFSET = 2047.0f;

    for (size_t i = 0; i < count; i++) {
        out[i][0] = (uint16_t)(SCALE * in[i][0] + OFFSET);
        out[i][1] = (uint16_t)(SCALE * in[i][1] + OFFSET);
    }
}

void synth_init_state(void)
{
    init_osc_tables();
    init_cfg(&gcfg);
    init_state(&gstate, &gcfg);
}

void synth_gen_samples(synth_frame *frames_out, size_t frame_count)
{
    static frame frames[SYNTH_BEST_FRAME_COUNT];

    if (frame_count == SYNTH_BEST_FRAME_COUNT && saved_frame_count == 0) {
        gen_chunk(&gstate, frames, SYNTH_BEST_FRAME_COUNT);
        convert_frames(frames, frames_out, frame_count);
    } else {
        assert(false && "untested code!");
    }
}

static void assign_voices(void)
{
    if (active_note_count) {
        int note = 0;
        for (size_t i = 0; i < VOICE_COUNT; i++) {
            while (!gstate.notes[note].is_on)
                note = (note + 1) % NOTE_COUNT;
            // Uniformly distribute voices
            // between C0 (MIDI note 12, ~16.4Hz)
            // and A5 (MIDI note 81, 880Hz).
            int octave = i * V_DEST_OCTAVES * NOTES_per_OCTAVE;
            octave += (23 - note) * VOICE_COUNT;
            octave /= NOTES_per_OCTAVE * VOICE_COUNT;
            note += NOTES_per_OCTAVE * octave;
            assert(12 <= note && note <= 81);
            float freq = midi_to_freq[note];
            voice_state *vstate = &gstate.voices[i];
            WITH_INTERRUPTS_MASKED {
                vstate->prev_note = vstate->dest_note;
                vstate->prev_freq = vstate->dest_freq;
                vstate->dest_note = note;
                vstate->dest_freq = freq;
            }
            note = (note + 1) % NOTE_COUNT;
        }
    } else {
        for (size_t i = 0; i < VOICE_COUNT; i++) {
            voice_state *vstate = &gstate.voices[i];
            WITH_INTERRUPTS_MASKED {
                vstate->prev_note = vstate->dest_note;
                vstate->prev_freq = vstate->dest_freq;
                vstate->dest_note = NO_NOTE;
            }
        }
    }
}

static void focus_voices(void)
{
    for (size_t i = 0; i < VOICE_COUNT; i++) {
        sweeper_state *sstate = &gstate.voices[i].sweep;
        const sweeper_cfg *scfg = &gstate.cfg->voices[i].sweep;
        WITH_INTERRUPTS_MASKED {
            sweeper_focus(sstate, scfg);
        }
    }
}

static void blur_voices(void)
{
    for (size_t i = 0; i < VOICE_COUNT; i++) {
        sweeper_state *sstate = &gstate.voices[i].sweep;
        const sweeper_cfg *scfg = &gstate.cfg->voices[i].sweep;
        WITH_INTERRUPTS_MASKED {
            sweeper_blur(sstate, scfg);
        }
    }
}

void synth_note_on(int note)
{
    note %= NOTE_COUNT;
    bool was_on = gstate.notes[note].is_on;
    gstate.notes[note].is_on = true;
    if (!was_on) {
        if (!active_note_count)
            focus_voices();
        active_note_count++;
        assign_voices();
    }
}

void synth_note_off(int note)
{
    note %= NOTE_COUNT;
    bool was_on = gstate.notes[note].is_on;
    gstate.notes[note].is_on = false;
    if (was_on) {
        active_note_count--;
        assign_voices();
        if (!active_note_count)
            blur_voices();
    }
}
