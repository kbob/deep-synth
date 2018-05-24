#ifndef SYNTH_included
#define SYNTH_included

#include <stddef.h>
#include <stdint.h>

#define SYNTH_SAMPLE_RATE      44100
#define SYNTH_BEST_FRAME_COUNT    88 // 500 Hz

typedef uint16_t synth_frame[2];

extern void synth_init_state  (void);
extern void synth_gen_samples (synth_frame *, size_t);

extern void synth_note_on     (int note);
extern void synth_note_off    (int note);

#endif /* !SYNTH_included */
