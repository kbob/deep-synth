#ifndef DEMO_included
#define DEMO_included

#include <stddef.h>
#include <stdint.h>

#define DEMO_SAMPLE_RATE      44100
#define DEMO_BEST_FRAME_COUNT    87 /* Trust me... */

typedef struct demo_state {
    void  *gstate;
    size_t frame_count;
} demo_state;

typedef uint16_t demo_frame[2];

extern void demo_init_state  (demo_state *);
extern void demo_gen_samples (demo_state *, demo_frame *, size_t);

#endif /* !DEMO_included */
