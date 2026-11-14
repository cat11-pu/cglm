/*
 * Randomized property tests for cglm.
 *
 * Each property is exercised with at least PROP_MIN_ITERS (100k)
 * pseudo-random inputs.  On failure the smallest (lowest iteration
 * index) counter-example is recorded and printed together with the
 * measured error.
 */

#ifndef cglm_test_property_common_h
#define cglm_test_property_common_h

#include <cglm/cglm.h>

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PROP_MIN_ITERS 100000

#ifndef PROP_SEED
#  define PROP_SEED 0xC61B5EEDu
#endif

/*
 * Deterministic xorshift32 PRNG: the same seed always reproduces the
 * same counter-example, which keeps failures debuggable.
 */
typedef struct prop_rng_t {
  uint32_t state;
} prop_rng_t;

static inline uint32_t
prop_rng_u32(prop_rng_t *rng) {
  uint32_t x = rng->state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  rng->state = x;
  return x;
}

/* uniform float in [-1, 1] */
static inline float
prop_rng_f(prop_rng_t *rng) {
  return (float)((int32_t)prop_rng_u32(rng)) / (float)INT32_MAX;
}

/* uniform float in [0, 1] */
static inline float
prop_rng_uf(prop_rng_t *rng) {
  return (float)(prop_rng_u32(rng) >> 1) / (float)(1u << 31);
}

/* uniform float in [lo, hi] */
static inline float
prop_rng_range(prop_rng_t *rng, float lo, float hi) {
  return lo + (hi - lo) * prop_rng_uf(rng);
}

/*
 * magnitude bucket: lets every test deliberately cover the edge
 * regimes requested by the task (tiny / normal / large scale).
 *
 *   bucket 0: subnormal-to-tiny,  [1e-38, 1e-30]
 *   bucket 1: very small,          [1e-12, 1e-6]
 *   bucket 2: ordinary,            [1e-3, 1e3]
 *   bucket 3: large,               [1e3, 1e6]
 */
static const float PROP_MAG_LO[4] = {1e-38f, 1e-12f, 1e-3f, 1e3f};
static const float PROP_MAG_HI[4] = {1e-30f, 1e-6f,  1e3f,  1e6f};

static inline float
prop_rng_magnitude(prop_rng_t *rng) {
  uint32_t bucket = prop_rng_u32(rng) & 3u;
  float lo = PROP_MAG_LO[bucket];
  float hi = PROP_MAG_HI[bucket];
  return lo * powf(hi / lo, prop_rng_uf(rng));
}

static inline int
prop_finite(float v) {
  return isfinite(v);
}

static inline float
prop_maxf(float a, float b) {
  return a > b ? a : b;
}

/*
 * Counter-example bookkeeping.  A property keeps one of these; the
 * first (smallest iteration index) failing input is retained so the
 * report is stable and minimal.
 */
typedef struct prop_fail_t {
  int      failed;
  uint32_t iter;
  float    error;
  float    tol;
} prop_fail_t;

static inline void
prop_fail_init(prop_fail_t *f) {
  f->failed = 0;
  f->iter   = 0;
  f->error  = 0.0f;
  f->tol    = 0.0f;
}

/* returns 1 while the property should keep running after a failure */
static inline int
prop_fail_record(prop_fail_t *f, uint32_t iter, float error, float tol) {
  if (error <= tol)
    return 1;

  if (!f->failed || iter < f->iter) {
    f->failed = 1;
    f->iter   = iter;
    f->error  = error;
    f->tol    = tol;
  }
  return 0;
}

static inline void
prop_print_vec(const char *name, const float *v, int n) {
  fprintf(stderr, "      %s = [", name);
  for (int i = 0; i < n; i++)
    fprintf(stderr, "%s%.9g", i ? ", " : "", (double)v[i]);
  fprintf(stderr, "]\n");
}

static inline void
prop_print_mat3(const char *name, mat3 m) {
  fprintf(stderr, "      %s =\n", name);
  for (int r = 0; r < 3; r++) {
    fprintf(stderr, "        [");
    for (int c = 0; c < 3; c++)
      fprintf(stderr, "%s%.9g", c ? ", " : "", (double)m[c][r]);
    fprintf(stderr, "]\n");
  }
}

static inline void
prop_print_mat4(const char *name, mat4 m) {
  fprintf(stderr, "      %s =\n", name);
  for (int r = 0; r < 4; r++) {
    fprintf(stderr, "        [");
    for (int c = 0; c < 4; c++)
      fprintf(stderr, "%s%.9g", c ? ", " : "", (double)m[c][r]);
    fprintf(stderr, "]\n");
  }
}

typedef struct prop_result_t {
  const char *name;
  int         passed;
  uint32_t    iters;
  float       max_error;
  prop_fail_t fail;
} prop_result_t;

#endif /* cglm_test_property_common_h */
