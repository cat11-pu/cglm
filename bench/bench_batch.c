/*
 * Benchmark + bit-exactness check for glm_mat4_mulv_batch / glm_mat3_mulv_batch
 * versus per-vector glm_mat4_mulv / glm_mat3_mulv.
 *
 * Build (from repo root):
 *   cc -O2 -Iinclude bench/bench_batch.c -o bench/bench_batch -lm
 *
 * All buffers are allocated once up front; nothing is allocated inside
 * timed regions or by the batch functions themselves.
 */

#include <cglm/cglm.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NRUNS 7

static double
now_ms(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec * 1.0e3 + (double)ts.tv_nsec * 1.0e-6;
}

/* deterministic xorshift32 so runs are reproducible */
static uint32_t rng_state = 0x9E3779B9u;

static float
frand(void) {
  uint32_t x = rng_state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  rng_state = x;
  return (float)(x >> 8) * (1.0f / 16777216.0f) * 2.0f - 1.0f;
}

/* reference: one matrix-vector multiply per point */
static void
per_point_mat4(CGLM_CONST mat4 m, const vec4 *src, vec4 *dst, size_t n) {
  size_t i;
  for (i = 0; i < n; i++) {
    glm_mat4_mulv(m, src[i], dst[i]);
  }
}

static void
per_point_mat3(CGLM_CONST mat3 m, const vec3 *src, vec3 *dst, size_t n) {
  size_t i;
  for (i = 0; i < n; i++) {
    glm_mat3_mulv(m, src[i], dst[i]);
  }
}

typedef void (*bench_fn)(void *ctx);

static double
best_of(bench_fn fn, void *ctx) {
  double best = 1.0e30;
  int    r;
  fn(ctx); /* warmup */
  for (r = 0; r < NRUNS; r++) {
    double t0 = now_ms();
    fn(ctx);
    double t1 = now_ms();
    if (t1 - t0 < best) {
      best = t1 - t0;
    }
  }
  return best;
}

typedef struct {
  mat4        m;
  const vec4 *src;
  vec4       *dst;
  size_t      n;
} ctx4_t;

typedef struct {
  mat3        m;
  const vec3 *src;
  vec3       *dst;
  size_t      n;
} ctx3_t;

static void run_pp4(void *c) { ctx4_t *x = c; per_point_mat4(x->m, x->src, x->dst, x->n); }
static void run_b4(void *c)  { ctx4_t *x = c; glm_mat4_mulv_batch(x->m, x->src, x->dst, x->n); }
static void run_pp3(void *c) { ctx3_t *x = c; per_point_mat3(x->m, x->src, x->dst, x->n); }
static void run_b3(void *c)  { ctx3_t *x = c; glm_mat3_mulv_batch(x->m, x->src, x->dst, x->n); }

int
main(void) {
  static const size_t sizes[] = {100000, 1000000, 4000000, 16000000};
  const size_t        nsize   = sizeof(sizes) / sizeof(sizes[0]);
  const size_t        maxn    = sizes[nsize - 1];
  size_t              si, i;
  int                 j, fails;

  vec4 *src4  = aligned_alloc(64, maxn * sizeof(vec4));
  vec4 *dst4  = aligned_alloc(64, maxn * sizeof(vec4));
  vec4 *ref4  = aligned_alloc(64, maxn * sizeof(vec4));
  vec3 *src3  = malloc(maxn * sizeof(vec3));
  vec3 *dst3  = malloc(maxn * sizeof(vec3));
  vec3 *ref3  = malloc(maxn * sizeof(vec3));

  if (!src4 || !dst4 || !ref4 || !src3 || !dst3 || !ref3) {
    fprintf(stderr, "allocation failed\n");
    return 1;
  }

  mat4 m4;
  mat3 m3;
  for (j = 0; j < 4; j++) {
    for (i = 0; i < 4; i++) {
      m4[j][i] = frand();
    }
  }
  glm_mat4_pick3(m4, m3);

  for (i = 0; i < maxn; i++) {
    for (j = 0; j < 4; j++) {
      src4[i][j] = frand();
    }
    for (j = 0; j < 3; j++) {
      src3[i][j] = frand();
    }
  }

  fails = 0;

  for (si = 0; si < nsize; si++) {
    size_t n = sizes[si];
    double t_pp4, t_b4, t_b4i, t_pp3, t_b3, t_b3i;
    int    ok4, ok4i, ok3, ok3i;

    ctx4_t c4 = {.src = src4, .dst = dst4, .n = n};
    ctx3_t c3 = {.src = src3, .dst = dst3, .n = n};
    glm_mat4_copy(m4, c4.m);
    glm_mat3_copy(m3, c3.m);

    /* bit-exactness: batch vs per-point, out-of-place and in-place */
    per_point_mat4(m4, src4, ref4, n);
    glm_mat4_mulv_batch(m4, src4, dst4, n);
    ok4 = memcmp(ref4, dst4, n * sizeof(vec4)) == 0;
    memcpy(dst4, src4, n * sizeof(vec4));
    glm_mat4_mulv_batch(m4, dst4, dst4, n);
    ok4i = memcmp(ref4, dst4, n * sizeof(vec4)) == 0;

    per_point_mat3(m3, src3, ref3, n);
    glm_mat3_mulv_batch(m3, src3, dst3, n);
    ok3 = memcmp(ref3, dst3, n * sizeof(vec3)) == 0;
    memcpy(dst3, src3, n * sizeof(vec3));
    glm_mat3_mulv_batch(m3, dst3, dst3, n);
    ok3i = memcmp(ref3, dst3, n * sizeof(vec3)) == 0;

    fails += !ok4 + !ok4i + !ok3 + !ok3i;

    t_pp4 = best_of(run_pp4, &c4);
    t_b4  = best_of(run_b4,  &c4);
    memcpy(dst4, src4, n * sizeof(vec4));
    c4.src = dst4; /* in-place: src == dst */
    t_b4i = best_of(run_b4,  &c4);
    c4.src = src4;

    t_pp3 = best_of(run_pp3, &c3);
    t_b3  = best_of(run_b3,  &c3);
    memcpy(dst3, src3, n * sizeof(vec3));
    c3.src = dst3;
    t_b3i = best_of(run_b3,  &c3);
    c3.src = src3;

    printf("\nn = %-9zu\n", n);
    printf("  mat4  per-point %8.3f ms | batch %8.3f ms | in-place %8.3f ms"
           " | speedup %5.2fx | per-1M %6.3f ms | bit-exact oop:%s inplace:%s\n",
           t_pp4, t_b4, t_b4i, t_pp4 / t_b4, t_b4 * 1.0e6 / (double)n,
           ok4 ? "PASS" : "FAIL", ok4i ? "PASS" : "FAIL");
    printf("  mat3  per-point %8.3f ms | batch %8.3f ms | in-place %8.3f ms"
           " | speedup %5.2fx | per-1M %6.3f ms | bit-exact oop:%s inplace:%s\n",
           t_pp3, t_b3, t_b3i, t_pp3 / t_b3, t_b3 * 1.0e6 / (double)n,
           ok3 ? "PASS" : "FAIL", ok3i ? "PASS" : "FAIL");
  }

  free(src4); free(dst4); free(ref4);
  free(src3); free(dst3); free(ref3);

  if (fails) {
    printf("\nBIT-EXACTNESS: %d check(s) FAILED\n", fails);
    return 1;
  }
  printf("\nBIT-EXACTNESS: all checks passed\n");
  return 0;
}
