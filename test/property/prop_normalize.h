#ifndef cglm_test_prop_normalize_h
#define cglm_test_prop_normalize_h

#include "property_common.h"

/*
 * Property 3: vector normalization stability at tiny magnitudes.
 *
 * For vec2 / vec3 / vec4 (both the _to and in-place variants):
 *
 *   - every output component is finite: no NaN / Inf may ever leak
 *   - the result is in exactly one of two valid states:
 *       (a) the zero vector -- the underflow collapse used for
 *           tiny/zero inputs, or
 *       (b) a unit vector parallel to the input (out == v / |v|)
 *   - the in-place variant and the _to variant always agree
 *
 * The exact underflow cut-off differs per backend (scalar paths compare
 * norm against FLT_EPSILON, the SIMD vec4 path compares the squared
 * norm), so the property asserts the observable two-state contract
 * rather than a platform-specific threshold.  Exact zero vectors are
 * injected explicitly and must always land in state (a).
 */

#define NORMALIZE_TOL 1.5e-6f

typedef enum prop_vec_kind_t {
  PROP_VEC2 = 0,
  PROP_VEC3 = 1,
  PROP_VEC4 = 2
} prop_vec_kind_t;

static inline float
prop_vec_norm(const float *v, int n) {
  float acc = 0.0f;
  for (int i = 0; i < n; i++)
    acc += v[i] * v[i];
  return sqrtf(acc);
}

static inline void
prop_random_vec(prop_rng_t *rng, float *v, int n, int force_zero) {
  if (force_zero) {
    for (int i = 0; i < n; i++)
      v[i] = 0.0f;
    return;
  }

  float mag = prop_rng_magnitude(rng);
  for (int i = 0; i < n; i++)
    v[i] = prop_rng_f(rng) * mag;

  /* sometimes make the vector point along a single axis */
  if ((prop_rng_u32(rng) & 7u) == 0) {
    uint32_t axis = prop_rng_u32(rng) % (uint32_t)n;
    for (int i = 0; i < n; i++)
      v[i] = (i == (int)axis)
               ? mag * (prop_rng_f(rng) >= 0.0f ? 1.0f : -1.0f)
               : 0.0f;
  }
}

/*
 * Returns the property error (0 means the contract holds).
 * INFINITY means a non-finite component escaped the library.
 */
static inline float
prop_check_normalize(const float *v, const float *out, int n,
                     int *finite_ok) {
  float norm, out_norm, err;
  int   i;

  *finite_ok = 1;
  for (i = 0; i < n; i++) {
    if (!prop_finite(out[i])) {
      *finite_ok = 0;
      return INFINITY;
    }
  }

  norm = prop_vec_norm(v, n);

  if (norm == 0.0f) {
    /* exact zero input: contract is an exact zero output */
    err = 0.0f;
    for (i = 0; i < n; i++)
      err = prop_maxf(err, fabsf(out[i]));
    return err;
  }

  out_norm = prop_vec_norm(out, n);

  if (out_norm == 0.0f) {
    /* state (a): legitimate underflow collapse, output is finite zero */
    return 0.0f;
  }

  /* state (b): unit length ... */
  err = fabsf(out_norm - 1.0f);

  /* ... and parallel to the input, out == v / norm component-wise */
  for (i = 0; i < n; i++)
    err = prop_maxf(err, fabsf(out[i] - v[i] / norm));

  return err;
}

static inline const char *
prop_vec_kind_name(prop_vec_kind_t kind) {
  switch (kind) {
    case PROP_VEC2: return "vec2";
    case PROP_VEC3: return "vec3";
    default:        return "vec4";
  }
}

static inline void
prop_normalize_one(prop_rng_t *rng, prop_vec_kind_t kind,
                   uint32_t iter, prop_fail_t *fail, float *max_err,
                   uint32_t *tiny_count, uint32_t *zero_count) {
  CGLM_ALIGN(16) vec4 v, out, inplace;
  int   n = kind == PROP_VEC2 ? 2 : (kind == PROP_VEC3 ? 3 : 4);
  int   force_zero = (prop_rng_u32(rng) % 25) == 0;
  float err, norm;
  int   finite_ok;

  prop_random_vec(rng, (float *)v, n, force_zero);
  norm = prop_vec_norm(v, n);
  if (norm == 0.0f)
    (*zero_count)++;
  else if (norm < 1e-4f)
    (*tiny_count)++;

  switch (kind) {
    case PROP_VEC2:
      glm_vec2_normalize_to(v, out);
      glm_vec2_copy(v, inplace);
      glm_vec2_normalize(inplace);
      break;
    case PROP_VEC3:
      glm_vec3_normalize_to(v, out);
      glm_vec3_copy(v, inplace);
      glm_vec3_normalize(inplace);
      break;
    default:
      glm_vec4_normalize_to(v, out);
      glm_vec4_copy(v, inplace);
      glm_vec4_normalize(inplace);
      break;
  }

  err = prop_check_normalize((float *)v, (float *)out, n, &finite_ok);

  /* in-place variant must agree with the _to variant */
  for (int i = 0; i < n; i++)
    err = prop_maxf(err, fabsf(inplace[i] - out[i]));

  if (err > *max_err) *max_err = err;

  if (!prop_fail_record(fail, iter, err, NORMALIZE_TOL) &&
      fail->iter == iter) {
    fprintf(stderr,
            "    [FAIL] %s normalize iter %u: error %.9g > tol %.9g, "
            "input norm %.9g, finite=%d\n",
            prop_vec_kind_name(kind), iter, (double)err,
            (double)NORMALIZE_TOL, (double)norm, finite_ok);
    prop_print_vec("v", (float *)v, n);
    prop_print_vec("normalize_to(v)", (float *)out, n);
  }
}

static inline prop_result_t
prop_normalize(uint32_t iters) {
  prop_rng_t    rng = {.state = PROP_SEED ^ 0xA1u};
  prop_result_t res;
  prop_fail_t   fail;
  float         max_err = 0.0f;
  uint32_t      tiny = 0, zeros = 0;

  res.name = "vec normalize stability at tiny / zero magnitudes";
  prop_fail_init(&fail);

  for (uint32_t i = 0; i < iters; i++) {
    prop_vec_kind_t kind = (prop_vec_kind_t)(prop_rng_u32(&rng) % 3u);
    prop_normalize_one(&rng, kind, i, &fail, &max_err, &tiny, &zeros);
  }

  res.passed    = !fail.failed;
  res.iters     = iters;
  res.max_error = max_err;
  res.fail      = fail;

  fprintf(stderr,
          "    coverage: %u sub-1e-4 vectors, %u exact zero vectors\n",
          tiny, zeros);
  return res;
}

#endif /* cglm_test_prop_normalize_h */
