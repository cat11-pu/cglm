#ifndef cglm_test_prop_interp_h
#define cglm_test_prop_interp_h

#include "property_common.h"

/*
 * Property 4: interpolation endpoint consistency.
 *
 * For every randomized pair of inputs:
 *
 *   interp(a, b, 0) == a      interp(a, b, 1) == b
 *
 * Checked for vec2/vec3/vec4 lerp, quaternion lerp/lerpc/slerp/nlerp.
 * Zero vectors and identity / near-identical / opposite quaternions
 * are explicitly covered.
 */

/*
 * Endpoint equality is limited by float round-off.
 *
 * t = 0 is exact: lerp(a, b, 0) evaluates a + 0*(b-a) and must equal a
 * bit-for-bit (covered by INTERP_TOL on the output).
 *
 * t = 1 goes through a + 1*(b-a); when |a| >> |b| the subtraction and
 * re-addition lose low bits of b (catastrophic cancellation at a mild
 * level), so its tolerance is ULP-scaled to the larger operand rather
 * than to b alone.
 */
#define INTERP_TOL      1.0e-6f
#define INTERP_TOL_REL  2.0e-4f
#define INTERP_TOL_QUAT 2.0e-4f

static inline float
prop_vec_maxdiff(const float *a, const float *b, int n) {
  float err = 0.0f;
  for (int i = 0; i < n; i++)
    err = prop_maxf(err, fabsf(a[i] - b[i]));
  return err;
}

static inline float
prop_endpoint_error(const float *expect, const float *got, int n) {
  float err = 0.0f;
  for (int i = 0; i < n; i++) {
    float scale = prop_maxf(1.0f, fabsf(expect[i]));
    err = prop_maxf(err, fabsf(expect[i] - got[i]) / scale);
  }
  return err;
}

/* t=1 error scaled to the larger operand's magnitude */
static inline float
prop_endpoint_error_scaled(const float *a, const float *b,
                           const float *got, int n) {
  float err = 0.0f;
  for (int i = 0; i < n; i++) {
    float scale = prop_maxf(1.0f, prop_maxf(fabsf(a[i]), fabsf(b[i])));
    err = prop_maxf(err, fabsf(b[i] - got[i]) / scale);
  }
  return err;
}

/*
 * Quaternions q and -q represent the same rotation.  nlerp flips the
 * target hemisphere internally, so endpoint equality is checked up to
 * a global sign.
 */
static inline float
prop_quat_eq_up_to_sign(const versor a, const versor b) {
  float d1 = 0.0f, d2 = 0.0f;
  for (int i = 0; i < 4; i++) {
    d1 = prop_maxf(d1, fabsf(a[i] - b[i]));
    d2 = prop_maxf(d2, fabsf(a[i] + b[i]));
  }
  return fminf(d1, d2);
}

static inline void
prop_random_pair(prop_rng_t *rng, float *a, float *b, int n) {
  float mag_a = prop_rng_magnitude(rng);
  float mag_b = prop_rng_magnitude(rng);

  for (int i = 0; i < n; i++) {
    a[i] = prop_rng_f(rng) * mag_a;
    b[i] = prop_rng_f(rng) * mag_b;
  }

  switch (prop_rng_u32(rng) % 12) {
    case 0: /* both zero */
      for (int i = 0; i < n; i++) { a[i] = 0.0f; b[i] = 0.0f; }
      break;
    case 1: /* a zero */
      for (int i = 0; i < n; i++) a[i] = 0.0f;
      break;
    case 2: /* b zero */
      for (int i = 0; i < n; i++) b[i] = 0.0f;
      break;
    case 3: /* identical endpoints */
      for (int i = 0; i < n; i++) b[i] = a[i];
      break;
    default:
      break;
  }
}

static inline void
prop_random_quat_pair(prop_rng_t *rng, versor a, versor b) {
  for (int i = 0; i < 4; i++) {
    a[i] = prop_rng_f(rng);
    b[i] = prop_rng_f(rng);
  }
  if (glm_quat_norm(a) == 0.0f) a[3] = 1.0f;
  if (glm_quat_norm(b) == 0.0f) b[3] = 1.0f;
  /*
   * slerp / nlerp are only defined on unit quaternions (the library
   * assumes |dot| <= 1); normalize first so the property exercises
   * valid rotations.
   */
  glm_quat_normalize(a);
  glm_quat_normalize(b);

  switch (prop_rng_u32(rng) % 12) {
    case 0: /* identity pair */
      glm_quat_identity(a);
      glm_quat_identity(b);
      break;
    case 1: /* identical quaternions (slerp short-circuit path) */
      glm_quat_copy(a, b);
      break;
    case 2: { /* opposite quaternions (dot < 0 hemisphere flip) */
      versor q;
      for (int i = 0; i < 4; i++)
        q[i] = prop_rng_f(rng);
      if (glm_quat_norm(q) == 0.0f) q[3] = 1.0f;
      glm_quat_normalize(q);
      glm_quat_copy(q, a);
      for (int i = 0; i < 4; i++) b[i] = -q[i];
      break;
    }
    case 3: /* nearly identical unit quaternions (lerp fallback path) */
      glm_quat_copy(a, b);
      for (int i = 0; i < 4; i++)
        b[i] += prop_rng_range(rng, -1e-4f, 1e-4f);
      glm_quat_normalize(b);
      break;
    default:
      break;
  }
}

static inline void
prop_interp_failure(prop_fail_t *fail, uint32_t iter, const char *op,
                    float err, const float *a, const float *b,
                    const float *at0, const float *at1, int n) {
  if (fail->iter != iter)
    return;
  fprintf(stderr,
          "    [FAIL] %s endpoint iter %u: error %.9g > tol %.9g\n",
          op, iter, (double)err, (double)INTERP_TOL);
  prop_print_vec("a (t=0 expected)", a, n);
  prop_print_vec("b (t=1 expected)", b, n);
  prop_print_vec("interp(...,0)", at0, n);
  prop_print_vec("interp(...,1)", at1, n);
}

static inline prop_result_t
prop_interp(uint32_t iters) {
  prop_rng_t    rng = {.state = PROP_SEED ^ 0x17u};
  prop_result_t res;
  prop_fail_t   fail;
  float         max_err = 0.0f;
  uint32_t      edge_cases = 0;

  res.name = "interpolation endpoint consistency (t=0 -> a, t=1 -> b)";
  prop_fail_init(&fail);

  for (uint32_t i = 0; i < iters; i++) {
    CGLM_ALIGN(16) vec4 a, b, at0, at1;
    versor qa, qb, q0, q1;
    float err;
    int   n;
    uint32_t mode = prop_rng_u32(&rng) % 7u;

    if (mode < 3) {
      /* ---------- vector lerp ---------- */
      n = mode == 0 ? 2 : (mode == 1 ? 3 : 4);
      prop_random_pair(&rng, (float *)a, (float *)b, n);
      if (prop_vec_maxdiff((float *)a, (float *)b, n) == 0.0f ||
          prop_vec_norm((float *)a, n) == 0.0f ||
          prop_vec_norm((float *)b, n) == 0.0f)
        edge_cases++;

      if (n == 2) {
        glm_vec2_lerp(a, b, 0.0f, at0);
        glm_vec2_lerp(a, b, 1.0f, at1);
      } else if (n == 3) {
        glm_vec3_lerp(a, b, 0.0f, at0);
        glm_vec3_lerp(a, b, 1.0f, at1);
      } else {
        glm_vec4_lerp(a, b, 0.0f, at0);
        glm_vec4_lerp(a, b, 1.0f, at1);
      }

      /* t=0 is exact; t=1 tolerates operand-magnitude-scaled round-off */
      err = prop_maxf(prop_endpoint_error((float *)a, (float *)at0, n),
                      prop_endpoint_error_scaled((float *)a, (float *)b,
                                                 (float *)at1, n));
      if (err > max_err) max_err = err;
      if (!prop_fail_record(&fail, i, err, INTERP_TOL_REL))
        prop_interp_failure(&fail, i, "vec lerp", err,
                            (float *)a, (float *)b,
                            (float *)at0, (float *)at1, n);
    } else {
      /* ---------- quaternion interpolation ---------- */
      n = 4;
      prop_random_quat_pair(&rng, qa, qb);
      if (glm_vec4_eqv(qa, qb) || fabsf(glm_quat_dot(qa, qb)) >= 0.999f)
        edge_cases++;

      switch (mode) {
        case 3:
          glm_quat_lerp(qa, qb, 0.0f, q0);
          glm_quat_lerp(qa, qb, 1.0f, q1);
          err = prop_maxf(prop_vec_maxdiff(qa, q0, 4),
                          prop_vec_maxdiff(qb, q1, 4));
          if (!prop_fail_record(&fail, i, err, INTERP_TOL_QUAT))
            prop_interp_failure(&fail, i, "quat lerp", err,
                                qa, qb, q0, q1, 4);
          break;
        case 4:
          glm_quat_lerpc(qa, qb, 0.0f, q0);
          glm_quat_lerpc(qa, qb, 1.0f, q1);
          err = prop_maxf(prop_vec_maxdiff(qa, q0, 4),
                          prop_vec_maxdiff(qb, q1, 4));
          if (!prop_fail_record(&fail, i, err, INTERP_TOL_QUAT))
            prop_interp_failure(&fail, i, "quat lerpc", err,
                                qa, qb, q0, q1, 4);
          break;
        case 5:
          glm_quat_slerp(qa, qb, 0.0f, q0);
          glm_quat_slerp(qa, qb, 1.0f, q1);
          /*
           * slerp flips the start hemisphere when dot(from,to) < 0,
           * so the t=0 endpoint may come back as -from.  The t=1
           * endpoint is always the unflipped `to`.  Compare both up
           * to sign since q and -q are the same rotation.
           */
          err = prop_maxf(prop_quat_eq_up_to_sign(qa, q0),
                          prop_quat_eq_up_to_sign(qb, q1));
          if (!prop_fail_record(&fail, i, err, INTERP_TOL_QUAT))
            prop_interp_failure(&fail, i, "quat slerp", err,
                                qa, qb, q0, q1, 4);
          break;
        default: {
          /* nlerp normalizes and may flip the target hemisphere */
          versor na, nb;
          glm_quat_normalize_to(qa, na);
          glm_quat_normalize_to(qb, nb);
          glm_quat_nlerp(qa, qb, 0.0f, q0);
          glm_quat_nlerp(qa, qb, 1.0f, q1);
          err = prop_maxf(prop_quat_eq_up_to_sign(na, q0),
                          prop_quat_eq_up_to_sign(nb, q1));
          if (!prop_fail_record(&fail, i, err, INTERP_TOL_QUAT))
            prop_interp_failure(&fail, i, "quat nlerp", err,
                                na, nb, q0, q1, 4);
          break;
        }
      }
      if (err > max_err) max_err = err;
    }
  }

  res.passed    = !fail.failed;
  res.iters     = iters;
  res.max_error = max_err;
  res.fail      = fail;

  fprintf(stderr,
          "    coverage: %u degenerate / edge-case endpoint pairs\n",
          edge_cases);
  return res;
}

#endif /* cglm_test_prop_interp_h */
