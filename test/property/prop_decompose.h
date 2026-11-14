#ifndef cglm_test_prop_decompose_h
#define cglm_test_prop_decompose_h

#include "property_common.h"

/*
 * Property 2: affine decomposition / reconstruction round-trip.
 *
 * Build an affine matrix  M = T * R * S  from a random rotation
 * (quaternion), per-axis scale and translation, then:
 *
 *   glm_decompose(M) -> t, R', s'
 *   M' = T(t) * R' * diag(s')
 *
 * and require ||M' - M||_inf <= tol * max(1, ||M||_inf).  Near-singular
 * matrices are covered by injecting scales down to 1e-20 and fully
 * zero columns; the zero-column case is checked for the contract that
 * actually holds there (translation recovered exactly, zero scales
 * preserved) rather than exact reconstruction, since a zero column
 * destroys the rotation information irreversibly.
 */

#define DECOMPOSE_TOL 5.0e-4f
#define DECOMPOSE_SCALE_MIN 1.0e-20f
#define DECOMPOSE_SCALE_MAX 1.0e15f

static inline float
prop_mat4_inf_diff(mat4 a, mat4 b) {
  float err = 0.0f;
  for (int c = 0; c < 4; c++)
    for (int r = 0; r < 4; r++)
      err = prop_maxf(err, fabsf(a[c][r] - b[c][r]));
  return err;
}

static inline float
prop_mat4_inf_norm(mat4 m) {
  float nrm = 0.0f;
  for (int c = 0; c < 4; c++)
    for (int r = 0; r < 4; r++)
      nrm = prop_maxf(nrm, fabsf(m[c][r]));
  return nrm;
}

static inline void
prop_make_affine(prop_rng_t *rng, mat4 m, int singular_kind) {
  versor q;
  mat4   rot, trans, scale;
  vec3   s;

  for (int i = 0; i < 4; i++)
    q[i] = prop_rng_f(rng);
  if (glm_quat_norm(q) == 0.0f)
    q[3] = 1.0f;

  for (int i = 0; i < 3; i++) {
    if (singular_kind == 2 && (uint32_t)i == (prop_rng_u32(rng) % 3)) {
      s[i] = 0.0f;                       /* exactly singular column */
    } else if (singular_kind == 1 && (prop_rng_u32(rng) % 3) == 0) {
      s[i] = prop_rng_range(rng, DECOMPOSE_SCALE_MIN, 1e-6f);
    } else {
      s[i] = prop_rng_range(rng, 1e-3f, DECOMPOSE_SCALE_MAX);
    }
  }

  glm_quat_mat4(q, rot);

  glm_mat4_identity(scale);
  scale[0][0] = s[0];
  scale[1][1] = s[1];
  scale[2][2] = s[2];

  glm_mat4_identity(trans);
  trans[3][0] = prop_rng_f(rng);
  trans[3][1] = prop_rng_f(rng);
  trans[3][2] = prop_rng_f(rng);

  /* M = T * R * S */
  glm_mat4_mul(rot, scale, m);
  glm_mat4_mul(trans, m, m);
}

static inline int
prop_mat4_all_finite(mat4 m) {
  for (int c = 0; c < 4; c++)
    for (int r = 0; r < 4; r++)
      if (!prop_finite(m[c][r])) return 0;
  return 1;
}

static inline prop_result_t
prop_decompose(uint32_t iters) {
  prop_rng_t    rng = {.state = PROP_SEED ^ 0xD3u};
  prop_result_t res;
  prop_fail_t   fail;
  float         max_err = 0.0f;
  uint32_t      near_singular = 0, fully_singular = 0;

  res.name = "mat4 decompose -> reconstruct round-trip";
  prop_fail_init(&fail);

  for (uint32_t i = 0; i < iters; i++) {
    mat4   m, r, rebuilt, trans, diag;
    vec4   t;
    vec3   s;
    vec3   col_norm;
    float  min_col_norm;
    float  err;
    int    singular_kind;

    /*
     * ~10% exactly singular (one zero scale), ~20% near-singular
     * (at least one scale below 1e-6), otherwise regular.
     */
    uint32_t pick = prop_rng_u32(&rng) % 10;
    singular_kind = (pick == 0) ? 2 : (pick < 3 ? 1 : 0);

    prop_make_affine(&rng, m, singular_kind);

    glm_decompose(m, t, r, s);

    col_norm[0] = glm_vec3_norm(m[0]);
    col_norm[1] = glm_vec3_norm(m[1]);
    col_norm[2] = glm_vec3_norm(m[2]);
    min_col_norm = col_norm[0];
    for (int k = 1; k < 3; k++)
      min_col_norm = fminf(min_col_norm, col_norm[k]);

    if (min_col_norm == 0.0f) {
      fully_singular++;
      /*
       * A zero scale makes exact reconstruction impossible (the
       * rotation column is unrecoverable, and the library's 1/scale
       * normalization leaves NaN/Inf in that rotation column).  The
       * enforceable contract is:
       *   1. the call returns without trapping;
       *   2. translation round-trips exactly (it is copied, not scaled);
       *   3. the zero scale is reported as exactly zero;
       *   4. every non-zero scale and the translation stay finite.
       */
      float t_err = prop_maxf(fabsf(t[0] - m[3][0]),
                     prop_maxf(fabsf(t[1] - m[3][1]),
                               fabsf(t[2] - m[3][2])));
      err = t_err;
      for (int k = 0; k < 3; k++) {
        if (col_norm[k] == 0.0f) {
          if (s[k] != 0.0f) err = INFINITY; /* zero scale must survive */
        } else if (!prop_finite(s[k])) {
          err = INFINITY;
        }
      }
      if (!prop_finite(t[0]) || !prop_finite(t[1]) ||
          !prop_finite(t[2]) || !prop_finite(t[3]))
        err = INFINITY;
    } else {
      /* rebuild M' = T(t) * R * diag(s); safe: every scale is non-zero */
      glm_mat4_identity(diag);
      diag[0][0] = s[0];
      diag[1][1] = s[1];
      diag[2][2] = s[2];

      glm_mat4_identity(trans);
      trans[3][0] = t[0];
      trans[3][1] = t[1];
      trans[3][2] = t[2];

      glm_mat4_mul(r, diag, rebuilt);
      glm_mat4_mul(trans, rebuilt, rebuilt);

      if (min_col_norm < 1e-6f)
        near_singular++;
      /* scale-aware tolerance: float round-off grows with magnitude */
      float scale = prop_maxf(1.0f, prop_mat4_inf_norm(m));
      err = prop_mat4_inf_diff(rebuilt, m) / scale;
    }

    if (err > max_err) max_err = err;

    if (!prop_fail_record(&fail, i, err, DECOMPOSE_TOL)) {
      if (fail.iter == i) {
        fprintf(stderr,
                "    [FAIL] decompose iter %u: error %.9g > tol %.9g "
                "(min col norm %.9g)\n",
                i, (double)err, (double)DECOMPOSE_TOL,
                (double)min_col_norm);
        prop_print_mat4("M (input)", m);
        if (min_col_norm > 0.0f)
          prop_print_mat4("M' (rebuilt)", rebuilt);
        prop_print_vec("t", t, 4);
        prop_print_vec("s", s, 3);
      }
    }
  }

  res.passed    = !fail.failed;
  res.iters     = iters;
  res.max_error = max_err;
  res.fail      = fail;

  fprintf(stderr,
          "    coverage: %u near-singular (col norm < 1e-6), "
          "%u fully singular matrices\n",
          near_singular, fully_singular);
  return res;
}

#endif /* cglm_test_prop_decompose_h */
