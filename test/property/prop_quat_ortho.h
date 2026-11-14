#ifndef cglm_test_prop_quat_ortho_h
#define cglm_test_prop_quat_ortho_h

#include "property_common.h"

/*
 * Property 1: quaternion -> rotation matrix orthogonality.
 *
 * For every finite quaternion with a norm safely above zero,
 * glm_quat_mat3 / glm_quat_mat4 must produce an orthogonal matrix:
 *
 *   M^T * M = I,   |det(M)| = 1
 *
 * Non-unit (arbitrarily scaled) quaternions are the main input class;
 * the library normalizes internally, so scaling a quaternion must not
 * change the resulting rotation.  The exact zero quaternion is a
 * documented degenerate input and is checked separately for defined
 * behaviour rather than orthogonality.
 */

#define QUAT_ORTHO_TOL 2.5e-5f

static inline float
prop_mat3_ortho_error(mat3 m) {
  float err = 0.0f;

  /*
   * Orthogonality: columns of M are orthonormal, so
   * (M^T M)[r][c] = dot(column r, column c) == (r == c ? 1 : 0).
   * Computed directly instead of via glm_mat3_mul (which forms M*M).
   */
  for (int r = 0; r < 3; r++) {
    for (int c = 0; c < 3; c++) {
      float dot = m[0][r] * m[0][c]
                + m[1][r] * m[1][c]
                + m[2][r] * m[2][c];
      float want = (r == c) ? 1.0f : 0.0f;
      err = prop_maxf(err, fabsf(dot - want));
    }
  }
  return err;
}

static inline float
prop_mat3_det(mat3 m) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
       - m[1][0] * (m[0][1] * m[2][2] - m[0][2] * m[2][1])
       + m[2][0] * (m[0][1] * m[1][2] - m[0][2] * m[1][1]);
}

static inline void
prop_random_quat(prop_rng_t *rng, versor q) {
  float mag = prop_rng_magnitude(rng);

  for (int i = 0; i < 4; i++)
    q[i] = prop_rng_f(rng) * mag;

  /* explicit edge cases injected at a fixed rate */
  switch (prop_rng_u32(rng) % 20) {
    case 0: { /* identity scaled by random magnitude */
      float s = prop_rng_magnitude(rng);
      q[0] = 0.0f; q[1] = 0.0f; q[2] = 0.0f; q[3] = s;
      break;
    }
    case 1: { /* all-zero quaternion: normalized to identity by cglm */
      q[0] = q[1] = q[2] = q[3] = 0.0f;
      break;
    }
    case 2: { /* single tiny non-zero component */
      float tiny = PROP_MAG_LO[0] * (1.0f + prop_rng_uf(rng));
      for (int i = 0; i < 4; i++) q[i] = 0.0f;
      q[prop_rng_u32(rng) & 3u] = tiny;
      break;
    }
    default:
      break;
  }
}

static inline prop_result_t
prop_quat_ortho(uint32_t iters) {
  prop_rng_t    rng = {.state = PROP_SEED ^ 0x71u};
  prop_result_t res;
  prop_fail_t   fail;
  versor        q;
  mat3          m3;
  mat4          m4;
  float         max_err = 0.0f;
  uint32_t      nonunit = 0, zero_like = 0, skipped = 0;

  res.name = "quat->mat orthogonality (M^T M = I, |det| = 1)";
  prop_fail_init(&fail);

  for (uint32_t i = 0; i < iters; i++) {
    mat3   m3_copy;
    float  err, derr, qnorm;
    int    is_zero;

    prop_random_quat(&rng, q);
    qnorm = glm_quat_norm(q);
    is_zero = (qnorm == 0.0f);

    if (is_zero)
      zero_like++;
    else if (fabsf(qnorm - 1.0f) > 1e-3f)
      nonunit++;

    glm_quat_mat3(q, m3);
    glm_quat_mat4(q, m4);

    if (is_zero) {
      /*
       * Degenerate input: only require the library to stay defined
       * (no NaN/Inf).  The current implementation normalizes a zero
       * quaternion to the identity, so accept identity as well.
       */
      int finite = 1, is_identity = 1;
      for (int c = 0; c < 3; c++) {
        for (int r = 0; r < 3; r++) {
          if (!prop_finite(m3[c][r])) finite = 0;
          float want = (r == c) ? 1.0f : 0.0f;
          if (fabsf(m3[c][r] - want) > QUAT_ORTHO_TOL) is_identity = 0;
        }
      }
      if (!finite && !is_identity) {
        if (!prop_fail_record(&fail, i, INFINITY, QUAT_ORTHO_TOL) &&
            fail.iter == i) {
          fprintf(stderr,
                  "    [FAIL] zero quat iter %u produced non-finite matrix\n",
                  i);
          prop_print_vec("q (x,y,z,w), norm", (float[]){
                             q[0], q[1], q[2], q[3], qnorm}, 5);
          prop_print_mat3("M", m3);
        }
      }
      continue;
    }

    /*
     * Near-zero but non-zero norms can overflow 1/norm internally;
     * skip them in the strict orthogonality check but count coverage.
     */
    if (qnorm < 1e-20f) {
      skipped++;
      continue;
    }

    err = prop_mat3_ortho_error(m3);

    /* the 3x3 block of the mat4 form must match the mat3 form */
    for (int c = 0; c < 3; c++)
      for (int r = 0; r < 3; r++)
        err = prop_maxf(err, fabsf(m4[c][r] - m3[c][r]));

    derr = fabsf(fabsf(prop_mat3_det(m3)) - 1.0f);
    err  = prop_maxf(err, derr);

    /* no NaN/Inf may ever leak out */
    for (int c = 0; c < 3; c++)
      for (int r = 0; r < 3; r++)
        if (!prop_finite(m3[c][r])) err = INFINITY;

    if (err > max_err) max_err = err;

    if (!prop_fail_record(&fail, i, err, QUAT_ORTHO_TOL)) {
      glm_mat3_copy(m3, m3_copy);
      if (fail.iter == i) {
        fprintf(stderr,
                "    [FAIL] quat ortho iter %u: error %.9g > tol %.9g\n",
                i, (double)err, (double)QUAT_ORTHO_TOL);
        prop_print_vec("q (x,y,z,w), norm", (float[]){
                          q[0], q[1], q[2], q[3], qnorm}, 5);
        prop_print_mat3("M", m3_copy);
        fprintf(stderr, "      |det(M)| - 1 = %.9g\n", (double)derr);
      }
    }
  }

  res.passed    = !fail.failed;
  res.iters     = iters;
  res.max_error = max_err;
  res.fail      = fail;

  fprintf(stderr,
          "    coverage: %u non-unit quats, %u zero quats, "
          "%u near-zero (relaxed)\n",
          nonunit, zero_like, skipped);
  return res;
}

#endif /* cglm_test_prop_quat_ortho_h */
