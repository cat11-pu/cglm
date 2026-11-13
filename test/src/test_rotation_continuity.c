/*
 * Copyright (c), Recep Aslantas.
 *
 * MIT License (MIT), http://opensource.org/licenses/MIT
 * Full license can be found in the LICENSE file
 */

/*
 * Regression tests for continuous rotation across angle branch cuts:
 *   - angle normalization / shortest-arc helpers
 *   - euler angle extraction near the gimbal-lock poles
 *   - slerp sign handling in the near-parallel fallback
 *   - euler <-> quaternion round trips with a 1e-4 rad budget
 */

#include "test_common.h"

#define ANGLE_TWO_PI (2.0f * GLM_PIf)

static float
quat_orient_dist(const versor a, const versor b) {
  float d = fabsf(glm_quat_dot(a, b));
  if (d > 1.0f) d = 1.0f;
  return 2.0f * acosf(d);
}

static float
mat4_max_diff(CGLM_CONST mat4 a, CGLM_CONST mat4 b) {
  float e = 0.0f;
  for (int r = 0; r < 4; r++)
    for (int c = 0; c < 4; c++)
      e = fmaxf(e, fabsf(a[r][c] - b[r][c]));
  return e;
}

TEST_IMPL(angle_normalize) {
  /* every representative angle maps into (-PI, PI] */
  float samples[] = {
    0.0f, 1.0f, -1.0f, GLM_PIf, -GLM_PIf,
    GLM_PIf + 0.5f, -GLM_PIf - 0.5f,
    3.0f * GLM_PIf, -3.0f * GLM_PIf,
    12.7f, -12.7f, 100.0f * GLM_PIf
  };

  for (size_t i = 0; i < sizeof(samples) / sizeof(samples[0]); i++) {
    float n = glm_angle_normalize(samples[i]);
    ASSERT(n > -GLM_PIf)
    ASSERT(n <= GLM_PIf)

    /* normalization must not change the represented orientation */
    versor q1, q2;
    glm_quatv(q1, samples[i], (vec3){0.0f, 1.0f, 0.0f});
    glm_quatv(q2, n, (vec3){0.0f, 1.0f, 0.0f});
    ASSERT(test_eq_th(quat_orient_dist(q1, q2), 0.0f, 1e-6f))
  }

  /* the half-open cut is (-PI, PI]: -PI normalizes to +PI, which is the
     same orientation, and +PI stays put */
  ASSERT(test_eq_th(glm_angle_normalize(GLM_PIf), GLM_PIf, 1e-6f))
  ASSERT(test_eq_th(glm_angle_normalize(-GLM_PIf), GLM_PIf, 1e-6f))
  ASSERT(test_eq_th(glm_angle_normalize(3.0f * GLM_PIf), GLM_PIf, 1e-6f))

  TEST_SUCCESS
}

TEST_IMPL(angle_shortest_and_lerp) {
  /* shortest distance always lands in (-PI, PI] and points the short way */
  ASSERT(test_eq_th(glm_angle_shortest(3.0f, 3.2f), 0.2f, 1e-6f))
  ASSERT(test_eq_th(glm_angle_shortest(-3.0f, 3.0f), 6.0f - ANGLE_TWO_PI, 1e-6f))
  ASSERT(test_eq_th(glm_angle_shortest(3.1f, -3.1f), -6.2f + ANGLE_TWO_PI, 1e-6f))

  /* interpolation across the seam follows the short arc instead of ~2PI */
  float mid = glm_angle_lerp(GLM_PIf - 0.1f, -GLM_PIf + 0.1f, 0.5f);
  ASSERT(fabsf(mid - GLM_PIf) < 1e-5f || fabsf(mid + GLM_PIf) < 1e-5f)

  /* endpoints are reproduced exactly */
  ASSERT(test_eq_th(glm_angle_lerp(0.7f, 2.3f, 0.0f), 0.7f, 1e-6f))
  ASSERT(test_eq_th(glm_angle_lerp(0.7f, 2.3f, 1.0f), 2.3f, 1e-6f))

  /* vec3 variant normalizes every component independently */
  vec3 a = {GLM_PIf + 0.2f, -GLM_PIf - 0.2f, 7.0f * GLM_PIf};
  vec3 b = {0.0f, 0.0f, 0.0f};
  vec3 out;
  glm_vec3_angle_shortest(b, a, out);
  for (int i = 0; i < 3; i++) {
    ASSERT(out[i] > -GLM_PIf)
    ASSERT(out[i] <= GLM_PIf)
  }

  glm_vec3_angle_lerp(b, a, 0.5f, out);
  vec3 full;
  glm_vec3_angle_lerp(b, a, 1.0f, full);
  /* t=1 reaches the target orientation; the scalar value may sit on the
     other side of the branch cut, so compare as rotations */
  for (int i = 0; i < 3; i++)
    ASSERT(test_eq_th(glm_angle_shortest(full[i], a[i]), 0.0f, 1e-5f))

  TEST_SUCCESS
}

TEST_IMPL(euler_continuous_through_poles) {
  /*
   * Sweep the pitch angle through both gimbal-lock poles. The orientation
   * rebuilt from the extracted euler angles must track the source matrix
   * with no branch-cut jump, including frames that hit the exact pole.
   */
  mat4 m, rebuilt;
  vec3 in, out;
  float worst = 0.0f;
  const int n = 4000;

  for (int i = 0; i <= n; i++) {
    in[0] = 0.25f;
    in[1] = -GLM_PIf + (float)i * (ANGLE_TWO_PI / (float)n);
    in[2] = -0.4f;

    glm_euler_xyz(in, m);
    glm_euler_angles(m, out);
    glm_euler_xyz(out, rebuilt);

    worst = fmaxf(worst, mat4_max_diff(m, rebuilt));
  }

  ASSERT(test_eq_th(worst, 0.0f, 2e-5f))

  /* exact locks must preserve the roll stored at the pole */
  for (int pole = -1; pole <= 1; pole += 2) {
    for (int k = -3; k <= 3; k++) {
      in[0] = (float)k * 0.37f;
      in[1] = pole > 0 ? GLM_PI_2f : -GLM_PI_2f;
      in[2] = 0.0f;

      glm_euler_xyz(in, m);
      glm_euler_angles(m, out);
      glm_euler_xyz(out, rebuilt);

      ASSERT(test_eq_th(mat4_max_diff(m, rebuilt), 0.0f, 1e-5f))
    }
  }

  TEST_SUCCESS
}

TEST_IMPL(quat_slerp_fallback_sign) {
  /*
   * Nearly antipodal, nearly parallel quaternions take the sinTheta < 0.001
   * fallback. The fallback must honor the short-arc sign flip instead of
   * lerping toward the opposite hemisphere.
   */
  versor q1 = {0.0f, 0.0f, 0.1f, 0.994987f};
  versor q2, r;

  glm_quat_normalize(q1);
  glm_vec4_negate_to(q1, q2);
  q2[2] -= 1e-4f;
  glm_quat_normalize(q2);

  ASSERT(glm_quat_dot(q1, q2) < 0.0f)

  glm_quat_slerp(q1, q2, 0.5f, r);
  ASSERT(test_eq_th(quat_orient_dist(r, q1), 0.0f, 1e-3f))

  /* endpoint invariants hold on both paths */
  glm_quat_slerp(q1, q2, 0.0f, r);
  ASSERT(test_eq_th(quat_orient_dist(r, q1), 0.0f, 1e-6f))
  glm_quat_slerp(q1, q2, 1.0f, r);
  ASSERT(test_eq_th(quat_orient_dist(r, q2), 0.0f, 1e-6f))
  glm_quat_slerp_longest(q1, q2, 0.0f, r);
  ASSERT(test_eq_th(quat_orient_dist(r, q1), 0.0f, 1e-6f))
  glm_quat_slerp_longest(q1, q2, 1.0f, r);
  ASSERT(test_eq_th(quat_orient_dist(r, q2), 0.0f, 1e-6f))

  /* longest path keeps the opposite convention */
  glm_quat_slerp_longest(q1, q2, 0.5f, r);
  ASSERT(quat_orient_dist(r, q1) > GLM_PI_2f - 1e-2f)

  TEST_SUCCESS
}

TEST_IMPL(euler_quat_roundtrip_grid) {
  /*
   * Dense grid over several turns, including the PI seams and the poles.
   * The euler -> quat -> matrix -> euler cycle must rebuild an equivalent
   * rotation within the 1e-4 radian error budget.
   */
  versor q;
  mat4 m, ref;
  vec3 in, out;
  float worst = 0.0f;

  for (int ix = -6; ix <= 6; ix++) {
    for (int iy = -6; iy <= 6; iy++) {
      for (int iz = -6; iz <= 6; iz++) {
        in[0] = (float)ix * 0.5f;
        in[1] = (float)iy * 0.5f;
        in[2] = (float)iz * 0.5f;

        glm_euler_xyz(in, ref);
        glm_euler_xyz_quat(in, q);
        glm_quat_mat4(q, m);
        glm_euler_angles(m, out);
        glm_euler_xyz(out, m);

        worst = fmaxf(worst, mat4_max_diff(ref, m));
      }
    }
  }

  /* matrix entries are unit-scale, so this is comparable to a radian bound */
  ASSERT(worst < 1e-4f)

  TEST_SUCCESS
}

TEST_IMPL(continuous_keyframe_slerp) {
  /*
   * Five full turns around a single axis with per-frame angle
   * normalization. Adjacent quaternions must always be connected by the
   * short arc and slerp midpoints must stay on that arc.
   */
  versor prev, cur, mid;
  vec3 z_axis = {0.0f, 0.0f, 1.0f};
  const int per_turn = 360;
  float max_step = 0.0f;
  int seams = 0;

  for (int i = 0; i <= per_turn * 5; i++) {
    float raw = (float)i * (ANGLE_TWO_PI / (float)per_turn);
    float ang = glm_angle_normalize(raw);

    glm_quatv(cur, ang, z_axis);

    if (i > 0) {
      float prev_raw = (float)(i - 1) * (ANGLE_TWO_PI / (float)per_turn);
      if (fabsf(ang - glm_angle_normalize(prev_raw)) > GLM_PIf)
        seams++;

      float step = quat_orient_dist(prev, cur);
      max_step = fmaxf(max_step, step);
      ASSERT(step < 0.2f)

      glm_quat_slerp(prev, cur, 0.5f, mid);

      /* slerp stays on the short arc: the midpoint never sits farther from
         either endpoint than the endpoints are from each other, and across
         a seam it must not snap back toward the previous orientation */
      float a_mid = glm_angle_lerp(prev_raw, raw, 0.5f);
      versor q_mid;
      glm_quatv(q_mid, a_mid, z_axis);
      ASSERT(quat_orient_dist(prev, mid) <= step + 1e-6f)
      ASSERT(quat_orient_dist(mid, cur) <= step + 1e-6f)
      ASSERT(quat_orient_dist(mid, q_mid) < 0.02f)
    }

    glm_quat_copy(cur, prev);
  }

  ASSERT(seams >= 4)

  TEST_SUCCESS
}
