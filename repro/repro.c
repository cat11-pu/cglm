/*
 * Reproduction for: single-axis continuous rotation flips past half a turn,
 * interpolation jumps, angle normalization and gimbal-branch discontinuities.
 *
 * Build:
 *   cc -O2 -I include repro/repro.c -L build -lcglm -lm -o build/repro
 */

#include <cglm/cglm.h>
#include <math.h>
#include <stdio.h>

#define PI2 6.283185307179586f
#define PI  3.141592653589793f

/* angular distance between two normalized orientations on SO(3) */
static float
orient_dist(const versor a, const versor b) {
  float d = fabsf(glm_quat_dot(a, b));
  if (d > 1.0f) d = 1.0f;
  return 2.0f * acosf(d);
}

static int failures = 0;

static void
check(int ok, const char *name) {
  printf("  [%s] %s\n", ok ? "PASS" : "FAIL", name);
  if (!ok) failures++;
}

/*
 * Case 1: a keyframe stream rotates around Z continuously for several turns.
 * Each frame normalizes its own angle, converts to quat, and adjacent frames
 * are slerped. Across the +-PI seam the quaternion sign flips and slerp
 * snaps back by almost a full turn.
 */
static void
case_wrapped_keyframes(void) {
  versor prev, cur, mid;
  float max_step = 0.0f;
  int n = 360; /* 5 deg steps, covers 5 full turns */
  int seam_crossings = 0;

  for (int i = 0; i <= n * 5; i++) {
    float raw = (float)i * (PI2 / n);
    float ang = glm_angle_normalize(raw);
    if (i > 0 && fabsf(ang - glm_angle_normalize((float)(i - 1) * (PI2 / n))) > PI)
      seam_crossings++;

    glm_quatv(cur, ang, (vec3){0.0f, 0.0f, 1.0f});

    if (i > 0) {
      float step = orient_dist(prev, cur);
      max_step = fmaxf(max_step, step);

      /* midpoint between adjacent keyframes must stay on the short arc */
      glm_quat_slerp(prev, cur, 0.5f, mid);
      float to_mid = orient_dist(prev, mid);
      if (to_mid > 0.5f) {
        check(0, "case1 slerp midpoint continuous across PI seam");
        printf("         frame=%d seam step=%.4f rad midpoint drift=%.4f rad\n",
               i, step, to_mid);
        return;
      }

      /* the angle-domain helper must pick the same short arc */
      float prev_raw = (float)(i - 1) * (PI2 / n);
      float a_mid = glm_angle_lerp(prev_raw, raw, 0.5f);
      versor q_mid;
      glm_quatv(q_mid, a_mid, (vec3){0.0f, 0.0f, 1.0f});
      vec3 probe = {1.0f, 0.3f, -0.7f};
      vec3 va, vb;
      glm_quat_rotatev(mid, probe, va);
      glm_quat_rotatev(q_mid, probe, vb);
      float dm = fabsf(va[0] - vb[0]) + fabsf(va[1] - vb[1]) + fabsf(va[2] - vb[2]);
      if (dm > 1e-4f) {
        printf("         frame=%d prev_raw=%.4f raw=%.4f a_mid=%.4f vec-diff=%.6f\n",
               i, prev_raw, raw, a_mid, dm);
        check(0, "case1 angle_lerp matches quat slerp across PI seam");
        return;
      }
    }
    glm_quat_copy(cur, prev);
  }

  printf("         %d seam crossings over 5 turns, max per-frame step %.5f rad\n",
         seam_crossings, max_step);
  check(seam_crossings >= 4 && max_step < 0.2f,
        "case1 slerp midpoint continuous across PI seam");
  check(1, "case1 angle_lerp matches quat slerp across PI seam");

  /* near-antipodal, near-parallel fallback must keep the short arc too */
  versor q1 = {0.0f, 0.0f, 0.1f, 0.994987f};
  versor q2, r;
  glm_quat_normalize(q1);
  glm_vec4_negate_to(q1, q2);
  q2[2] -= 1e-4f;
  glm_quat_normalize(q2);
  glm_quat_slerp(q1, q2, 0.5f, r);
  check(orient_dist(r, q1) < 1e-3f, "case1 slerp fallback sign near antipode");
}

/*
 * Case 2: euler -> matrix -> euler round trip must give the same rotation,
 * including angles outside the principal range (multi-turn input).
 */
static void
case_euler_roundtrip(void) {
  mat4 m1, m2;
  vec3 in, out;
  float worst = 0.0f;

  in[0] = 0.3f; in[1] = 8.7f; in[2] = -6.4f; /* > one turn on Y and Z */
  glm_euler_xyz(in, m1);
  glm_euler_angles(m1, out);
  glm_euler_xyz(out, m2);

  for (int r = 0; r < 4; r++)
    for (int c = 0; c < 3; c++)
      worst = fmaxf(worst, fabsf(m1[r][c] - m2[r][c]));

  printf("         multi-turn euler matrix round-trip err=%.3e\n", worst);
  check(worst < 1e-4f, "case2 multi-turn euler round trip preserves matrix");
}

/*
 * Case 3: sweep Y (the pitch axis for XYZ) through +-PI/2. The extracted
 * angles may switch between equivalent euler representations at the pole,
 * but the orientation rebuilt from those angles must not jump: any branch
 * switch is only allowed if both representations describe the same matrix.
 */
static void
case_gimbal_sweep(void) {
  mat4 m, rebuilt, prev_rebuilt;
  vec3 out;
  float worst_orient = 0.0f;
  int jump_at = -1;
  int have_prev = 0;

  for (int i = 0; i <= 4000; i++) {
    float y = -PI + (float)i * (PI2 / 4000.0f);
    vec3 in = {0.25f, y, -0.4f};
    glm_euler_xyz(in, m);
    glm_euler_angles(m, out);
    glm_euler_xyz(out, rebuilt);

    if (have_prev) {
      float j = 0.0f;
      for (int r = 0; r < 4; r++)
        for (int c = 0; c < 3; c++)
          j = fmaxf(j, fabsf(rebuilt[r][c] - prev_rebuilt[r][c]));
      if (j > worst_orient) { worst_orient = j; jump_at = i; }
    }
    glm_mat4_copy(rebuilt, prev_rebuilt);
    have_prev = 1;
  }

  float step = PI2 / 4000.0f;
  printf("         worst rebuilt-orientation jump=%.3e (input step %.5f, frame %d)\n",
         worst_orient, step, jump_at);
  check(worst_orient < 4.0f * step, "case3 orientation continuous near gimbal lock");

  /* exact gimbal-lock matrices must round-trip exactly (roll info lives in
     m01/m11 at the pole and must not be discarded) */
  float lock_err = 0.0f;
  for (int pole = -1; pole <= 1; pole += 2) {
    for (int k = -3; k <= 3; k++) {
      vec3 in = {(float)k * 0.37f, pole > 0 ? GLM_PI_2f : -GLM_PI_2f, 0.0f};
      glm_euler_xyz(in, m);
      glm_euler_angles(m, out);
      glm_euler_xyz(out, rebuilt);
      for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++)
          lock_err = fmaxf(lock_err, fabsf(m[r][c] - rebuilt[r][c]));
    }
  }
  printf("         exact-lock round-trip matrix err=%.3e\n", lock_err);
  check(lock_err < 1e-5f, "case3 exact gimbal lock round trip preserves roll");
}

/*
 * Case 4: euler -> quat -> matrix -> euler round trip across a dense grid
 * including the PI seam and gimbal neighborhoods; bound the angle error.
 */
static void
case_grid_roundtrip(void) {
  versor q;
  mat4 m;
  vec3 in, out;
  float worst = 0.0f;
  int total = 0;

  for (int ix = -4; ix <= 4; ix++) {
    for (int iy = -4; iy <= 4; iy++) {
      for (int iz = -4; iz <= 4; iz++) {
        in[0] = ix * 0.5f;
        in[1] = iy * 0.5f;
        in[2] = iz * 0.5f;

        glm_euler_xyz_quat(in, q);
        glm_quat_mat4(q, m);
        glm_euler_angles(m, out);

        /* compare rotations as matrices built from both angle sets */
        mat4 ma, mb;
        glm_euler_xyz(in, ma);
        glm_euler_xyz(out, mb);
        float e = 0.0f;
        for (int r = 0; r < 4; r++)
          for (int c = 0; c < 4; c++)
            e = fmaxf(e, fabsf(ma[r][c] - mb[r][c]));
        worst = fmaxf(worst, e);
        total++;
      }
    }
  }

  printf("         grid %d points, worst matrix err=%.3e\n", total, worst);
  check(worst < 1e-4f, "case4 euler<->quat round trip within 1e-4");
}

int
main(void) {
  printf("rotation continuity reproduction\n");
  case_wrapped_keyframes();
  case_euler_roundtrip();
  case_gimbal_sweep();
  case_grid_roundtrip();
  printf("\n%s: %d failing case(s)\n", failures ? "FAILED" : "OK", failures);
  return failures ? 1 : 0;
}
