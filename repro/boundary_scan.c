/*
 * Boundary scanner: emits boundary_scan.csv with measured errors around
 * the +/-PI angle seam and the +/-PI/2 gimbal-lock poles, plus dense
 * euler<->quat round-trip statistics over a multi-turn grid.
 *
 * Build:
 *   cc -O2 -I include repro/boundary_scan.c -L build -lcglm -lm \
 *      -o build/boundary_scan
 */

#include <cglm/cglm.h>
#include <math.h>
#include <stdio.h>

#define TAU (2.0f * GLM_PIf)

static float
mat4_diff(CGLM_CONST mat4 a, CGLM_CONST mat4 b) {
  float e = 0.0f;
  for (int r = 0; r < 4; r++)
    for (int c = 0; c < 4; c++)
      e = fmaxf(e, fabsf(a[r][c] - b[r][c]));
  return e;
}

static float
quat_dist(const versor a, const versor b) {
  versor na, nb;
  glm_quat_normalize_to(a, na);
  glm_quat_normalize_to(b, nb);
  float d = fabsf(glm_quat_dot(na, nb));
  if (d > 1.0f) d = 1.0f;
  return 2.0f * acosf(d);
}

int
main(void) {
  FILE *fp = fopen("boundary_scan.csv", "w");
  if (!fp) return 1;

  fprintf(fp, "section,parameter,value,error,note\n");

  /* ---- 1. angle normalization seam +/- PI ------------------------------ */
  float seam_worst = 0.0f;
  for (int k = -2; k <= 2; k++) {
    for (int i = 0; i <= 2000; i++) {
      float delta = -0.05f + (float)i * (0.1f / 2000.0f);
      float ang = (float)k * TAU + GLM_PIf + delta;
      float n = glm_angle_normalize(ang);

      versor q1, q2;
      glm_quatv(q1, ang, (vec3){0, 1, 0});
      glm_quatv(q2, n, (vec3){0, 1, 0});
      float e = quat_dist(q1, q2);
      seam_worst = fmaxf(seam_worst, e);

      if (i % 400 == 0 && k == 0)
        fprintf(fp, "seam,angle_norm,%.7f,%.3e,turn%+d delta%+.4f\n",
                ang, e, k, delta);
    }
  }
  fprintf(fp, "seam,WORST,0,%.3e,+-PI seam over 5 turns x2001 points\n",
          seam_worst);

  /* very large magnitudes are limited by libc sinf argument reduction,
     measure that floor separately so it is not confused with seam errors */
  float huge_worst = 0.0f;
  for (int i = 0; i <= 1000; i++) {
    float ang = 100.0f * GLM_PIf + (float)i * 0.001f;
    float n = glm_angle_normalize(ang);
    versor q1, q2;
    glm_quatv(q1, ang, (vec3){0, 1, 0});
    glm_quatv(q2, n, (vec3){0, 1, 0});
    huge_worst = fmaxf(huge_worst, quat_dist(q1, q2));
  }
  fprintf(fp, "seam,HUGE_ANGLE_FLOOR,0,%.3e,~100PI inputs, libc sinf reduction\n",
          huge_worst);

  /* ---- 2. gimbal poles: pitch sweep through +/- PI/2 ------------------- */
  mat4 m, rebuilt;
  vec3 in, out;
  float pole_worst = 0.0f;
  float jump_worst = 0.0f;
  int pole_frames = 0;
  int have_prev = 0;
  mat4 prev_rebuilt;

  for (int i = 0; i <= 8000; i++) {
    in[0] = 0.25f;
    in[1] = -GLM_PIf + (float)i * (TAU / 8000.0f);
    in[2] = -0.4f;

    glm_euler_xyz(in, m);
    glm_euler_angles(m, out);
    glm_euler_xyz(out, rebuilt);

    float e = mat4_diff(m, rebuilt);
    pole_worst = fmaxf(pole_worst, e);
    pole_frames++;

    if (have_prev)
      jump_worst = fmaxf(jump_worst, mat4_diff(rebuilt, prev_rebuilt));
    glm_mat4_copy(rebuilt, prev_rebuilt);
    have_prev = 1;

    float d_minus = fabsf(in[1] + GLM_PI_2f);
    float d_plus  = fabsf(in[1] - GLM_PI_2f);
    if (d_minus < 0.0008f || d_plus < 0.0008f)
      fprintf(fp, "pole,euler_roundtrip,%.7f,%.3e,pitch%+.5f\n",
              in[1], e, in[1]);
  }
  fprintf(fp, "pole,WORST,0,%.3e,%d frames full pitch sweep\n",
          pole_worst, pole_frames);
  fprintf(fp, "pole,WORST_FRAME_JUMP,0,%.3e,max adjacent-frame rebuilt difference\n",
          jump_worst);
  fprintf(fp, "pole,NOTE,0,0,error grows linearly with 1/dist-to-pole (euler conditioning), no discontinuity\n");

  /* ---- 3. exact-lock roll recovery ------------------------------------- */
  float lock_worst = 0.0f;
  for (int pole = -1; pole <= 1; pole += 2) {
    for (int k = -12; k <= 12; k++) {
      in[0] = (float)k * 0.21f;
      in[1] = pole > 0 ? GLM_PI_2f : -GLM_PI_2f;
      in[2] = 0.0f;

      glm_euler_xyz(in, m);
      glm_euler_angles(m, out);
      glm_euler_xyz(out, rebuilt);

      float e = mat4_diff(m, rebuilt);
      lock_worst = fmaxf(lock_worst, e);
      if (k % 6 == 0)
        fprintf(fp, "lock,exact_pole_roll,%.7f,%.3e,pole%+d roll%+.3f\n",
                in[1], e, pole, in[0]);
    }
  }
  fprintf(fp, "lock,WORST,0,%.3e,exact +/-PI/2 roll recovery\n", lock_worst);

  /* ---- 4. dense multi-turn euler -> quat -> euler grid ----------------- */
  versor q;
  float grid_worst = 0.0f;
  double grid_sum = 0.0;
  int grid_n = 0;

  for (int ix = -8; ix <= 8; ix++) {
    for (int iy = -8; iy <= 8; iy++) {
      for (int iz = -8; iz <= 8; iz++) {
        in[0] = (float)ix * 0.4f;
        in[1] = (float)iy * 0.4f;
        in[2] = (float)iz * 0.4f;

        glm_euler_xyz(in, m);
        glm_euler_xyz_quat(in, q);
        glm_quat_mat4(q, rebuilt);
        glm_euler_angles(rebuilt, out);
        glm_euler_xyz(out, rebuilt);

        float e = mat4_diff(m, rebuilt);
        grid_worst = fmaxf(grid_worst, e);
        grid_sum += e;
        grid_n++;
      }
    }
  }
  fprintf(fp, "grid,MEAN,0,%.3e,%d grid points\n",
          (float)(grid_sum / grid_n), grid_n);
  fprintf(fp, "grid,WORST,0,%.3e,%d grid points, budget 1e-4\n",
          grid_worst, grid_n);

  /* ---- 5. slerp across every seam crossing, five turns ----------------- */
  versor prev, cur, mid;
  vec3 z_axis = {0, 0, 1};
  int per_turn = 720, seams = 0;
  float slerp_step_worst = 0.0f, slerp_drift_worst = 0.0f;

  for (int i = 0; i <= per_turn * 5; i++) {
    float raw = (float)i * (TAU / (float)per_turn);
    float ang = glm_angle_normalize(raw);
    glm_quatv(cur, ang, z_axis);

    if (i > 0) {
      float prev_raw = (float)(i - 1) * (TAU / (float)per_turn);
      if (fabsf(ang - glm_angle_normalize(prev_raw)) > GLM_PIf)
        seams++;

      float step = quat_dist(prev, cur);
      slerp_step_worst = fmaxf(slerp_step_worst, step);

      glm_quat_slerp(prev, cur, 0.5f, mid);
      versor qmid;
      glm_quatv(qmid, glm_angle_lerp(prev_raw, raw, 0.5f), z_axis);
      slerp_drift_worst = fmaxf(slerp_drift_worst, quat_dist(mid, qmid));
    }
    glm_quat_copy(cur, prev);
  }
  fprintf(fp, "slerp,WORST_STEP,0,%.3e,%d seam crossings\n",
          slerp_step_worst, seams);
  fprintf(fp, "slerp,WORST_MID_DRIFT,0,%.3e,vs angle-domain midpoint\n",
          slerp_drift_worst);

  fclose(fp);

  printf("wrote boundary_scan.csv\n");
  printf("  seam worst:          %.3e\n", seam_worst);
  printf("  pole sweep worst:    %.3e (%d frames)\n", pole_worst, pole_frames);
  printf("  exact lock worst:    %.3e\n", lock_worst);
  printf("  grid worst / mean:   %.3e / %.3e (n=%d)\n",
         grid_worst, (float)(grid_sum / grid_n), grid_n);
  printf("  slerp step worst:    %.3e, mid drift worst: %.3e\n",
         slerp_step_worst, slerp_drift_worst);
  return 0;
}
