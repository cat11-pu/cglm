#define _GNU_SOURCE
#include <cglm/cglm.h>
#include <math.h>
#include <stdio.h>

/* double-precision XYZ euler -> matrix -> euler reference */
static void
euler_xyz_d(const double e[3], double m[3][3]) {
  double cx = cos(e[0]), sx = sin(e[0]);
  double cy = cos(e[1]), sy = sin(e[1]);
  double cz = cos(e[2]), sz = sin(e[2]);
  m[0][0] = cy * cz; m[0][1] = cz * sx * sy + cx * sz; m[0][2] = -cx * cz * sy + sx * sz;
  m[1][0] = -cy * sz; m[1][1] = cx * cz - sx * sy * sz; m[1][2] = cz * sx + cx * sy * sz;
  m[2][0] = sy; m[2][1] = -cy * sx; m[2][2] = cx * cy;
}

static void
angles_d(const double m[3][3], double out[3]) {
  if (m[2][0] <= -1.0) {
    out[1] = -M_PI_2;
    out[0] = atan2(-m[1][0], m[0][0]);
    out[2] = 0.0;
  } else if (m[2][0] >= 1.0) {
    out[1] = M_PI_2;
    out[0] = atan2(-m[1][0], m[0][0]);
    out[2] = 0.0;
  } else {
    out[1] = asin(m[2][0]);
    out[0] = atan2(-m[2][1], m[2][2]);
    out[2] = atan2(-m[1][0], m[0][0]);
  }
}

int
main(void) {
  printf(" dist-to-pole   float-rt-err   double-rt-err   angle-extraction-err(float vs double)\n");
  for (int j = 6; j >= 0; j--) {
    float dist = powf(10.0f, -(float)j);
    float in_f[3] = {0.25f, GLM_PI_2f - dist, -0.4f};

    mat4 mf, rbf;
    vec3 of;
    glm_euler_xyz(in_f, mf);
    glm_euler_angles(mf, of);
    glm_euler_xyz(of, rbf);
    float ferr = 0.0f;
    for (int r = 0; r < 3; r++)
      for (int c = 0; c < 3; c++)
        ferr = fmaxf(ferr, fabsf(mf[r][c] - rbf[r][c]));

    double in_d[3] = {0.25, M_PI_2 - (double)dist, -0.4};
    double md[3][3], rbd[3][3], od[3];
    euler_xyz_d(in_d, md);
    angles_d(md, od);
    euler_xyz_d(od, rbd);
    double derr = 0.0;
    for (int r = 0; r < 3; r++)
      for (int c = 0; c < 3; c++)
        derr = fmax(derr, fabs(md[r][c] - rbd[r][c]));

    printf(" %+.1e       %.3e       %.3e        x=%.2e y=%.2e z=%.2e\n",
           (double)dist, ferr, derr,
           fabs(of[0] - od[0]), fabs(of[1] - od[1]), fabs(of[2] - od[2]));
  }
  return 0;
}
