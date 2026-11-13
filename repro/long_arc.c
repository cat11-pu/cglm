#include <cglm/cglm.h>
#include <math.h>
#include <stdio.h>

static float dist(const versor a, const versor b) {
  float d = fabsf(glm_quat_dot(a, b));
  return 2.0f * acosf(d > 1 ? 1 : d);
}

int main(void) {
  versor q1, q2, r;
  vec3 z = {0, 0, 1};
  glm_quatv(q1, 0.3f, z);
  glm_quatv(q2, 1.0f, z);
  float sep = dist(q1, q2);

  glm_quat_slerp(q1, q2, 0.5f, r);
  printf("short arc mid: dist from q1=%.4f (sep=%.4f)\n", dist(q1, r), sep);

  glm_quat_slerp_longest(q1, q2, 0.5f, r);
  printf("long  arc mid: dist from q1=%.4f (expect ~PI-%.4f=%.4f)\n",
         dist(q1, r), sep / 2, GLM_PIf - sep / 2);

  glm_quat_slerp_longest(q1, q2, 0.0f, r);
  printf("long t=0 endpoint: %.3e (expect 0)\n", dist(q1, r));
  glm_quat_slerp_longest(q1, q2, 1.0f, r);
  printf("long t=1 endpoint: %.3e (expect 0)\n", dist(q2, r));
  return 0;
}
