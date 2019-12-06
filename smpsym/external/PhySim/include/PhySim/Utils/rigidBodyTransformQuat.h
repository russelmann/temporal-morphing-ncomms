#include <math.h>

void rigidBodyTransformQuat (
  const double *vec,
  const double *cen,
  const double *rot,
  double *pos)
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t13;
  double t14;
  double t17;
  double t2;
  double t21;
  double t22;
  double t23;
  double t24;
  double t27;
  double t29;
  double t3;
  double t31;
  double t34;
  double t36;
  double t39;
  double t4;
  double t40;
  double t42;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  t9 = rot[0];
  t1 = t9 * t9;
  t11 = rot[1];
  t2 = t11 * t11;
  t13 = rot[2];
  t3 = t13 * t13;
  t4 = t1 + t2 + t3;
  t5 = sqrt(t4);
  t6 = t5 / 0.2e1;
  t7 = sin(t6);
  t8 = t7 * t7;
  t10 = t8 / t4;
  t12 = 0.2e1 * t10 * t2;
  t14 = 0.2e1 * t10 * t3;
  t17 = cos(t6);
  t22 = t17 * t7 / t5;
  t21 = t22 * t13;
  t24 = t10 * t9;
  t23 = t24 * t11;
  t27 = t22 * t11;
  t29 = t24 * t13;
  t36 = 0.2e1 * t1 * t10;
  t40 = t22 * t9;
  t42 = t10 * t11 * t13;
  t31 = vec[0];
  t34 = vec[1];
  t39 = vec[2];
  pos[0] = cen[0] + (0.1e1 - t12 - t14) * t31 + 0.2e1 * (t21 + t23) * t34 + 0.2e1 * (-t27 + t29) * t39;
  pos[1] = cen[1] + 0.2e1 * (-t21 + t23) * t31 + (0.1e1 - t36 - t14) * t34 + 0.2e1 * (t40 + t42) * t39;
  pos[2] = cen[2] + 0.2e1 * (t27 + t29) * t31 + 0.2e1 * (-t40 + t42) * t34 + (0.1e1 - t36 - t12) * t39;
}
