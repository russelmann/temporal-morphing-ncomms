#include <math.h>

void rigidBodyTransform (
  const double *trans,
  const double *center,
  const double *eulerR,
  double *position)
{
  double t1;
  double t10;
  double t13;
  double t18;
  double t2;
  double t21;
  double t24;
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  t3 = eulerR[2];
  t1 = cos(t3);
  t4 = eulerR[1];
  t2 = cos(t4);
  t5 = sin(t4);
  t6 = t1 * t5;
  t8 = eulerR[0];
  t7 = sin(t8);
  t9 = sin(t3);
  t10 = cos(t8);
  t21 = t9 * t5;
  t13 = trans[0];
  t18 = trans[1];
  t24 = trans[2];
  position[0] = center[0] + t1 * t2 * t13 + (-t10 * t9 + t6 * t7) * t18 + (t10 * t6 + t7 * t9) * t24;
  position[1] = center[1] + t9 * t2 * t13 + (t1 * t10 + t21 * t7) * t18 + (-t1 * t7 + t10 * t21) * t24;
  position[2] = t10 * t2 * t24 + t18 * t2 * t7 - t13 * t5 + center[2];
}
