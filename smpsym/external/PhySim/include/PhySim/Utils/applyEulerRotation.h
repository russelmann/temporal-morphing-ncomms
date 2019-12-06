#include <math.h>

inline void applyEulerRotation (double *eulerAngles, double *eulerFromVector, double cgret[3])
{
  double t1;
  double t10;
  double t11;
  double t13;
  double t14;
  double t17;
  double t2;
  double t22;
  double t26;
  double t3;
  double t4;
  double t6;
  double t8;
  double t9;
  double eulerRotationApplied[3];
  eulerRotationApplied[0] = 0;
  eulerRotationApplied[1] = 0;
  eulerRotationApplied[2] = 0;
  t1 = eulerAngles[2];
  t2 = cos(t1);
  t3 = eulerAngles[1];
  t4 = cos(t3);
  t6 = eulerFromVector[0];
  t8 = sin(t3);
  t9 = t2 * t8;
  t10 = eulerAngles[0];
  t11 = sin(t10);
  t13 = sin(t1);
  t14 = cos(t10);
  t17 = eulerFromVector[1];
  t22 = eulerFromVector[2];
  eulerRotationApplied[0] = t2 * t4 * t6 + (t11 * t9 - t13 * t14) * t17 + (t11 * t13 + t14 * t9) * t22;
  t26 = t13 * t8;
  eulerRotationApplied[1] = t13 * t4 * t6 + (t11 * t26 + t14 * t2) * t17 + (-t11 * t2 + t14 * t26) * t22;
  eulerRotationApplied[2] = t11 * t17 * t4 + t14 * t22 * t4 - t6 * t8;
  cgret[0] = eulerRotationApplied[0];
  cgret[1] = eulerRotationApplied[1];
  cgret[2] = eulerRotationApplied[2];
}
double functionname (double eulerAngles, double eulerFromVector)
{
  return(eulerRotationApplied);
}
