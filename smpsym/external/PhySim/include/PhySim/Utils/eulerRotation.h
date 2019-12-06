#include <math.h>

void eulerRotation (const double *eulerAngles, const double *vectorOriginal, double cgret[3])
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
  double vectorRotated[3];
  vectorRotated[0] = 0;
  vectorRotated[1] = 0;
  vectorRotated[2] = 0;
  t1 = eulerAngles[2];
  t2 = cos(t1);
  t3 = eulerAngles[1];
  t4 = cos(t3);
  t6 = vectorOriginal[0];
  t8 = sin(t3);
  t9 = t2 * t8;
  t10 = eulerAngles[0];
  t11 = sin(t10);
  t13 = sin(t1);
  t14 = cos(t10);
  t17 = vectorOriginal[1];
  t22 = vectorOriginal[2];
  vectorRotated[0] = t2 * t4 * t6 + (t9 * t11 - t13 * t14) * t17 + (t9 * t14 + t13 * t11) * t22;
  t26 = t13 * t8;
  vectorRotated[1] = t13 * t4 * t6 + (t26 * t11 + t2 * t14) * t17 + (t26 * t14 - t2 * t11) * t22;
  vectorRotated[2] = -t8 * t6 + t4 * t11 * t17 + t4 * t14 * t22;
  cgret[0] = vectorRotated[0];
  cgret[1] = vectorRotated[1];
  cgret[2] = vectorRotated[2];
}