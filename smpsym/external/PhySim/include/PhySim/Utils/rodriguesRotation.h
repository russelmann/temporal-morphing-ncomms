#include <math.h>

inline void rodriguesRotation (
  double *k,
  double theta,
  double *v,
  double cgret[3])
{
  double r[3];
  double t1;
  double t13;
  double t19;
  double t2;
  double t4;
  double t5;
  double t6;
  double t8;
  double t9;
  r[0] = 0;
  r[1] = 0;
  r[2] = 0;
  t1 = cos(theta);
  t2 = v[0];
  t4 = sin(theta);
  t5 = k[1];
  t6 = v[2];
  t8 = k[2];
  t9 = v[1];
  t13 = k[0];
  t19 = (t13 * t2 + t5 * t9 + t6 * t8) * (0.1e1 - t1);
  r[0] = t1 * t2 + t4 * (t5 * t6 - t8 * t9) + t19 * t13;
  r[1] = t1 * t9 + t4 * (-t13 * t6 + t2 * t8) + t19 * t5;
  r[2] = t1 * t6 + t4 * (t13 * t9 - t2 * t5) + t19 * t8;
  cgret[0] = r[0];
  cgret[1] = r[1];
  cgret[2] = r[2];
}
