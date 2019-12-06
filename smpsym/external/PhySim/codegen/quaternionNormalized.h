t1 = rot[0] * rot[0];
t2 = rot[1] * rot[1];
t3 = rot[2] * rot[2];
t5 = sqrt(t1 + t2 + t3);
t6 = t5 / 0.2e1;
t7 = cos(t6);
t8 = sin(t6);
t10 = t8 / t5;
cos(sqrt(rot[0] * rot[0] + rot[1] * rot[1] + rot[2] * rot[2]) / 0.2e1) = t7;
sin(sqrt(rot[0] * rot[0] + rot[1] * rot[1] + rot[2] * rot[2]) / 0.2e1) * pow(rot[0] * rot[0] + rot[1] * rot[1] + rot[2] * rot[2], -0.1e1 / 0.2e1) * rot[0] = t10 * rot[0];
sin(sqrt(rot[0] * rot[0] + rot[1] * rot[1] + rot[2] * rot[2]) / 0.2e1) * pow(rot[0] * rot[0] + rot[1] * rot[1] + rot[2] * rot[2], -0.1e1 / 0.2e1) * rot[1] = t10 * rot[1];
sin(sqrt(rot[0] * rot[0] + rot[1] * rot[1] + rot[2] * rot[2]) / 0.2e1) * pow(rot[0] * rot[0] + rot[1] * rot[1] + rot[2] * rot[2], -0.1e1 / 0.2e1) * rot[2] = t10 * rot[2];
quatout[0] = undefined;
quatout[1] = undefined;
quatout[2] = undefined;
quatout[3] = undefined;
t1 = rot[0] * rot[0];
t2 = rot[1] * rot[1];
t3 = rot[2] * rot[2];
t5 = sqrt(t1 + t2 + t3);
t6 = t5 / 0.2e1;
t7 = cos(t6);
t8 = sin(t6);
t10 = t8 / t5;
quat[0] = t7;
quat[1] = t10 * rot[0];
quat[2] = t10 * rot[1];
quat[3] = t10 * rot[2];
#include <math.h>

void quaternionNormalized (double *rot, double *quatOut)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  t4 = rot[0];
  t1 = t4 * t4;
  t9 = rot[1];
  t2 = t9 * t9;
  t11 = rot[2];
  t3 = t11 * t11;
  t5 = sqrt(t1 + t2 + t3);
  t6 = t5 / 0.2e1;
  t7 = cos(t6);
  t8 = sin(t6);
  t10 = t8 / t5;
  quatOut[0] = t7;
  quatOut[1] = t10 * t4;
  quatOut[2] = t10 * t9;
  quatOut[3] = t10 * t11;
}
#include <math.h>

void quaternionNormalized (double *rot, double *quatOut)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  t4 = rot[0];
  t1 = t4 * t4;
  t9 = rot[1];
  t2 = t9 * t9;
  t11 = rot[2];
  t3 = t11 * t11;
  t5 = sqrt(t1 + t2 + t3);
  t6 = t5 / 0.2e1;
  t7 = cos(t6);
  t8 = sin(t6);
  t10 = t8 / t5;
  quatOut[0] = t7;
  quatOut[1] = t10 * t4;
  quatOut[2] = t10 * t9;
  quatOut[3] = t10 * t11;
}
#include <math.h>

void quaternionNormalized (double *rot, double *quatOut)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  t4 = rot[0];
  t1 = t4 * t4;
  t9 = rot[1];
  t2 = t9 * t9;
  t11 = rot[2];
  t3 = t11 * t11;
  t5 = sqrt(t1 + t2 + t3);
  t6 = t5 / 0.2e1;
  t7 = cos(t6);
  t8 = sin(t6);
  t10 = t8 / t5;
  quatOut[0] = t7;
  quatOut[1] = t10 * t4;
  quatOut[2] = t10 * t9;
  quatOut[3] = t10 * t11;
}
#include <math.h>

void quaternionNormalized (double *rot, double *quatOut)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  t4 = rot[0];
  t1 = t4 * t4;
  t9 = rot[1];
  t2 = t9 * t9;
  t11 = rot[2];
  t3 = t11 * t11;
  t5 = sqrt(t1 + t2 + t3);
  t6 = t5 / 0.2e1;
  t7 = cos(t6);
  t8 = sin(t6);
  t10 = t8 / t5;
  quatOut[0] = t7;
  quatOut[1] = t10 * t4;
  quatOut[2] = t10 * t9;
  quatOut[3] = t10 * t11;
}
