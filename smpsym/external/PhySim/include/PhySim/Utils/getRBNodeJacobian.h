#include <math.h>

void getRBNodeJacobian (
  const double *trans,
  const double *eulerR,
  double *nodeJacobian)
{
  double t1;
  double t10;
  double t13;
  double t16;
  double t17;
  double t18;
  double t19;
  double t2;
  double t20;
  double t21;
  double t24;
  double t26;
  double t29;
  double t3;
  double t33;
  double t34;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  t5 = eulerR[2];
  t1 = cos(t5);
  t8 = eulerR[1];
  t2 = sin(t8);
  t3 = t1 * t2;
  t10 = eulerR[0];
  t4 = cos(t10);
  t6 = sin(t5);
  t7 = sin(t10);
  t9 = t3 * t4 + t6 * t7;
  t13 = -t3 * t7 + t4 * t6;
  t17 = cos(t8);
  t18 = t1 * t17;
  t16 = trans[1];
  t19 = t7 * t16;
  t20 = trans[2];
  t21 = t4 * t20;
  t24 = t6 * t17;
  t26 = t6 * t2;
  t29 = -t1 * t4 - t26 * t7;
  t33 = t1 * t7 - t26 * t4;
  nodeJacobian[0] = 1;
  nodeJacobian[1] = 0;
  nodeJacobian[2] = 0;
  nodeJacobian[3] = (t13 * t20 + t16 * t9);
  t34 = trans[0];
  nodeJacobian[4] = (t19 * t18 + t18 * t21 - t3 * t34);
  nodeJacobian[5] = (t16 * t29 + t20 * t33 - t24 * t34);
  nodeJacobian[6] = 0;
  nodeJacobian[7] = 1;
  nodeJacobian[8] = 0;
  nodeJacobian[9] = (-t16 * t33 + t20 * t29);
  nodeJacobian[10] = (t19 * t24 + t21 * t24 - t26 * t34);
  nodeJacobian[11] = (-t13 * t16 + t18 * t34 + t20 * t9);
  nodeJacobian[12] = 0;
  nodeJacobian[13] = 0;
  nodeJacobian[14] = 1;
  nodeJacobian[15] = (t16 * t17 * t4 - t17 * t20 * t7);
  nodeJacobian[16] = (-t16 * t2 * t7 - t2 * t20 * t4 - t17 * t34);
  nodeJacobian[17] = 0;
}
