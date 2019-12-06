#pragma once
inline void parallelTransportNormalized (
  double fromV[3],
  double toV[3],
  double x[3],
  double parallelTransportResult[3]) {
  double t1;
  double t10;
  double t11;
  double t15;
  double t16;
  double t2;
  double t20;
  double t21;
  double t25;
  double t3;
  double t32;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  t1 = fromV[0];
  t2 = toV[0];
  t3 = t1 * t2;
  t4 = fromV[1];
  t5 = toV[1];
  t6 = t4 * t5;
  t7 = fromV[2];
  t8 = toV[2];
  t9 = t7 * t8;
  t10 = t3 + t6 + t9;
  t11 = x[0];
  t15 = -t1 * t8 + t2 * t7;
  t16 = x[2];
  t20 = t1 * t5 - t2 * t4;
  t21 = x[1];
  t25 = t4 * t8 - t5 * t7;
  t32 = (t11 * t25 + t15 * t21 + t16 * t20) / (t3 + t6 + t9 + 0.1e1);
  parallelTransportResult[0] = t10 * t11 + t15 * t16 - t20 * t21 + t25 * t32;
  parallelTransportResult[1] = t10 * t21 + t11 * t20 + t15 * t32 - t16 * t25;
  parallelTransportResult[2] = t10 * t16 - t11 * t15 + t20 * t32 + t21 * t25;
  return;
}
