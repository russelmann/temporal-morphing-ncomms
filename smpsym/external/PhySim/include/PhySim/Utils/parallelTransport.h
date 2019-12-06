void parallelTransport (
  double *v0,
  double *v1,
  double *vectorOriginal,
  double cgret[3])
{
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
  double vectorTransported[3];
  vectorTransported[0] = 0;
  vectorTransported[1] = 0;
  vectorTransported[2] = 0;
  t1 = v0[0];
  t2 = v1[0];
  t3 = t1 * t2;
  t4 = v0[1];
  t5 = v1[1];
  t6 = t5 * t4;
  t7 = v0[2];
  t8 = v1[2];
  t9 = t7 * t8;
  t10 = t3 + t6 + t9;
  t11 = vectorOriginal[0];
  t15 = t7 * t2 - t1 * t8;
  t16 = vectorOriginal[2];
  t20 = t1 * t5 - t4 * t2;
  t21 = vectorOriginal[1];
  t25 = t4 * t8 - t7 * t5;
  t32 = (t25 * t11 + t15 * t21 + t20 * t16) / (0.1e1 + t3 + t6 + t9);
  vectorTransported[0] = t10 * t11 + t15 * t16 - t20 * t21 + t32 * t25;
  vectorTransported[1] = t10 * t21 + t20 * t11 - t25 * t16 + t32 * t15;
  vectorTransported[2] = t10 * t16 + t25 * t21 - t15 * t11 + t32 * t20;
  cgret[0] = vectorTransported[0];
  cgret[1] = vectorTransported[1];
  cgret[2] = vectorTransported[2];
}
