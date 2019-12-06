#include <math.h>

void getRBNodeJacobianQuat (
  const double *vec,
  const double *cen,
  const double *rot,
  double *nodeJacobian)
{
  double t1;
  double t10;
  double t102;
  double t106;
  double t109;
  double t11;
  double t111;
  double t113;
  double t115;
  double t116;
  double t118;
  double t119;
  double t12;
  double t121;
  double t122;
  double t124;
  double t125;
  double t128;
  double t129;
  double t13;
  double t132;
  double t14;
  double t140;
  double t143;
  double t145;
  double t146;
  double t149;
  double t15;
  double t150;
  double t153;
  double t16;
  double t17;
  double t18;
  double t19;
  double t2;
  double t20;
  double t21;
  double t22;
  double t25;
  double t27;
  double t3;
  double t30;
  double t31;
  double t32;
  double t33;
  double t34;
  double t35;
  double t36;
  double t37;
  double t4;
  double t41;
  double t42;
  double t43;
  double t45;
  double t47;
  double t48;
  double t49;
  double t5;
  double t52;
  double t53;
  double t54;
  double t55;
  double t58;
  double t59;
  double t6;
  double t62;
  double t64;
  double t65;
  double t66;
  double t7;
  double t70;
  double t73;
  double t75;
  double t76;
  double t77;
  double t8;
  double t80;
  double t82;
  double t85;
  double t86;
  double t87;
  double t9;
  double t90;
  double t91;
  double t92;
  double t95;
  double t96;
  double t99;
  t8 = rot[0];
  t1 = t8 * t8;
  t13 = rot[1];
  t2 = t13 * t13;
  t14 = rot[2];
  t3 = t14 * t14;
  t4 = t1 + t2 + t3;
  t5 = sqrt(t4);
  t6 = t5 / 0.2e1;
  t7 = sin(t6);
  t18 = 0.1e1 / t5;
  t20 = 0.1e1 / t4;
  t9 = t18 * t20;
  t10 = t7 * t9;
  t11 = t2 * t8;
  t12 = cos(t6);
  t15 = 0.2e1 * t10 * t11 * t12;
  t16 = t7 * t7;
  t17 = t4 * t4;
  t19 = t16 / t17;
  t21 = 0.4e1 * t19 * t11;
  t22 = t3 * t8;
  t25 = 0.2e1 * t10 * t22 * t12;
  t27 = 0.4e1 * t19 * t22;
  t30 = t20;
  t31 = t16 * t30;
  t32 = t8 * t14;
  t33 = t31 * t32;
  t34 = t12 * t12;
  t35 = t34 * t30;
  t36 = t35 * t32;
  t37 = t12 * t7;
  t43 = t37 * t9 * t14;
  t41 = 0.2e1 * t43 * t8;
  t42 = t1 * t13;
  t45 = 0.2e1 * t10 * t42 * t12;
  t47 = 0.4e1 * t19 * t42;
  t48 = t31 * t13;
  t49 = 0.2e1 * t48;
  t52 = t8 * t13;
  t53 = t31 * t52;
  t54 = t35 * t52;
  t55 = t37 * t9;
  t58 = 0.2e1 * t55 * t52;
  t59 = t1 * t14;
  t62 = 0.2e1 * t10 * t59 * t12;
  t64 = 0.4e1 * t19 * t59;
  t65 = t31 * t14;
  t66 = 0.2e1 * t65;
  t70 = t2 * t13;
  t73 = 0.2e1 * t10 * t70 * t12;
  t75 = 0.4e1 * t19 * t70;
  t76 = 0.4e1 * t48;
  t77 = t3 * t13;
  t80 = 0.2e1 * t10 * t77 * t12;
  t82 = 0.4e1 * t19 * t77;
  t85 = t13 * t14;
  t86 = t31 * t85;
  t87 = t35 * t85;
  t90 = 0.2e1 * t43 * t13;
  t91 = t31 * t8;
  t92 = 0.2e1 * t91;
  t95 = t31 * t2;
  t96 = t35 * t2;
  t99 = 0.2e1 * t55 * t2;
  t102 = 0.2e1 * t37 * t18;
  t106 = 0.2e1 * t10 * t8 * t85 * t12;
  t109 = 0.4e1 * t19 * t32 * t13;
  t113 = t2 * t14;
  t116 = 0.2e1 * t10 * t113 * t12;
  t118 = 0.4e1 * t19 * t113;
  t119 = t3 * t14;
  t122 = 0.2e1 * t10 * t119 * t12;
  t124 = 0.4e1 * t19 * t119;
  t125 = 0.4e1 * t65;
  t128 = t31 * t3;
  t129 = t35 * t3;
  t132 = 0.2e1 * t55 * t3;
  t140 = t1 * t8;
  t143 = 0.2e1 * t10 * t140 * t12;
  t145 = 0.4e1 * t19 * t140;
  t146 = 0.4e1 * t91;
  t149 = t31 * t1;
  t150 = t35 * t1;
  t153 = 0.2e1 * t55 * t1;
  nodeJacobian[0] = 1;
  nodeJacobian[1] = 0;
  nodeJacobian[2] = 0;
  t111 = vec[0];
  t115 = vec[1];
  t121 = vec[2];
  nodeJacobian[3] = ((-t15 + t21 - t25 + t27) * t111 + (-t33 + t36 - t41 + t45 - t47 + t49) * t115 + (t53 - t54 + t58 + t62 - t64 + t66) * t121);
  nodeJacobian[4] = ((-t73 + t75 - t76 - t80 + t82) * t111 + (-t86 + t87 - t90 + t15 - t21 + t92) * t115 + (t95 - t96 + t99 - t102 + t106 - t109) * t121);
  nodeJacobian[5] = ((-t116 + t118 - t122 + t124 - t125) * t111 + (-t128 + t129 - t132 + t102 + t106 - t109) * t115 + (t86 - t87 + t90 + t25 - t27 + t92) * t121);
  nodeJacobian[6] = 0;
  nodeJacobian[7] = 1;
  nodeJacobian[8] = 0;
  nodeJacobian[9] = ((t33 - t36 + t41 + t45 - t47 + t49) * t111 + (-t143 + t145 - t146 - t25 + t27) * t115 + (-t149 + t150 - t153 + t102 + t106 - t109) * t121);
  nodeJacobian[10] = ((t86 - t87 + t90 + t15 - t21 + t92) * t111 + (-t45 + t47 - t80 + t82) * t115 + (-t53 + t54 - t58 + t116 - t118 + t66) * t121);
  nodeJacobian[11] = ((t128 - t129 + t132 - t102 + t106 - t109) * t111 + (-t62 + t64 - t122 + t124 - t125) * t115 + (-t33 + t36 - t41 + t80 - t82 + t49) * t121);
  nodeJacobian[12] = 0;
  nodeJacobian[13] = 0;
  nodeJacobian[14] = 1;
  nodeJacobian[15] = ((-t53 + t54 - t58 + t62 - t64 + t66) * t111 + (t149 - t150 + t153 - t102 + t106 - t109) * t115 + (-t143 + t145 - t146 - t15 + t21) * t121);
  nodeJacobian[16] = ((-t95 + t96 - t99 + t102 + t106 - t109) * t111 + (t53 - t54 + t58 + t116 - t118 + t66) * t115 + (-t45 + t47 - t73 + t75 - t76) * t121);
  nodeJacobian[17] = ((-t86 + t87 - t90 + t25 - t27 + t92) * t111 + (t33 - t36 + t41 + t80 - t82 + t49) * t115 + (-t62 + t64 - t116 + t118) * t121);
}
