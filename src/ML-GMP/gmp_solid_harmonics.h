#include "gmp_helper.h"


typedef void (*SolidGMPFunction) (double, double, double, double, double, double, double, double, double *, double *);
SolidGMPFunction get_solid_gmp_function(int mcsh_order, int group_num);

void calc_solid_GMP_n1_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);

void calc_solid_GMP_0_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);

void calc_solid_GMP_1_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);

void calc_solid_GMP_2_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_2_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);

void calc_solid_GMP_3_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_3_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_3_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);

void calc_solid_GMP_4_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_4_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_4_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_4_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);

void calc_solid_GMP_5_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_5_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_5_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_5_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_5_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);

void calc_solid_GMP_6_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_6_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_6_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_6_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_6_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_6_6(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_6_7(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);

void calc_solid_GMP_7_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_7_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_7_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_7_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_7_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_7_6(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_7_7(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_7_8(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);

void calc_solid_GMP_8_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_8_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_8_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_8_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_8_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_8_6(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_8_7(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_8_8(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_8_9(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_8_10(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);

void calc_solid_GMP_9_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_9_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_9_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_9_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_9_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_9_6(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_9_7(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_9_8(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_9_9(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_9_10(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_9_11(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_solid_GMP_9_12(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
