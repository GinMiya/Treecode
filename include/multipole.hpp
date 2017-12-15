#pragma once

#include "simulation_option.hpp"  // シミュレーションの内容を変えるための define 群
#include "bhnode.hpp"             // "particle.hpp" "myvector.hpp" "constant.hpp"
#include "complex.hpp"

double Factorial1(int n);
double Factorial2(int n);
double Legendre(int m, int l, double x);
double Legendre_Cos(int l, int m, double x, double theta);

double NormalizedLegendre(int m, int l, double x);
double ScalingFactor(int l, int m);
complex ComputeShericalHarmonics(int l, int m, double theta, double phi);

void Y_term_ij(complex *Y, double theta, double phi);
vector calc_r_theta_phi_gradient(double r_i, double costh, double sinth, double sin_square, double phi);
vector th_ph_acc(double r_i, double theta, double phi, complex alpha[], complex Y_ij[]);
void get_rthph_acc(vector add_temp, double r_i, double theta, double phi, complex alpha[], complex Y_ij[]);
vector get_xyz_acc(vector acc_sph, double theta, double phi);
