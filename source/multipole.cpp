// 面倒なのでstdのnamespaceをつかう
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <assert.h>

#include "simulation_option.hpp"  // シミュレーションの内容を変えるための define 群
#include "bhnode.hpp"             // "particle.hpp" "myvector.hpp" "constant.hpp"
#include "complex.hpp"
#include "multipole.hpp"

using namespace std;
#define real double

double Factorial1(int n){
  if(n<=0) return 1.0;
  double F = 1.0;
  for(int i=n; i>=2; i--){
    F *= double(i);
  }
  return F;
}
double Factorial2(int n){
  if(n<=0) return 1.0;
  double F = 1.0;
  for(int i=n; i>=2; i=i-2){
    F *= double(i);
  }
  return F;
}
double Legendre(int l, int m, double x){
  int mm = fabs(m);
  if(mm>l) return 0;
  double r0, r1, r2;
  r0 = 0.0;
  r1 = pow(1.0-x*x, double(mm)/2.0) * Factorial2(2*mm-1);

  if(mm==l && m>=0) return r1;
  if(mm==l && m<0)  return r1 * pow(-1.0, mm) * Factorial1(l-mm)/Factorial1(l+mm);
  for(int ll=mm+1; ll<=l; ll++){
    r2 = ((2.0*ll-1.0)*x*r1 - (ll+mm-1.0)*r0)/(ll-mm);
    r0 = r1;
    r1 = r2;
  }
  if(m>=0) return r2;
  else return r2 * pow(-1.0,mm) * Factorial1(l-mm)/Factorial1(l+mm);
}
double Legendre_Cos(int l, int m, double x, double theta){
  int mm = fabs(m);
  if(mm>l) return 0;
  double r0, r1, r2;
  r0 = 0.0;
  r1 = sqrt((1.0+x)*(1.0-x));
  if(theta<0.0) { r1 *= -1.0;}
  r1 = pow(r1, mm);
  r1 *= Factorial2(2*mm-1);

  if(mm==l && m>=0) return r1;
  if(mm==l && m<0)  return r1*pow(-1.0, mm)*Factorial1(l-mm)/Factorial1(l+mm);
  for(int ll=mm+1; ll<=l; ll++){
    r2 = ((2.0*ll-1.0)*x*r1-(ll+mm-1.0)*r0)/(ll-mm);
    r0 = r1;
    r1 = r2;
  }
  if(m>=0) return r2;
  else return r2*pow(-1.0, mm)*Factorial1(l-mm)/Factorial1(l+mm);
}
double NormalizedLegendre(int l, int m, double x){
  int mm = abs(m);
  if(mm>l) return 0;
  double r0, r1, r2;
    r0 = 0.0;
    r1 = pow(1.0-x*x, double(mm)/2.0) * sqrt(1.0/2.0 * Factorial2(2*mm+1)/Factorial2(2*mm));
    if(mm==l && m>=0) return r1;
    if(mm==l && m<0)  return r1 * pow(-1.0,mm);
    for(int ll = mm+1; ll<=l; ll++){
      r2 = sqrt((2.0*ll-1.0)*(2.0*ll+1.0)/((ll+mm)*(ll-mm))) *x*r1 - sqrt((ll+mm-1.0)*(ll-mm-1.0)/((ll+mm)*(ll-mm)) * (2.0*ll+1.0)/(2.0*ll-3.0) ) *r0;
      r0 = r1;
      r1 = r2;
    }
  if(m>=0) return r2;
  else return r2 * pow(-1.0,mm);
}
double ScalingFactor(int l, int m){
  // (-1)^mはいれてないので注意
  double lpm = 1.0, lnm = 1.0;
  int mm = abs(m);
  lnm = Factorial1(l-fabs(m));
  lpm = Factorial1(l+fabs(m));
  // return sqrt(((2.0*l+1.0)*lnm)/(4.0*M_PI*lpm));
  return sqrt(lnm/lpm);
}
complex ComputeShericalHarmonics(int l, int m, double theta, double phi){
  assert(l >= 0);
  assert((-l<=m) && (m<=l));
  int mm = fabs(m);
  double m_phi = mm*phi;
  double factor = pow(-1.0, (mm+m)/2) * ScalingFactor(l, mm);
  complex temp;
  if(m>0){
    temp.set_real(factor*cos(m_phi)*Legendre_Cos(l, mm, cos(theta), theta));
    temp.set_imag(factor*sin(m_phi)*Legendre_Cos(l, mm, cos(theta), theta));
  }else if(m<0){
    temp.set_real( factor*cos(m_phi)*Legendre_Cos(l, mm, cos(theta), theta));
    temp.set_imag(-factor*sin(m_phi)*Legendre_Cos(l, mm, cos(theta), theta));
  }else{
    // m=0のときφ依存性は消えてθ依存する実部のみが現れる
    temp.set_real(factor*Legendre_Cos(l, mm, cos(theta), theta));
    // temp.set_imag(factor*Legendre_Cos(l, mm, cos(theta), theta));
    temp.set_imag(0.0);
  }
  return temp;
}
void Y_term_ij(complex *Y, double theta, double phi){
  complex temp;
  for(int l=0; l<=p_max; l++){
    for(int m=-l, j=0; m<=l; m++, j++){
      temp = ComputeShericalHarmonics(l, m, theta, phi);
      Y[l*l+j].set_real(temp.get_real());
      Y[l*l+j].set_imag(temp.get_imag());
    }
  }
}

// vector th_ph_grad(double r, double theta, double phi, complex alpha[], complex Y_ij[]){
//   // θ微分
//   vector acc_sph = 0.0;
//   double br, coeff, coeff1, coeff2;
//   double common_term, legendre_1, legendre_2;
//   int mm;
//   for(int l=0; l<=p_max; l++){
//     // br = 4.0*M_PI*G0*pow(1.0/r, l+2)/(2.0*l+1.0);
//     br = pow(1.0/r, l+2);
//     for(int m=-l, j=0; m<=l; m++, j++){
//       if(m==0) {coeff = 1.0;} else{coeff=2.0;}
//       mm = fabs(m);
//       common_term = pow(-1.0, (mm+m)/2) * ScalingFactor(l, mm);
//       legendre_1 = common_term * Legendre_Cos(l, fabs(m), cos(theta), theta);
//       legendre_2 = common_term * Legendre_Cos(l, fabs(m)+1, cos(theta), theta);
//       coeff  = ((cos(theta)*fabs(m)*legendre_1)/sin(theta)-legendre_2);
//       // 実部虚部成分
//       coeff1 = coeff*cos(m*phi);
//       coeff2 = coeff*sin(m*phi);
//       acc_sph[1] -= br * (alpha[l*l+j].get_real()*coeff1
//                         + alpha[l*l+j].get_imag()*coeff2);
//       acc_sph[2] -= br * (m/sin(theta)) * (-alpha[l*l+j].get_real()*Y_ij[l*l+j].get_imag()
//                         + alpha[l*l+j].get_imag()*Y_ij[l*l+j].get_real());
//     }
//   }
//   return acc_sph;
// }
// void get_rthph_acc(vector add_temp, double r_i, double theta, double phi, complex alpha[], complex Y_ij[]){
//   // 極座標での加速度を得る関数
//   add_temp += th_ph_grad(r_i, theta, phi, alpha, Y_ij);
// }
// vector get_xyz_acc(vector acc_sph, double theta, double phi){
//   // カーテシアンでの加速度を得る関数
//   myvector acc_xyz = 0.0;
//   acc_xyz[0] = sin(theta)*cos(phi)*acc_sph[0]+cos(theta)*cos(phi)*acc_sph[1]-sin(phi)*acc_sph[2];
//   acc_xyz[1] = sin(theta)*sin(phi)*acc_sph[0]+cos(theta)*sin(phi)*acc_sph[1]+cos(phi)*acc_sph[2];
//   acc_xyz[2] = cos(theta)*acc_sph[0]-sin(theta)*acc_sph[1];
//   return acc_xyz;
// }
