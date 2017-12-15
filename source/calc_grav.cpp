// Calcurate gravity using tree structre
// Using point to point
// In the future, I wanna use multi expansion

// include packages
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sstream>

// define namespace
using namespace std;
#define real double

// include my header files
#include "simulation_option.hpp"
#include "bhnode.hpp" // BHtree.hpp include "particle.hpp" and "myvector.hpp" and "constant.hpp"
#include "calc_grav.hpp"

void clear_acc_and_phi(real_particle * bp, int n){
  // 加速度と重力ポテンシャルの初期化
  real_particle * p_temp = bp;
  for(int i=0; i<n; i++){
    (p_temp+i)->clear_acc_phi_gravity();
  }
}
// main.cpp に移植中　なぜかここで定義すると遅くなる
// void accumulate_force_from_point(vector dx, real r2, real eps2,
//                         				 vector & acc, real & phi,real jmass){
//   // 質点からうける重力
//   double r2inv = 1.0/(r2+eps2);
//   double rinv  = sqrt(r2inv);
//   double r3inv = r2inv*rinv;
//   phi -= jmass*rinv;
//   acc += jmass*r3inv*dx;
// }
void accumulate_mutual_gravity(real_particle & p1, real_particle & p2, real eps2){
  // 直接計算法
  vector dx = p1.pos-p2.pos;
  double r2inv = 1.0/(dx*dx+eps2);
  double rinv  = sqrt(r2inv);
  double r3inv = r2inv*rinv;
  p1.acc_direct -= p2.mass*r3inv*dx;
  p2.acc_direct += p1.mass*r3inv*dx;
}
void calculate_gravity_direct(real_particle * pb, int n, double eps2){
  int i, j;
  real_particle * pi;
  real_particle * pj;
  for(i=0, pi=pb; i<n-1; i++,pi++){
    for(j=i+1, pj=pi+1; j<n; j++,pj++){
      accumulate_mutual_gravity(*pi, *pj, eps2);
    }
  }
}
