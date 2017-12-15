// Calcurate energy Potential,Kinetic,Momentum,Relative error
// ライブラリ読み込み
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>

// define namespace
using namespace std;
#define real double

// ヘッダーファイルの読み込み
#include "simulation_option.hpp"
#include "bhnode.hpp" // BHtree.hpp include "particle.hpp" and "myvector.hpp" and "constant.hpp"

real real_system::kinetic_energy(int n, real r0){
  real ke = 0.0;
  #if defined(COLD_COLLAPSE) || defined(PLUMMER_ONLY) || defined(INCLUDE_DATASET)
    for(int i=0; i<n; i++){
      ke += (pb+i)->get_ke();
    }
  #else
    for(int i=0; i<n; i++){
      real r2 = (pb+i)->get_radius_cm();
      if(r2 <= r0){
        ke += (pb+i)->get_ke_rcut();
      }
    }
  #endif
  return ke;
}
real real_system::potential_energy(int n, double r0){
  real pe = 0.0;
  #if defined(COLD_COLLAPSE) || defined(PLUMMER_ONLY) || defined(INCLUDE_DATASET)
    for(int i=0; i<n; i++){
      pe += (pb+i)->get_pe();
    }
  #else
    for(int i=0; i<n; i++){
      real r2 = (pb+i)->get_radius_cm();
      if(r2 <= r0){
        pe += (pb+i)->get_pe();
      }
    }
  #endif
  return pe*0.5;
}
real real_system::momentum(int n, real r0){
  real mv = 0.0;
  #if defined(COLD_COLLAPSE) || defined(PLUMMER_ONLY) || defined(INCLUDE_DATASET)
    for(int i=0; i<n; i++){
      mv += (pb+i)->get_momentum();
    }
  #else
    for(int i=0; i<n; i++){
      real r2 = (pb+i)->get_radius_cm();
      if(r2 <= r0){
        mv += (pb+i)->get_momentum();
      }
    }
  #endif
  return mv;
}
real real_system::energy_error(int n, int check, real ke, real pe){
  static real e_tot0;
  real de;
  real e_tot = pe+ke;
  if(check == 0){
    e_tot0 = e_tot;
    return e_tot;
  }
  if(check){
    de = fabs((e_tot-e_tot0)/e_tot0);
    cout << "relative error of energy : " << de << endl;
    return de;
  }
}
