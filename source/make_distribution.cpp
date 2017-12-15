// Make particle distribution
// You can use Plummer/Cold Collapse/Include from File

// ライブラリ読み込み
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <omp.h>

// define namespace
using namespace std;
#define real double

// ヘッダーファイルの読み込み
#include "simulation_option.hpp"
#include "bhnode.hpp" // BHtree.hpp include "particle.hpp" and "myvector.hpp" and "constant.hpp"
#include "morton.hpp"
#include "calc_grav.hpp"
#include "make_distribution.hpp"

// =============COLD_COLLAPSE===============
const vector uniform_random_position_in_sphere(double power_index){
  vector x;
  do{
	  for(int i=0; i<3; i++) x[i] = drand48()*2-1;
  }while(x*x >= 1);
  x *=  pow(x*x, 3.0/(power_index+3)-1);
  return x;
}
void nbody_system::create_uniform_sphere(int nbody, real power_index, real r0){
  n = nbody;
  pb = new real_particle[n];
  nsize = n;
  real_particle * p = pb;
  for(int i=0; i<n; i++){
    p->set_pos(uniform_random_position_in_sphere(power_index)*r0);
    p->set_vel(0.0);
    p->set_mass(1.0/n);
    p->set_index(i);
    p++;
  }
}
real gaussian(void){
	double x,y,z,r2;
	do{
		x=2.0*drand48()-1.0;
		y=2.0*drand48()-1.0;
		r2=(x*x)+(y*y);
	}while (r2>=1.0 || r2==0.0);
	z=sqrt(-2.0*log(r2)/r2)*x;
	return (z);
}
const vector gaussian_velocity(real v_rms){
  return v_rms*gaussian();
}
void nbody_system::create_coldcollapse(int nbody, real power_index, real r0, real vir0, real eps2){
  n = nbody;
  pb = new real_particle[n];
  nsize = n;
  real_particle * p = pb;
	for(int i=0; i<n; i++){
    (p+i)->set_pos(uniform_random_position_in_sphere(power_index)*r0);
    (p+i)->set_vel(0.0);
    (p+i)->set_mass(1.0/n);
    (p+i)->set_index(i);
	}
  real imass = M0/n;
  real pot = 0.0;
  vector r_ij, pos1, pos2;
  for(int i=0; i<(n-1); i++){
    pos1 = (p+i)->get_pos();
		for(int j=i+1; j<n; j++){
      pos2 = (p+i)->get_pos();
			r_ij = pos2 - pos1;
			pot += -(imass*imass)/sqrt(r_ij*r_ij+eps2);
		}
	}

	real v_rms = sqrt((2.0*vir0*fabs(pot))/(3.0*imass));
	for(int i=0; i<n; i++){
    (p+i)->set_vel(gaussian_velocity(v_rms));
	}
}

// =============INCLUDE_DATASET===============
void real_system::INCLUDE_PARTICLE_DATASET(int nbody){
  int i = 0;
  n = nbody;
  pb = new real_particle[n];
  nsize = n;
  real_particle * p = pb;
  vector position = 0.0;
  vector velocity = 0.0;
  vector dposition = 0.0;
  cerr << "Input initial particle dataset from filename 'initial_dataset.dat'\n";
  ifstream ifs_ini("./include/initial_dataset.dat");
  string str;
  // 3列のデータセットを再利用するとき（位置、速度、重心からの位置）
  // while(getline(ifs_ini, str)){
  //   sscanf(str.data(), "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
  //                       &position[0], &position[1], &position[2],
  //                       &velocity[0], &velocity[1], &velocity[2],
  //                       &dposition[0], &dposition[1], &dposition[2]);
  //   (p+i)->pos  = position;
  //   (p+i)->vel  = velocity;
  //   (p+i)->dpos = dposition;
  //   (p+i)->index = i;
  //   (p+i)->mass = M0/nbody;
  //   i++;
  // }

  // ２列のデータセットを再利用するとき（位置、速度）
  while(getline(ifs_ini, str)){
    sscanf(str.data(), "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                        &position[0], &position[1], &position[2],
                        &velocity[0], &velocity[1], &velocity[2]);
    (p+i)->pos  = position;
    (p+i)->vel  = velocity;
    (p+i)->index = i;
    (p+i)->mass = M0/nbody;
    i++;
  }
}


// =============PLUMMER_SPHERE===============
// Plummer球のHill半径からずらす角度を計算する関数
void shift_theta_hill(int nplummer, double shift_theta[], double shift_radius[]){
  double r_hill = fabs(shift_radius[1]-shift_radius[0])*pow((1.0/3.0),(1.0/3.0));
  shift_theta[0] = 0.0;
  for(int I=1; I<nplummer; I++){
    shift_theta[I] = acos(((shift_radius[I-1]*shift_radius[I-1]+shift_radius[I]*shift_radius[I])-(2.0*r_hill)*(2.0*r_hill))
                     /(2.0*shift_radius[I-1]*shift_radius[I]))*(180.0/M_PI);
  }
  cerr << "hill radius[kpc]=" << r_hill <<endl;
}
//全ての粒子で計算した重心だけずらす関数
void shift_pos_cm(real_particle * bp, int nbody){
  real_particle * p = bp;
  vector cmpos = 0.0;
  real cmmass = 0.0;
  for(int i=0; i<nbody; i++){
    cmpos  += (p+i)->get_mass()*(p+i)->get_pos();
    cmmass += (p+i)->get_mass();
  }
  cmpos /= cmmass;
  for(int i=0; i<nbody; i++){
    (p+i)->get_pos() -= cmpos;
  }
}
//中心からr_cut以下の粒子だけずらす
void shift_pos_cm_radius(real_particle * bp, int I, int nbody, double r_cut){
  real_particle * p = bp+I*nbody;
  vector cmpos = 0.0;
  real cmmass = 0.0;
  for(int i=0; i<nbody; i++){
    real r = (p+i)->get_radius();
    if(r<=r_cut){
      cmpos  += (p+i)->get_mass()*(p+i)->get_pos();
      cmmass += (p+i)->get_mass();
    }
  }
  cmpos /= cmmass;
  for(int i=0; i<nbody; i++){
    (p+i)->get_pos() -= cmpos;
  }
}
//shift_theta/shift_radiusだけずらす関数
void shift_pos_theta_radius(real_particle * bp, int I, int nbody,
                            double shift_theta[], double shift_radius[]){
  real_particle * p = bp+I*nbody;
  vector dpos;
  dpos[0] = shift_radius[I]*cos(shift_theta[I]*M_PI/180.0);
  dpos[1] = shift_radius[I]*sin(shift_theta[I]*M_PI/180.0);
  dpos[2] = 0.0;
  for(int i=0; i<nbody; i++){
    (p+i)->inc_pos(dpos);
  }
}
// ランダム値を返す関数
real randunit(int seed){
  const real MAXN = 2147483647;  // the maximum value which rand() can return
  static int randx;
  if (seed){
    randx = seed;
    return(0);        // to make the compiler happy, we return a real value
  }else{
    return((real)((randx = randx * 1103515245 + 12345) & 0x7fffffff)/MAXN);
  }
}
real randinter(real a, real b){
  return(a + (b-a)*randunit(0));
}
// Plummer sphereの位置を返す関数
const vector uniform_position_in_plummer(){
  vector x;
  double bx[3];
  for(int k=0; k<3; k++){
    bx[k]=0.0; bx[k]=randinter(0,1);
  }
  double r_plummer = 1.0/sqrt(1.0/pow(bx[0],(2.0/3.0))-1.0);
  // double r_plummer = fabs(1.0/((bx[0]*bx[0])-1.0));
  x[2] = (1.0-2.0*bx[1])*r_plummer;
  x[0] = sqrt(pow(r_plummer,2.0)-pow(x[2],2.0))*cos(2.0*M_PI*bx[2]);
  x[1] = sqrt(pow(r_plummer,2.0)-pow(x[2],2.0))*sin(2.0*M_PI*bx[2]);
  return x;
}
// Plummer sphereの速度を返す関数
const vector uniform_velocity_in_plummer(vector pos){
  vector v;
  double bv[4];
  bv[0] = 0.0; bv[1] = 0.1;
  while((pow(bv[0],2.0)*pow((1.0-pow(bv[0],2.0)),(7.0/2.0)))<0.1*bv[1]){
    bv[0] = 0.0;
    bv[1] = 0.0;
    bv[0] = randinter(0,1);
    bv[1] = randinter(0,1);
  }
  double V = bv[0]*sqrt(2.0)/pow((1.0+(pow(pos[0],2.0)+pow(pos[1],2.0)+pow(pos[2],2.0))),(0.25));
  for(int k=2; k<4; k++){
    bv[k] = 0.0; bv[k] = randinter(0,1);
  }
  v[2] = (1.0-2.0*bv[2])*V;
  v[0] = sqrt(pow(V,2.0)-pow(v[2],2.0))*sin(2.0*M_PI*bv[3]);
  v[1] = sqrt(pow(V,2.0)-pow(v[2],2.0))*cos(2.0*M_PI*bv[3]);
  return v;
}
// Plummer sphereを作る関数
void nbody_system::create_plummer_sphere(int nplummer, int nbody, double shift_theta[],
                                         double shift_radius[]){
  // srand48((int)time(NULL));
  srand48(1223L);
  n = nplummer*nbody;
  pb = new real_particle[n];
  nsize = n;
  real_particle * p = pb;
  // I個のPlummer球を作成
  for(int I=0; I<nplummer; I++){
    for(int i=0; i<nbody; i++){
      (p+I*nbody+i)->pos  = uniform_position_in_plummer();
      (p+I*nbody+i)->mass = M0/nbody;
      (p+I*nbody+i)->set_index(i);
    }
    // 1kpc以内の粒子で重心をとって移動
    shift_pos_cm_radius(p, I, nbody, 1.0);
    // 重心移動したあとの位置でもって速度の計算
    for(int i=0; i<nbody; i++){
      (p+I*nbody+i)->vel = uniform_velocity_in_plummer((p+I*nbody+i)->pos);
    }
    // offset分だけ動径方向に移動させる
    shift_pos_theta_radius(p, I, nbody, shift_theta, shift_radius);
  }
}


// NFW profile の DMH についての Circular velocity と加速度
void NFW_circular_velocity(real_particle * pb, int nbody){
  real_particle * p = pb;
  for(int i=0; i<nbody; i++){
    vector pos_temp = (p+i)->get_pos();
    vector vel_nfw;
    double r2 = (p+i)->get_radius();
    double V_c = sqrt(4.0*M_PI*rho_nfw*(r2*rs_nfw)*((log(1.0+r2/rs_nfw)/((r2/rs_nfw)*(r2/rs_nfw)))-(1.0/((r2/rs_nfw)*(1.0+r2/rs_nfw)))));
    vel_nfw[0] += V_c*(sqrt(pow(pos_temp[0],2.0)+pow(pos_temp[1],2.0))/r2)
                    *(-pos_temp[1]/(sqrt(pow(pos_temp[0],2.0)+pow(pos_temp[1],2.0))));
    vel_nfw[1] += V_c*(sqrt(pow(pos_temp[0],2.0)+pow(pos_temp[1],2.0))/r2)
                    *( pos_temp[0]/(sqrt(pow(pos_temp[0],2.0)+pow(pos_temp[1],2.0))));
    vel_nfw[2] += V_c*(pos_temp[2]/r2);
    (p+i)->vel_cycle = vel_nfw;
    (p+i)->inc_vel(vel_nfw);
  }
}
void NFW_acc(real_particle * pb, int nbody){
  real_particle * p = pb;
  for(int i=0; i<nbody; i++){
    vector pos_temp = (p+i)->get_pos();
    vector acc_nfw = 0.0;
    double r2 = (p+i)->get_radius();
    double nfw_pot = -(4.0*M_PI*rho_nfw*pow(rs_nfw,3.0))
                      *(log(1.0+r2/rs_nfw)-(r2/rs_nfw)/(1.0+r2/rs_nfw))/(r2*r2);
    acc_nfw[0] += nfw_pot*((sqrt(pow(pos_temp[0],2.0)+pow(pos_temp[1],2.0)))/r2)
                  *(pos_temp[0]/(sqrt(pow(pos_temp[0],2.0)+pow(pos_temp[1],2.0))));;
    acc_nfw[1] += nfw_pot*((sqrt(pow(pos_temp[0],2.0)+pow(pos_temp[1],2.0)))/r2)
                  *(pos_temp[1]/(sqrt(pow(pos_temp[0],2.0)+pow(pos_temp[1],2.0))));
    acc_nfw[2] += nfw_pot*(pos_temp[2]/r2);
    (p+i)->acc_external = acc_nfw;
    (p+i)->inc_acc(acc_nfw);
  }
}


// Plummer 球の重心計算
void Calculate_center_plummers(real_particle * pb, int nbody, double r_cut, vector cmdx){
  real_particle * p = pb;
  double cmmass = 0.0;
  for(int i=0; i<nbody; i++){
    double r2 = (p+i)->get_radius_cm();
    if(r2<=r_cut){
      cmdx   += (p+i)->get_mass()*(p+i)->get_dpos();
      cmmass += (p+i)->get_mass();
    }
  }
  cmdx /= cmmass;
  if(cmmass==0.0){cmdx = 0.0;}
  for(int i=0; i<nbody; i++){
    (p+i)->dpos -= cmdx;
  }
}
void rotating_center_of_mass(real_particle * pb, int nplummer, int nbody, int check,
                              double t, double shift_radius[], double shift_theta[],
                              double center_pos[3]){
  real_particle * p = pb;
  static double shift_radius_ave;
  static double shift_theta_ave;
  static double omega_center;
  vector dpos_temp = 0.0;
  vector cmdx = 0.0;
  vector pos_temp;
  if(check==0){
    for(int I=0; I<nplummer; I++){
      shift_theta_ave += shift_theta[I];
      shift_radius_ave += shift_radius[I];
    }
    shift_theta_ave = shift_theta_ave/nplummer;
    shift_radius_ave = shift_radius_ave/nplummer;
    omega_center = sqrt(4.0*M_PI*rho_nfw*(shift_radius_ave*rs_nfw)*
                        ((log(1.0+shift_radius_ave/rs_nfw)
                          /((shift_radius_ave/rs_nfw)*(shift_radius_ave/rs_nfw)))
                          -(1.0/((shift_radius_ave/rs_nfw)*(1.0+shift_radius_ave/rs_nfw)))))
                          /shift_radius_ave;
  }//end t=0
  center_pos[0] = shift_radius_ave*cos(omega_center*t+shift_theta_ave*M_PI/180.0);
  center_pos[1] = shift_radius_ave*sin(omega_center*t+shift_theta_ave*M_PI/180.0);
  center_pos[2] = 0.0;
  for(int i=0; i<(int)nplummer*nbody; i++){
    pos_temp = (p+i)->get_pos();
    for(int k=0; k<3; k++){
      dpos_temp[k] = pos_temp[k]-center_pos[k];
    }
    (p+i)->set_dpos(dpos_temp);
  }
  Calculate_center_plummers(p,(int)nplummer*nbody,5.0,cmdx);
  for(int i=0; i<3; i++){center_pos[i] -= cmdx[i];}
  for(int i=0; i<(int)nplummer*nbody; i++){
    double r2 = sqrt(center_pos[0]*center_pos[0]+center_pos[1]*center_pos[1]+center_pos[2]*center_pos[2]);
    double V_c = sqrt(4.0*M_PI*rho_nfw*(r2*rs_nfw)*((log(1.0+r2/rs_nfw)/((r2/rs_nfw)*(r2/rs_nfw)))-(1.0/((r2/rs_nfw)*(1.0+r2/rs_nfw)))));
    (p+i)->vel_cycle[0] = V_c*(sqrt(pow(center_pos[0],2.0)+pow(center_pos[1],2.0))/r2)
                           *(-center_pos[1]/(sqrt(pow(center_pos[0],2.0)+pow(center_pos[1],2.0))));
    (p+i)->vel_cycle[1] = V_c*(sqrt(pow(center_pos[0],2.0)+pow(center_pos[1],2.0))/r2)
                           *( center_pos[0]/(sqrt(pow(center_pos[0],2.0)+pow(center_pos[1],2.0))));
    (p+i)->vel_cycle[2] = V_c*(center_pos[2]/r2);
  }
}


void set_initial_distribution(real_system pb, int N, int n, int n_all, double eps2,
                              double shift_theta[], double shift_radius[], double center_pos[]){
  #ifdef COLD_COLLAPSE
    pb.create_coldcollapse(n_all, 0, 50.0, 0.5, eps2);
  #elif defined PLUMMER_ONLY
    pb.create_plummer_sphere(N, n, shift_theta, shift_radius);
  #elif defined INCLUDE_DATASET
    pb.INCLUDE_PARTICLE_DATASET(n_all);
    // pb.create_uniform_sphere(n_all, 0.0, 10.0);
  #else
    cerr << "create_some_plummer_sphere_in_DMH" << endl;
    pb.create_plummer_sphere(N, n, shift_theta, shift_radius);
    NFW_circular_velocity(pb.get_particle_pointer(), n_all);
    rotating_center_of_mass(pb.get_particle_pointer(), N, n, 0, 0.0, shift_radius, shift_theta, center_pos);
  #endif
}
