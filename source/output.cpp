// other usefull functions
// include packages
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <sys/stat.h> // ファイルの状態を得る（フォルダ作成時のchmod用）

using namespace std;
#define real double

// include my header files
#include "simulation_option.hpp"
#include "bhnode.hpp" // BHtree.hpp include "particle.hpp" and "myvector.hpp" and "constant.hpp"
#include "output.hpp"

void integrate(int nbody, real_system pb, bhnode * bn, bhparticle * bp,
                real eps2, real theta2, real dt, int check){
  real_particle * p_temp = pb.get_particle_pointer();
  for(int i=0; i<nbody; i++){
    (p_temp+i)->predict(dt);
  }
  calculate_gravity(nbody, pb, bn, bp, eps2, theta2);
  #if defined(COLD_COLLAPSE) || defined(PLUMMER_ONLY) || defined(INCLUDE_DATASET)
  #else
    NFW_acc(p_temp, nbody);
  #endif
  for(int i=0; i<nbody; i++){
    (p_temp+i)->correct(dt);
  }
}

int compare_double(const void *a, const void *b){
	int ret;
	const double *value1 = (const double *)a;
	const double *value2 = (const double *)b;
	if(*value1<*value2){
		ret =-1;
	}else if(*value1>*value2){
		ret=1;
	}
	return ret;
}
void compare_acc(int nbody, real_particle * pp, real Sacc[]){
  vector acc_error;
  real acc_direct2 = 0.0;
  for(int i=0; i<nbody; i++){
    acc_error   = (pp+i)->get_acc_gravity() - (pp+i)->get_acc_direct();
    acc_direct2 = sqrt((pp+i)->get_acc_direct() * (pp+i)->get_acc_direct());
    Sacc[i] = sqrt(acc_error*acc_error)/(acc_direct2);
  }
  qsort(Sacc, nbody, sizeof(double), compare_double);
}
void compare_time(int nbody, real_particle * pp, real Stime[]){
  for(int i=0; i<nbody; i++){
    Stime[i] = (pp+i)->get_time();
  }
  qsort(Stime, nbody, sizeof(double), compare_double);
}

void print_quantities(int nplummer, int n, int n_all, double eps2, double theta,
                      double dt, double tend, int iout, double M, double R,
                      double shift_theta[], double shift_radius[]){
  cerr << "A number of Plummer=" << nplummer
       << "\nOne plummer n=" << n << "\nTotal n=" << n_all
       << "\nsquared eps=" << eps2
       << "\nsquared theta=" <<theta <<"\ndt=" << dt << "\ttend=" << tend
       << "\tiout=" << iout << "\nM=" << M << "\tR=" << R
       << "\ndelta_interaction=" << delta_int << endl;
  #if !defined(INCLUDE_DATASET)
    for(int I=0; I<nplummer; I++){
      cerr << "shift_theta=" << shift_theta[I]
           << "\tshift_radius=" << shift_radius[I] << "\n";
    }
  #endif
}

void density_particle(real_particle * pp, int nbody, double SR[], double dens[]){
  real_particle * p_test = pp;
  #if defined(COLD_COLLAPSE) || defined(PLUMMER_ONLY) || defined(INCLUDE_DATASET)
  // Plane of Satelliteのシミュレーションでないとき
    for(int i=0; i<nbody; i++){
      // 単純に原点からの距離を半径とする
      SR[i] = (p_test+i)->get_radius();
    }
    qsort(SR, nbody, sizeof(double), compare_double);
    for(int i=0; i<nbody; i++){
      if(i%nbin==0 && i!=0){
        // bin内粒子数の中央値を用いて密度計算
        double RAVE1 = 0.5*(SR[(i-(int)(nbin*0.5))]+SR[(i-(int)(nbin*0.5-1))]);
        double RAVE2 = 0.5*(SR[(i+(int)(nbin*0.5))]+SR[(i+(int)(nbin*0.5+1))]);
        dens[i] = ((double)nbin/nbody)/fabs((4.0*M_PI*pow(RAVE2,3.0)/3.0)-(4.0*M_PI*pow(RAVE1,3.0)/3.0));
      }
    }
  #else
  // Plane of Satelliteのシミュレーションのとき
    for(int i=0; i<nbody; i++){
      // Plummer球同士の重心からの距離を半径とする
      SR[i] = (p_test+i)->get_radius_cm();
    }
    qsort(SR, nbody, sizeof(double), compare_double);
    for(int i=0; i<nbody; i++){
      if(i%nbin==0 && i!=0){
        double RAVE1 = 0.5*(SR[(i-(int)(nbin*0.5))]+SR[(i-(int)(nbin*0.5-1))]);
        double RAVE2 = 0.5*(SR[(i+(int)(nbin*0.5))]+SR[(i+(int)(nbin*0.5+1))]);
        dens[i] = ((double)nbin/nbody)/fabs((4.0*M_PI*pow(RAVE2,3.0)/3.0)-(4.0*M_PI*pow(RAVE1,3.0)/3.0));
      }
    }
  #endif
}

void input_number_of_particles(int *n_p, int *N_pl, int *n_all, real shift_radius[]){
  cerr << "Number of particles::";
  cin  >> *n_p;
  *N_pl = 1;
  *n_all = (*N_pl)*(*n_p);
  for(int I=0; I<(*N_pl); I++){
    cerr << "Shift Radius of No." << I+1 << " sphere[kpc]:";
    for( ;!(cin >> shift_radius[I]);){
      cout << "入力が間違っています。" << endl
      << "Shift Radius of No." << I+1 << "plummer";
    }
  }
}

void include_quantities(char *arg[], int *N, int *n, int *n_all, double *eps2, double *theta2,
                        double *dt, double *tend,int *iout, double *M, double *R){
  cerr << "Input quantities from filename " << arg[1] << endl;
  ifstream ifs(arg[1]);
  string str;
  while(getline(ifs, str)){
    sscanf(str.data(),
            "%d\t%d\t%lf\t%lf\t%lf\t%lf\t%d\t%lf\t%lf",
            N, n, eps2, theta2, dt, tend, iout, M, R);
  }
  *n_all = (*N) * (*n);
  *eps2 *= *eps2;
  *theta2 *= *theta2;
  return;
}
void print_simulation_mode(int *N, int *n, int *n_all, double shift_theta[3], double shift_radius[3]){
  #ifdef COLD_COLLAPSE
    cerr << "\n===Cold Collapse Simulation===" << endl;
    input_number_of_particles(n, N, n_all, shift_radius);
  #elif defined PLUMMER_ONLY
    cerr << "\n===Only Plummer Simulation===" << endl;
    input_number_of_particles(n, N, n_all, shift_radius);
  #elif defined INCLUDE_DATASET
    cerr << "\n===INCLUDE_DATASET Simulation===" << endl;
  #else
    cerr << "\n===Plane of Satellites Simulation===" << endl;
    shift_theta_hill(&N, shift_theta, shift_radius);
  #endif
}

stringstream make_directory(double shift_radius[]){
  //保存用フォルダの作成
  stringstream dirname;
  dirname << "data_" << (int)shift_radius[0] << "_" << (int)shift_radius[1];
  mkdir(dirname.str().c_str(), 0777);
  chmod(dirname.str().c_str(), 0777);
  return dirname;
}

void output_acc_error(const char * dirname, real_system pb, int n_all, double Sacc[], double Stime[]){
  // 出力ファイル名
  stringstream filename;
  filename << dirname << "/acc_error.dat";

  // 粒子データがほしいときはここを有効にしてアクセスする
  // real_particle * p_temp = pb.get_particle_pointer();

  ofstream output (filename.str().c_str());
  output << Sacc[n_all-1] << "\t" << Stime[n_all-1] << "\n\n\n";
  for(int i=0; i<n_all; i++){
    output << 100.0*i/n_all << "\t" << Sacc[i] << "\t" << Stime[i] << endl;
  }
  output.close();
}
void output_x_v(const char * dirname, real_system pb, int n_all, int id){
  // 出力ファイル名
  stringstream filename;
  filename << dirname << "/data_" << id << ".dat";

  // 粒子データアクセス用のポインタ
  real_particle * p_temp = pb.get_particle_pointer();

  ofstream output (filename.str().c_str());
  for(int i=0; i<n_all; i++){
    output << (p_temp+i)->get_pos()  << "\t"
           << (p_temp+i)->get_vel()  << "\t"
           << (p_temp+i)->get_dpos() << endl;
  }
  output .close();
}
void output_acc(const char * dirname, real_system pb, int n_all, int id){
  // 出力ファイル名
  stringstream filename;
  filename << dirname << "/acc_data_" << id << ".dat";

  // 粒子データアクセス用のポインタ
  real_particle * p_temp = pb.get_particle_pointer();

  ofstream output (filename.str().c_str());
  for(int i=0; i<n_all; i++){
    output << (p_temp+i)->get_acc_gravity()  << "\t"
           << (p_temp+i)->get_acc_direct()   << "\t"
           << (p_temp+i)->get_acc_external() << endl;
  }
  output .close();
}
void output_energy(const char * dirname, real_system pb, int n_all, double t, double center_pos[3]){
  // 出力ファイル名
  stringstream filename;
  filename << dirname << "/energy_data.dat";

  double PE=0.0, KE=0.0, DE=0.0, MV=0.0;
  KE = pb.kinetic_energy(n_all, 5.0);     // 運動エネルギーの計算
  PE = pb.potential_energy(n_all, 5.0);   // ポテンシャルエネルギーの計算
  MV = pb.momentum(n_all, 5.0);           // 運動量の計算

  ofstream output;
  if(t==0.0){
    DE = pb.energy_error(n_all, 0, KE, PE); // 初期エネルギーとの誤差計算
    output .open(filename.str().c_str());
  }else{
    DE = pb.energy_error(n_all, 1, KE, PE); // 初期エネルギーとの誤差計算 2回目以降
    // 2回目以降の出力はファイルに上書きしないといけない
    output .open(filename.str().c_str(), ios_base::app);
  }
  output << scientific << t << "\t"
         << scientific << KE << "\t"
         << scientific << PE << "\t"
         << scientific << 2.0*KE/fabs(PE) << "\t"
         << scientific << DE << "\t"
         << scientific << MV << "\t"
         << scientific << center_pos[0] << "\t"
         << scientific << center_pos[1] << "\t"
         << scientific << center_pos[2] << "\t" << endl;
  output .close();
}
void output_dens(const char * dirname, real_system pb, int n_all, int id, double SR[], double dens[]){
  // 出力ファイル名
  stringstream filename;
  filename << dirname << "/density_" << id << ".dat";

  // 粒子データアクセス用のポインタ
  // real_particle * p_temp = pb.get_particle_pointer();

  ofstream output(filename.str().c_str());
  for(int i=0; i<n_all; i++){
    if(i%nbin==0 && i!=0){
      output << scientific << SR[i] << "\t"
             << scientific << dens[i] << endl;
    }
  }
  output .close();
}

void gnuplot_plot(){
  FILE * fpg;
  #ifdef COMP_ACC
    fpg = popen("gnuplot", "w");
    // gnuplot で morton key をカラーマップで画像出力
    fputs("load 'gnuplot/gnu_acc_error_percentile_snap.gp'\n", fpg);
    fflush(fpg);
    pclose(fpg);
  #endif
  #ifdef MORTON_OUTPUT
    fpg = popen("gnuplot", "w");
    // gnuplot で morton key をカラーマップで画像出力
    fputs("load './gnuplot/gnu_morton_particle.gp'\n", fpg);
    fflush(fpg);
    pclose(fpg);
  #endif
  #ifdef ENERGY_OUTPUT
    fpg = popen("gnuplot", "w");
    fputs("load './gnuplot/gnu_virial_ratio.gp'\n", fpg);
    fputs("load './gnuplot/gnu_energy_error.gp'\n", fpg);
    fflush(fpg);
    pclose(fpg);
  #endif
  #ifdef GNUPLOT_POS
    fpg = popen("gnuplot", "w");
    // fputs("load 'gnu_position2d_gc_gif.gp'\n", fpg);
    fputs("load './gnuplot/gnu_position2d_gif.gp'\n", fpg);
    fflush(fpg);
    pclose(fpg);
  #endif
  #ifdef DENSITY_OUTPUT
    fpg = popen("gnuplot", "w");
    fputs("load './gnuplot/gnu_density_snap.gp'\n", fpg);
    fflush(fpg);
    pclose(fpg);
  #endif
  #ifdef TREE_WALK_OUTPUT
    fpg = popen("gnuplot", "w");
    //計算時間をpercentile-timeの図を出力
    fputs("load './gnuplot/gnu_time_percentile_snap.gp'\n", fpg);
    //半径-計算時間の図を出力
    fputs("load './gnuplot/gnu_time_position.gp'\n", fpg);
    //半径-計算時間の図を密度分布とともに出力
    fputs("load './gnuplot/gnu_time_position_another.gp'\n", fpg);
    // //2次元座標と計算時間の図を出力
    // fputs("load 'gnu_time_position_2d.gp'\n", fpg);
    //セル内の粒子数-半径の図を出力
    fputs("load './gnuplot/gnu_nbody_radius.gp'\n",fpg);
    fflush(fpg);
    pclose(fpg);
  #endif
}
