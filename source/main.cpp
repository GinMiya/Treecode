// N-body simulation [Main stream code]
// 重力計算法：Tree / Multipole
// Morton keyによる粒子データのsort有り
// 相互作用判定は rθ > l もしくは、 Warren & Salmon(1993) を用いることが可能

#include <stdlib.h>   // C言語の標準ライブラリ
#include <stdio.h>    // ファイルアクセス関数・書式付き入出力関数
#include <math.h>     // C言語の数学関数ライブラリ
#include <iostream>   // 基本的なstream入出力機能
#include <sstream>    // stringstream利用機能
#include <fstream>    // filestream ファイルの入出力関数
#include <string>     // 文字列クラス
#include <limits>     // min,maxなどの算術関数
#include <sys/stat.h> // ファイルの状態を得る（フォルダ作成時のchmod用）
#include <omp.h>      // OpenMPを使うためのライブラリ

using namespace std;  // std::入力不要

#define real double   // real型はdouble
// #define real float // real型はfloat

// 出力用マクロ関数
#define PR(x)  cerr << #x << " = " << x << " "    // 引数の値を返すだけ
#define PRC(x) cerr << #x << " = " << x << ",  "  // 引数の値を返して , で区切る
#define PRL(x) cerr << #x << " = " << x << "\n"   // 引数の値を返して改行する

// 使うヘッダーファイル
#include "simulation_option.hpp"  // シミュレーションの内容を変えるための define 群
#include "bhnode.hpp"             // "particle.hpp" "myvector.hpp" "constant.hpp"
#include "debug.hpp"              // debug用関数
#include "morton.hpp"             // morton keyによるsort用関数
#include "make_distribution.hpp"  // 初期粒子分布作成用関数
#include "calc_grav.hpp"          // 相互作用計算用関数
#include "multipole.hpp"          // 多重極展開の関数
#include "output.hpp"             // ファイル入出力用関数


// 重力計算パートはここで定義するとなぜか早くなる現象が見られる 原因不明
static real total_interactions[16];  // 総相互作用計算回数
static int nisum[16];                // 重力計算済み粒子数
static real nisleaf[16];             // Treeの末端までいたった回数
void clear_tree_counters(){
  for(int i=0; i<16; i++){
    total_interactions[i] = 0;
    nisum[i] = 0;
    nisleaf[i] = 0;
  }
}
void print_tree_counters(){
  real avg[16];
  real average[4];
  #ifdef OPENMP
    for(int i=0; i<n_threads; i++){
      avg[i] = total_interactions[i]/nisum[i];
      PRC(nisum[i]); PRC(total_interactions[i]); PRC(nisleaf[i]);
      PRL(avg[i]);
      cout << "n_interact_average:" << i << " " << avg[i] << endl;
      average[0] += nisum[i];
      average[1] += total_interactions[i];
      average[2] += nisleaf[i];
    }
    // average[0] /= n_threads;
    // average[1] /= n_threads;
    // average[2] /= n_threads;
    average[3] = average[1]/average[0];
    cout << "n_total = " << average[0] << " total_interactions = " << average[1]
    << " nisleaf_total = " << average[2] << endl;
    cout << "n_interact_average:" << average[3] << endl;
  #else
    avg[0] = total_interactions[0]/nisum[0];
    PRC(nisum[0]); PRC(total_interactions[0]); PRC(nisleaf[0]);
    PRL(avg[0]);
    cout << "n_interact_average:" << " " << avg[0] << endl;
  #endif
}
void accumulate_force_from_point(vector dx, real r2, real eps2,
                        				 vector & acc, real & phi,real jmass){
  // 質点からうける重力 もともとは calc_grav.cpp で定義されている関数
  double r2inv = 1.0/(r2+eps2);
  double rinv  = sqrt(r2inv);
  double r3inv = r2inv*rinv;
  phi -= jmass*rinv;
  acc += jmass*r3inv*dx;
}
void bhnode::accumulate_force_from_tree(vector & ipos, real eps2, real theta2,
					                              vector & acc, real & phi, int n_th){
	// Treeを使った重力の計算
  vector dx = cmpos - ipos;
  real r2 = dx * dx;
	#ifdef MULTIPOLE
    // r_crit = l*l;
		// if((r2*theta2 > r_crit)){
    if((r2 > r_crit)){
    // node and position is well separated
    // or the number of particles in node is Only 1
	    accumulate_force_multipole(dx, r2, b_max, acc, phi);
    	total_interactions[n_th] += 1;
		}else if((nparticle == 1)){
	    accumulate_force_from_point(dx,r2,eps2,acc,phi,cmmass);
	    total_interactions[n_th] += 1;
		}else{
	    int i;
	    if(isleaf){
	      //セルの末端まで来てしまったときは中の粒子をすべてdirectで解く
	      bhparticle * bp = bpfirst;
	      for(i=0; i<nparticle; i++){
	        nisleaf[n_th] ++;
	        dx = ((bp+i)->get_rp())->get_pos()-ipos;
	        r2 = dx*dx;
	        accumulate_force_from_point(dx,r2,eps2,acc,phi,(bp+i)->get_rp()->get_mass());
	        total_interactions[n_th] += 1;
	      }
	    }else{
	      for(i=0; i<8; ++i){
	        if(child[i] != NULL){
	          //近いけどまだ下にセルがあるとき
	          child[i]->accumulate_force_from_tree(ipos,eps2,theta2,acc,phi,n_th);
	        }
	      }
	    }
		}
  #else
		if((r2*theta2 > l*l) || (nparticle == 1)){
    // node and position is well separated
    // or the number of particles in node is Only 1
    accumulate_force_from_point(dx,r2,eps2,acc,phi,cmmass);
    total_interactions[n_th] += 1;
  } else {
    int i;
    if(isleaf){
      //セルの末端まで来てしまったときは中の粒子をすべてdirectで解く
      bhparticle * bp = bpfirst;
      for(i=0; i<nparticle; i++){
        nisleaf[n_th] ++;
        dx = ((bp+i)->get_rp())->get_pos()-ipos;
        r2 = dx*dx;
        accumulate_force_from_point(dx,r2,eps2,acc,phi,(bp+i)->get_rp()->get_mass());
        total_interactions[n_th] += 1;
      }
    }else{
      for(i=0; i<8; ++i){
        if(child[i] != NULL){
          //近いけどまだ下にセルがあるとき
          child[i]->accumulate_force_from_tree(ipos,eps2,theta2,acc,phi,n_th);
        }
      }
    }
  }
	#endif
}
void calculate_gravity_using_tree(real_particle * rp, bhnode * bn, real eps2, real theta2, int n_th){
	// 値渡しの関数
  rp->acc_gravity = 0;
  rp->phi_gravity = rp->mass/sqrt(eps2);
  bn->accumulate_force_from_tree(rp->pos,eps2,theta2,rp->acc_gravity, rp->phi_gravity, n_th);
  nisum[n_th] += 1;
}
void calculate_gravity(int nbody, real_system pb, bhnode * bn, bhparticle * bp, real eps2, real theta2){
											  double start_key,      end_key;
											  double start_tree,     end_tree;
											  double start_particle, end_particle;
											  double start_test,     end_test;
	// morton keyを設定　Treeの作成　重力計算まで行う関数

  // Morton key計算 Part
  int nkey = 0;
  start_key = omp_get_wtime();
    real rsize = initialize_key(nbody, pb.get_particle_pointer(), nkey, bp);
    PRL(rsize);
  end_key = omp_get_wtime();
  cerr << "initialize_key time:" << (end_key-start_key) << "[sec]" << endl;

  // Tree作成のための準備
  for(int j=0; j<nbody; j++) (bn+j)->clear();
  //ここの初期サイズrsizeを大きくすると精度が増す
  bn->assign_root(vector(0.0), 2*rsize, bp, nbody);
  bhnode * btmp = bn+1;
  int heap_remainder = nbody*2;
  BHlong key = 0;

  // Tree作成 Part
  start_tree = omp_get_wtime();
    bn->create_tree_recursive(btmp, heap_remainder, key, default_key_length, 1);
  end_tree = omp_get_wtime();
  cerr << "create tree constructure:" << (end_tree-start_tree) << "[sec]" << endl;

  // 重心情報計算 Part
  start_tree = omp_get_wtime();
    bn->set_cm_quantities();
    // bn->sanity_check();
    // bn->dump();
  end_tree = omp_get_wtime();
  cerr << "calc. cm quantities:" << (end_tree-start_tree) << "[sec]" << endl;

  // 重力計算のための準備
  clear_acc_and_phi(pb.get_particle_pointer(),nbody);
  // PRL(bn->sanity_check());
  // bn->dump();
  clear_tree_counters();
  bhparticle * bp_test = bn->get_bpfirst();
  real_particle * p_test = NULL;

  // 重力計算 Part
  int n_th = 0;
  start_test = omp_get_wtime();
    #ifdef OPENMP
      #pragma omp parallel private(n_th)
      {
        n_th = omp_get_thread_num();
        #pragma omp for
    #endif
        for(int i=0; i<nbody; i++){
          start_particle = omp_get_wtime();
            p_test = (bp_test+i)->get_rp();
            // p_test = pb.get_particle_pointer()+i;
            calculate_gravity_using_tree(p_test, bn, eps2, theta2, n_th);
          end_particle = omp_get_wtime();
          p_test->set_time(end_particle-start_particle);
        }
    #ifdef OPENMP
      }
    #endif
  end_test = omp_get_wtime();
  cerr << "time of calc. gravity:" << (end_test-start_test) << "[sec.]" << endl;

  // 結果の表示
  print_tree_counters();
}

vector th_ph_grad(double r, double theta, double phi, complex alpha[], complex Y_ij[]){
  // θ微分
  vector acc_sph = 0.0;
  double br, coeff, coeff1, coeff2;
  double common_term, legendre_1, legendre_2;
  int mm;
  for(int l=0; l<=p_max; l++){
    // br = 4.0*M_PI*G0*pow(1.0/r, l+2)/(2.0*l+1.0);
    br = pow(1.0/r, l+2);
    for(int m=-l, j=0; m<=l; m++, j++){
      if(m==0) {coeff = 1.0;} else{coeff=2.0;}
      mm = fabs(m);
      common_term = pow(-1.0, (mm+m)/2) * ScalingFactor(l, mm);
      legendre_1 = common_term * Legendre_Cos(l, fabs(m), cos(theta), theta);
      legendre_2 = common_term * Legendre_Cos(l, fabs(m)+1, cos(theta), theta);
      coeff  = ((cos(theta)*fabs(m)*legendre_1)/sin(theta)-legendre_2);
      // 実部虚部成分
      coeff1 = coeff*cos(m*phi);
      coeff2 = coeff*sin(m*phi);
      acc_sph[1] -= br * (alpha[l*l+j].get_real()*coeff1
                        + alpha[l*l+j].get_imag()*coeff2);
      acc_sph[2] -= br * (m/sin(theta)) * (-alpha[l*l+j].get_real()*Y_ij[l*l+j].get_imag()
                        + alpha[l*l+j].get_imag()*Y_ij[l*l+j].get_real());
    }
  }
  return acc_sph;
}
void get_rthph_acc(vector add_temp, double r_i, double theta, double phi, complex alpha[], complex Y_ij[]){
  // 極座標での加速度を得る関数
  add_temp += th_ph_grad(r_i, theta, phi, alpha, Y_ij);
}
vector get_xyz_acc(vector acc_sph, double theta, double phi){
  // カーテシアンでの加速度を得る関数
  myvector acc_xyz = 0.0;
  acc_xyz[0] = sin(theta)*cos(phi)*acc_sph[0]+cos(theta)*cos(phi)*acc_sph[1]-sin(phi)*acc_sph[2];
  acc_xyz[1] = sin(theta)*sin(phi)*acc_sph[0]+cos(theta)*sin(phi)*acc_sph[1]+cos(phi)*acc_sph[2];
  acc_xyz[2] = cos(theta)*acc_sph[0]-sin(theta)*acc_sph[1];
  return acc_xyz;
}

void bhnode::accumulate_force_multipole(vector dx, real r2, real b, vector & acc, real & phi){
  // accumulate_force_from_tree の multipole バージョン
  double r = sqrt(r2);
  double th = acos(dx[2]/r);
  double ph = atan2(dx[1],dx[0]);
  complex Y_ij[(2*p_max+1)+(p_max*p_max)];
  double add_phi = 0.0;
  vector add_acc = 0.0;
  Y_term_ij(Y_ij, th, ph);
  double br;
  for(int l=0; l<=p_max; l++){
    // br = 4.0*M_PI*G0*pow(1.0/r, l+1)/(2.0*l+1.0);
    br = pow(1.0/r, l+1);
    for(int m=-l, j=0; m<=l; m++, j++){
      add_phi    -= br * (alpha[l*l+j].get_real()*Y_ij[l*l+j].get_real()
                       +  alpha[l*l+j].get_imag()*Y_ij[l*l+j].get_imag());
      add_acc[0] -= -(l+1)*br/r * (alpha[l*l+j].get_real()*Y_ij[l*l+j].get_real()
                    + alpha[l*l+j].get_imag()*Y_ij[l*l+j].get_imag());
    }
  }
  phi += add_phi;
  get_rthph_acc(add_acc, r, th, ph, alpha, Y_ij);
  acc += get_xyz_acc(add_acc, th, ph);
}

// 引数に物理量をまとめたファイル "quantities.dat" を入力する
int main(int argc, char *arg[]){
  // 'quantities.dat'からの物理量の読み込み
  int N, n, n_all, iout, istep;
  real eps2, theta2, dt, tend, t=0.0, M, R;
  include_quantities(arg, &N, &n, &n_all, &eps2, &theta2, &dt, &tend, &iout, &M, &R);

  // Hill半径の計算とoffset角度の計算
  real shift_theta[3], shift_radius[3], center_pos[3];
  print_simulation_mode(&N, &n, &n_all, shift_theta, shift_radius);

  // 物理量の確認用出力
  print_quantities(N, n, n_all, eps2, theta2, dt, tend, iout, M, R, shift_theta, shift_radius);

  // 保存用ディレクトリの作成
  stringstream dirname;
  dirname = make_directory(shift_radius);

  // nbody_system の作成→シミュレーションの土台
  static nbody_system pb;
  // Tree法のノードアクセス用ポインタ
  bhnode * bn;

  // nbody_system クラスにおいて　初期粒子分布の作成
  #ifdef COLD_COLLAPSE
    pb.create_coldcollapse(n_all, 0, 50.0, 0.5, eps2);
  #elif defined PLUMMER_ONLY
    pb.create_plummer_sphere(N, n, shift_theta, shift_radius);
  #elif defined INCLUDE_DATASET
    pb.INCLUDE_PARTICLE_DATASET(n_all);
  #else
    cerr << "create_some_plummer_sphere_in_DMH" << endl;
    pb.create_plummer_sphere(N, n, shift_theta, shift_radius);
    NFW_circular_velocity(pb.get_particle_pointer(), n_all);
    rotating_center_of_mass(pb.get_particle_pointer(), N, n, 0, 0.0, shift_radius, shift_theta, center_pos);
  #endif

  // Tree作成の準備　bhparticle と Treeのnode を確保しておく
  bhparticle * bp = NULL;
  bp = new bhparticle[n_all];
  // ここなぞの処理
  real_particle * p_temp = pb.get_particle_pointer();
  for(int i=0; i<n_all; i++) {(bp+i)->set_rp(p_temp+i);}
  // ノード数（セル数）はだいたい粒子数の２倍程度用意
  int nnode = 2 * n_all + 100;
  bn = new bhnode[nnode];

  // Tree法での加速度計算　Tree make もここで行う
  calculate_gravity(n_all, pb, bn, bp, eps2, theta2);

  #ifdef TREE_WALK_OUTPUT
    // Tree walk がうまくいっているのかを確認するためのデバッグ関数
    search_tree_walk(n_all,pb,bn,bp,eps2,theta2);
  #endif

  #ifdef COMP_ACC
    // Tree法 と Direct法 の加速度相対誤差を評価する場合に有効
    double start, end;
    start = omp_get_wtime();
    // ダイレクトで相互作用計算 / 計算時間の測定
      calculate_gravity_direct(pb.get_particle_pointer(), n_all, eps2);
    end = omp_get_wtime();
    cerr << "One timestep in Direct:" << (end-start) << "[sec.]" << endl;

    real Sacc[n_all];   // 加速度の大きい順に並べるための配列
    real Stime[n_all];  // 計算時間が長い順に並べるための配列
    compare_acc(n_all, pb.get_particle_pointer(), Sacc);  // 加速度のクイックソート
    compare_time(n_all,pb.get_particle_pointer(), Stime); // 計算時間のクイックソート
    output_acc_error(dirname.str().c_str(), pb, n_all, Sacc, Stime); // 加速度errorの出力
  #endif

  #if !defined(COLD_COLLAPSE) && !defined(PLUMMER_ONLY) && !defined(INCLUDE_DATASET)
    // NFW profile の Dark matter halo からうける重力を加算する必要がある
    NFW_acc(pb.get_particle_pointer(),n_all);
  #endif

  // 位置と速度の出力
  output_x_v(dirname.str().c_str(), pb, n_all, 0);

  #ifdef MORTON_OUTPUT
    // Morton順序に沿って粒子が並べられているかの確認用デバッグ関数
    bn->output_data_in_node(0, 3); // 第１引数は出力時間　第２引数は出力するレイヤーレベル
  #endif

  #ifdef ACC_OUTPUT
    // 粒子の加速度を出力　第４引数はファイル番号
    output_acc(dirname.str().c_str(), pb, n_all, 0);
  #endif

  #ifdef ENERGY_OUTPUT
    // 系のエネルギーの出力　第４引数は時間　第５引数はplummer球の重心
    output_energy(dirname.str().c_str(), pb, n_all, 0.0, center_pos);
  #endif

  #ifdef DENSITY_OUTPUT
    real dens[n_all], SR[n_all]; // 密度と半径　小さい順に並べる
    // 密度の計算　半径をある bin サイズで区切ってその中での中央値における体積比を用いる
    density_particle(pb.get_particle_pointer(), n_all, SR, dens);
    // 系の密度の出力　第４引数はファイル番号
    output_dens(dirname.str().c_str(), pb, n_all, 0, SR, dens);
  #endif

  #ifdef INTEGRATE
    // iteration にかかる時間の計測
    double start_roop, end_roop;
    start_roop=omp_get_wtime();
      cerr << "starttime[μsec.]:" << start_roop << endl;

      // 時間積分
      while(t<=tend){
        // leap-frog 法の実行と NFW 加速度の加算
        #if !defined(COLD_COLLAPSE) && !defined(PLUMMER_ONLY) && !defined(INCLUDE_DATASET)
          // NFW profile の Dark matter halo からうける重力を加算する必要がある
          integrate(n_all, pb, bn, bp, eps2, theta2, dt, 0);
        #else
          //NFWなしのPlummer球が１個の場合は第８引数1にしておこう
          integrate(n_all, pb, bn, bp, eps2, theta2, dt, 1);
        #endif
        istep++;
        cerr << "Time=" << t << endl;

        if((istep%iout==0)){
          // ファイル出力する場合の if 文
          cerr << "OUTPUT FILE" << endl;

          #if !defined(COLD_COLLAPSE) && !defined(PLUMMER_ONLY) && !defined(INCLUDE_DATASET)
            // 回転するいくつかの衛生銀河の重心をもとめる必要がある
            rotating_center_of_mass(pb.get_particle_pointer(), N, n, 1, t, shift_radius, shift_theta, center_pos);
          #endif

          // 粒子の位置と速度の出力
          output_x_v(dirname.str().c_str(), pb, n_all, istep);

          #ifdef ACC_OUTPUT
            // 粒子の加速度の出力
            output_acc(dirname.str().c_str(), pb, n_all, istep);
          #endif

          #ifdef ENERGY_OUTPUT
            // 系のエネルギーの出力
            output_energy(dirname.str().c_str(), pb, n_all, t, center_pos);
          #endif

          #ifdef DENSITY_OUTPUT
            // 系の密度の出力
            density_particle(pb.get_particle_pointer(), n_all, SR, dens);
            output_dens(dirname.str().c_str(), pb, n_all, istep, SR, dens);
          #endif
        }

        t += dt;

      }
    // 終了時間の計測とかかった時間の計算
    end_roop=omp_get_wtime();
    cerr << "uptime:" << (end_roop-start_roop)/1000000.0 << "[sec.]" << endl;
  #endif

  #ifdef GNUPLOT
    // gnuplot による画像 / gif の出力を行う
    gnuplot_plot();
  #endif

  return 0;
}
