#include <stdlib.h>   // C言語の標準ライブラリ
#include <stdio.h>    // ファイルアクセス関数・書式付き入出力関数
#include <math.h>     // C言語の数学関数ライブラリ
#include <iostream>   // 基本的なstream入出力機能
#include <sstream>    // stringstream利用機能
#include <fstream>    // filestream ファイルの入出力関数
#include <string>     // 文字列クラス

using namespace std;  // std::入力不要
#define real double   // real型はdouble
// #define real float // real型はfloat

#include "simulation_option.hpp"  // シミュレーションの内容を変えるための define 群
#include "bhnode.hpp"             // "particle.hpp" "myvector.hpp" "constant.hpp"
#include "debug.hpp"              // debug用の関数

int n_intaract = 0.0;
void bhnode::output_data_in_node(int tstep, int output_level){
  // セル情報を集める関数
  static int check = 0;
  char nodename[50];
  if(check==0){
    sprintf(nodename,"data_node_%d.dat",tstep);
    ofstream nodefile(nodename);
    nodefile .close();
    check = 1;
    cerr << nodename << " is cleared." << endl;
  }
  real_particle * p_out = NULL;
  sprintf(nodename,"data_node_%d.dat",tstep);
  ofstream nodefile(nodename, ios_base::app);
  for(int i=0; i<nparticle; i++){
    p_out = (bpfirst+i)->get_rp();
    nodefile << (p_out)->get_pos()  << "\t" << (p_out)->get_vel() << "\t"
             << (p_out)->get_dpos() << "\t" << (p_out)->get_index() << "\t"
             << i << endl;
  }
  nodefile .close();
}

void tree_walk_counter(real_particle * rp, bhparticle * bp, int nbody, vector cmpos, real length, int level){
  // ある粒子のTreeWalk情報を集める関数
  real_particle * p_test = rp;
  bhparticle * bp_out = bp;
  ofstream outputtree("./data_0_0/tree_walk.dat",ios_base::app);
  outputtree << cmpos << " " << nbody << " " << length << endl;
  // // もろもろの木構造を出力、gnuplotには利用できない構造
  //   for(int i=0; i<level; i++) outputtree << "  ";
  //   outputtree << "level:" << level
  //              << " cm:"   << cmpos
  //              << " nparticle:" << nbody
  //              << " index:";
  //   for(int i=0; i<nbody; i++){
  //     outputtree << ((bp_out+i)->get_rp())->get_index() << " ";
  //   }
  // outputtree << endl;
  outputtree .close();
}

void bhnode::accumulate_tree_walk(real_particle * rp, real theta2, int n_th){
  real_particle * p_test = rp;
  vector dx = cmpos-p_test->get_pos();
  real r2 = dx*dx;
  bhparticle * bp = NULL;
  if((r2*theta2 > l*l) || (nparticle == 1)){
    bp = bpfirst;
    n_intaract += get_nparticle();
    tree_walk_counter(p_test,bp,get_nparticle(),get_cmpos(),get_length(),get_level());
  }else{
    for(int i=0; i<8; ++i){
      if(child[i] != NULL){
        //近いけどまだ下にセルがあるとき
        child[i]->accumulate_tree_walk(rp,theta2,n_th);
      }
    }
  }
}
void search_tree_walk(int nbody, real_system pb, bhnode * bn, bhparticle * bp, real eps2, real theta2){
  // この関数でTree Walkを探査することができる
  // clear_tree_counters();
  bhparticle * bp_test = bn->get_bpfirst();
  real_particle * p_test = NULL;
  real t_crit = 60.0;
  int n_th = 0, nt_max = 0;
  float t_max = 0.0, t_temp = 0.0;
  for(int i=0; i<nbody; i++){
    //ほしい粒子の条件を入れてtree walkをどうやってるかを出力
    real r_temp = (pb.get_particle_pointer()+i)->get_radius();
    if((r_temp>0.99) && (r_temp<1.02))
      t_temp = (pb.get_particle_pointer()+i)->get_time();
    if(t_temp>t_max){
      nt_max = i;
      t_max = t_temp;
    }
    // この場合一番時間のかかっている粒子をひっぱってこれる
  }
  // 条件にあう粒子をp_testに確保
  p_test = pb.get_particle_pointer()+nt_max;

  // 出力準備
  ofstream outputtree("./data_0_0/tree_walk.dat");
  outputtree << p_test->get_pos() << " " << p_test->get_radius()
             << "\t" << p_test->get_time() << "\n\n\n";
  outputtree .close();

  // Tree Walkの探査開始
  bn->accumulate_tree_walk(p_test, theta2, n_th);
  // for(int i=0; i<nbody; i++){
  //   p_test = pb.get_particle_pointer()+i;
  //   if(p_test->get_time()>t_crit){
  //     ofstream outputtree("./data_0_0/tree_walk.dat",ios_base::app);
  //     outputtree << p_test->get_pos() << "\n\n" << endl;
  //       // outputtree << "I'm No." << p_test->get_index()
  //       //            << " time:"  << p_test->get_time() << endl;
  //     outputtree .close();
  //     bn->accumulate_tree_walk(p_test, theta2, n_th);
  //   }
  // }

  cerr << "All intaraction number: " << n_intaract << endl;
  // print_tree_counters();
}
