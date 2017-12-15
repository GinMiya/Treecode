// bhnode.hpp　で定義した関数
// bhnode クラスの メンバ関数 が主な構成員

// ライブラリ読み込み
#include <stdlib.h>   // C言語の標準ライブラリ
#include <stdio.h>    // ファイルアクセス関数・書式付き入出力関数
#include <math.h>     // C言語の数学関数ライブラリ
#include <iostream>   // 基本的なstream入出力機能
#include <fstream>    // filestream ファイルの入出力関数
#include <omp.h>      // OpenMPを使うためのライブラリ

using namespace std;  // std:: 入力不要
#define real double   // double と float を使い分ける用

// マクロ関数定義
#define PR(x)  cerr << #x << " = " << x << " "    // 引数の値を返すだけ
#define PRC(x) cerr << #x << " = " << x << ",  "  // 引数の値を返して , で区切る
#define PRL(x) cerr << #x << " = " << x << "\n"   // 引数の値を返して改行する

#include "simulation_option.hpp"  // シミュレーションの内容を変えるための define 群
#include "bhnode.hpp"             // "particle.hpp" "myvector.hpp" "constant.hpp"
#include "debug.hpp"              // debug用 関数
#include "morton.hpp"             // morton key による sort用 関数
#include "calc_grav.hpp"          // 相互作用計算用 関数


// Tree construction 作成関数　の　メイン部分
inline BHlong construct_key(const vector & pos, real rscale, int ioffset, int keybits){
  // 第一引数posに入れた３次元座標を整数値で番地を割り当てる関数
  BHlong ix[3];
  vector p = pos;
  for(int i=0; i<3; i++){
    //セルの1辺の長さがわかっているので ioffset 分のセルが使われていると思って計算する
  	ix[i] = (BHlong) (p[i]*rscale+ioffset);
  }
  // morton keyを返す
  return getMorton3D(ix[0],ix[1],ix[2]);
}
void bhnode::assign_root(vector root_pos, real length, bhparticle * bp, int np){
  //BHtreeを作るにあたって最初の bhポインタ をつなぐ関数
  pos = root_pos;
  l = length;
  bpfirst = bp;
  nparticle = np;
}
void bhnode::create_tree_recursive(bhnode * & heap_top, int & heap_remainder,
                        				   BHlong current_key,  int current_level,
                        				   int n_critical){
  // 実際にTree構造を作る関数
  // cerr << "create tree called ";
  // PRC(nparticle); PRL(n_critical);PRL(heap_remainder);
  if (heap_remainder <= 0){
    // ヒープ領域にストックがなくなったらエラー吐いて終了
  	cerr << "create_tree: no more free node... exit\n";
  	exit(1);
  }

  // n_criticalよりもすくなくなった、もしくは一番上の階層まで来たら終了
  if (nparticle <= n_critical) return;
  if (current_level == 0) return;
  // dump();

  // 今の階層におけるmorton keyの大きさ
  BHlong keyscale = ((BHlong) 1)<<((current_level-1)*3);

  // 最初の bhparticle ポインタの入手
  bhparticle * bptmp = bpfirst;
  int npremain = nparticle;
  // 子ノードの作成
  for(int i=0; i<8; i++) child[i] = NULL;
  isleaf = 1;

  for(int i=0; i<8; i++){
    // 子ノードでの morton key の特徴的な大きさ
  	BHlong new_key = current_key + keyscale * i;
    // 子ノードの中心位置の計算
  	vector new_pos = pos + vector(((i&4)*0.5-1)*l/4,
                  				        ((i&2)    -1)*l/4,
                  				        ((i&1)*2  -1)*l/4);
  	if(bptmp->get_key() - new_key < keyscale){
      // current bptmp is actually in the current subnode
      // search for the end location
      int p0 = 0;
      int p1 = npremain-1;
      if ((bptmp+p1)->get_key() - new_key >= keyscale){
    		while (p1-p0 > 1){
    	    int pnew = (p0+p1)/2;
    	    if ((bptmp+pnew)->get_key() - new_key < keyscale){
        		p0 = pnew;
    	    }else{
        		p1 = pnew;
    	    }
      	}
    		p1 = p0;
      }
      p1 ++;
      isleaf = 0;
      child[i] = heap_top;
      heap_top ++;
      heap_remainder --;
      child[i]->bpfirst = bptmp;
      child[i]->pos = new_pos;
      child[i]->l = l*0.5;
      child[i]->nparticle = p1;
      child[i]->isleaf = 1;
      child[i]->level = default_key_length-current_level+1;
      child[i]->create_tree_recursive(heap_top, heap_remainder,
                          				    new_key, current_level-1, n_critical);
      bptmp += p1;
      npremain -= p1;
      if(npremain<=0) return;
  	}
  }
  // dump();
}


// Tree construction デバッグ用関数
void spc(int indent){
  // 引数の数だけスペースを空ける関数
  for(int i=0;i<indent;i++)cerr << " ";
}
void bhnode::dump(int indent){
  // 木構造がちゃんとつくられているかを確認するための関数
  int i;
  // spc(indent); cerr << "node pos " << pos ;
  cerr << endl;
  spc(indent);
  cerr << "node cm  " << cmpos << " m " << cmmass << " bmax " << b_max ;
  if (isleaf){
  	cerr << " IS LEAF" ;PRL(nparticle);
  	bhparticle * bp = bpfirst;
  	for(i=0; i<nparticle; i++){
  	    for(int j=0; j<indent+2; j++) cerr << " ";
  	    real_particle * p = (bp+i)->get_rp();
  	    PR(p->get_index());
        // PR(to_binString((bp+i)->get_key()));
        PRL(p->get_pos());
  	}
  }else{
  	cerr << " IS _not_ LEAF "; PRL(nparticle);
  	for(i=0; i<8; i++){
	    if (child[i] != NULL){
    		child[i]->dump(indent + 2);
	    }
    }
  }
}
int inbox(vector  & cpos, vector  & pos, real l){
  // ボックスに粒子が入っていれば0を返す関数
  for(int i=0; i<ndim; i++){
  	if(fabs(pos[i]-cpos[i]) > l*0.5) return 1;
  }
  return 0;
}
int bhnode::sanity_check(){
  //最終的に粒子がはみ出したりしていないかを確認できる関数
  int i, iret = 0;
  if (isleaf){
  	// this is the lowest level node. Things to check:
  	// all particles are in the cell
  	bhparticle * bp = bpfirst;
  	cout << "Leaf np="<< nparticle <<endl;
  	for(i = 0; i < nparticle; i++){
      real_particle * p = (bp+i)->get_rp();
      vector ppos = p->get_pos();
      if(inbox(pos,ppos,l)){
    		cerr << "Error, particle out of box ... \n";
    		dump();
    		return 1;
      }
  	}
  }else{
  	// This is the non-leaf node. Check the position and side
  	// length of the child cells and then check recursively..
  	cout << "Non Leaf " << pos  <<endl;
  	for(i=0;i<8;i++){
	    if (child[i] != NULL){
    		int err = 0;
  	        err = child[i]->sanity_check();
    		if (l*0.5 != child[i]->get_length()) err += 2;
    		vector relpos = pos-child[i]->get_pos();
    		for (int k = 0 ; k<ndim;k++){
  		    if (fabs(relpos[k]) != l*0.25) err += 4;
    		}
    		if (err){
  		    cerr << "Child " << i << " Error type = " << err << endl;
  		    dump();
    		}
    		iret += err;
	    }
  	}
  }
  return iret;
}

double set_rcritical(double B_2, real b_max){
  return (b_max/2.0)+sqrt((b_max*b_max/4.0)+sqrt(3.0*B_2/delta_int));
}
void bhnode::set_cm_quantities(){
  // セル内の重心情報を計算する関数
  int i;
  cmpos = 0.0;
  cmmass = 0.0;
  bhparticle * bp = bpfirst;
  bhparticle * bp_temp;
  if(isleaf){
    //ノードが末端であればそのままモーメント計算
  	for(i=0; i<nparticle; i++){
      bp_temp = bp + i;
  	    real mchild = bp_temp->get_rp()->get_mass();
  	    cmpos += mchild*bp_temp->get_rp()->get_pos();
  	    cmmass += mchild;
	  }
  }else{
    //ノードが末端でなければ再帰的に子ノードについてモーメント計算（重心計算）
  	for(i=0; i<8; i++){
	    if(child[i] != NULL){
    		child[i]->set_cm_quantities();
    		real mchild = child[i]->cmmass;
    		cmpos += mchild*child[i]->cmpos;
    		cmmass += mchild;
	    }
	  }
  }
  //重心位置は重心まわりのモーメントを総質量で割ればいい
  cmpos /= cmmass;

  //今のセル内で一番遠い粒子をさがす
  vector dx = 0.0;
  real r = 0.0, r_max = 0.0;
  if(nparticle == 1){
    b_max = 0.0;
  }else{
  	for(i=0; i<nparticle; i++){
      bp_temp = bp + i;
      dx = bp_temp->get_rp()->get_pos() - cmpos;
      r = sqrt(dx*dx);
      if(r > r_max){
        r_max = r;
      }
    }
  }
  b_max = r_max;

  #ifdef MULTIPOLE
    // ポテンシャルの多重極展開をするならばセル内部で多重極展開係数を求めておく必要がある
    if(nparticle != 1){
      // 球座標での角度と 半径に関する common term
      double theta = 0.0, phi = 0.0, br = 0.0, B_2 = 0.0;
      complex Y_term;
      for(i=0; i<nparticle; i++){
        bp_temp = bp + i;
        dx = bp_temp->get_rp()->get_pos() - cmpos;
        r = sqrt(dx*dx);
        real imass = bp_temp->get_rp()->get_mass();
        // セルの中心から見た時の角度を計算する 球座標系であることに注意
        if(r==0.0){theta = 0.0;}else{theta = acos(dx[2]/r);}
        if((dx[0]==0.0)&&(dx[1]==0.0)){phi = 0.0;}else{phi = atan2(dx[1],dx[0]);}
        // 多重極展開係数の計算
        for(int l=0; l<=p_max; l++){
          br = imass*pow(r, l);
          for(int m=-l, j=0; m<=l; m++, j++){
            Y_term = ComputeShericalHarmonics(l, -m, theta, phi);
            alpha[l*l+j].complex_add(br, Y_term);
          }
        }
        // 2次のモーメント計算
        B_2 += imass*r*r;
      }
      r_crit = set_rcritical(B_2, b_max);
      // cout << "r critical: " << r_crit << endl;
    }
  #endif
  // cout << "b_max:" << b_max << " l*l:" << l*l << endl;
}

// Tree法　での　重力計算関数　mainに移植中　なぜか main.cpp で定義すると早くなる

// static real total_interactions[16];  // 総相互作用計算回数
// static int nisum[16];                // 重力計算済み粒子数
// static real nisleaf[16];             // Treeの末端までいたった回数
// void clear_tree_counters(){
//   for(int i=0; i<16; i++){
//     total_interactions[i] = 0;
//     nisum[i] = 0;
//     nisleaf[i] = 0;
//   }
// }
// void print_tree_counters(){
//   real avg[16];
//   real average[4];
//   #ifdef OPENMP
//     for(int i=0; i<n_threads; i++){
//       avg[i] = total_interactions[i]/nisum[i];
//       PRC(nisum[i]); PRC(total_interactions[i]); PRC(nisleaf[i]);
//       PRL(avg[i]);
//       cout << "n_interact_average:" << i << " " << avg[i] << endl;
//       average[0] += nisum[i];
//       average[1] += total_interactions[i];
//       average[2] += nisleaf[i];
//     }
//     // average[0] /= n_threads;
//     // average[1] /= n_threads;
//     // average[2] /= n_threads;
//     average[3] = average[1]/average[0];
//     cout << "n_total = " << average[0] << " total_interactions = " << average[1]
//     << " nisleaf_total = " << average[2] << endl;
//     cout << "n_interact_average:" << average[3] << endl;
//   #else
//     avg[0] = total_interactions[0]/nisum[0];
//     PRC(nisum[0]); PRC(total_interactions[0]); PRC(nisleaf[0]);
//     PRL(avg[0]);
//     cout << "n_interact_average:" << " " << avg[0] << endl;
//   #endif
// }
// void bhnode::accumulate_force_from_tree(vector & ipos, real eps2, real theta2,
// 					                              vector & acc, real & phi, int n_th){
//   // Tree法 を使った重力の計算
//   vector dx = cmpos-ipos;
//   real r2 = dx*dx;
//   real r_crit = l*l;
//   int i;
//   if((r2*theta2 > r_crit)){
//     // node and position is well separated
//     // or the number of particles in node is Only 1
//
//     accumulate_force_from_point(dx,r2,eps2,acc,phi,cmmass);
//     // accumulate_force_multipole(dx, r2, b_max, acc, phi);
//     // total_interactions[n_th] += 1;
//   } else if((nparticle == 1)) {
//     accumulate_force_from_point(dx,r2,eps2,acc,phi,cmmass);
//     // total_interactions[n_th] += 1;
//   } else {
//     int i;
//     if(isleaf){
//       //セルの末端まで来てしまったときは中の粒子をすべてdirectで解く
//       bhparticle * bp = bpfirst;
//       for(i=0; i<nparticle; i++){
//         // nisleaf[n_th] ++;
//         dx = ((bp+i)->get_rp())->get_pos()-ipos;
//         r2 = dx*dx;
//         accumulate_force_from_point(dx,r2,eps2,acc,phi,(bp+i)->get_rp()->get_mass());
//         // total_interactions[n_th] += 1;
//       }
//     }else{
//       for(i=0; i<8; ++i){
//         if(child[i] != NULL){
//           //近いけどまだ下にセルがあるとき
//           child[i]->accumulate_force_from_tree(ipos,eps2,theta2,acc,phi,n_th);
//         }
//       }
//     }
//   }
//   // if 文を入れ替えたver. isleafがカウントされる（ノード内に１粒子しか入っていないもの）
//   // if(isleaf){
//   //   //セルの末端まで来てしまったときは中の粒子をすべてdirectで解く
//   //   bhparticle * bp = bpfirst;
//   //   for(i=0; i<nparticle; i++){
//   //     nisleaf[n_th] ++;
//   //     dx = ((bp+i)->get_rp())->get_pos()-ipos;
//   //     r2 = dx*dx;
//   //     accumulate_force_from_point(dx,r2,eps2,acc,phi,(bp+i)->get_rp()->get_mass());
//   //     total_interactions[n_th] ++;
//   //   }
//   // }else{
//   //   if((r2*theta2 > l*l) || (nparticle == 1)){
//   //     // node and position is well separated
//   //     // or the number of particles in node is Only 1
//   //     accumulate_force_from_point(dx,r2,eps2,acc,phi,cmmass);
//   //     total_interactions[n_th] ++;
//   //   }else{
//   //     for(i=0; i<8; i++){
//   //       if(child[i] != NULL){
//   //         //近いけどまだ下にセルがあるとき
//   //         child[i]->accumulate_force_from_tree(ipos,eps2,theta2,acc,phi,n_th);
//   //       }
//   //     }
//   //   }
//   // }
// }
// void calculate_gravity_using_tree(real_particle * rp, bhnode * bn, real eps2, real theta2, int n_th){
//   // 値渡しの関数
//   rp->acc_gravity = 0;
//   rp->phi_gravity = rp->mass/sqrt(eps2);
//   bn->accumulate_force_from_tree(rp->pos,eps2,theta2,rp->acc_gravity, rp->phi_gravity, n_th);
//   // nisum[n_th] += 1;
// }
// void calculate_gravity(int nbody, real_system pb, bhnode * bn, bhparticle * bp, real eps2, real theta2){
//   // morton keyを設定　Treeの作成　重力計算まで行う関数
//   double start_key,      end_key;
//   double start_tree,     end_tree;
//   double start_particle, end_particle;
//   double start_test,     end_test;
//
//   // Morton key計算 Part
//   int nkey = 0;
//   start_key = omp_get_wtime();
//     real rsize = initialize_key(nbody, pb.get_particle_pointer(), nkey, bp);
//     PRL(rsize);
//   end_key = omp_get_wtime();
//   cerr << "initialize_key time:" << (end_key-start_key) << "[sec]" << endl;
//
//   // Tree作成のための準備
//   for(int j=0; j<nbody; j++) (bn+j)->clear();
//   //ここの初期サイズrsizeを大きくすると精度が増す
//   bn->assign_root(vector(0.0), 2*rsize, bp, nbody);
//   bhnode * btmp = bn+1;
//   int heap_remainder = nbody*2;
//   BHlong key = 0;
//
//   // Tree作成 Part
//   start_tree = omp_get_wtime();
//     bn->create_tree_recursive(btmp, heap_remainder, key, default_key_length, 1);
//   end_tree = omp_get_wtime();
//   cerr << "create tree constructure:" << (end_tree-start_tree) << "[sec]" << endl;
//
//   // 重心情報計算 Part
//   start_tree = omp_get_wtime();
//     bn->set_cm_quantities();
//   end_tree = omp_get_wtime();
//   cerr << "calc. cm quantities:" << (end_tree-start_tree) << "[sec]" << endl;
//
//   // 重力計算のための準備
//   clear_acc_and_phi(pb.get_particle_pointer(),nbody);
//   // PRL(bn->sanity_check());
//   // bn->dump();
//   clear_tree_counters();
//   bhparticle * bp_test = bn->get_bpfirst();
//   real_particle * p_test = NULL;
//
//   // 重力計算 Part
//   int n_th = 0;
//   start_test = omp_get_wtime();
//     #ifdef OPENMP
//       #pragma omp parallel private(n_th)
//       {
//         // n_th = omp_get_thread_num();
//         #pragma omp for
//     #endif
//         for(int i=0; i<nbody; i++){
//           // start_particle = omp_get_wtime();
//             p_test = (bp_test+i)->get_rp();
//             // p_test = pb.get_particle_pointer()+i;
//             calculate_gravity_using_tree(p_test, bn, eps2, theta2, n_th);
//           // end_particle = omp_get_wtime();
//           // p_test->set_time(end_particle-start_particle);
//         }
//     #ifdef OPENMP
//       }
//     #endif
//   end_test = omp_get_wtime();
//   cerr << "time of calc. gravity:" << (end_test-start_test) << "[sec.]" << endl;
//
//   // 結果の表示
//   print_tree_counters();
// }

// Multipole expansion による重力計算 part
// #ifdef MULTIPOLE
//   void bhnode::accumulate_force_multipole(vector dx, real r2, real b, vector & acc, real & phi){
//     // accumulate_force_from_tree の multipole バージョン
//     double r = sqrt(r2);
//     double th = acos(dx[2]/r);
//     double ph = atan2(dx[1],dx[0]);
//     complex Y_ij[(2*p_max+1)+(p_max*p_max)];
//     double add_phi = 0.0;
//     vector add_acc = 0.0;
//     Y_term_ij(Y_ij, th, ph);
//     double br;
//     for(int l=0; l<=p_max; l++){
//       // br = 4.0*M_PI*G0*pow(1.0/r, l+1)/(2.0*l+1.0);
//       br = pow(1.0/r, l+1);
//       for(int m=-l, j=0; m<=l; m++, j++){
//         add_phi    -= br * (alpha[l*l+j].get_real()*Y_ij[l*l+j].get_real()
//                       + alpha[l*l+j].get_imag()*Y_ij[l*l+j].get_imag());
//         add_acc[0] -= -(l+1)*br/r * (alpha[l*l+j].get_real()*Y_ij[l*l+j].get_real()
//                       + alpha[l*l+j].get_imag()*Y_ij[l*l+j].get_imag());
//       }
//     }
//     phi += add_phi;
//     get_rthph_acc(add_acc, r, th, ph, alpha, Y_ij);
//     acc += get_xyz_acc(add_acc, th, ph);
//   }
// #endif
