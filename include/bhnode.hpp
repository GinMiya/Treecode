#pragma once
/*-----------------------------------------------------------------------------
 *  BHtree : basic class for C++ implementation of BH treecode
 *-----------------------------------------------------------------------------*/
//default_key_lengthは最深のlevelに対応
const int default_key_length = 20;
//default_ix_offsetは最深セルでの1辺の長さのファクター2^levelに対応
const unsigned long int default_ix_offset = 1 << (default_key_length-1);

using namespace std;
#define real double

#include "myvector.hpp"
#include "particle.hpp"
#include "complex.hpp"
#include "multipole.hpp"
#include "constant.hpp"

// 名前の変更
typedef nbody_particle real_particle;
typedef nbody_system real_system;
typedef long BHlong;
// typedef  long long BHlong;

//node内での粒子情報（bhparticle）のクラス定義
class bhparticle{
 private:
    //粒子へのポインタ
    real_particle * rp;
    //morton key
    BHlong key;
 public:
    bhparticle(){
      //初期化関数
     	rp = NULL;
     	key = 0;
    }

    void set_key(real rscale, int ioffset, int keybits);
    // void set_key(real rscale, int ioffset);
    void set_rp(real_particle * p){rp = p;}

    BHlong get_key(){return key;}
    real_particle * get_rp(){return rp ;}

    int friend compare_key( bhparticle * p1,  bhparticle * p2);
    void friend sort_bh_array( bhparticle * r, int lo, int up );
};

//nodeのクラス定義
class bhnode{
  private:
    vector pos;           //nodeの重心
    real l;               //nodeの1辺の長さ
    bhnode * child[8];    //子nodeへのポインタ
    bhparticle * bpfirst; //node内の最初の粒子へのポインタ
    int nparticle;        //node内に入っている粒子の個数
    int isleaf;           //このノードが末端なら1をもつような定数
    int level;            //ノードのレベル
    vector cmpos;         //nodeの中心
    real cmmass;          //nodeの総質量
    real b_max;           //ノード内で一番遠い粒子までの半径　精度いまひとつ
    real r_crit;          //ノードの持つ臨界半径（Warren Salmon）精度低
    complex alpha[(2*p_max+1)+(p_max*p_max)];  //多重極展開係数
  public:
    bhnode(){
      //初期化関数
  	  pos = 0.0;
  	  l = 0.0;
  	  for(int i=0; i<8; i++)child[i] = NULL;
  	  bpfirst = NULL;
    	nparticle = 0;
    	isleaf = 1;
      level = 0;
    	cmpos = 0.0;
    	cmmass = 0.0;
      b_max = 0.0;
      r_crit = 0.0;
      for(int l=0; l<p_max; l++){
        for(int m=-l, j=0; m<=l; m++, j++){
          alpha[l*l+j].set_real(0.0);
          alpha[l*l+j].set_imag(0.0);
        }
      }
    }
    void clear(){
      //初期化関数その2
    	pos = 0.0;
    	l = 0.0;
    	for(int i=0; i<8; i++)child[i] = NULL;
    	bpfirst = NULL;
    	nparticle = 0;
    	isleaf = 1;
    	cmpos = 0.0;
    	cmmass = 0.0;
      b_max = 0.0;
      r_crit = 0.0;
      for(int l=0; l<p_max; l++){
        for(int m=-l, j=0; m<=l; m++, j++){
          alpha[l*l+j].set_real(0.0);
          alpha[l*l+j].set_imag(0.0);
        }
      }
    }
    // ほしい物理量を返す
    bhparticle * get_bpfirst() {return bpfirst;}
    vector  get_pos()    {return pos;}
  	vector* get_posp()   {return & pos;}
    vector  get_cmpos()  {return cmpos;}
    vector* get_cmposp() {return & cmpos;}
    real    get_length() {return l;}
    real    get_bmax()   {return b_max;}
    int     get_nparticle()  {return   nparticle;}
    int*    get_nparticlep() {return & nparticle;}
    int     get_level()  {return level;}
    int     get_isleaf() {return isleaf;}

    // 設定したい物理用に引数を代入
    void set_pos(vector newpos) {pos = newpos;}
    void set_length(real newl)  {l   = newl;}

    // 木構造を作る関数
    void assign_root(vector root_pos, real length, bhparticle * bp, int nparticle);
    void create_tree_recursive(bhnode * & heap_top, int & heap_remainder,
	                              BHlong current_key, int current_level, int n_critical);

    // 確認用関数
    void dump(int indent = 0);
    int sanity_check();

    // 重心計算関数
    void set_cm_quantities();
    // 多重極展開係数の計算
    void multipole_exp_coeff();

    // セル情報やTree Walkに関する情報を集める関数
    void output_data_in_node(int tstep, int level);
    void accumulate_force_from_tree(vector & ipos, real eps2, real theta2,
                          				   vector & acc, real & phi, int n_th);
    void accumulate_tree_walk(real_particle * rp, real theta2, int n_th);

    void accumulate_force_multipole(vector dx, real r2, real b, vector & acc, real & phi);
};

void clear_tree_counters();
void print_tree_counters();

void calculate_gravity_using_tree(real_particle * rp, bhnode * bn, real eps2, real theta2, int n_th);
void calculate_gravity(int nbody, real_system pb, bhnode * bn, bhparticle * bp, real eps2, real theta2);

void tree_walk_counter(real_particle * rp, bhparticle * bp, int nbody, vector cmpos, real length, int level);
void search_tree_walk(int nbody, real_system pb, bhnode * bn, bhparticle * bp, real eps2, real theta2);
