#pragma once
/*-----------------------------------------------------------------------------
 *  nbody-particle : basic class for simple nbody implementation
 *-----------------------------------------------------------------------------
 */
#include "myvector.hpp"
#ifndef ONED
#define THREED
#endif
#define REAL_GRAVITY

//particleクラスの定義
class nbody_particle{
  private:
  public:
    vector pos;                //粒子の位置
    vector dpos;               //粒子の重心からの位置
  	vector vel;                //粒子の速度
    vector vel_cycle;          //回転速度
  	vector acc_gravity;        //粒子の加速度
    vector acc_direct;         //粒子の加速度（直接計算）
    vector acc_external;       //粒子が外場から受ける加速度
  	real phi_gravity;          //粒子のポテンシャル
  	real phi_gravity_external; //外場から受けるポテンシャル
  	real mass;                 //粒子の質量
  	int index;                 //粒子の番号
    real time_gravity;         //1粒子の計算にかかる時間

  	nbody_particle(){
      //初期化関数
	    pos = 0.0;
      dpos = 0.0;
	    vel = 0.0;
      vel_cycle = 0.0;
	    acc_gravity = 0.0;
      acc_direct = 0.0;
      acc_external = 0.0;
	    phi_gravity =  mass =  0.0;
	    index = 0;
      time_gravity = 0.0;
  	}
    //引数にいれたvectorを粒子のpos/vec/acc/phiに入れる関数
    void set_pos(const vector& new_pos) {pos = new_pos;}
    void set_dpos(const vector & new_dpos) {dpos = new_dpos;}
    void set_vel(const vector& new_vel) {vel = new_vel;}
    void set_acc_gravity(const vector& new_acc) {acc_gravity = new_acc;}
    void set_phi_gravity(real new_phi)  {phi_gravity = new_phi;}

    //pos/vel/acc/phiの初期化関数
  	void clear_pos() {pos = 0.0;}
  	void clear_vel() {vel = 0.0;}
  	void clear_acc_phi_gravity(){
      acc_gravity = 0.0;
      acc_direct = 0.0;
      acc_external = 0.0;
      phi_gravity = 0.0;
    }
    //無限遠を0とするようなポテンシャルではr->∞でmass/epsとなるので最初にいれておく
  	void correct_phi_self_gravity(real epsinv) {phi_gravity += mass*epsinv;}

    //引数に入れただけpos/velに加える関数
  	void inc_pos(vector d_pos) {pos += d_pos;}
    void inc_dpos(vector d_dpos) {dpos += d_dpos;}
  	void inc_vel(vector d_vel) {vel += d_vel;}
    void inc_acc(vector d_acc) {acc_gravity += d_acc;}
  	void update_vel(real dt) {vel += dt*acc_gravity;}
  	void update_pos(real dt) {pos = (pos+dt*vel).readjust();}

    //スケールファクターをpos/velにかける関数
  	void scale_pos(const real scale_factor) {pos *= scale_factor; }
  	void scale_vel(const real scale_factor) {vel *= scale_factor; }

    //粒子のposとposのポインタを返す関数
  	vector  get_pos()  {return   pos;}
  	vector* get_posp() {return & pos;}
    vector  get_dpos() {return   dpos;}
    real get_radius()    {return sqrt(pos*pos);}
    real get_radius_cm() {return sqrt(dpos*dpos);}

    //粒子のvel/phi/accを返す関数
  	vector get_vel() {return vel;}
  	real   get_phi_gravity() {return phi_gravity;}
  	real   get_phi_gravity_external() {return phi_gravity_external;}
  	vector get_acc_gravity() {return acc_gravity;}
    vector get_acc_direct() {return acc_direct;}
    vector get_acc_external() {return acc_external;}


    //引数を粒子のmassに入れる関数 //粒子の質量を返す関数
  	void set_mass(real m)	{mass = m;}
    real get_mass() {return mass;}

    //引数を粒子のindexに入れる関数 //粒子のindexを返す関数
  	void set_index(int i){index = i;}
  	int get_index(){return index;}

    //1粒子あたりにかかる計算時間の出力
    void set_time(real new_time) {time_gravity = new_time;}
    real get_time(){return time_gravity;}

    //時間積分（leap-frog法）
    void predict(real dt){
      real dt2 = dt*dt*0.5;
      pos += dt*vel+dt2*acc_gravity;
      vel += (dt*0.5)*acc_gravity;
    }//leap-frog法*/
    void correct(real dt){
      vel += 0.5*dt*acc_gravity;
    }//leap-frog後の速度を１/2タイムステップ進めて位置の時間と一致させる
    void predict2(real dt){
      pos += dt*vel;
    }//leap-frog法その２*/

    //print用関数
  	void read(istream & );
  	void write(ostream & );
    //木構造がうまく作れているかを確認するための関数
  	void dump();

    //エネルギー計算関係の関数
  	real kinetic_energy();
  	real energy();
  	real get_ke(){return 0.5*mass*vel*vel;}
    real get_ke_rcut(){
      double ke = ((vel-vel_cycle)*(vel-vel_cycle))*mass*0.5;
      return ke;
    }
    real get_pe(){return phi_gravity*mass;}
    real get_momentum(){return sqrt((vel-vel_cycle)*(vel-vel_cycle))*mass;}
};

//nbody_systemクラスの定義（粒子全体の情報をもつクラス）
class nbody_system{
  private:
  public:
    int n;               //粒子数
    int nsize;           //粒子数
    int nplummer;        //Plummer球の個数
    nbody_particle * pb; //粒子へのポインタ

	  nbody_system(){
      //初期化関数
	    n = 0;
	    nsize = 0;
	    pb = NULL;
	  }

    //エネルギー関係の関数
  	real kinetic_energy(int nbody, real r0);
    real potential_energy(int nbody, real r0);
    real momentum(int nbody, real r0);
    real energy_error(int n, int check, real ke, real pe);

    //一様球をつくる関数
    void create_uniform_sphere(int nbody, real power_index, real r0);
    // COLD_COLLAPSE
    void create_coldcollapse(int nbody, real power_index, real r0, real vir0, real eps2);
    // INCLUDE_DATASET
    void INCLUDE_PARTICLE_DATASET(int nbody);
    // PLUMMER
    void create_plummer_sphere(int nplummer, int nbody, real shift_theta[], real shift_radius[]);

    //時間積分の関数
  	void integrate( real dt);

  	nbody_particle * get_particle_pointer() {return pb;}
};
