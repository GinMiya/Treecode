// morton keyを計算する関数

// ライブラリ読み込み
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sstream>

// 面倒なのでstdのnamespaceをつかう
using namespace std;
#define real double

// ヘッダーファイルの読み込み
#include "simulation_option.hpp"
#include "bhnode.hpp"         //  "particle.hpp" "myvector.hpp" "constant.hpp"

// 値を一番近い偶数にする関数
double ConvertEvenValue(double val){
  int val_int = (int)val;
  if(val_int%2){
    val_int += 1;
  }
  return (double)val_int;
}

// morton keyを２進数で出力するための関数
string to_binString(unsigned int val){
  if(!val) return string("0");
  string str;
  while(val != 0){
    if((val & 1) == 0){ //valは偶数かどうか
      str.insert(str.begin(), '0'); //偶数なら0を入れる
    }else{
      str.insert(str.begin(), '1'); //奇数なら1を入れる
    }
    val >>= 1;
  }
  return str;
}

// morton key作成関数　newとoldでかわらず
BHlong conv_to_morton_old(int ix, int keybits){
  BHlong dum = 0;
  //    char st[256];
  //    cerr << "conv_to_morton "; PR(ix); PRL(keybits);
  //    sprintf(st,"%lo",ix);
  //    cerr << "ix = " << st << endl;
  int i, j;
  for(i=j=0; i<keybits; i++,j+=3){
    //3次元なのでjは3bitずつ増やす
  	if (ix & (1<<i)){
      //引数ixが、1からibit左シフトさせた値に対して余り0のとき
      //dumには１bitとjbitの和を代入する
	    dum |= ((BHlong) 1)<<j;
	  }
  }
  //sprintf(st,"%lo",dum);
  //    cerr << "dum = " << st << endl;
  return dum;
}
BHlong conv_to_morton(int ix, int keybits){
  //newタイプのモートンキーを作成
  BHlong dum = 0;
  //    char st[256];
  //    cerr << "conv_to_morton "; PR(ix); PRL(keybits);
  //    sprintf(st,"%lo",ix);
  //    cerr << "ix = " << st << endl;
  int i, j;
  for(i= j= 0; i<keybits; i++,j+=3){
  	if ((ix>>i) & 1){
      //ixをibit右シフトさせたのが1に対して余り0のとき
	    dum |= ((BHlong) 1)<<j;
  	}
  }
  //sprintf(st,"%lo",dum);
  //    cerr << "dum = " << st << endl;
  return dum;
}

// 3次元空間でmorton keyを算出するための関数
BHlong dilate3D(BHlong val){
  #ifdef using_30bits
    val = (val * 0x000010001) & 0xff0000ff;
    val = (val * 0x000000101) & 0x0f00f00f;
    val = (val * 0x000000011) & 0xc30c30c3;
    val = (val * 0x000000005) & 0x49249249;
  #else
    val = (val * 0x100000001) & 0x7fff00000000ffff;
    val = (val * 0x000010001) & 0x00ff0000ff0000ff;
    val = (val * 0x000000101) & 0x700f00f00f00f00f;
    val = (val * 0x000000011) & 0x30c30c30c30c30c3;
    val = (val * 0x000000005) & 0x1249249249249249;
  #endif
  return (val);
}
BHlong getMorton3D(const BHlong ix, const BHlong iy, const BHlong iz){
  return ((dilate3D(ix) << 2) | (dilate3D(iy) << 1) | (dilate3D(iz)));
}
// 第一引数posに入れた３次元座標を整数値で番地を割り当てる関数
inline BHlong construct_key(const vector & pos, real rscale, int ioffset, int keybits){
  BHlong ix[3];
  vector p = pos;
  for(int i=0; i<3; i++){
    //セルの1辺の長さがわかっているのでioffset分のセルが使われていると思って、
  	ix[i] = (BHlong) (p[i]*rscale+ioffset);
  }
  return getMorton3D(ix[0],ix[1],ix[2]);
}

//qsortのためのmorton keyの比較を行う関数
int compare_key(const void * p1, const void * p2){
  //２粒子を引数にして粒子間のkeyを比較
  bhparticle * pp = (bhparticle*)p1;
  bhparticle * qq = (bhparticle*)p2;
  BHlong comp = ((BHlong) pp->get_key()) - ((BHlong) qq->get_key());
  if(comp > 0L){
    //keyの差が0よりも大きいなら1
  	return 1;
  }else if (comp == 0L){
    //keyが同じであれば0
  	return 0;
  }else{
    //keyの差が0よりも小さければ-1を返す
  	return -1;
  }
}

//construct_keyを実行するための引数を与える関数
void bhparticle::set_key(real rscale, int ioffset, int keybits){
  key = construct_key(rp->get_pos(), rscale, ioffset, keybits);
  //    PRL(key);
}
//morton key を初期化して系のサイズを返す関数
real initialize_key(int nbody, real_particle * rp, int & nkeysize, bhparticle * bhp){
  //全粒子に対して系の大きさが足りなければ系の大きさを２倍にしていく
  real rmax = 1;
  for(int i=0; i<nbody; i++){
  	vector p = (rp+i)->get_pos();
  	for(int k=0; k<3; k++){
	    // if(fabs(p[k])>=rmax) rmax *= 2;
      if(fabs(p[k])>=rmax) rmax = p[k];
  	}
  }
  cout << "clear" << endl;
  rmax = ConvertEvenValue(rmax);
  //最深のセルの一辺の長さ
  real rscale = default_ix_offset/rmax;
  //全粒子に対してmorton keyを作る部分
  #ifdef OPENMP
    #pragma omp parallel for
  #endif
  for(int i=0; i<nbody; i++){
    bhparticle * p = bhp + i;
    p->set_rp(rp+i);
    p->set_key(rscale, default_ix_offset, default_key_length);
    // PR(i); PRL(p->get_key());
  }
  // morton keyの大きさでソートする
  qsort(bhp, nbody, sizeof(bhparticle), compare_key);
  return rmax;
}
