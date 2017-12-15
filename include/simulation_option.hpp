#pragma once

// 並列化を行う
// #define OPENMP

// 多重極展開して計算する
#define MULTIPOLE

// morton keyの有効桁数を３０までにする
// #define using_30bits

// 同じディレクトリ内にあるファイルを初期条件とする
// #define INCLUDE_DATASET

// 一様球のCold Collapseを行う
// #define COLD_COLLAPSE

// Plummer球１つのみを用いる
#define PLUMMER_ONLY

// 時間積分を行う
// #define INTEGRATE

// 加速度の相対誤差の計算をおこなう
#define COMP_ACC

// ファイル出力まで行う それぞれ　加速度 / エネルギー / 密度
// #define ACC_OUTPUT
// #define ENERGY_OUTPUT
// #define DENSITY_OUTPUT

// デバッグ用関数の結果出力　それぞれ　Mortonナンバー / Treewalk
// #define MORTON_OUTPUT
// #define TREE_WALK_OUTPUT

// gnuplotを用いた出力を行うかどうか
#define GNUPLOT
// #define GNUPLOT_POS
