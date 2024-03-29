/**
 * @file experimental_model.hpp
 * @brief 運送モデルについての情報を保存します
 */
#pragma once
#include <iostream>
#include <vector>

#include "carrier.hpp"
#include "directed_graph.hpp"
#include "state.hpp"

struct MODEL {                                 // 運送モデル
  int64_t V;                                   // 頂点数
  int64_t E;                                   // 辺の数
  DirectedGraph G;                             // グラフ
  int64_t T;                                   // 最大シミュレート時間
  std::vector<int64_t> M;                      // 頂点ごとの運送数
  int64_t N;                                   // クエリ数
  std::vector<std::pair<int64_t, int64_t>> Q;  // クエリ
  std::vector<std::vector<int64_t>> DIST;      // 頂点間距離
  std::vector<std::vector<int64_t>> NEXT;      // ある頂点に行くため次どこに行けばいいか
  STATE generateState();
  friend std::istream &operator>>(std::istream &is, MODEL &model);
};