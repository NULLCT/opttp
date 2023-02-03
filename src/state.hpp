/**
 * @file simulated_annealing.hpp
 * @brief 焼きなまし法での状態を保存します
 */
#pragma once
#include <set>
#include <tuple>
#include <unordered_set>
#include <vector>

struct STATE {                                // 状態構造体
  std::vector<std::vector<int64_t>> acthist;  // 頂点移動履歴
  std::vector<int64_t> pathlen;               // パス長
};