#pragma once
#include <vector>

struct STATE {                                // 状態構造体
  std::vector<std::vector<int64_t>> acthist;  // 頂点移動履歴
  std::vector<int64_t> pathlen;               // パス長
};