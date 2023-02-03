#pragma once

#include <unordered_set>

class Carrier {
 public:
  int pos;                                      // Carrierの頂点番号
  std::unordered_multiset<int64_t> passengers;  // 頂点passengers[i]に行きたい人たち

  Carrier(int _pos) : pos(_pos) {}
};