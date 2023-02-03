#pragma once
#include <vector>

class UnionFind {
 private:
  int n;
  std::vector<int> par;

 public:
  UnionFind(int _n);
  bool merge(int a, int b);
  bool isSame(int a, int b);
  int root(int a);
  int size(int a);
  std::vector<std::vector<int>> groups();
};