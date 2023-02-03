#include "unionfind.hpp"

UnionFind::UnionFind(int _n) : n(_n), par(_n, -1) {}
bool UnionFind::merge(int a, int b) {
  a = root(a), b = root(b);
  if (a == b)
    return false;
  if (par[a] > par[b])
    std::swap(a, b);
  par[a] += par[b];
  par[b] = a;
  return true;
}
bool UnionFind::isSame(int a, int b) {
  return root(a) == root(b);
}
int UnionFind::root(int a) {
  if (par[a] < 0)
    return a;
  return par[a] = root(par[a]);
}
int UnionFind::size(int a) {
  return -par[root(a)];
}
std::vector<std::vector<int>> UnionFind::groups() {
  std::vector<int> leader_buf(n), group_size(n);
  for (int i = 0; i < n; i++) {
    leader_buf[i] = root(i);
    group_size[leader_buf[i]]++;
  }
  std::vector<std::vector<int>> result(n);
  for (int i = 0; i < n; i++)
    result[i].reserve(group_size[i]);
  for (int i = 0; i < n; i++)
    result[leader_buf[i]].push_back(i);
  result.erase(remove_if(result.begin(), result.end(), [&](const std::vector<int> &v) { return v.empty(); }), result.end());
  return result;
}