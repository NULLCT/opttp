#pragma once

#include "util.hpp"

using namespace std;

constexpr int V = 100;  // 頂点数
constexpr int T = 1000;  // 最大シミュレート時間
constexpr int M = 10;    // 運送者数
constexpr int Q = 100;   // 荷物数
constexpr int L = 5;     // パス最大長

void generateTestCase(MODEL &model) {
  random_device rand;
  UnionFind uf(V);

  model.G.init(V);
  for (int i = 0; i < V; i++) {
    int to = rand() % V;
    while (i == to)
      to = rand() % V;
    if (uf.merge(i, to)) {
      int l = (rand() % L) + 1;
      model.G.add(i, to, l);
      model.G.add(to, i, l);
    }
  }
  for (int i = 0; i < V - 1; i++) {
    if (uf.merge(i, i + 1)) {
      int l = (rand() % L) + 1;
      model.G.add(i, i + 1, l);
      model.G.add(i + 1, i, l);
    }
  }

  model.V = V;
  model.E = [&model]() {
    int sum = 0;
    for (auto &i : model.G.g)
      sum += i.size();
    return sum;
  }();

  model.T = T;

  {
    model.M.resize(V);

    map<int, int> mp;
    for (int i = 0; i < M; i++)
      mp[rand() % V]++;

    for (int i = 0; i < V; i++)
      model.M[i] = mp[i];
  }

  for (int i = 0; i < Q; i++) {
    int a = rand() % V, b = rand() % V;
    while (a == b)
      b = rand() % V;

    model.Q.push_back({a, b});
  }

  auto res = model.G.warshallfloyd();
  model.DIST = res.first;
  model.NEXT = res.second;
}
