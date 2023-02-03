#pragma once

#include <queue>
#include <vector>

#include "consts.hpp"

class DirectedGraph {
  int64_t N;
  template <class T>
  bool chmin(T &a, const T &b);

 public:
  struct Edge {
    int64_t to, cost;
  };
  std::vector<std::vector<Edge>> g;
  void init(int64_t _n);
  DirectedGraph(int64_t _n = 0);
  void add(int64_t s, int64_t t, int64_t cost);
  struct warshallfloyd_return {
    std::vector<std::vector<int64_t>> dist, next;
  };
  warshallfloyd_return warshallfloyd();      // O(V^3)
  std::vector<int64_t> dijkstra(int64_t s);  // O(E+VlogV)
  struct bellmanford_return {
    std::vector<int64_t> path, distances;
    bool hascycle;
  };
  bellmanford_return bellmanford_SPFA(int64_t st, int64_t ed);  // O(V*E)
  std::vector<int64_t> topologicalSort();                       // O(V+E)
};