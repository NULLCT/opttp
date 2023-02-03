#include "directed_graph.hpp"

#include <bits/stdint-intn.h>

template <class T>
bool DirectedGraph::chmin(T &a, const T &b) {
  if (a > b) {
    a = b;
    return true;
  }
  return false;
}
void DirectedGraph::init(int64_t _n) {
  DirectedGraph::N = _n;
  g.resize(_n);
}
DirectedGraph::DirectedGraph(int64_t _n) {
  init(_n);
}
void DirectedGraph::add(int64_t s, int64_t t, int64_t cost) {
  g[s].push_back({t, cost});
}
// O(V^3)
DirectedGraph::warshallfloyd_return DirectedGraph::warshallfloyd() {
  std::vector<std::vector<int64_t>> dist(N, std::vector<int64_t>(N, INF)), next(N, std::vector<int64_t>(N, INF));
  for (int64_t i = 0; i < N; i++)
    dist[i][i] = 0;
  for (int64_t i = 0; i < N; i++)
    for (int64_t j = 0; j < N; j++)
      next[i][j] = j;
  for (int64_t i = 0; i < N; i++)
    for (auto &j : g[i])
      chmin(dist[i][j.to], j.cost);
  for (int64_t k = 0; k < N; k++)
    for (int64_t i = 0; i < N; i++)
      for (int64_t j = 0; j < N; j++)
        if (chmin(dist[i][j], dist[i][k] + dist[k][j]))
          next[i][j] = next[i][k];
  return {dist, next};
}
// O(E+VlogV)
std::vector<int64_t> DirectedGraph::dijkstra(int64_t s) {
  std::vector<int64_t> d(N, INF);
  d[s] = 0;
  std::priority_queue<std::pair<int64_t, int64_t>, std::vector<std::pair<int64_t, int64_t>>, std::greater<std::pair<int64_t, int64_t>>> que;
  que.push({0, s});
  while (!que.empty()) {
    std::pair<int64_t, int64_t> p = que.top();
    que.pop();
    int64_t v = p.second;
    if (d[v] < p.first)
      continue;
    for (auto e : g[v]) {
      if (d[e.to] > d[v] + e.cost) {
        d[e.to] = d[v] + e.cost;
        que.push({d[e.to], e.to});
      }
    }
  }
  return d;
}
// O(V*E)
// https://zenn.dev/reputeless/books/standard-cpp-for-competitive-programming/viewer/bellman-ford
DirectedGraph::bellmanford_return DirectedGraph::bellmanford_SPFA(int64_t st, int64_t ed) {
  std::vector<int64_t> counts(g.size()), distances(g.size(), INF), path;
  std::vector<bool> inqueue(g.size());
  std::queue<int64_t> q;
  std::vector<int64_t> p(g.size(), -1);
  distances[st] = 0;
  q.push(st);
  inqueue[st] = true;
  while (!q.empty()) {
    const int64_t from = q.front();
    q.pop();
    inqueue[from] = false;
    for (const auto &edge : g[from]) {
      const long long d = (distances[from] + edge.cost);
      if (d < distances[edge.to]) {
        distances[edge.to] = d;
        p[edge.to] = from;
        if (!inqueue[edge.to]) {
          q.push(edge.to);
          inqueue[edge.to] = true;
          ++counts[edge.to];
          if (g.size() < counts[edge.to])
            return {std::vector<int64_t>(), std::vector<int64_t>(), true};
        }
      }
    }
  }
  if (distances[ed] != INF) {
    for (int64_t i = ed; i != -1; i = p[i])
      path.push_back(i);
    reverse(path.begin(), path.end());
  }
  return {path, distances, false};
}
// O(V+E)
std::vector<int64_t> DirectedGraph::topologicalSort() {
  std::vector<int64_t> d, ind(N);
  for (int64_t i = 0; i < N; i++)
    for (auto e : g[i])
      ind[e.to]++;
  std::queue<int64_t> que;
  for (int64_t i = 0; i < N; i++)
    if (ind[i] == 0)
      que.push(i);
  while (!que.empty()) {
    int64_t now = que.front();
    d.push_back(now);
    que.pop();
    for (auto e : g[now]) {
      ind[e.to]--;
      if (ind[e.to] == 0)
        que.push(e.to);
    }
  }
  return d;
}