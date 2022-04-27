#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cctype>
#include <chrono>
#include <cinttypes>
#include <climits>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <deque>
#include <functional>
#include <iomanip>
#include <ios>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <numeric>
#include <ostream>
#include <queue>
#include <set>
#include <stack>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#define int long long
#define double long double
#define rep(i, e) for (int i = 0; i < (e); i++)
#define ALL(var) ((var).begin()), ((var).end())
#define D(var)                                                                 \
  cerr << "(\e[34m" << #var << "\e[0m \e[2m@" << __LINE__                      \
       << "\e[0m): \e[36m\e[1m" << var << "\e[0m" << endl;
#define DD(var) for_each(var.begin(), var.end(), [](const auto &i) { D(i); })

using namespace std;

constexpr array<int, 4> dx{1, 1, -1, -1}, dy{-1, 1, 1, -1};

template <class T, size_t S>
ostream &operator<<(ostream &_ostr, const array<T, S> &_v);
template <class T> ostream &operator<<(ostream &_ostr, const vector<T> &_v);
template <class T> ostream &operator<<(ostream &_ostr, const deque<T> &_v);
template <class T> ostream &operator<<(ostream &_ostr, const list<T> &_v);
template <class T, class Y>
ostream &operator<<(ostream &_ostr, const pair<T, Y> &_v);
template <class... Ts>
ostream &operator<<(ostream &_ostr, const tuple<Ts...> &t);
template <class T, class Y>
ostream &operator<<(ostream &_ostr, const map<T, Y> &_v);
template <class T, class Y>
ostream &operator<<(ostream &_ostr, const multimap<T, Y> &_v);
template <class T, class Y>
ostream &operator<<(ostream &_ostr, const unordered_map<T, Y> &_v);
template <class T, class Y>
ostream &operator<<(ostream &_ostr, const unordered_multimap<T, Y> &_v);
template <class T> ostream &operator<<(ostream &_ostr, const set<T> &_v);
template <class T> ostream &operator<<(ostream &_ostr, const multiset<T> &_v);
template <class T>
ostream &operator<<(ostream &_ostr, const unordered_set<T> &_v);
template <class T>
ostream &operator<<(ostream &_ostr, const unordered_multiset<T> &_v);

template <class T> void orange(ostream &_ostr, const T &_v) {
  if (_v.size() == 0)
    return;
  _ostr << *_v.begin();
  for (auto itr = next(_v.begin()); itr != _v.end(); itr++)
    _ostr << " " << *itr;
}
template <class T, size_t S>
ostream &operator<<(ostream &_ostr, const array<T, S> &_v) {
  orange(_ostr, _v);
  return _ostr;
}
template <class T> ostream &operator<<(ostream &_ostr, const vector<T> &_v) {
  orange(_ostr, _v);
  return _ostr;
}
template <class T> ostream &operator<<(ostream &_ostr, const deque<T> &_v) {
  orange(_ostr, _v);
  return _ostr;
}
template <class T> ostream &operator<<(ostream &_ostr, const list<T> &_v) {
  orange(_ostr, _v);
  return _ostr;
}
template <class T, class Y>
ostream &operator<<(ostream &_ostr, const pair<T, Y> &_v) {
  _ostr << _v.first << " " << _v.second;
  return _ostr;
}
template <class... Ts>
ostream &operator<<(ostream &_ostr, const tuple<Ts...> &_v) {
  bool first = true;
  apply(
      [&_ostr, &first](auto &&...args) {
        auto print = [&](auto &&val) {
          if (!first)
            _ostr << " ";
          (_ostr << val);
          first = false;
        };
        (print(args), ...);
      },
      _v);
  return _ostr;
}
template <class T, class Y>
ostream &operator<<(ostream &_ostr, const map<T, Y> &_v) {
  orange(_ostr, _v);
  return _ostr;
}
template <class T, class Y>
ostream &operator<<(ostream &_ostr, const multimap<T, Y> &_v) {
  orange(_ostr, _v);
  return _ostr;
}
template <class T, class Y>
ostream &operator<<(ostream &_ostr, const unordered_map<T, Y> &_v) {
  orange(_ostr, _v);
  return _ostr;
}
template <class T> ostream &operator<<(ostream &_ostr, const set<T> &_v) {
  orange(_ostr, _v);
  return _ostr;
}
template <class T> ostream &operator<<(ostream &_ostr, const multiset<T> &_v) {
  orange(_ostr, _v);
  return _ostr;
}
template <class T>
ostream &operator<<(ostream &_ostr, const unordered_set<T> &_v) {
  orange(_ostr, _v);
  return _ostr;
}
template <class T>
ostream &operator<<(ostream &_ostr, const unordered_multiset<T> &_v) {
  orange(_ostr, _v);
  return _ostr;
}

template <class T, size_t S>
istream &operator>>(istream &_istr, array<T, S> &_v);
template <class T> istream &operator>>(istream &_istr, vector<T> &_v);
template <class T> istream &operator>>(istream &_istr, deque<T> &_v);
template <class T, class Y> istream &operator>>(istream &_istr, pair<T, Y> &_v);

template <class T, size_t S>
istream &operator>>(istream &_istr, array<T, S> &_v) {
  for (auto &i : _v)
    _istr >> i;
  return _istr;
}
template <class T> istream &operator>>(istream &_istr, vector<T> &_v) {
  for (auto &i : _v)
    _istr >> i;
  return _istr;
}
template <class T> istream &operator>>(istream &_istr, deque<T> &_v) {
  for (auto &i : _v)
    _istr >> i;
  return _istr;
}
template <class T, class Y>
istream &operator>>(istream &_istr, pair<T, Y> &_v) {
  _istr >> _v.first >> _v.second;
  return _istr;
}

class DirectedGraph {
public:
  struct Edge {
    int to, cost;
  };
  int n;
  vector<vector<Edge>> g;

  DirectedGraph(int _n) : g(_n) { n = _n; }
  void add(int s, int t, int cost) {
    Edge e;
    e.to = t, e.cost = cost;
    g[s].push_back(e);
  }
  // O(V^3)
  vector<vector<int>> warshallfloyd() {
    vector<vector<int>> d(n, vector<int>(n, INT_MAX));
    for (int i = 0; i < n; i++)
      d[i][i] = 0;
    for (int i = 0; i < n; i++)
      for (Edge &j : g[i])
        d[i][j.to] = min(d[i][j.to], j.cost);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        for (int k = 0; k < n; k++)
          d[j][k] = min(d[j][k], d[j][i] + d[i][k]);
    return d;
  }
  // O(E+VlogV)
  vector<int> dijkstra(int s) {
    vector<int> d(n, INT_MAX);
    d[s] = 0;
    priority_queue<pair<int, int>, vector<pair<int, int>>,
                   greater<pair<int, int>>>
        que;
    que.push(make_pair(0, s));
    while (!que.empty()) {
      pair<int, int> p = que.top();
      que.pop();
      int v = p.second;
      if (d[v] < p.first)
        continue;
      for (auto e : g[v]) {
        if (d[e.to] > d[v] + e.cost) {
          d[e.to] = d[v] + e.cost;
          que.push(make_pair(d[e.to], e.to));
        }
      }
    }
    return d;
  }
  // O(V*E)
  vector<int> bellmanford(int s) {
    vector<int> d(n, INT_MAX);
    d[s] = 0;
    for (int _ = 0; _ < n; _++) {
      bool upd = false;
      for (int u = 0; u < n; u++)
        if (d[u] < INT_MAX)
          for (const auto &e : g[u]) {
            int v = e.to;
            if (d[v] > d[u] + e.cost)
              d[v] = d[u] + e.cost, upd = true;
          }
      if (!upd)
        return d;
    }
    queue<int> Q;
    for (int u = 0; u < n; u++)
      if (d[u] < INT_MAX)
        Q.emplace(u);
    while (!Q.empty()) {
      int u = Q.front();
      Q.pop();
      for (const auto &e : g[u]) {
        int v = e.to;
        if (d[v] > INT_MIN && (d[u] == INT_MIN || d[v] > d[u] + e.cost))
          d[v] = INT_MIN, Q.emplace(v);
      }
    }
    return d;
  }
  // O(V+E)
  vector<int> topologicalSort() {
    vector<int> d, ind(n);
    for (int i = 0; i < n; i++)
      for (auto e : g[i])
        ind[e.to]++;
    queue<int> que;
    for (int i = 0; i < n; i++)
      if (ind[i] == 0)
        que.push(i);
    while (!que.empty()) {
      int now = que.front();
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
};

signed main() {
  // Initializing vals
  int v, e;
  cin >> v >> e;
  DirectedGraph g(v);
  for (int i = 0; i < e; i++) {
    int a, b, d;
    cin >> a >> b >> d;
    g.add(a, b, d);
  }
  vector<int> m(v);
  cin >> m;
  int t;
  cin >> t;
  vector<vector<int>> x(v, vector<int>(t));
  cin >> x;

  auto res = g.warshallfloyd();
  DD(res);
}
