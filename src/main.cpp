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

#define ALL(var) (var).begin(), (var).end()

using namespace std;

template <class T, size_t S>
ostream &operator<<(ostream &_ostr, const array<T, S> &_v);
template <class T>
ostream &operator<<(ostream &_ostr, const vector<T> &_v);
template <class T>
ostream &operator<<(ostream &_ostr, const deque<T> &_v);
template <class T>
ostream &operator<<(ostream &_ostr, const list<T> &_v);
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
template <class T>
ostream &operator<<(ostream &_ostr, const set<T> &_v);
template <class T>
ostream &operator<<(ostream &_ostr, const multiset<T> &_v);
template <class T>
ostream &operator<<(ostream &_ostr, const unordered_set<T> &_v);
template <class T>
ostream &operator<<(ostream &_ostr, const unordered_multiset<T> &_v);

template <class T>
void orange(ostream &_ostr, const T &_v) {
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
template <class T>
ostream &operator<<(ostream &_ostr, const vector<T> &_v) {
  orange(_ostr, _v);
  return _ostr;
}
template <class T>
ostream &operator<<(ostream &_ostr, const deque<T> &_v) {
  orange(_ostr, _v);
  return _ostr;
}
template <class T>
ostream &operator<<(ostream &_ostr, const list<T> &_v) {
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
template <class T>
ostream &operator<<(ostream &_ostr, const set<T> &_v) {
  orange(_ostr, _v);
  return _ostr;
}
template <class T>
ostream &operator<<(ostream &_ostr, const multiset<T> &_v) {
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
template <class T>
istream &operator>>(istream &_istr, vector<T> &_v);
template <class T>
istream &operator>>(istream &_istr, deque<T> &_v);
template <class T, class Y>
istream &operator>>(istream &_istr, pair<T, Y> &_v);

template <class T, size_t S>
istream &operator>>(istream &_istr, array<T, S> &_v) {
  for (auto &i : _v)
    _istr >> i;
  return _istr;
}
template <class T>
istream &operator>>(istream &_istr, vector<T> &_v) {
  for (auto &i : _v)
    _istr >> i;
  return _istr;
}
template <class T>
istream &operator>>(istream &_istr, deque<T> &_v) {
  for (auto &i : _v)
    _istr >> i;
  return _istr;
}
template <class T, class Y>
istream &operator>>(istream &_istr, pair<T, Y> &_v) {
  _istr >> _v.first >> _v.second;
  return _istr;
}

template <class T>
bool chmax(T &a, const T &b) {
  if (a < b) {
    a = b;
    return true;
  }
  return false;
}
template <class T>
bool chmin(T &a, const T &b) {
  if (a > b) {
    a = b;
    return true;
  }
  return false;
}

class DirectedGraph {
public:
  struct Edge {
    int64_t to, cost;
  };
  int64_t n;
  vector<vector<Edge>> g;

  DirectedGraph(int64_t _n) : g(_n) { n = _n; }
  void add(int64_t s, int64_t t, int64_t cost) {
    Edge e;
    e.to = t, e.cost = cost;
    g[s].push_back(e);
  }
  // O(V^3)
  pair<vector<vector<int64_t>>, vector<vector<int64_t>>> warshallfloyd() {
    vector<vector<int64_t>> DIST(n, vector<int64_t>(n, LLONG_MAX));
    vector<vector<int64_t>> NEXT(n, vector<int64_t>(n, LLONG_MAX));
    for (int64_t i = 0; i < n; i++)
      DIST[i][i] = 0;
    for (int64_t i = 0; i < n; i++)
      for (int64_t j = 0; j < n; j++)
        NEXT[i][j] = j;
    for (int64_t i = 0; i < n; i++)
      for (Edge &j : g[i])
        chmin(DIST[i][j.to], j.cost);
    for (int64_t k = 0; k < n; k++)
      for (int64_t i = 0; i < n; i++)
        for (int64_t j = 0; j < n; j++)
          if (DIST[i][k] != LLONG_MAX and DIST[k][j] != LLONG_MAX)
            if (chmin(DIST[i][j], DIST[i][k] + DIST[k][j]))
              NEXT[i][j] = NEXT[i][k];
    return {DIST, NEXT};
  }
  // O(E+VlogV)
  vector<int64_t> dijkstra(int64_t s) {
    vector<int64_t> d(n, INT_MAX);
    d[s] = 0;
    priority_queue<pair<int64_t, int64_t>, vector<pair<int64_t, int64_t>>,
                   greater<pair<int64_t, int64_t>>>
        que;
    que.push(make_pair(0, s));
    while (!que.empty()) {
      pair<int64_t, int64_t> p = que.top();
      que.pop();
      int64_t v = p.second;
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
  vector<int64_t> bellmanford(int64_t s) {
    vector<int64_t> d(n, INT_MAX);
    d[s] = 0;
    for (int64_t _ = 0; _ < n; _++) {
      bool upd = false;
      for (int64_t u = 0; u < n; u++)
        if (d[u] < INT_MAX)
          for (const auto &e : g[u]) {
            int64_t v = e.to;
            if (d[v] > d[u] + e.cost)
              d[v] = d[u] + e.cost, upd = true;
          }
      if (!upd)
        return d;
    }
    queue<int64_t> Q;
    for (int64_t u = 0; u < n; u++)
      if (d[u] < INT_MAX)
        Q.emplace(u);
    while (!Q.empty()) {
      int64_t u = Q.front();
      Q.pop();
      for (const auto &e : g[u]) {
        int64_t v = e.to;
        if (d[v] > INT_MIN && (d[u] == INT_MIN || d[v] > d[u] + e.cost))
          d[v] = INT_MIN, Q.emplace(v);
      }
    }
    return d;
  }
  // O(V+E)
  vector<int64_t> topologicalSort() {
    vector<int64_t> d, ind(n);
    for (int64_t i = 0; i < n; i++)
      for (auto e : g[i])
        ind[e.to]++;
    queue<int64_t> que;
    for (int64_t i = 0; i < n; i++)
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
};

void putLogo() {
  cout << "              _   _          \n";
  cout << "   ___  _ __ | |_| |_ _ __   \n";
  cout << "  / _ \\| '_ \\| __| __| '_ \\  \n";
  cout << " | (_) | |_) | |_| |_| |_) | \n";
  cout << "  \\___/| .__/ \\__|\\__| .__/  \n";
  cout << "       |_|           |_|     \n";
  cout << "----------------------------\n";
}

void checkArgs(const int64_t V, const int64_t E, const DirectedGraph &G, const int64_t T, const vector<int64_t> &M, const int64_t N, const vector<pair<int64_t, int64_t>> &Q, const vector<vector<int64_t>> &DIST, const vector<vector<int64_t>> &NEXT) {
  cout << "G = (" << V << ", " << E << ")\n";
  cout << "頂点間距離:\n";
  for (auto &i : DIST)
    cout << i << "\n";
  cout << "全点対最短経路:\n";
  for (auto &i : NEXT)
    cout << i << "\n";
  cout << "最大シミュレート時間T: " << T << "\n";
  cout << "運送者M: " << M << "\n";
  cout << "クエリ数N: " << N << "\n";
  for (auto &i : Q)
    cout << i << "\n";
  cout << "頂点出次: "
       << "\n";
  for (auto &i : G.g) {
    for (auto &j : i) {
      cout << j.to << " ";
    }
    cout << "\n";
  }
  cout << "----------------------------\n";
}

class Carrier {
public:
  int pos;                                // Carrierの頂点番号
  unordered_multiset<int64_t> passengers; // 頂点passengers[i]に行きたい人たち

  Carrier(int _pos) : pos(_pos) {}
};

void byGreedy(const int64_t V, const int64_t E, DirectedGraph &G, const int64_t T, const vector<int64_t> &M, const int64_t N, const vector<pair<int64_t, int64_t>> &Q, const vector<vector<int64_t>> &DIST, const vector<vector<int64_t>> &NEXT) {
  vector<Carrier> carriers; // 運送者たち id:i
  for (size_t i = 0; i < M.size(); i++)
    for (int64_t j = 0; j < M[i]; j++)
      carriers.emplace_back(Carrier(i)); // 運送者情報をcarriersに詰める

  vector<unordered_multiset<int64_t>> m(V); // 配達するもの(頂点m[i]に行きたい)
  for (auto &[x, y] : Q)                    // 頂点xからyに行きたい
    m[x].insert(y);                         // 運送物登録

  set<pair<int64_t, int64_t>> que; // pair<時間, 運送者番号>
  for (size_t i = 0; i < carriers.size(); i++)
    que.insert({0, i}); // queは運送者がいつ空いてるかを管理する setなので時間が先な方が優先

  int t = 0; // 現在シュミレート時間
  cout << "que: " << que << "\n";
  while (not que.empty()) {
    if ((*que.begin()).first == t) {
      int64_t num = que.begin()->second; // 運送者番号 -> num
      que.erase(que.begin());            // queの上から一つとる

      cout << "pos: " << carriers[num].pos << "\n";

      while (carriers[num].passengers.find(carriers[num].pos) != carriers[num].passengers.end()) { // 荷物下ろせるやつは全て下ろす
        cout << "dropped at: " << carriers[num].pos << "\n";
        carriers[num].passengers.erase(carriers[num].pos);
      }

      for (auto &passenger : m[carriers[num].pos]) // 頂点[carriers[num].pos]の荷物を全て移す
        carriers[num].passengers.insert(passenger);
      m[carriers[num].pos].clear();

      /*
      int mode = -1; // modeに遷移する
      {              // 最頻値を取得し、最頻の頂点番号をmodeに
        vector<int64_t> cnts(V, 0);
        for (auto &i : G.g[carriers[num].pos])
          cnts[i.to]++;
        mode = max_element(cnts.begin(), cnts.end()) - cnts.begin();
      }
      assert(mode != -1);
      */

      if(carriers[num].passengers.size() == 0)
        continue;

      int next = -1; // nextに遷移する
      {
        vector<int64_t> cnts(V, 0); // とりあえず頂点iに行きたい荷物の数
        for (auto &i : carriers[num].passengers) {
          cnts[NEXT[carriers[num].pos][i]]++;
        }
        cout<<"cnts: "<<cnts<<"\n";
        next = max_element(cnts.begin(), cnts.end()) - cnts.begin();
      }
      cout << "next: " << next << "\n";
      que.insert({t+DIST[carriers[num].pos][next],num});
      carriers[num].pos = next;

      //どちらに進むのか決定する(TODO: グラフが一整列している場合のみ)
      //出次数は高々二つなので、いちいち探索を行わない
      /*
      int back = 0;  //後ろに下がると良い乗客の数
      int front = 0; //前に進む良い乗客の数

      cout << carriers[num].passengers << "\n";

      for (auto &passenger : carriers[num].passengers) {
        (carriers[num].pos < passenger ? front : back)++;
      }
      cout << back << " " << front << "\n";
      if (back < front) {
        que.insert({t + DIST[carriers[num].pos][carriers[num].pos + 1], num});
        carriers[num].pos++;
      } else if (back > front) {
        que.insert({t + DIST[carriers[num].pos][carriers[num].pos - 1], num});
        carriers[num].pos--;
      } else {
        cout << "dotti ittara iika wakaranaiyo~~\n";
        exit(0);
      }
      cout << flush;
      cout << carriers[num].pos << "\n";
      */

    } else {
      t = (*que.begin()).first; // 時刻tを追従させる
    }
  }
}

int main() {
  int64_t V; // 頂点数
  cin >> V;
  int64_t E; // 辺の数
  cin >> E;
  DirectedGraph G(V); // 有向重み付きグラフ
  for (int64_t i = 0; i < E; i++) {
    int64_t a, b, d;
    cin >> a >> b >> d;
    G.add(a, b, d);
  }
  int64_t T; // 最大シミュレート時間
  cin >> T;
  vector<int64_t> M(V); // 頂点ごとの運送者数
  cin >> M;
  int64_t N; // クエリ数
  cin >> N;
  vector<pair<int64_t, int64_t>> Q(N);
  cin >> Q;
  const auto [DIST, NEXT] = G.warshallfloyd();

  putLogo();
  checkArgs(V, E, G, T, M, N, Q, DIST, NEXT);
  byGreedy(V, E, G, T, M, N, Q, DIST,NEXT);
}
