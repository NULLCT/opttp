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

struct MODEL {                      // 運送モデル
  int64_t V;                        // 頂点数
  int64_t E;                        // 辺の数
  DirectedGraph G;                  // グラフ
  int64_t T;                        // 最大シミュレート時間
  vector<int64_t> M;                // 頂点ごとの運送数
  int64_t N;                        // クエリ数
  vector<pair<int64_t, int64_t>> Q; // クエリ
  vector<vector<int64_t>> DIST;     // 頂点間距離
  vector<vector<int64_t>> NEXT;     // ある頂点に行くため次どこに行けばいいか
};

void checkArgs(const MODEL &model) {
  cout << "G = (" << model.V << ", " << model.E << ")\n";
  cout << "頂点間距離(dist):\n";
  for (auto &i : model.DIST)
    cout << i << "\n";
  cout << "全点対最短経路(index):\n";
  for (auto &i : model.NEXT)
    cout << i << "\n";
  cout << "最大シミュレート時間T: " << model.T << "\n";
  cout << "運送者M: " << model.M << "\n";
  cout << "クエリ数N: " << model.N << "\n";
  for (auto &i : model.Q)
    cout << i << "\n";
  cout << "頂点出次: "
       << "\n";
  for (auto &i : model.G.g) {
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

// https://jetbead.hatenablog.com/entry/20120623/1340419446
struct STATE {           // 状態構造体
  vector<vector<int>> x; // 頂点移動履歴
  vector<int> l;         // パス長
};

STATE fake(const MODEL &model) {
  vector<Carrier> carriers; // 運送者たち id:i
  for (size_t i = 0; i < model.M.size(); i++)
    for (int64_t j = 0; j < model.M[i]; j++)
      carriers.emplace_back(Carrier(i)); // 運送者情報をcarriersに詰める

  STATE res; // SA法用状態保持
  res.x.resize(carriers.size());
  res.l.resize(carriers.size());

  for (int i = 0; i < carriers.size(); i++)
    res.x[i].push_back(carriers[i].pos);

  for (int i = 0; i < res.x.size(); i++) {
    int sum = 0;
    for (int j = 0; j < res.x[i].size() - 1; j++)
      sum += model.DIST[res.x[i][j]][res.x[i][j + 1]];
    res.l[i] = sum;
  }

  return res;
}

STATE byGreedy(const MODEL &model) {
  vector<Carrier> carriers; // 運送者たち id:i
  for (size_t i = 0; i < model.M.size(); i++)
    for (int64_t j = 0; j < model.M[i]; j++)
      carriers.emplace_back(Carrier(i)); // 運送者情報をcarriersに詰める

  STATE res; // SA法用状態保持
  res.x.resize(carriers.size());
  res.l.resize(carriers.size());

  vector<unordered_multiset<int64_t>> m(model.V); // 配達するもの(頂点m[i]に行きたい)
  for (auto &[x, y] : model.Q)                    // 頂点xからyに行きたい
    m[x].insert(y);                               // 運送物登録

  set<pair<int64_t, int64_t>> que; // pair<時間, 運送者番号>
  for (size_t i = 0; i < carriers.size(); i++)
    que.insert({0, i}); // queは運送者がいつ空いてるかを管理する setなので時間が先な方が優先

  int t = 0; // 現在シミュレート時間
  while (not que.empty()) {
    if ((*que.begin()).first == t) {
      int64_t num = que.begin()->second; // 運送者番号 -> num
      que.erase(que.begin());            // queの上から一つとる
      res.x[num].push_back(carriers[num].pos);

      while (carriers[num].passengers.find(carriers[num].pos) != carriers[num].passengers.end()) { // 荷物下ろせるやつは全て下ろす
        carriers[num].passengers.erase(carriers[num].pos);
      }

      for (auto &passenger : m[carriers[num].pos]) // 頂点[carriers[num].pos]の荷物を全て移す
        carriers[num].passengers.insert(passenger);
      m[carriers[num].pos].clear();

      if (carriers[num].passengers.size() == 0) { // 荷物を載せていなかった場合
        int next = 0;
        for (int i = 0; i < model.V; i++)
          if (m[next].size() < m[i].size())
            next = i;

        if (m[next].size() == 0) {
          continue;
        }
        next = model.NEXT[carriers[num].pos][next];
        if (model.DIST[carriers[num].pos][next] == LLONG_MAX) {
          continue;
        }
        que.insert({t + model.DIST[carriers[num].pos][next], num});
        carriers[num].pos = next;
      } else {         // 荷物を載せていた場合
        int next = -1; // nextに遷移する
        {
          vector<int64_t> cnts(model.V, 0); // とりあえず頂点iに行きたい荷物の数
          for (auto &i : carriers[num].passengers)
            if (model.DIST[carriers[num].pos][i] != LLONG_MAX)
              cnts[model.NEXT[carriers[num].pos][i]]++;
          next = max_element(cnts.begin(), cnts.end()) - cnts.begin();
        }
        if (model.DIST[carriers[num].pos][next] == LLONG_MAX) {
          continue;
        }
        que.insert({t + model.DIST[carriers[num].pos][next], num});
        carriers[num].pos = next;
      }
    } else {
      t = (*que.begin()).first; // 時刻tを追従させる
    }
  }

  for (int i = 0; i < res.x.size(); i++) {
    int sum = 0;
    for (int j = 0; j < res.x[i].size() - 1; j++)
      sum += model.DIST[res.x[i][j]][res.x[i][j + 1]];
    res.l[i] = sum;
  }

  return res;
}

// 焼きなまし法
class SA {
  STATE state;              // 現在の状態
  double temp;              // 温度
  const int repnum;         // 反復回数
  const double coolingcoef; // 冷却係数
  MODEL model;              // モデル
  STATE beststate;          // 暫定最適状態
  double beststate_score;   // 暫定最適状態ansを評価関数に通したスコア

  // [0,1]の乱数を返す
  double frand() {
    return ((double)rand() / (RAND_MAX));
  }

  // 評価関数
  double evalScore(const STATE &state) {
    double res = 0;
    set<tuple<int, int, int>> que; // 時刻 番号 頂点
    for (int i = 0; i < state.x.size(); i++) {
      int time = 0;
      for (int j = 0; j < state.x[i].size(); j++) {
        que.insert({time, i, state.x[i][j]});
        if (j + 1 != state.x[i].size())
          time += model.DIST[state.x[i][j]][state.x[i][j + 1]];
      }
    }

    vector<set<int64_t>> m(model.V); // 配達するもの(頂点m[i]に行きたい)
    for (auto &[x, y] : model.Q)     // 頂点xからyに行きたい
      m[x].insert(y);                // 運送物登録

    vector<set<int64_t>> carring(state.x.size()); // 各運送者が何を保持しているか

    for (auto &q : que) {
      auto &[t, n, v] = q; // 時刻 番号 頂点
      for (auto &i : m[v])
        carring[n].insert(i);
      m[v].clear();

      while (carring[n].find(v) != carring[n].end()) {
        res += model.T - t;
        carring[n].erase(v);
      }
    }
    return res;
  }

  // s-opt法 近傍に行き帰りするパスを追加
  void modify(STATE &state) {
    int carriernum = (rand() % state.x.size());                             // 何番目の運送者か
    int posnum = (rand() % state.x[carriernum].size());                     // 運送経路履歴の何番目頂点を選ぶか
    int outdnum = (rand() % model.G.g[state.x[carriernum][posnum]].size()); // 何本目の辺を選ぶか
    int pos = state.x[carriernum][posnum];                                  // 運送者がいる頂点番号

    state.x[carriernum].insert(state.x[carriernum].begin() + posnum + 1, model.G.g[pos][outdnum].to);
    state.x[carriernum].insert(state.x[carriernum].begin() + posnum + 2, pos);

    state.l[carriernum] += model.G.g[pos][outdnum].cost * 2; // TODO: 偏重の場合対応

    while (model.T < state.l[carriernum]) {
      state.l[carriernum] -= model.DIST[state.x[carriernum].end()[-1]][state.x[carriernum].end()[-2]];
      state.x[carriernum].pop_back();
    }
  }

  // 温度の更新
  double next_T(double t) {
    return t * coolingcoef;
  }

public:
  SA(STATE &_state, double _temp, int _repnum, double _coolingcoef, MODEL &_model) : state(_state),             // 状態
                                                                                     temp(_temp),               // 温度
                                                                                     repnum(_repnum),           // 試行回数
                                                                                     coolingcoef(_coolingcoef), // 冷却係数
                                                                                     model(_model),             // 運送モデル
                                                                                     beststate(_state),
                                                                                     beststate_score(evalScore(_state)) {}

  // 焼きなまし法
  STATE simulated_annealing() {
    double initscore = evalScore(state);
    while (temp > 1.0) { // 十分冷えるまで
      for (int i = 0; i < repnum; i++) {
        STATE newstate = state;
        modify(newstate);
        double statescore = evalScore(state);
        double newstatescore = evalScore(newstate);
        double delta = newstatescore - statescore;

        if (0 < delta)
          state = newstate;
        else if (frand() < exp(delta / temp))
          state = newstate;

        if (beststate_score < newstatescore) {
          beststate_score = newstatescore;
          beststate = newstate;
        }
      }
      temp = next_T(temp);
      cout << temp << "\n";
    }
    cout << initscore << "," << beststate_score << "\n";
    return beststate;
  }
};

int main() {
  MODEL model;
  cin >> model.V;
  cin >> model.E;
  model.G.n = model.V;
  model.G.g.resize(model.V);
  for (int64_t i = 0; i < model.E; i++) {
    int64_t a, b, d;
    cin >> a >> b >> d;
    model.G.add(a, b, d);
  }
  cin >> model.T;
  model.M.resize(model.V);
  cin >> model.M;
  cin >> model.N;
  model.Q.resize(model.N);
  cin >> model.Q;
  model.DIST = model.G.warshallfloyd().first;
  model.NEXT = model.G.warshallfloyd().second;

  STATE res;
  res = fake(model);

  STATE state = res;
  SA sa(state, 100000, 10, 0.9, model);

  cout << fixed << setprecision(32);
  auto ans = sa.simulated_annealing();
  cout << ans.l << ": " << ans.x << "\n";
}
