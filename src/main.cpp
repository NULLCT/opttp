#include <algorithm>
#include <climits>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <unordered_set>
#include <vector>

#include "gen.hpp"
#include "util.hpp"

using namespace std;

STATE init(const MODEL &model) {
  vector<Carrier> carriers;  // 運送者たち id:i
  for (size_t i = 0; i < model.M.size(); i++)
    for (int64_t j = 0; j < model.M[i]; j++)
      carriers.emplace_back(Carrier(i));  // 運送者情報をcarriersに詰める

  STATE res;  // SA法用状態保持
  res.acthist.resize(carriers.size());
  res.pathlen.resize(carriers.size());

  for (int i = 0; i < carriers.size(); i++) {
    res.acthist[i].push_back(carriers[i].pos);
    while (true) {
      int nextpos = rand() % model.G.g[res.acthist[i].back()].size();
      if (res.pathlen[i] + model.G.g[res.acthist[i].back()][nextpos].cost <= model.T) {
        res.pathlen[i] += model.G.g[res.acthist[i].back()][nextpos].cost;
        res.acthist[i].push_back(model.G.g[res.acthist[i].back()][nextpos].to);
      } else {
        break;
      }
    }
  }

  return res;
}

STATE byGreedy(const MODEL &model) {
  vector<Carrier> carriers;  // 運送者たち id:i
  for (size_t i = 0; i < model.M.size(); i++)
    for (int64_t j = 0; j < model.M[i]; j++)
      carriers.emplace_back(Carrier(i));  // 運送者情報をcarriersに詰める

  STATE res;  // SA法用状態保持
  res.acthist.resize(carriers.size());
  res.pathlen.resize(carriers.size());

  vector<unordered_multiset<int64_t>> m(model.V);  // 配達するもの(頂点m[i]に行きたい)
  for (auto &[x, y] : model.Q)                     // 頂点xからyに行きたい
    m[x].insert(y);                                // 運送物登録

  set<pair<int64_t, int64_t>> que;  // pair<時間, 運送者番号>
  for (size_t i = 0; i < carriers.size(); i++)
    que.insert({0, i});  // queは運送者がいつ空いてるかを管理する setなので時間が先な方が優先

  int t = 0;  // 現在シミュレート時間
  while (not que.empty()) {
    if ((*que.begin()).first == t) {
      int64_t num = que.begin()->second;  // 運送者番号 -> num
      que.erase(que.begin());             // queの上から一つとる
      res.acthist[num].push_back(carriers[num].pos);

      while (carriers[num].passengers.find(carriers[num].pos) != carriers[num].passengers.end()) {  // 荷物下ろせるやつは全て下ろす
        carriers[num].passengers.erase(carriers[num].pos);
      }

      for (auto &passenger : m[carriers[num].pos])  // 頂点[carriers[num].pos]の荷物を全て移す
        carriers[num].passengers.insert(passenger);
      m[carriers[num].pos].clear();

      if (carriers[num].passengers.size() == 0) {  // 荷物を載せていなかった場合
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
      } else {          // 荷物を載せていた場合
        int next = -1;  // nextに遷移する
        {
          vector<int64_t> cnts(model.V, 0);  // とりあえず頂点iに行きたい荷物の数
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
      t = (*que.begin()).first;  // 時刻tを追従させる
    }
  }

  for (int i = 0; i < res.acthist.size(); i++) {
    int sum = 0;
    for (int j = 0; j < res.acthist[i].size() - 1; j++)
      sum += model.DIST[res.acthist[i][j]][res.acthist[i][j + 1]];
    res.pathlen[i] = sum;
  }

  return res;
}

// 焼きなまし法
class SA {
  STATE state;               // 現在の状態
  double temp;               // 温度
  const int repnum;          // 反復回数
  const double coolingcoef;  // 冷却係数
  MODEL model;               // モデル
  STATE beststate;           // 暫定最適状態
  int64_t beststate_score;   // 暫定最適状態ansを評価関数に通したスコア

  // [0,1]の乱数を返す
  double frand() {
    return ((double)rand() / (RAND_MAX));
  }

  // 評価関数
  int64_t evalScore(const STATE &state) {
    int64_t res = 0;                // 評価値
    set<tuple<int, int, int>> que;  // 時刻 番号 頂点
    for (int i = 0; i < state.acthist.size(); i++) {
      int time = 0;
      for (int j = 0; j < state.acthist[i].size(); j++) {
        que.insert({time, i, state.acthist[i][j]});
        if (j + 1 != state.acthist[i].size())
          time += model.DIST[state.acthist[i][j]][state.acthist[i][j + 1]];
      }
    }

    vector<multiset<int64_t>> m(model.V);  // m[x] = y -> 頂点xからyに行きたい
    for (auto &[x, y] : model.Q)           // 頂点xからyに行きたい
      m[x].insert(y);                      // 運送物登録

    vector<set<int64_t>> carring(state.acthist.size());  // 各運送者が何を保持しているか

    for (auto &q : que) {
      auto &[t, n, v] = q;  // 時刻 番号 頂点
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
    int carriernum = (rand() % state.acthist.size());                              // 何番目の運送者か
    int posnum = (rand() % state.acthist[carriernum].size());                      // 運送経路履歴の何番目頂点を選ぶか
    int outdnum = (rand() % model.G.g[state.acthist[carriernum][posnum]].size());  // 何本目の辺を選ぶか
    int pos = state.acthist[carriernum][posnum];                                   // 運送者がいる頂点番号
    int nextpos = model.G.g[pos][outdnum].to;                                      // これから行き来する頂点

    state.acthist[carriernum].insert(state.acthist[carriernum].begin() + posnum + 1, nextpos);
    state.acthist[carriernum].insert(state.acthist[carriernum].begin() + posnum + 2, pos);

    state.pathlen[carriernum] += model.DIST[pos][nextpos] + model.DIST[nextpos][pos];

    while (model.T < state.pathlen[carriernum]) {
      state.pathlen[carriernum] -= model.DIST[state.acthist[carriernum].end()[-1]][state.acthist[carriernum].end()[-2]];
      state.acthist[carriernum].pop_back();
    }
  }

 public:
  SA(STATE &_state, double _temp, int _repnum, double _coolingcoef, MODEL &_model) : state(_state),              // 状態
                                                                                     temp(_temp),                // 温度
                                                                                     repnum(_repnum),            // 試行回数
                                                                                     coolingcoef(_coolingcoef),  // 冷却係数
                                                                                     model(_model),              // 運送モデル
                                                                                     beststate(_state),
                                                                                     beststate_score(evalScore(_state)) {}

  // 焼きなまし法
  STATE simulated_annealing() {
    int64_t initscore = evalScore(state);
    int nonexped = 0;
    int exped = 0;
    while (temp > 1.0) {  // 十分冷えるまで
      for (int i = 0; i < repnum; i++) {
        STATE newstate = state;
        modify(newstate);
        int64_t statescore = evalScore(state);
        int64_t newstatescore = evalScore(newstate);
        int64_t delta = newstatescore - statescore;

        if (0 < delta) {
          state = newstate;
          nonexped++;
        } else if (frand() < exp(-delta / temp)) {
          state = newstate;
          exped++;
        }

        if (beststate_score < newstatescore) {
          beststate_score = newstatescore;
          beststate = newstate;
        }
      }
      temp *= coolingcoef;
      cerr << temp << " ("
           << "exp/nexp = " << exped << "/" << nonexped << ") "<<beststate_score<<"\n";
      nonexped = exped = 0;
    }
    cout << initscore << "," << beststate_score << "\n";
    return beststate;
  }
};

int main() {
  MODEL model;
  cerr << "generating model... " << flush;
  generateTestCase(model);
  cerr << "done" << endl;

  STATE state;
  cerr << "initializing model... " << flush;
  state = init(model);
  cerr << "done" << endl;

  cerr<<"=====DIRECTED WEIGHTED GRAPH CODE BEGIN=====\n";
  cerr<<model.V<<" "<<model.E<<"\n";
  for(int i=0;i<model.V;i++){
    for(int j=0;j<model.G.g[i].size();j++){
      cerr<<i<<" "<<model.G.g[i][j].to<<" "<<model.G.g[i][j].cost<<" ";
    }
  }
  cerr<<"\n";
  cerr<<"=====DIRECTED WEIGHTED GRAPH CODE END=====\n";
  cerr<<"=====QUERY BEGIN=====\n";
  for(auto &i:model.Q){
    cerr<<i.first<<" -> "<<i.second<<"\n";
  }
  cerr<<"=====QUERY END=====\n";

  SA sa(state, 100, 1000, 0.99, model);

  cout << fixed << setprecision(32);
  auto ans = sa.simulated_annealing();
  for (int i = 0; i < ans.pathlen.size(); i++) {
    cerr << "carrier#" << i << ":\n";
    cerr << "pathlen: " << ans.pathlen[i] << "\n";
    cerr << "acthist: " << ans.acthist[i] << "\n";
  }
}
