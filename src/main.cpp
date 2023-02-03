#include <algorithm>
#include <climits>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <unordered_set>
#include <vector>

#include "experimental_model.hpp"
#include "gen.hpp"
#include "simulated_annealing.hpp"

using namespace std;

class Carrier {
 public:
  int pos;                                 // Carrierの頂点番号
  unordered_multiset<int64_t> passengers;  // 頂点passengers[i]に行きたい人たち

  Carrier(int _pos) : pos(_pos) {}
};

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

int main() {
  MODEL model;
  cerr << "generating model... " << flush;
  generateTestCase(model);
  cerr << "done" << endl;

  STATE state;
  cerr << "initializing model... " << flush;
  state = init(model);
  cerr << "done" << endl;

  cerr << "=====DIRECTED WEIGHTED GRAPH CODE BEGIN=====\n";
  cerr << model.V << " " << model.E << "\n";
  for (int i = 0; i < model.V; i++) {
    for (int j = 0; j < model.G.g[i].size(); j++) {
      cerr << i << " " << model.G.g[i][j].to << " " << model.G.g[i][j].cost << " ";
    }
  }
  cerr << "\n";
  cerr << "=====DIRECTED WEIGHTED GRAPH CODE END=====\n";
  cerr << "=====QUERY BEGIN(size: " << model.Q.size() << ")=====\n";
  for (auto &i : model.Q)
    cout << i.first << " -> " << i.second << "\n";
  cerr << "=====QUERY END=====\n";
  cerr << "=====CARRIER BEGIN(size: " << model.M.size() << ")=====\n";
  for (auto &i : model.M)
    cout << i << " ";
  cerr << "\n=====QUERY END=====\n";
  cerr << "Limit: " << model.T << "\n";

  SA sa(state, 100, 100, 0.99, model);

  cout << fixed << setprecision(32);
  auto ans = sa.simulated_annealing();
  for (int i = 0; i < ans.pathlen.size(); i++) {
    cerr << "carrier#" << i << ":\n";
    cerr << "pathlen: " << ans.pathlen[i] << "\n";
    cerr << "acthist: ";
    for(auto &act:ans.acthist[i])
      cout<<act<<" ";
    cout<<"\n";
  }
}
