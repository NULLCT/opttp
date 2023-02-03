#include "simulated_annealing.hpp"

double SA::frand() {
  return ((double)rand() / (RAND_MAX));
}

int64_t SA::evalScore(const STATE &state) {
  int64_t res = 0;                          // 評価値
  std::set<std::tuple<int, int, int>> que;  // 時刻 番号 頂点
  for (int i = 0; i < state.acthist.size(); i++) {
    int time = 0;
    for (int j = 0; j < state.acthist[i].size(); j++) {
      que.insert({time, i, state.acthist[i][j]});
      if (j + 1 != state.acthist[i].size())
        time += model.DIST[state.acthist[i][j]][state.acthist[i][j + 1]];
    }
  }

  std::vector<std::multiset<int64_t>> m(model.V);  // m[x] = y -> 頂点xからyに行きたい
  for (auto &[x, y] : model.Q)                     // 頂点xからyに行きたい
    m[x].insert(y);                                // 運送物登録

  std::vector<std::set<int64_t>> carring(state.acthist.size());  // 各運送者が何を保持しているか

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

void SA::modify(STATE &state) {
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

STATE SA::simulated_annealing() {
  int64_t initscore = evalScore(state);
  std::cerr << "init: " << initscore << "\n";
  int nonexped = 0;
  int exped = 0;
  while (temp > 1) {  // 十分冷えるまで
    for (int i = 0; i < repnum; i++) {
      STATE newstate = state;
      modify(newstate);
      int64_t statescore = evalScore(state);
      int64_t newstatescore = evalScore(newstate);
      int64_t delta = newstatescore - statescore;

      if (0 < delta) {
        state = newstate;
        nonexped++;
      } else if (frand() < std::exp(delta / temp)) {
        state = newstate;
        exped++;
      }

      if (beststate_score < newstatescore) {
        beststate_score = newstatescore;
        beststate = newstate;
      }
    }
    temp *= coolingcoef;
    std::cerr << "\e[2K\e[0E" << std::left << std::setw(8) << temp << " ("
              << "exp/nexp = " << exped << "/" << nonexped << ") " << beststate_score;
    nonexped = exped = 0;
  }
  std::cout << "\n"
            << initscore << "," << beststate_score << "\n";
  return beststate;
}