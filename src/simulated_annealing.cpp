#include "simulated_annealing.hpp"

double SA::frand() {
  return ((double)rand() / (RAND_MAX));
}

int64_t SA::evalScore(const STATE &_state) {
  int64_t res = 0;                          // 評価値
  std::set<std::tuple<int, int, int>> que;  // 時刻 番号 頂点
  for (int i = 0; i < _state.acthist.size(); i++) {
    int time = 0;
    for (int j = 0; j < _state.acthist[i].size(); j++) {
      que.insert({time, i, _state.acthist[i][j]});
      if (j + 1 != _state.acthist[i].size())
        time += model.DIST[_state.acthist[i][j]][_state.acthist[i][j + 1]];
    }
  }

  std::vector<std::multiset<int64_t>> m(model.V);  // m[x] = y -> 頂点xからyに行きたい
  for (auto &[x, y] : model.Q)                     // 頂点xからyに行きたい
    m[x].insert(y);                                // 運送物登録

  std::vector<std::set<int64_t>> carring(_state.acthist.size());  // 各運送者が何を保持しているか

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

int64_t evalScoreFast(const STATE &_state_before_adaptation, const STATE &_state_after_adaptation) {
}

void SA::modify(STATE &_state) {
  int carriernum = (rand() % _state.acthist.size());                              // 何番目の運送者か
  int posnum = (rand() % _state.acthist[carriernum].size());                      // 運送経路履歴の何番目頂点を選ぶか
  int outdnum = (rand() % model.G.g[_state.acthist[carriernum][posnum]].size());  // 何本目の辺を選ぶか
  int pos = _state.acthist[carriernum][posnum];                                   // 運送者がいる頂点番号
  int nextpos = model.G.g[pos][outdnum].to;                                       // これから行き来する頂点

  _state.acthist[carriernum].insert(_state.acthist[carriernum].begin() + posnum + 1, nextpos);
  _state.acthist[carriernum].insert(_state.acthist[carriernum].begin() + posnum + 2, pos);

  _state.pathlen[carriernum] += model.DIST[pos][nextpos] + model.DIST[nextpos][pos];

  while (model.T < _state.pathlen[carriernum]) {
    _state.pathlen[carriernum] -= model.DIST[_state.acthist[carriernum].end()[-1]][_state.acthist[carriernum].end()[-2]];
    _state.acthist[carriernum].pop_back();
  }
}

STATE SA::simulated_annealing() {
  int64_t initscore = evalScore(state);
  std::cerr << "init: " << initscore << "\n";
  int nonexped = 0;
  int exped = 0;
  while (temperature > 0.1) {
    for (int i = 0; i < repnum; i++) {
      STATE newstate = state;
      modify(newstate);
      int64_t statescore = evalScore(state);
      int64_t newstatescore = evalScore(newstate);
      int64_t delta = newstatescore - statescore;

      if (0 < delta) {
        state = newstate;
        nonexped++;
      } else if (frand() < std::exp(delta / temperature)) {
        state = newstate;
        exped++;
      }

      if (beststate_score < newstatescore) {
        beststate_score = newstatescore;
        beststate = newstate;
      }
    }
    temperature *= coolingcoef;
    std::cerr << "\e[2K\e[0E" << std::left << std::setw(8) << temperature << " ("
              << "exp/nexp = " << exped << "/" << nonexped << ") " << beststate_score;
    nonexped = exped = 0;
  }
  std::cout << "\n"
            << initscore << "," << beststate_score << "\n";
  return beststate;
}
