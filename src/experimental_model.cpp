#include "experimental_model.hpp"

STATE MODEL::generateState() {
  std::vector<Carrier> carriers;  // 運送者たち id:i
  for (size_t i = 0; i < M.size(); i++)
    for (int64_t j = 0; j < M[i]; j++)
      carriers.emplace_back(Carrier(i));  // 運送者情報をcarriersに詰める

  STATE res;  // SA法用状態保持
  res.acthist.resize(carriers.size());
  res.pathlen.resize(carriers.size());

  for (int i = 0; i < carriers.size(); i++) {
    res.acthist[i].push_back(carriers[i].pos);
    while (true) {
      int nextpos = rand() % G.g[res.acthist[i].back()].size();
      if (res.pathlen[i] + G.g[res.acthist[i].back()][nextpos].cost <= T) {
        res.pathlen[i] += G.g[res.acthist[i].back()][nextpos].cost;
        res.acthist[i].push_back(G.g[res.acthist[i].back()][nextpos].to);
      } else {
        break;
      }
    }
  }
  return res;
}

std::istream &operator>>(std::istream &is, MODEL &model) {
  is >> model.V >> model.E;
  model.G.init(model.V);
  model.M.resize(model.V);
  for (int i = 0; i < model.E; i++) {
    int64_t s, t, d;
    is >> s >> t >> d;
    model.G.add(s, t, d);
  }
  is >> model.T;
  for (auto &i : model.M) {
    is >> i;
  }
  int64_t q;
  is >> q;
  model.Q.resize(q);
  for (auto &i : model.Q) {
    is >> i.first >> i.second;
  }
  auto res = model.G.warshallfloyd();
  model.DIST = res.dist;
  model.NEXT = res.next;
  return is;
}