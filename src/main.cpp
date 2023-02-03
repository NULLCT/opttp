#include <algorithm>
#include <climits>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <unordered_set>
#include <vector>

#include "carrier.hpp"
#include "experimental_model.hpp"
#include "generate_testcase.hpp"
#include "simulated_annealing.hpp"
#include "state.hpp"

using namespace std;

STATE init(const MODEL &model) {
  std::vector<Carrier> carriers;  // 運送者たち id:i
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
    for (auto &act : ans.acthist[i])
      cout << act << " ";
    cout << "\n";
  }
}
