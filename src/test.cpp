#include <sstream>
#include <stdexcept>

#include "acutest.h"
#include "directed_graph.hpp"
#include "experimental_model.hpp"
#include "simulated_annealing.hpp"
#include "state.hpp"
#include "unionfind.hpp"

void test_UnionFind() {
  for (int i = 10; i < 1000; i++) {
    UnionFind uf(i);
    TEST_CHECK(uf.groups().size() == i);
    TEST_CHECK(uf.merge(0, 1));
    TEST_CHECK(uf.groups().size() == i - 1);
    TEST_CHECK(uf.merge(1, 2));
    TEST_CHECK(uf.groups().size() == i - 2);
    TEST_CHECK(not uf.merge(0, 2));
    TEST_CHECK(uf.isSame(0, 2));
    TEST_CHECK(uf.size(1) == 3);
  }
}

void test_DirectedGraph() {
  for (int i = 10; i < 100; i++) {
    DirectedGraph g(i);
    g.add(0, 1, 10);
    g.add(1, 2, 100);
    g.add(2, 3, 100);
    auto dij = g.dijkstra(0);
    TEST_CHECK(dij[1] == 10);
    TEST_CHECK(dij[2] == 110);
    TEST_CHECK(dij[3] == 210);
    auto war = g.warshallfloyd();
    TEST_CHECK(war.dist[0][1] == 10);
    TEST_CHECK(war.dist[0][2] == 110);
    TEST_CHECK(war.dist[0][3] == 210);
    auto bel = g.bellmanford_SPFA(0, 3);
    TEST_CHECK(not bel.hascycle);
    TEST_CHECK(bel.distances[1] == 10);
    TEST_CHECK(bel.distances[2] == 110);
    TEST_CHECK(bel.distances[3] == 210);
    auto top = g.topologicalSort();
  }
}

void test_SimulatedAnnealing() {
  MODEL model;
  std::string s =
      "5 8\n\
0 1 1\n\
1 2 1\n\
2 3 1\n\
3 4 1\n\
4 3 1\n\
3 2 1\n\
2 1 1\n\
1 0 1\n\
100\n\
1 0 0 0 0\n\
1\n\
0 4";
  std::stringbuf strbuf(s.c_str());
  std::istream ist(&strbuf);
  ist >> model;
  STATE state = model.generateState();
  SA sa(state,100,100,0.99,model);
  auto res = sa.simulated_annealing();
  TEST_CHECK(sa.evalScore(res) == 96);
  for(int i=0;i<res.acthist[0].size()-1;i++)
    TEST_CHECK(abs(res.acthist[0][i] - res.acthist[0][i+1]) == 1);
}

TEST_LIST = {
    {"unionfind", test_UnionFind},
    {"directedgraph", test_DirectedGraph},
    {"simulatedannealing", test_SimulatedAnnealing},
    {nullptr, nullptr}};
