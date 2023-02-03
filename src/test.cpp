#include <stdexcept>

#include "acutest.h"
#include "directed_graph.hpp"
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

TEST_LIST = {
  {"unionfind test", test_UnionFind},
  {"directedgraph test", test_DirectedGraph},
  {nullptr, nullptr}
};