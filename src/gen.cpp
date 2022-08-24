#include <algorithm>
#include <ctime>
#include <iostream>
#include <list>
#include <random>
#include <set>
#include <tuple>
#include <vector>
using namespace std;

constexpr int V = 100;
constexpr int T = 1000;
constexpr int M = 10;
constexpr int Q = 50;

class UnionFind {
public:
  int n;
  vector<int> par;
  UnionFind(int _n) : n(_n), par(_n, -1) {}
  // O(α(n))
  bool merge(int a, int b) {
    a = root(a), b = root(b);
    if (a == b)
      return false;
    if (par[a] > par[b])
      swap(a, b);
    par[a] += par[b];
    par[b] = a;
    return true;
  }
  // O(α(n))
  bool isSame(int a, int b) {
    return root(a) == root(b);
  }
  int root(int a) {
    if (par[a] < 0)
      return a;
    return par[a] = root(par[a]);
  }
  int size(int a) {
    return -par[root(a)];
  }
  vector<vector<int>> groups() {
    vector<int> leader_buf(n), group_size(n);
    for (int i = 0; i < n; i++) {
      leader_buf[i] = root(i);
      group_size[leader_buf[i]]++;
    }
    vector<vector<int>> result(n);
    for (int i = 0; i < n; i++)
      result[i].reserve(group_size[i]);
    for (int i = 0; i < n; i++)
      result[leader_buf[i]].push_back(i);
    result.erase(remove_if(result.begin(), result.end(), [&](const vector<int> &v) { return v.empty(); }), result.end());
    return result;
  }
};

int main() {
  srand(time(nullptr));
  vector<vector<pair<int, int>>> g(V);
  UnionFind uf(V);

  while (true) {
  BG:
    int a = rand() % V;
    int b = rand() % V;
    int l = (rand() % 10) + 1;
    while (a == b)
      b = rand() % V;

    if (uf.isSame(a, b)) {
      if (uf.groups().size() == 1)
        break;
      goto BG;
    } else {
      uf.merge(a, b);
      g[a].push_back({b, l});
      g[b].push_back({a, l});
    }
  }

  cout << V << " " << [&g]() {
    int sum = 0;
    for (auto &i : g)
      sum += i.size();
    return sum;
  }() << "\n";

  for (int i = 0; i < g.size(); i++)
    for (int j = 0; j < g[i].size(); j++)
      cout << i << " " << g[i][j].first << " " << g[i][j].second << "\n";

  cout << T << "\n";
  {
    list<int> ls;
    for (int i = 0; i < V; i++)
      ls.push_back(i);

    while (ls.size() > M)
      ls.erase(next(ls.begin(), rand() % ls.size()));

    for (int i = 0; i < V; i++)
      cout << (find(ls.begin(), ls.end(), i) == ls.end() ? 0 : 1) << " ";
    cout << "\n";
  }

  cout << Q << "\n";
  {
    set<pair<int, int>> st;
    for (int i = 0; i < Q; i++) {
      pair<int, int> m;
      do {
        m.first = rand() % V;
        m.second = rand() % V;
      } while (st.find(m) != st.end() and m.first != m.second);
      st.insert(m);
    }

    for(auto &i:st){
      cout<<i.first<<" "<<i.second<<"\n";
    }
  }
}
