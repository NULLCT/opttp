#include<iostream>
#include<vector>
#include<queue>
using namespace std;
class DirectedGraph {
  public:
    struct Edge {
      int to, cost;
    };
    int n;
    vector<vector<Edge>> g;

    DirectedGraph(int _n) : g(_n) {
      n = _n;
    }
    void add(int s, int t, int cost) {
      Edge e;
      e.to = t, e.cost = cost;
      g[s].push_back(e);
    }
    // O(V^3)
    vector<vector<int>> warshallfloyd() {
      vector<vector<int>> d(n, vector<int>(n, INT_MAX));
      for (int i = 0; i < n; i++)
        d[i][i] = 0;
      for (int i = 0; i < n; i++)
        for (Edge &j : g[i])
          d[i][j.to] = min(d[i][j.to], j.cost);
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          for (int k = 0; k < n; k++)
            d[j][k] = min(d[j][k], d[j][i] + d[i][k]);
      return d;
    }
    // O(E+VlogV)
    vector<int> dijkstra(int s) {
      vector<int> d(n, INT_MAX);
      d[s] = 0;
      priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> que;
      que.push(make_pair(0, s));
      while (!que.empty()) {
        pair<int, int> p = que.top();
        que.pop();
        int v = p.second;
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
    vector<int> bellmanford(int s) {
      vector<int> d(n, INT_MAX);
      d[s] = 0;
      for (int _ = 0; _ < n; _++) {
        bool upd = false;
        for (int u = 0; u < n; u++)
          if (d[u] < INT_MAX)
            for (const auto &e : g[u]) {
              int v = e.to;
              if (d[v] > d[u] + e.cost)
                d[v] = d[u] + e.cost, upd = true;
            }
        if (!upd)
          return d;
      }
      queue<int> Q;
      for (int u = 0; u < n; u++)
        if (d[u] < INT_MAX)
          Q.emplace(u);
      while (!Q.empty()) {
        int u = Q.front();
        Q.pop();
        for (const auto &e : g[u]) {
          int v = e.to;
          if (d[v] > INT_MIN && (d[u] == INT_MIN || d[v] > d[u] + e.cost))
            d[v] = INT_MIN, Q.emplace(v);
        }
      }
      return d;
    }
    // O(V+E)
    vector<int> topologicalSort() {
      vector<int> d, ind(n);
      for (int i = 0; i < n; i++)
        for (auto e : g[i])
          ind[e.to]++;
      queue<int> que;
      for (int i = 0; i < n; i++)
        if (ind[i] == 0)
          que.push(i);
      while (!que.empty()) {
        int now = que.front();
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

int main(){
  int v,e;cin>>v>>e;
  DirectedGraph g(v);
  for(int i=0;i<e;i++){
    int a,b,d;cin>>a>>b>>d;
    g.add(a,b,d);
    g.add(b,a,d);
  }
  for(auto &i:g.g){
    for(auto &j:i){
      cout<<j.to<<" ";
    }
    cout<<"\n";
  }
}
