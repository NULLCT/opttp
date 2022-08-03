// https://jetbead.hatenablog.com/entry/20120623/1340419446
#include <climits>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>

using namespace std;

struct STATE { // 状態構造体
  double x;
};

//焼きなまし
class SA {
  STATE state;  // 現在の状態
  STATE ans;    // 暫定最適状態
  double score; // 暫定最適状態ansを評価関数に通したスコア
  double T;     // 温度
  const int R;  // 反復回数

  double frand() {
    return ((double) rand() / (RAND_MAX));
  }

  //評価関数
  double calc_score(STATE &state) {
    const double x = state.x;
    return -0.7*pow(x,4)+1.3*pow(x,3)+4*pow(x,2)-3.5*x+1.9;
  }

  //近傍からランダムに選ぶ
  void modify(STATE &state) {
    if (frand() < 0.5)
      state.x += 0.01;
    else
      state.x -= 0.01;
  }

  //温度の更新
  double next_T(double t) {
    return t * 0.995;
  }

public:
  SA(STATE &_state, double _t, int _r) : T(_t), R(_r) { // 温度の初期値、反復回数の初期値
    state = _state;
    ans = _state;
    score = calc_score(ans);
  }

  //探索
  STATE simulated_annealing() {
    while (T > 1.0) { //十分冷えるまで
      for (int i = 0; i < R; i++) {
        // xの近傍からランダムに選ぶ
        STATE new_state = state;
        modify(new_state);
        //変化量
        double delta = calc_score(state) - calc_score(new_state);

        if (delta < 0.0) { //よりよい解が見つかった場合
          state = new_state;
        } else if (exp(-delta / T) > frand()) { //ある程度の確率で探索を許す
          state = new_state;
        }

        //最適解の更新
        if (calc_score(state) > score) {
          score = calc_score(state);
          ans = state;
        }
      }

      //温度の更新
      T = next_T(T);
    }

    return ans;
  }
};

int main() {
  srand(0);
  STATE state;
  state.x = -0.5;
  SA sa(state, 1000, 1000);

  cout<<fixed<<setprecision(32);
  cout << sa.simulated_annealing().x << "\n";
}
