// https://gasin.hatenadiary.jp/entry/2019/09/03/162613
#pragma once
#include <vector>
#include <cmath>
#include <tuple>

using namespace std;

constexpr double START_TEMP = 0;
constexpr double END_TEMP = 0;
constexpr double START_TIME = 0;
constexpr double TIME_LIMIT = 0;
constexpr int INF = 0;

struct STATE { // 状態
  vector<vector<pair<int,int>>> transition; // transition[i][j]{t,p} はi番目の運送者が0..j..で時刻tに座標pに進むこと
};

void init(STATE &state) { // 初期化
}

void modify(STATE &state) { // 状態遷移
}

int calc_score(STATE &state) { // 評価関数計算
}

void simulated_annealing() { // 焼きなまし法
  STATE state;
  init(state);
  while (true) {
    double now_time;
    if (now_time - START_TIME > TIME_LIMIT)
      break;

    STATE new_state = state;
    modify(new_state);
    int new_score = calc_score(new_state);
    int pre_score = calc_score(state);

    double temp = START_TEMP + (END_TEMP - START_TEMP) * (now_time - START_TIME) / TIME_LIMIT;
    double prob = exp((new_score - pre_score) / temp);

    if(prob > (rand()%INF)/(double)INF) {
      state = new_state;
    }
  }
}
