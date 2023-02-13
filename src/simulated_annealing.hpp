/**
 * @file simulated_annealing.hpp
 * @brief 焼きなまし法を適応します
 */
#pragma once
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <tuple>

#include "experimental_model.hpp"
#include "state.hpp"

class SA {
  STATE state;               // 現在の状態
  double temperature;               // 温度
  const int repnum;          // 反復回数
  const double coolingcoef;  // 冷却係数
  MODEL model;               // モデル
  STATE beststate;           // 暫定最適状態
  int64_t beststate_score;   // 暫定最適状態ansを評価関数に通したスコア

  double frand();                         // [0,1]の乱数を返す

 public:
  int64_t evalScore(const STATE &_state);  // 評価関数
  int64_t evalScoreFast(const STATE &_state_before_adaptation, const STATE &_state_after_adaptation);  // 評価関数
  void modify(STATE &_state);              // s-opt法 近傍に行き帰りするパスを追加
  SA(STATE &_state, double _temperature, int _repnum, double _coolingcoef, MODEL &_model) : state(_state),                      // 状態
                                                                                     temperature(_temperature),                        // 温度
                                                                                     repnum(_repnum),                    // 試行回数
                                                                                     coolingcoef(_coolingcoef),          // 冷却係数
                                                                                     model(_model),                      // 運送モデル
                                                                                     beststate(_state),                  // 最良状態
                                                                                     beststate_score(evalScore(_state))  // 最良状態の時のスコア
  {}

  // 焼きなまし法
  STATE simulated_annealing();
};
