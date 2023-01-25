![opttp](https://user-images.githubusercontent.com/45321757/214467985-86846b51-1af3-4303-a4b3-57b075626a75.png)
# About
大域的に利用可能な運送経路最適化手法について、簡易的に表現した運送モデルの構築後、ヒューリスティック的最適化手法である焼きなまし法に独自考案した評価関数、状態遷移手法を実装し、焼きなましパラメーターを調整して理想的な運送経路の導出を高速に行います

# Build
- C++20対応コンパイラ
- CMake
が必要です
### CLI
```bash
git clone https://github.com/NULLCT/opttp
cd opttp
mkdir build
cd build
cmake ..
cmake --build .
```
### Features
- [x] CUIでの最適化&焼きなまし温度出力
- [ ] 焼きなましのログを吐く
- [ ] GUIラッパーの作成
