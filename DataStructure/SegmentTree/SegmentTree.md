# SegmentTree(セグメント木)
　モノイドに対する区間演算の結果を高速に返すことのできるデータ構造．

# コンストラクタ
```cpp
(1) SegmentTree(LambdaMonoid(Op op, E e), int n)
(2) SegmentTree(LambdaMonoid(Op op, E e), vector<S> v)
```
1. 型 `S` ，演算 `op` ，単位元 `e` ，サイズ `n` で初期化する．各要素には単位元が格納されている．
2. 型 `S` ，演算 `op` ，単位元 `e` ，配列 `v` で初期化する．

　演算 `op` ，単位元 `e` はラムダ式の形で構築する．　

　例えば，int型の区間minを求めるサイズ `n` のSegmentTreeを構築したい時，
```cpp
auto op(int a, int b){ return min(a, b); };
auto e(){ return 2147483647; };
SegmentTree seg(LambdaMonoid(op, e), n);
```
のように記述する．

## 制約
- 型 `S` ，演算 `op` ，単位元 `e` でモノイドをなすこと．

## 計算量
- $O(n)$

# build
```cpp
void build(vector<S> v)
```
　配列 `v` で初期化

## 制約
- コンストラクタの `n` と `v` の長さが一致する．

## 計算量
- $O(n)$

# set
```cpp
void set(int k, S x)
```
　`k` 番目の要素を `x` に変更．

## 制約
- $0 \leq k \lt n$

## 計算量
- $O(\log n)$

# get
```cpp
S get(int k)
```
　`k` 番目の要素を返す．`seg[k]` の形でもできる．

## 制約
- $0 \leq k \lt n$

##
- $O(1)$

# apply
```cpp
void apply(int k, S x)
```
　`k` 番目の要素を `x` と二項演算したものに変更．

## 制約
- $0 \leq k \lt n$

## 計算量
- $O(\log n)$

# query
```cpp
S query(int l, int r)
```
　区間 $[l, r)$ に対して二項演算した結果を返す。

## 制約
- $0 \leq l \leq r \leq n$

## 計算量
- $O(\log n)$

# all_query
```cpp
S all_query()
```
　すべての要素に対して二項演算した結果を返す．

## 計算量
- $O(1)$

# find_first
```cpp
int find_first(int l, C check)
```
　$[l, x)$ が `check` を満たす最初の要素位置 $x$ を返す．存在しないとき $n$ を返す．`check` はラムダ式で与える．
 
## 制約
- $0 \leq l \leq n$

## 計算量
- $O(\log n)$

# find_last
```cpp
int find_last(int r, C check)
```
　$[x, r)$ が `check` を満たす最後の要素位置 $x$ を返す．存在しないとき $-1$ を返す．`check` はラムダ式で与える．

## 制約
- $0 \leq l \leq n$

## 計算量
- $O(\log n)$
