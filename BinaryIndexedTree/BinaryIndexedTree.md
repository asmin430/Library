# BinaryIndexedTree

# コンストラクタ

```cpp
(1) BinaryIndexedTree(int n)
(2) BinaryIndexedTree(const vector<T> &v)
```

1. 数列を長さ `n` で初期化．各要素は $0$ ．
2. 数列の各要素を 配列 `v` で初期化．

## 計算量

- $O(n)$

# build

```cpp
void build(const vector<T> &v)
```

数列の各要素を配列 `v` で初期化．

## 制約

- `n` と `v` の長さが一致する．

## 計算量

- $O(n)$

# add

```cpp
void add(int k, const T &x)
```

`k` 番目の要素に値 `x` を加える．

## 制約

- $0 \leq k \lt n$

## 計算量

- $O(\log n)$

# sum

```cpp
(1) T sum(int r) const 
(2) T sum(int l, int r) const
```

1. 数列の区間 $[0, r)$ の要素の総和を返す．
2. 数列の区間 $[l, r)$ の要素の総和を返す．

## 制約

- $0 \leq l \leq r \leq n$

## 計算量

- $O(\log n)$

# lower_bound

```cpp
int lower_bound(T x) const
```

数列の区間 $[0,k]$ の要素の総和が $x$ 以上となる最小の $k$ を返す．

## 制約

- 数列は広義単調増加

## 計算量

- $O(\log n)$

## upper_bound

```cpp
int upper_bound(T x) const
```

数列の区間 $[0,k]$ の要素の総和が $x$ を上回る最小の $k$ を返す．

## 制約

- 数列は広義単調増加

## 計算量

* $O(\log n)$
