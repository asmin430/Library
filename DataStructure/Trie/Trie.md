# Trie木(PrefixTree)
　効率的な検索のために使われるデータ構造．

# コンストラクタ
```cpp
Trie<int char_size, int base>
```
char_size := 文字の種類  
base      := 開始文字

## 計算量
- $O(1)$

# insert
```cpp
void insert(string S)
```
文字列 `S` を辞書に追加．

## 計算量

- $O(|S|)$

# search
```cpp
bool search(string S)
```
辞書に文字列 `S` が存在するかどうかを出力．

## 計算量
- $O(|S|)$

# start_with
```cpp
bool start_with(string S)
```
文字列 `S` を接頭辞に持つ文字列が辞書に存在するかどうかを出力．

## 計算量
- $O(|S|)$

# count
```cpp
int count()
```
追加した文字列の総数を出力．

## 計算量
- $O(1)$

# size
```cpp
int size()
```
Trie木の頂点数を出力．

## 計算量
- $O(1)$
