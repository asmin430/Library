# Levenshtein距離(編集距離)

Levenshtein距離(編集距離)とは，できる操作を，「文字の挿入」，「削除」，「置換」とした時，文字列Sを文字列Tに変換するのに必要な最小の操作回数のことである．

# LevenshteinDistance
```cpp
LevenshteinDistance(string S, string T)
```
文字列 `S`, `T` のLevenshtein距離を返す．

## 計算量
- $O(|S||T|)$
