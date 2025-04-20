#pragma region header
#include<algorithm>    // ソート, 二分探索, 最大・最小
#include<array>        // 固定長配列（C++11~）
#include<bitset>       // ビット演算・フラグ管理
#include<cassert>      // デバッグ用の assert
#include<chrono>       // 時間計測（C++11~）
#include<cinttypes>    // int64_t, PRIu64 などのフォーマット指定
#include<climits>      // INT_MAX, INT_MIN など
#include<cmath>        // sqrt, pow, sin, cos など
#include<complex>      // FFT（高速フーリエ変換）など
#include<cstdio>       // printf, scanf（Cスタイル）
#include<cstring>      // memset, memcpy
#include<deque>        // 両端キュー（deque）
#include<functional>   // 比較関数 (greater, less)
#include<iomanip>      // 小数点表示の調整（setprecision など）
#include<iostream>     // 標準入出力
#include<iterator>     // イテレータ操作
#include<limits>       // numeric_limits
#include<map>         // std::map（連想配列）
#include<numeric>      // gcd, lcm, accumulate など
#include<queue>        // 優先度付きキュー (priority_queue)
#include<random>       // 乱数生成（乱択アルゴリズムなど）
#include<set>         // std::set（集合）
#include<sstream>      // 文字列ストリーム（stringstream）
#include<stack>        // スタック（LIFO）
#include<string>       // 文字列操作
#include<tuple>        // tuple, tie
#include<type_traits>  // 型特性（C++11~）
#include<unordered_map> // ハッシュマップ（unordered_map）
#include<unordered_set> // ハッシュセット（unordered_set）
#include<utility>      // pair, swap, move など
#include<vector>       // 動的配列（最重要）
using namespace std;
struct Init{Init(){std::cin.tie(0); ios::sync_with_stdio(false); cout << setprecision(20) << fixed;}}init;
using ll = long long;
using ull = unsigned long long;
using ld = long double;
#define all(x) begin((x)), end((x))
#define pb push_back
#define mp make_pair
#define mt make_tuple
#define uq(v) v.erase(unique(begin(v), end(v)), end(v))
#define _overload4(_1,_2,_3,_4,name,...) name
#define _overload3(_1,_2,_3,name,...) name
#define _rep1(n) for(int i=0;i<n;++i)
#define _rep2(i,n) for(int i=0;i<n;++i)
#define _rep3(i,a,b) for(int i=a;i<b;++i)
#define _rep4(i,a,b,c) for(int i=a;i<b;i+=c)
#define rep(...) _overload4(__VA_ARGS__,_rep4,_rep3,_rep2,_rep1)(__VA_ARGS__)
#define _rrep1(n) for(int i=(n)-1;i>=0;i--)
#define _rrep2(i,n) for(int i=(n)-1;i>=0;i--)
#define _rrep3(i,a,b) for(int i=(b)-1;i>=(a);i--)
#define _rrep4(i,a,b,c) for(int i=a+(b-a-1)/c*c;i>=a;i-=c)
#define rrep(...) _overload4(__VA_ARGS__,_rrep4,_rrep3,_rrep2,_rrep1)(__VA_ARGS__)
template<class T> using pq = priority_queue<T>;
template<class T> using pq_g = priority_queue<T, vector<T>, greater<T>>;
template<class T> bool chmax(T &a, const T &b){if(a < b){a = b; return 1; } return 0;}
template<class T> bool chmin(T &a, const T &b){if(a > b){a = b; return 1; } return 0;}
template<class T> auto min(const T& a){ return *min_element(all(a)); }
template<class T> auto max(const T& a){ return *max_element(all(a)); }
constexpr ull INF = (1ULL << 61) + (1ULL << 30);
constexpr int inf = (1 << 30);
constexpr ld EPS = 1e-9;
constexpr ld PI = std::acos(ld(-1));
constexpr int dx[] = {1, 0, -1, 0, 1, 1, -1, -1};
constexpr int dy[] = {0, 1 ,0, -1, 1, -1, 1, -1};

#pragma endregion header



int main(){
    
}
