template <uint32_t mod_, bool fast = false>
struct MontgomeryModInt {
private:
    using mint = MontgomeryModInt;
    using i32 = int32_t;
    using i64 = int64_t;
    using u32 = uint32_t;
    using u64 = uint64_t;

    static constexpr u32 get_r() {
        u32 ret = mod_;
        for (i32 i = 0; i < 4; i++) ret *= 2 - mod_ * ret;
        return ret;
    }

    static constexpr u32 r = get_r();

    static constexpr u32 n2 = -u64(mod_) % mod_;

    static_assert(r * mod_ == 1, "invalid, r * mod != 1");
    static_assert(mod_ < (1 << 30), "invalid, mod >= 2 ^ 30");
    static_assert((mod_ & 1) == 1, "invalid, mod % 2 == 0");

    u32 x;

    public:
    MontgomeryModInt() : x{} {}

    MontgomeryModInt(const i64 &a): x(reduce(u64(fast ? a : (a % mod() + mod())) * n2)) {}

    static constexpr u32 reduce(const u64 &b) {
        return u32(b >> 32) + mod() - u32((u64(u32(b) * r) * mod()) >> 32);
    }

    mint &operator+=(const mint &p) {
        if (i32(x += p.x - 2 * mod()) < 0) x += 2 * mod();
        return *this;
    }

    mint &operator-=(const mint &p) {
        if (i32(x -= p.x) < 0) x += 2 * mod();
        return *this;
    }

    mint &operator*=(const mint &p) {
        x = reduce(u64(x) * p.x);
        return *this;
    }

    mint &operator/=(const mint &p) {
        *this *= p.inv();
        return *this;
    }

    mint operator-() const { return mint() - *this; }

    mint operator+(const mint &p) const { return mint(*this) += p; }

    mint operator-(const mint &p) const { return mint(*this) -= p; }

    mint operator*(const mint &p) const { return mint(*this) *= p; }

    mint operator/(const mint &p) const { return mint(*this) /= p; }

    bool operator==(const mint &p) const {
        return (x >= mod() ? x - mod() : x) == (p.x >= mod() ? p.x - mod() : p.x);
    }

    bool operator!=(const mint &p) const {
        return (x >= mod() ? x - mod() : x) != (p.x >= mod() ? p.x - mod() : p.x);
    }

    u32 val() const {
        u32 ret = reduce(x);
        return ret >= mod() ? ret - mod() : ret;
    }

    mint pow(u64 n) const {
        mint ret(1), mul(*this);
        while (n > 0) {
            if (n & 1) ret *= mul;
            mul *= mul;
            n >>= 1;
        }
        return ret;
    }

    mint inv() const { return pow(mod() - 2); }

    friend ostream &operator<<(ostream &os, const mint &p) {
        return os << p.val();
    }

    friend istream &operator>>(istream &is, mint &a) {
        i64 t;
        is >> t;
        a = mint(t);
        return is;
    }

    static constexpr u32 mod() { return mod_; }

    static constexpr bool isNTT(){
        int p = mod() - 1;
        int cnt = 0;
        while(!(p & 1)){
            p >>= 1;
            ++cnt;
        }
        if(cnt >= 23) return true;
        else return false;
    }
};

template<class T_VAL, class T_MOD>
constexpr T_VAL mod_pow(T_VAL a, T_VAL n, T_MOD m) {
    T_VAL res = 1;
    while (n > 0) {
        if (n % 2 == 1) res = res * a % m;
        a = a * a % m;
        n >>= 1;
    }
    return res;
}

constexpr int calc_primitive_root(long long m){
    if(m == 1) return -1;
    if(m == 2) return 1;
    if(m == 998244353) return 3;
    if(m == 167772161) return 3;
    if(m == 469762049) return 3;
    if(m == 754974721) return 11;
    if(m == 645922817) return 3;
    if(m == 897581057) return 3;
    long long divs[20] = {};
    divs[0] = 2;
    long long cnt = 1;
    long long x = (m - 1) / 2;
    while(x % 2 == 0) x /= 2;
    for(long long i = 3; i * i <= x; i += 2){
        if(x % i == 0){
            divs[cnt++] = i;
            while(x % i == 0) x /= i;
        }
    }
    if(x > 1) divs[cnt++] = x;
    for(long long g = 2; ; g++){
        bool ok = true;
        for(int i = 0; i < cnt; i++){
            if(mod_pow(g, (m - 1) / divs[i], m) == 1){
                ok = false;
                break;
            }
        }
        if (ok) return g;
    }
}

template <uint32_t mod>
using modint = MontgomeryModInt<mod>;
using mint = modint<998244353>;

