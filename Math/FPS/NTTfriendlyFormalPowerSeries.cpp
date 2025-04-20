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
};

template <uint32_t mod>
using modint = MontgomeryModInt<mod>;
using mint = modint<998244353>;

template<typename Mint>
struct NTTFriendlyConvolution{
    static vector<Mint> roots, iroots, rate3, irate3;
    static int max_base;

    NTTFriendlyConvolution() = default;

    static void init(){
        if(roots.empty()){
            const unsigned mod = Mint::mod();
            assert(mod >= 3 && mod % 2 == 1);
            auto tmp = mod - 1;
            max_base = 0;
            while(tmp % 2 == 0) tmp >>= 1, ++max_base;
            Mint root = 2;
            while(root.pow((mod - 1) >> 1) == 1){
                root += 1;
            }
            assert(root.pow(mod - 1) == 1);

            roots.resize(max_base + 1);
            iroots.resize(max_base + 1);
            rate3.resize(max_base + 1);
            irate3.resize(max_base + 1);

            roots[max_base] = root.pow((mod - 1) >> max_base);
            iroots[max_base] = (Mint)1 / roots[max_base];
            for(int i = max_base - 1; i >= 0; --i){
                roots[i] = roots[i + 1] * roots[i + 1];
                iroots[i] = iroots[i + 1] * iroots[i + 1];
            }
            {
                Mint prod = 1, iprod = 1;
                for(int i = 0; i <= max_base - 3; ++i){
                    rate3[i] = roots[i + 3] * prod;
                    irate3[i] = iroots[i + 3] * iprod;
                    prod *= iroots[i + 3];
                    iprod *= roots[i + 3];
                }
            }
        }
    }

    static void ntt(vector<Mint> &a){
        init();
        const int n = (int)a.size();
        assert((n & (n - 1)) == 0);
        int h = __builtin_ctz(n);
        assert(h <= max_base);
        int len = 0;
        Mint imag = roots[2];
        if(h & 1){
            int p = 1 << (h - 1);
            Mint rot = 1;
            for(int i = 0; i < p; ++i){
                auto r = a[i + p];
                a[i + p] = a[i] - r;
                a[i] += r;
            }
            ++len;
        }
        for(;len + 1 < h; len += 2){
            int p = 1 << (h - len - 2);
            {
                for(int i = 0; i < p; ++i){
                    auto a0 = a[i];
                    auto a1 = a[i + p];
                    auto a2 = a[i + 2 * p];
                    auto a3 = a[i + 3 * p];
                    auto a1na3imag = (a1 - a3) * imag;
                    auto a0a2 = a0 + a2;
                    auto a1a3 = a1 + a3;
                    auto a0na2 = a0 - a2;
                    a[i] = a0a2 + a1a3;
                    a[i + p] = a0a2 - a1a3;
                    a[i + 2 * p] = a0na2 + a1na3imag;
                    a[i + 3 * p] = a0na2 - a1na3imag;
                }
            }
            Mint rot = rate3[0];
            for(int s = 1; s < (1 << len); ++s){
                int offset = s << (h - len);
                Mint rot2 = rot * rot;
                Mint rot3 = rot2 * rot;
                for(int i = 0; i < p; ++i){
                    auto a0 = a[i + offset];
                    auto a1 = a[i + offset + p] * rot;
                    auto a2 = a[i + offset + 2 * p] * rot2;
                    auto a3 = a[i + offset + 3 * p] * rot3;
                    auto a1na3imag = (a1 - a3) * imag;
                    auto a0a2 = a0 + a2;
                    auto a1a3 = a1 + a3;
                    auto a0na2 = a0 - a2;
                    a[i + offset] = a0a2 + a1a3;
                    a[i + offset + p] = a0a2 - a1a3;
                    a[i + offset + 2 * p] = a0na2 + a1na3imag;
                    a[i + offset + 3 * p] = a0na2 - a1na3imag;
                }
                rot *= rate3[__builtin_ctz(~s)];
            }
        }
    }

    static void intt(vector<Mint> &a, bool f = true){
        init();
        const int n = (int)a.size();
        assert((n & (n - 1)) == 0);
        int h = __builtin_ctz(n);
        assert(h <= max_base);
        int len = h;
        Mint iimag = iroots[2];
        for(; len > 1; len -= 2){
            int p = 1 << (h - len);
            {
                for(int i = 0; i < p; i++){
                    auto a0 = a[i];
                    auto a1 = a[i + 1 * p];
                    auto a2 = a[i + 2 * p];
                    auto a3 = a[i + 3 * p];
                    auto a2na3iimag = (a2 - a3) * iimag;
                    auto a0na1 = a0 - a1;
                    auto a0a1 = a0 + a1;
                    auto a2a3 = a2 + a3;
                    a[i] = a0a1 + a2a3;
                    a[i + 1 * p] = (a0na1 + a2na3iimag);
                    a[i + 2 * p] = (a0a1 - a2a3);
                    a[i + 3 * p] = (a0na1 - a2na3iimag);
                }
            }
            Mint irot = irate3[0];
            for(int s = 1; s < (1 << (len - 2)); s++){
                int offset = s << (h - len + 2);
                Mint irot2 = irot * irot;
                Mint irot3 = irot2 * irot;
                for(int i = 0; i < p; i++){
                    auto a0 = a[i + offset];
                    auto a1 = a[i + offset + 1 * p];
                    auto a2 = a[i + offset + 2 * p];
                    auto a3 = a[i + offset + 3 * p];
                    auto a2na3iimag = (a2 - a3) * iimag;
                    auto a0na1 = a0 - a1;
                    auto a0a1 = a0 + a1;
                    auto a2a3 = a2 + a3;
                    a[i + offset] = a0a1 + a2a3;
                    a[i + offset + 1 * p] = (a0na1 + a2na3iimag) * irot;
                    a[i + offset + 2 * p] = (a0a1 - a2a3) * irot2;
                    a[i + offset + 3 * p] = (a0na1 - a2na3iimag) * irot3;
                }
                irot *= irate3[__builtin_ctz(~s)];
            }
        }
        if(len >= 1){
            int p = 1 << (h - 1);
            for(int i = 0; i < p; i++){
                auto ajp = a[i] - a[i + p];
                a[i] += a[i + p];
                a[i + p] = ajp;
            }
        }
        if(f){
            Mint inv_sz = Mint(1) / n;
            for(int i = 0; i < n; i++) a[i] *= inv_sz;
        }
    }

    static vector<Mint> multiply(vector<Mint> a, vector<Mint> b){
        int need = a.size() + b.size() - 1;
        int nbase = 1;
        while((1 << nbase) < need) nbase++;
        int sz = 1 << nbase;
        a.resize(sz, 0);
        b.resize(sz, 0);
        ntt(a);
        ntt(b);
        Mint inv_sz = Mint(1) / sz;
        for(int i = 0; i < sz; i++) a[i] *= b[i] * inv_sz;
        intt(a, false);
        a.resize(need);
        return a;
    }
};

template <typename Mint>
vector<Mint> NTTFriendlyConvolution<Mint>::roots = vector<Mint>();
template <typename Mint>
vector<Mint> NTTFriendlyConvolution<Mint>::iroots = vector<Mint>();
template <typename Mint>
vector<Mint> NTTFriendlyConvolution<Mint>::rate3 = vector<Mint>();
template <typename Mint>
vector<Mint> NTTFriendlyConvolution<Mint>::irate3 = vector<Mint>();
template <typename Mint>
int NTTFriendlyConvolution<Mint>::max_base = 0;

template <typename T>
T mod_pow(T x, int64_t n, const T &p){
    T ret = 1;
    while(n > 0){
        if(n & 1)(ret *= x) %= p;
        (x *= x) %= p;
        n >>= 1;
    }
    return ret % p;
}

template<typename T>
T mod_sqrt(const T &a, const T &p){
    if(a == 0) return 0;
    if(p == 2) return a;
    if(mod_pow(a, (p - 1) >> 1, p) != 1) return -1;
    T b = 1;
    while(mod_pow(b, (p - 1) >> 1, p) == 1) ++b;
    T e = 0, m = p - 1;
    while(m % 2 == 0) m >>= 1, ++e;
    T x = mod_pow(a, (m - 1) >> 1, p);
    T y = a * (x * x % p) % p;
    (x *= a) %= p;
    T z = mod_pow(b, m, p);
    while(y != 1){
        T j = 0, t = y;
        while(t != 1){
            j += 1;
            (t *= t) %= p;
        }
        z = mod_pow(z, T(1) << (e - j - 1), p);
        (x *= z) %= p;
        (z *= z) %= p;
        (y *= z) %= p;
        e = j;
    }
    return x;
}


template<typename T>
struct NTTFriendlyFormalPowerSeries: vector<T>{
    using vector<T>::vector;
    using P = NTTFriendlyFormalPowerSeries;
    using NTT = NTTFriendlyConvolution<T>;

    P pre(int deg) const {
        return P(begin(*this), begin(*this) + min((int)this->size(), deg));
    }

    P rev(int deg = -1) const {
        P ret(*this);
        if(deg != -1) ret.resize(deg, T(0));
        reverse(begin(ret), end(ret));
        return ret;
    }

    void shrink(){
        while(this->size() && this->back() == T(0)) this->pop_back();
    }

    P operator+(const P &r) const{return P(*this) += r;}
    P operator+(const T &v) const{return P(*this) += v;}
    P operator-(const P &r) const{return P(*this) -= r;}
    P operator-(const T &v) const{return P(*this) -= v;}
    P operator*(const P &r) const{return P(*this) *= r;}
    P operator*(const T &v) const{return P(*this) *= v;}
    P operator/(const P &r) const{return P(*this) /= r;}
    P operator%(const P &r) const{return P(*this) %= r;}

    P &operator+=(const P &r){
        if(r.size() > this->size()) this->resize(r.size());
        for(int i = 0; i < (int)r.size(); ++i) (*this)[i] += r[i];
        return *this;
    }

    P &operator-=(const P &r){
        if(r.size() > this->size()) this->resize(r.size());
        for(int i = 0; i < (int)r.size(); ++i) (*this)[i] -= r[i];
        return *this;
    }

    P &operator*=(const P &r){
        if(this->empty() || r.empty()){
            this->clear();
            return *this;
        }
        auto ret = NTT::multiply(*this, r);
        return *this = {begin(ret), end(ret)};
    }

    P &operator/=(const P &r){
        if(this->size() < r.size()){
            this->clear();
            return *this;
        }
        int n = this->size() - r.size() + 1;
        return *this = (rev().pre(n) * r.rev().inv(n)).pre(n).rev(n);
    }

    P &operator%=(const P &r){
        *this -= *this / r * r;
        shrink();
        return *this;
    }

    pair<P, P> div_mod(const P &r){
        P q = *this / r;
        P x = *this - (q * r);
        x.shrink();
        return make_pair(q, x);
    }

    P operator-() const {
        P ret(this->size());
        for(int i = 0; i < (int)this->size(); ++i) ret[i] = -(*this)[i];
        return ret;
    }

    P &operator+=(const T &r){
        if(this->empty()) this->resize(1);
        (*this)[0] += r;
        return *this;
    }

    P &operator-=(const T &r){
        if(this->empty()) this->resize(1);
        (*this)[0] -= r;
        return *this;
    }

    P &operator*=(const T &v){
        for(int i = 0; i < (int)this->size(); ++i) (*this)[i] *= v;
        return *this;
    }

    P dot(P r) const{
        P ret(min(this->size(), r.size()));
        for(int i = 0; i < (int)ret.size(); ++i) ret[i] = (*this)[i] * r[i];
        return ret;
    }

    P operator>>(int sz) const{
        if((int)this->size() <= sz) return {};
        P ret(*this);
        ret.erase(ret.begin(), ret.begin() + sz);
        return ret;
    }

    P operator<<(int sz) const{
        P ret(*this);
        ret.insert(ret.begin(), sz, T(0));
        return ret;
    }

    T operator()(T x) const{
        T r = 0, w = 1;
        for(auto &v: *this){
            r += w * v;
            w *= x;
        }
        return r;
    }

    P diff() const{
        const int n = (int)this->size();
        P ret(max(0, n - 1));
        for(int i = 1; i < n; ++i) ret[i - 1] = (*this)[i] * T(i);
        return ret;
    }

    P integral() const{
        const int n = (int)this->size();
        P ret(n + 1);
        ret[0] = T(0);
        for(int i = 0; i < n; ++i) ret[i + 1] = (*this)[i] / T(i + 1);
        return ret;
    }

    // F(0) must not be 0
    P inv(int deg = -1) const{
        assert(((*this)[0]) != T(0));
        const int n = (int)this->size();
        if(deg == -1) deg = n;
        P res(deg);
        res[0] = {T(1) / (*this)[0]};
        for(int d = 1; d < deg; d <<= 1){
            P f(2 * d), g(2 * d);
            for(int j = 0; j < min(n, 2 * d); ++j) f[j] = (*this)[j];
            for(int j = 0; j < d; ++j) g[j] = res[j];
            NTT::ntt(f);
            NTT::ntt(g);
            f = f.dot(g);
            NTT::intt(f);
            for(int j = 0; j < d; ++j) f[j] = 0;
            NTT::ntt(f);
            for(int j = 0; j < 2 * d; ++j) f[j] *= g[j];
            NTT::intt(f);
            for(int j = d; j < min(2 * d, deg); ++j) res[j] = -f[j];
        }
        return res;
    }

    // F(0) must be 1
    P log(int deg = -1) const{
        assert((*this)[0] == T(1));
        const int n = (int)this->size();
        if(deg == -1) deg = n;
        return (this->diff() * this->inv(deg)).pre(deg - 1).integral();
    }

    P sqrt(int deg = -1, const function<T(T)> &get_sqrt = [](T){return T(1);}) const {
        const int n = (int)this->size();
        if(deg == -1) deg = n;
        if((*this)[0] == T(0)){
            for(int i = 1; i < n; i++){
                if((*this)[i] != T(0)){
                    if(i & 1) return {};
                    if(deg - i / 2 <= 0) break;
                    auto ret = (*this >> i).sqrt(deg - i / 2, get_sqrt);
                    if(ret.empty()) return {};
                    ret = ret << (i / 2);
                    if((int)ret.size() < deg) ret.resize(deg, T(0));
                    return ret;
                }
            }
            return P(deg, 0);
        }
        auto sqr = T(get_sqrt((*this)[0]));
        if(sqr * sqr != (*this)[0]) return {};
        P ret{sqr};
        T inv2 = T(1) / T(2);
        for(int i = 1; i < deg; i <<= 1){
            ret = (ret + pre(i << 1) * ret.inv(i << 1)) * inv2;
        }
        return ret.pre(deg);
    }

    P sqrt(const function<T(T)> &get_sqrt, int deg = -1) const {
        return sqrt(deg, get_sqrt);
    }

    // F(0) must be 0
    P exp(int deg = -1) const{
        if(deg == -1) deg = this->size();
        assert((*this)[0] == T(0));

        P inv;
        inv.reserve(deg + 1);
        inv.push_back(T(0));
        inv.push_back(T(1));

        auto inplace_integral = [&](P &F) -> void {
            const int n = (int)F.size();
            auto mod = T::mod();
            while((int)inv.size() <= n){
                int i = inv.size();
                inv.push_back((-inv[mod % i]) * (mod / i));
            }
            F.insert(begin(F), T(0));
            for(int i = 1; i <= n; ++i) F[i] *= inv[i];
        };

        auto inplace_diff = [](P &F) -> void {
            if(F.empty()) return;
            F.erase(begin(F));
            T coeff = 1, one = 1;
            for(int i = 0; i < (int)F.size(); i++){
                F[i] *= coeff;
                coeff += one;
            }
        };

        P b{1, 1 < (int)this->size() ? (*this)[1] : 0}, c{1}, z1, z2{1, 1};
        for(int m = 2; m < deg; m *= 2){
            auto y = b;
            y.resize(2 * m);
            NTT::ntt(y);
            z1 = z2;
            P z(m);
            for(int i = 0; i < m; ++i) z[i] = y[i] * z1[i];
            NTT::intt(z);
            fill(begin(z), begin(z) + m / 2, T(0));
            NTT::ntt(z);
            for(int i = 0; i < m; ++i) z[i] *= -z1[i];
            NTT::intt(z);
            c.insert(end(c), begin(z) + m / 2, end(z));
            z2 = c;
            z2.resize(2 * m);
            NTT::ntt(z2);
            P x(begin(*this), begin(*this) + min<int>(this->size(), m));
            inplace_diff(x);
            x.push_back(T(0));
            NTT::ntt(x);
            for(int i = 0; i < m; ++i) x[i] *= y[i];
            NTT::intt(x);
            x -= b.diff();
            x.resize(2 * m);
            for(int i = 0; i < m - 1; ++i) x[m + i] = x[i], x[i] = T(0);
            NTT::ntt(x);
            for(int i = 0; i < 2 * m; ++i) x[i] *= z2[i];
            NTT::intt(x);
            x.pop_back();
            inplace_integral(x);
            for(int i = m; i < min<int>(this->size(), 2 * m); ++i) x[i] += (*this)[i];
            fill(begin(x), begin(x) + m, T(0));
            NTT::ntt(x);
            for(int i = 0; i < 2 * m; ++i) x[i] *= y[i];
            NTT::intt(x);
            b.insert(end(b), begin(x) + m, end(x));
        }
        return P{begin(b), begin(b) + deg};
    }

    P pow(int64_t k, int deg = -1) const{
        const int n = (int)this->size();
        if(deg == -1) deg = n;
        if(k == 0){
            P ret(deg, T(0));
            ret[0] = T(1);
            return ret;
        }
        for(int i = 0; i < n; i++){
            if(i * k > deg) return P(deg, T(0));
            if((*this)[i] != T(0)){
                T rev = T(1) / (*this)[i];
                P ret = (((*this * rev) >> i).log() * k).exp() * ((*this)[i].pow(k));
                ret = (ret << (i * k)).pre(deg);
                if ((int)ret.size() < deg) ret.resize(deg, T(0));
                return ret;
            }
        }
        return *this;
    }

    P mod_pow(int64_t k, P g) const{
        P modinv = g.rev().inv();
        auto get_div = [&](P base){
            if(base.size() < g.size()){
                base.clear();
                return base;
            }
            int n = base.size() - g.size() + 1;
            return (base.rev().pre(n) * modinv.pre(n)).pre(n).rev(n);
        };
        P x(*this), ret{1};
        while(k > 0){
            if(k & 1){
                ret *= x;
                ret -= get_div(ret) * g;
                ret.shrink();
            }
            x *= x;
            x -= get_div(x) * g;
            x.shrink();
            k >>= 1;
        }
        return ret;
    }

    P taylor_shift(T c) const{
        int n = (int)this->size();
        vector<T> fact(n), rfact(n);
        fact[0] = rfact[0] = T(1);
        for(int i = 1; i < n; i++) fact[i] = fact[i - 1] * T(i);
        rfact[n - 1] = T(1) / fact[n - 1];
        for(int i = n - 1; i > 1; i--) rfact[i - 1] = rfact[i] * T(i);
        P p(*this);
        for(int i = 0; i < n; i++) p[i] *= fact[i];
        p = p.rev();
        P bs(n, T(1));
        for(int i = 1; i < n; i++) bs[i] = bs[i - 1] * c * rfact[i] * fact[i - 1];
        p = (p * bs).pre(n);
        p = p.rev();
        for(int i = 0; i < n; i++) p[i] *= rfact[i];
        return p;
    }
};

template<class Mint>
using FPS = NTTFriendlyFormalPowerSeries<Mint>;
