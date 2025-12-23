// depend on Convolution.cpp
template<class mint> struct BiCoef {
    vector<mint> fact_, inv_, finv_;
    constexpr BiCoef() {}
    constexpr BiCoef(int n) : fact_(n, 1), inv_(n, 1), finv_(n, 1) {
        init(n);
    }
    constexpr void init(int n) {
        fact_.assign(n, 1), inv_.assign(n, 1), finv_.assign(n, 1);
        int MOD = fact_[0].mod();
        for(int i = 2; i < n; i++){
            fact_[i] = fact_[i-1] * i;
            inv_[i] = -inv_[MOD%i] * (MOD/i);
            finv_[i] = finv_[i-1] * inv_[i];
        }
    }
    constexpr mint com(int n, int k) const {
        if (n < k || n < 0 || k < 0) return 0;
        return fact_[n] * finv_[k] * finv_[n-k];
    }
    constexpr mint fact(int n) const {
        if (n < 0) return 0;
        return fact_[n];
    }
    constexpr mint inv(int n) const {
        if (n < 0) return 0;
        return inv_[n];
    }
    constexpr mint finv(int n) const {
        if (n < 0) return 0;
        return finv_[n];
    }
};

template<typename T> 
struct FormalPowerSeries: vector<T> {
    using vector<T>::vector;
    using P = FormalPowerSeries;
    using NTT = Convolution_base<T>;
    static const int SPARSE_BOARDER = 50;
    P pre(int deg) const {
        return P(begin(*this), begin(*this) + min((int)this->size(), deg));
    }
    P rev(int deg = -1) const {
        P ret(*this);
        if(deg != -1) ret.resize(deg);
        reverse(begin(ret), end(ret));
        return ret;
    }
    void shrink(){
        while(this->size() && this->back() == T(0)) this->pop_back();
    }
    P operator+(const P &r) const {return P(*this) += r; }
    P operator+(const T &v) const {return P(*this) += v; }
    P operator-(const P &r) const {return P(*this) -= r; }
    P operator-(const T &v) const {return P(*this) -= v; }
    P operator*(const P &r) const {return P(*this) *= r; }
    P operator*(const T &v) const {return P(*this) *= v; }
    P operator/(const P &r) const {return P(*this) /= r; }
    P operator%(const P &r) const {return P(*this) %= r; }
    P operator<<(int x) const {return P(*this) <<= x; }
    P operator>>(int x) const {return P(*this) >>= x; }
    P& operator+=(const P &r){
        if(r.size() > this->size()) this->resize(r.size());
        for(int i = 0; i < (int)r.size(); ++i) (*this)[i] += r[i];
        return *this;
    }
    P& operator-=(const P &r){
        if(r.size() > this->size()) this->resize(r.size());
        for(int i = 0; i < (int)r.size(); ++i) (*this)[i] -= r[i];
        return *this;
    }
    P& operator*=(const P &r){
        if(this->empty() || r.empty()){
            this->clear();
            return *this;
        }
        auto ret = ArbitraryConvolution(*this, r);
        return *this = {begin(ret), end(ret)};
    }
    P& operator/=(const P &r){
        if(this->size() < r.size()){
            this->clear();
            return *this;
        }
        int n = this->size() - r.size() + 1;
        return *this = (rev().pre(n) * r.rev().inv(n)).pre(n).rev(n);
    }
    P& operator%=(const P &r){
        *this -= *this / r * r;
        shrink();
        return *this;
    }
    pair<P, P> div_mod(const P &r){
        P q = *this / r;
        P x = *this - q * r;
        x.shrink();
        return make_pair(q, x);
    }
    P operator-() const {
        P ret(this->size());
        for(int i = 0; i < (int)this->size(); ++i) ret[i] = -(*this)[i];
        return ret;
    }
    P& operator+=(const T &r){
        if(this->empty()) this->resize(1);
        (*this)[0] += r;
        return *this;
    }
    P& operator-=(const T &r){
        if(this->empty()) this->resize(1);
        (*this)[0] -= r;
        return *this;
    }
    P& operator*=(const T &r){
        for(int i = 0; i < (int)this->size(); ++i) (*this)[i] *= r;
        return *this;
    }
    P& operator<<=(int x){
        P res(x, 0);
        res.insert(res.end(), begin(*this), end(*this));
        return *this = res;
    }
    P& operator>>=(int x){
        P res;
        res.insert(res.end(), begin(*this) + x, end(*this));
        return *this = res;
    }
    T operator()(T x) const {
        T r = 0, w = 1;
        for(auto &v: *this){
            r += w * v;
            w *= x;
        }
        return r;
    }
    P diff() const {
        const int n = (int)this->size();
        if(n <= 0) return {};
        P ret(n - 1);
        for(int i = 1; i < n; ++i) ret[i - 1] = (*this)[i] * T(i);
        return ret;
    }
    P integral() const {
        const int n = (int)this->size();
        P ret(n + 1, 0);
        for(int i = 0; i < n; ++i) ret[i + 1] = (*this)[i] / T(i + 1);
        return ret;
    }
    constexpr int count_terms() const {
        int res = 0;
        for(int i = 0; i < (int)this->size(); ++i) if((*this)[i] != T(0)) ++res;
        return res;
    }
    P inv_ntt(int deg = -1) const {
        assert(this->size() >= 1 && (*this)[0] != 0);
        if(deg < 0) deg = (int)this->size();
        P res(deg, 0);
        res[0] = T(1) / (*this)[0];
        for(int d = 1; d < deg; d <<= 1){
            P g(d * 2), h(d * 2);
            for(int i = 0; i < min((int)this->size(), d * 2); ++i) g[i] = (*this)[i];
            for(int i = 0; i < d; ++i) h[i] = res[i];
            NTT::ntt(g);
            NTT::ntt(h);
            for(int i = 0; i < d * 2; ++i) g[i] *= h[i];
            NTT::intt(g);
            for(int i = 0; i < d; ++i) g[i] = 0;
            NTT::ntt(g);
            for(int i = 0; i < d * 2; ++i) g[i] *= h[i];
            NTT::intt(g);
            for(int i = d; i < min(deg, d * 2); ++i) res[i] = -g[i];
        }
        return res.pre(deg); 
    }
    P inv_sparse(int deg = -1) const {
        assert(this->size() >= 1 && (*this)[0] != 0);
        if(deg < 0) deg = (int)this->size();
        vector<pair<int, T>> dat;
        for(int i = 1; i < (int)this->size(); ++i) if((*this)[i] != T(0)){
            dat.emplace_back(i, (*this)[i]);
        }
        P res(deg);
        res[0] = (*this)[0].inv();
        for(int i = 1; i < deg; ++i){
            T r = 0;
            for(auto &&[k, val]: dat){
                if(k > i) break;
                r -= val * res[i - k];
            }
            res[i] = r * res[0];
        }
        return res;

    }
    P inv(int deg = -1) const {
        if(count_terms() <= SPARSE_BOARDER) return inv_sparse(deg);
        if(T::isNTT()) return inv_ntt(deg);
        assert(this->size() >= 1 && (*this)[0] != 0);
        if(deg < 0) deg = (int)this->size();
        P res({T(1) / (*this)[0]}); 
        for(int d = 1; d < deg; d <<= 1){
            res = (res + res - res * res * pre(d << 1)).pre(d << 1);
        }
        res.resize(deg);
        return res;
    }
    P log(int deg = -1) const {
        assert(this->size() >= 1 && (*this)[0] == 1);
        if(count_terms() <= SPARSE_BOARDER) return log_sparse(deg);
        if(deg < 0) deg = (int)this->size();
        return ((diff() * inv(deg)).pre(deg - 1)).integral();
    }
    P log_sparse(int deg = -1) const {
        assert(this->size() >= 1 && (*this)[0] == 1);
        if(deg < 0) deg = (int)this->size();
        vector<pair<int, T>> dat;
        for(int i = 1; i < (int)this->size(); ++i) if((*this)[i] != T(0)){
            dat.emplace_back(i, (*this)[i]);
        }
        BiCoef<T> bc(deg);
        P res(deg), tmp(deg);
        for(int i = 0; i < deg - 1; ++i){
            T r = T(i + 1) * (*this)[i + 1];
            for(auto &&[k, val]: dat){
                if(k > i) break;
                r -= val * tmp[i - k];
            }
            tmp[i] = r;
            res[i + 1] = r * bc.inv(i + 1);
        }
        return res;
    }
    P exp_ntt(int deg = -1) const {
        if((int)this->size() == 0) return {T(1)};
        assert((*this)[0] == 0);
        if(deg < 0) deg = (int)this->size();
        P fiv;
        fiv.reserve(deg + 1);
        fiv.emplace_back(T(0));
        fiv.emplace_back(T(1));
        auto integral_in = [&](P &F) -> void {
            const int n = (int)F.size();
            auto mod = T::mod();
            while((int)fiv.size() <= n){
                int i = fiv.size();
                fiv.emplace_back((-fiv[mod % i]) * (mod / i));
            }
            F.insert(begin(F), T(0));
            for(int i = 1; i <= n; ++i) F[i] *= fiv[i];
        };
        auto diff_in = [&](P &F) -> void {
            if(F.empty()) return;
            F.erase(begin(F));
            T coef = 1;
            for(int i = 0; i < (int)F.size(); ++i){
                F[i] *= coef;
                coef += 1;
            }
        };
        P b{1, (1 < (int)this->size() ? (*this)[1] : 0)}, c{1}, z1, z2{1, 1};
        for(int m = 2; m < deg; m <<= 1){
            auto y = b;
            y.resize(m * 2);
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
            z2.resize(m * 2);
            NTT::ntt(z2);
            P x(begin(*this), begin(*this) + min((int)this->size(), m));
            diff_in(x);
            x.emplace_back(T(0));
            NTT::ntt(x);
            for(int i = 0; i < m; ++i) x[i] *= y[i];
            NTT::intt(x);
            x -= b.diff();
            x.resize(m * 2);
            for(int i = 0; i < m - 1; ++i) x[m + i] = x[i], x[i] = T(0);
            NTT::ntt(x);
            for(int i = 0; i < m * 2; ++i) x[i] *= z2[i];
            NTT::intt(x);
            x.pop_back();
            integral_in(x);
            for(int i = m; i < min((int)this->size(), m * 2); ++i) x[i] += (*this)[i];
            fill(begin(x), begin(x) + m, T(0));
            NTT::ntt(x);
            for(int i = 0; i < m * 2; ++i) x[i] *= y[i];
            NTT::intt(x);
            b.insert(end(b), begin(x) + m, end(x));
        }
        return P(begin(b), begin(b) + deg);
    }
    P exp_sparse(int deg = -1) const {
        if((int)this->size() == 0) return {T(1)};
        assert((*this)[0] == 0);
        if(deg < 0) deg = (int)this->size();
        vector<pair<int, T>> dat;
        for(int i = 1; i < (int)this->size(); ++i) if((*this)[i] != T(0)){
            dat.emplace_back(i - 1, (*this)[i] * i);
        }
        BiCoef<T> bc(deg);
        P res(deg);
        res[0] = 1;
        for(int i = 1; i < deg; ++i){
            T r = 0;
            for(auto &&[k, val]: dat){
                if(k > i - 1) break;
                r += val * res[i - k - 1];
            }
            res[i] = r * bc.inv(i);
        }
        return res;
    }
    P exp(int deg = -1) const {
        if((int)this->size() == 0) return {T(1)};
        if(T::isNTT()) return exp_ntt(deg);
        assert((*this)[0] == 0);
        if(deg < 0) deg = (int)this->size();
        P res(1, 1);
        for(int d = 1; d < deg; d <<= 1){
            res = res * (pre(d << 1) - res.log(d << 1) + 1).pre(d << 1);
        }
        res.resize(deg);
        return res;
    }
    P pow_sparse(long long e, int deg = -1) const {
        assert(e >= 0);
        if(deg < 0) deg = (int)this->size();
        if(deg == 0) return {};
        if(e == 0){
            P res(deg, 0);
            res[0] = 1;
            return res;
        }
        long long ord = 0;
        while(ord < (int)this->size() && (*this)[ord] == 0) ++ord;
        if(ord == (int)this->size() || ord > (deg - 1) / e) return P(deg, 0);
        if((*this)[0] == 1) return pow_sparse_constant1(e, deg);
        auto f = (*this);
        rotate(f.begin(), f.begin() + ord, f.end());
        T con = f[0], icon = f[0].inv();
        for(int i = 0; i < deg; ++i) f[i] *= icon;
        auto res = f.pow_sparse_constant1(e, deg);
        int ord2 = e * ord;
        rotate(res.begin(), res.begin() + (deg - ord2), res.end());
        fill(res.begin(), res.begin() + ord2, T(0));
        T pw = con.pow(e);
        for(int i = ord2; i < deg; ++i) res[i] *= pw;
        return res;
    }
    P pow_sparse_constant1(T e, int deg = -1) const {
        assert((int)this->size() > 0 && (*this)[0] == 1);
        if(deg < 0) deg = (int)this->size();
        vector<pair<int, T>> dat;
        for(int i = 1; i < (int)this->size(); ++i) if((*this)[i] != T(0)){
            dat.emplace_back(i, (*this)[i]);
        }
        BiCoef<T> bc(deg);
        P res(deg);
        res[0] = 1;
        for(int i = 0; i < deg - 1; ++i){
            T &r = res[i + 1];
            for(auto &&[k, val]: dat){
                if(k > i + 1) break;
                T t = val * res[i - k + 1];
                r += t * (T(k) * e - T(i - k + 1));
            }
            r *= bc.inv(i + 1);
        }
        return res;
    }
    P pow(long long e, int deg = -1) const {
        if(count_terms() <= SPARSE_BOARDER) return pow_sparse(e, deg);
        assert(e >= 0);
        if(deg < 0) deg = (int)this->size();
        if(deg == 0) return {};
        if(e == 0){
            P res(deg, 0);
            res[0] = 1;
            return res;
        }
        long long ord = 0;
        while(ord < (int)this->size() && (*this)[ord] == 0) ++ord;
        if(ord == (int)this->size() || ord > (deg - 1) / e) return P(deg, 0);
        T k = (*this)[ord];
        P res = ((((*this) >> ord) * k.inv()).log(deg) * e).exp(deg) * T(k).pow(e) << (e * ord);
        res.resize(deg);
        return res;
    }
    // if not exist, out {}.
    P sqrt(int deg = -1) const {
        if(count_terms() <= SPARSE_BOARDER) return sqrt_sparse(deg);
        if(deg < 0) deg = (int)this->size();
        if((int)this->size() == 0) return P(deg, 0);
        if((*this)[0] == T(0)) {
            for(int i = 1; i < (int)this->size(); ++i){
                if((*this)[i] != T(0)){
                    if(i & 1) return {};
                    if(deg - i / 2 <= 0) return P(deg, 0);
                    auto res = ((*this) >> i).sqrt(deg - i / 2);
                    if (res.empty()) return {};
                    res = res << (i / 2);
                    if ((int)res.size() < deg) res.resize(deg, T(0));
                    return res;
                }
            }
            return P(deg, 0);
        }
        long long sqr;
        while(1){
            long long a = (*this)[0].val() % T::mod();
            if(a < 0) a += T::mod();
            if(a <= 1){
                sqr = a;
                break;
            }
            int p = T::mod();
            if(T(a).pow((p - 1) >> 1) != 1){
                sqr = -1;
                break;
            }
            T b = 1, one = 1;
            while(b.pow((p - 1) >> 1) == 1) b += 1;
            long long m = p - 1, e = 0;
            while(!(m & 1)) m >>= 1, ++e;
            T x = T(a).pow((m - 1) >> 1);
            T y = T(a) * x * x;
            x *= a;
            T z = T(b).pow(m);
            while(y != 1){
                long long j = 0;
                T t = y;
                while(t != one){
                    ++j;
                    t *= t;
                }
                z = z.pow((long long)(1) << (e - j - 1));
                x *= z, z *= z, y *= z;
                e = j;
            }
            long long res = x.val();
            if(res * 2 > p) res = p - res;
            sqr = res;
            break;
        }
        if(sqr == -1) return {};
        assert((*this)[0].val() == sqr * sqr % T::mod());
        P res = {T(sqr)};
        T iv2 = T(2).inv();
        for(int d = 1; d < deg; d <<= 1){
            res = (res + pre(d << 1) * res.inv(d << 1)).pre(d << 1) * iv2;
        }
        res.resize(deg);
        return res;
    }
    P sqrt_sparse(int deg) const {
        if(deg < 0) deg = (int)this->size();
        if((int)this->size() == 0) return P(deg, 0);
        if((*this)[0] == T(0)){
            for(int i = 1; i < (int)this->size(); ++i){
                if((*this)[i] != T(0)){
                    if(i & 1) return P();
                    if(deg - i / 2 <= 0) return P(deg, 0);
                    auto res = ((*this) >> i).sqrt_sparse(deg - i / 2);
                    if(res.empty()) return P();
                    res = res << (i / 2);
                    if((int)res.size() < deg) res.resize(deg);
                    return res;
                }
            }
            return P(deg, 0);
        }
        T con = (*this)[0], icon = con.inv();
        long long sqr;
        while(1){
            long long a = (*this)[0].val() % T::mod();
            if(a < 0) a += T::mod();
            if(a <= 1){
                sqr = a;
                break;
            }
            int p = T::mod();
            if(T(a).pow((p - 1) >> 1) != 1){
                sqr = -1;
                break;
            }
            T b = 1, one = 1;
            while(b.pow((p - 1) >> 1) == 1) b += 1;
            long long m = p - 1, e = 0;
            while(!(m & 1)) m >>= 1, ++e;
            T x = T(a).pow((m - 1) >> 1);
            T y = T(a) * x * x;
            x *= a;
            T z = T(b).pow(m);
            while(y != 1){
                long long j = 0;
                T t = y;
                while(t != one){
                    ++j;
                    t *= t;
                }
                z = z.pow((long long)(1) << (e - j - 1));
                x *= z, z *= z, y *= z;
                e = j;
            }
            long long res = x.val();
            if(res * 2 > p) res = p - res;
            sqr = res;
            break;
        }
        if(sqr == -1) return P();
        assert(con.val() == sqr * sqr % T::mod());
        auto res = (*this) * icon;
        return res.sqrt_sparse_constant1(deg) * sqr;
    } 
    P sqrt_sparse_constant1(int deg) const {
        return pow_sparse_constant1(T(2).inv(), deg);
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

template<typename T>
using FPS = FormalPowerSeries<T>;

template<typename T>
FPS<T> composition(FPS<T> g, FPS<T> f, int deg){
    using NTT = Convolution_base<T>;
    auto rec = [&](auto &&rec, FPS<T> Q, int n, int h, int k) -> FPS<T> {
        if(n == 0){
            FPS<T> nT(begin(Q), begin(Q) + k);
            nT.emplace_back(T(1));
            FPS<T> u = g * nT.rev().inv().rev();
            FPS<T> P(h * k);
            for(int i = 0; i < (int)g.size(); ++i) P[k - i - 1] = u[i + k];
            return P;
        }
        FPS<T> nQ(h * k * 4), nR(h * k * 2);
        for(int i = 0; i < k; ++i){
            copy(begin(Q) + i * h, begin(Q) + i * h + n + 1, begin(nQ) + i * h * 2);
        }
        nQ[h * k * 2] += 1;
        NTT::ntt(nQ);
        for(int i = 0; i < h * k * 4; i += 2) swap(nQ[i], nQ[i + 1]);
        for(int i = 0; i < h * k * 2; ++i) nR[i] = nQ[i * 2] * nQ[i * 2 + 1];
        NTT::intt(nR);
        nR[0] -= 1;
        Q.assign(h * k, 0);
        for(int i = 0; i < k * 2; ++i) for(int j = 0; j <= n / 2; ++j){
            Q[i * h / 2 + j] = nR[i * h + j];
        }
        auto P = rec(rec, Q, n / 2, h / 2, k * 2);
        FPS<T> nP(h * k * 4);
        for(int i = 0; i < k * 2; ++i) for(int j = 0; j <= n / 2; ++j){
            nP[i * h * 2 + j * 2 + n % 2] = P[i * h / 2 + j];
        }
        NTT::ntt(nP);
        for(int i = 1; i < h * k * 4; i <<= 1) reverse(begin(nQ) + i, begin(nQ) + i * 2);
        for(int i = 0; i < h * k * 4; ++i) nP[i] *= nQ[i];
        NTT::intt(nP);
        P.assign(h * k, 0);
        for(int i = 0; i < k; ++i){
            copy(begin(nP) + i * h * 2, begin(nP) + i * h * 2 + n + 1, begin(P) + i * h);
        }
        return P;
    };
    if (deg == -1) deg = max((int)f.size(), (int)g.size());
    f.resize(deg), g.resize(deg);
    int n = (int)f.size() - 1, h = 1, k = 1;
    while (h < n + 1) h *= 2;
    FPS<T> Q(h * k);
    for (int i = 0; i <= n; i++) Q[i] = -f[i];
    FPS<T> P = rec(rec, Q, n, h, k);
    return P.pre(n + 1).rev();
}
template<typename T>
FPS<T> composition(FPS<T> g, FPS<T> f){
    return composition(g, f, max(f.size(), g.size()));
}

template<class T, int MOD = T::mod(), int pr = calc_primitive_root(MOD)>
FPS<T> power_projection(FPS<T> f, FPS<T> g = {1}, int m = -1) {
    using NTT = Convolution_base<T>;
    int n = (int)f.size() - 1, k = 1, h = 1;
    g.resize(n + 1);
    if (m < 0) m = n;
    while (h < n + 1) h <<= 1;
    FPS<T> P((n + 1) * k), Q((n + 1) * k), nP, nQ, buf, buf2;
    for (int i = 0; i <= n; i++) P[i * k] = g[i];
    for (int i = 0; i <= n; i++) Q[i * k] = -f[i];
    Q[0] += 1;
    mint iv2 = T(2).inv();
    while (n) {
        mint w = T(pr).pow((MOD - 1) / (2 * k)), iw = w.inv();
        buf2.resize(k);
        auto ntt_doubling = [&]() {
            copy(begin(buf), end(buf), begin(buf2));
            NTT::intt(buf2);
            mint c = 1;
            for (int i = 0; i < k; i++) buf2[i] *= c, c *= w;
            NTT::ntt(buf2);
            copy(begin(buf2), end(buf2), back_inserter(buf));
        };
        nP.clear(), nQ.clear();
        for (int i = 0; i <= n; i++) {
            buf.resize(k);
            copy(begin(P) + i * k, begin(P) + (i + 1) * k, begin(buf));
            ntt_doubling();
            copy(begin(buf), end(buf), back_inserter(nP));
            buf.resize(k);
            copy(begin(Q) + i * k, begin(Q) + (i + 1) * k, begin(buf));
            if (i == 0) {
                for (int j = 0; j < k; j++) buf[j] -= 1;
                ntt_doubling();
                for (int j = 0; j < k; j++) buf[j] += 1;
                for (int j = 0; j < k; j++) buf[k + j] -= 1;
            } else {
                ntt_doubling();
            }
            copy(begin(buf), end(buf), back_inserter(nQ));
        }
        nP.resize(h * 2 * k * 2), nQ.resize(h * 2 * k * 2);
        FPS<T> p(h * 2), q(h * 2);
        w = T(pr).pow((MOD - 1) / (h * 2)), iw = w.inv();
        vector<int> btr;
        if (n % 2) {
            btr.resize(h);
            for (int i = 0, lg = __builtin_ctz(h); i < h; i++) {
                btr[i] = (btr[i >> 1] >> 1) + ((i & 1) << (lg - 1));
            }
        }
        for (int j = 0; j < k * 2; j++) {
            p.assign(h * 2, 0), q.assign(h * 2, 0);
            for (int i = 0; i < h; i++) p[i] = nP[i * k * 2 + j], q[i] = nQ[i * k * 2 + j];
            NTT::ntt(p), NTT::ntt(q);
            for (int i = 0; i < h * 2; i += 2) swap(q[i], q[i + 1]);
            for (int i = 0; i < h * 2; i++) p[i] *= q[i];
            for (int i = 0; i < h; i++) q[i] = q[i * 2] * q[i * 2 + 1];
            if (n & 1) {
                T c = iv2;
                buf.resize(h);
                for (int i : btr) buf[i] = (p[i * 2] - p[i * 2 + 1]) * c, c *= iw;
                swap(p, buf);
            } else {
                for (int i = 0; i < h; i++) p[i] = (p[i * 2] + p[i * 2 + 1]) * iv2;
            }
            p.resize(h), q.resize(h);
            NTT::intt(p), NTT::intt(q);
            for (int i = 0; i < h; i++) nP[i * k * 2 + j] = p[i];
            for (int i = 0; i < h; i++) nQ[i * k * 2 + j] = q[i];
        }
        nP.resize((n / 2 + 1) * k * 2), nQ.resize((n / 2 + 1) * k * 2);
        swap(P, nP), swap(Q, nQ);
        n /= 2, h /= 2, k *= 2;
    }
    FPS<mint> S{begin(P), begin(P) + k}, ST{begin(Q), begin(Q) + k};
    NTT::intt(S), NTT::intt(ST);
    ST[0] -= 1;
    if (ST[0] == 0) return S.rev().pre(m + 1);
    else return (S.rev() * (ST + (FPS<mint>{1} << k)).rev().inv(m + 1)).pre(m + 1);
}

template<class T>
FPS<T> composition_inverse(FPS<T> f, int deg = -1) {
    assert((int)f.size() >= 2 && f[1] != 0);
    if (deg == -1) deg = (int)f.size();
    if (deg < 2) return FPS<T>{0, f[1].inv()}.pre(deg);
    int n = deg - 1;
    FPS<T> h = power_projection(f) * n;
    for (int k = 1; k <= n; k++) h[k] /= k;
    h = h.rev(), h *= h[0].inv();
    FPS<T> g = (h.log() * mint(-n).inv()).exp();
    g *= f[1].inv();
    return (g << 1).pre(deg);
}

template<typename T>
FPS<T> middle_product(const FPS<T> &a, const FPS<T> &b){
    assert(a.size() >= b.size());
    if(b.empty()) return FPS<T>((int)a.size() - (int)b.size() + 1);
    int N = 1;
    while(N < (int)a.size()) N <<= 1;
    FPS<T> fa(N), fb(N);
    copy(a.begin(), a.end(), fa.begin());
    copy(b.rbegin(), b.rend(), fb.begin());
    fa *= fb;
    fa.resize(a.size());
    fa.erase(fa.begin(), fa.begin() + (int)b.size() - 1);
    return fa;
}

template<typename T>
struct SubProductTree{
    int num_points, siz;
    vector<FPS<T>> tree;
    SubProductTree(){}
    SubProductTree(const vector<T> &x){
        num_points = (int)x.size();
        siz = 1;
        while(siz < num_points) siz <<= 1;
        tree.resize(siz * 2);
        for(int i = 0; i < siz; ++i) tree[siz + i] = {1, (i < num_points ? -x[i]: 0)};
        for(int i = siz - 1; i >= 1; --i) tree[i] = tree[i * 2] * tree[i * 2 + 1];
    }
    vector<T> eval(FPS<T> f){
        int N = (int)f.size();
        if(N == 0) return vector<T>(num_points, T(0));
        f.resize(N * 2 - 1);
        vector<FPS<T>> g(siz * 2);
        g[1] = tree[1];
        g[1].resize(N);
        g[1] = g[1].inv();
        g[1] = middle_product(f, g[1]);
        g[1].resize(siz);
        for(int i = 1; i < siz; ++i){
            g[i * 2] = middle_product(g[i], tree[i * 2 + 1]);
            g[i * 2 + 1] = middle_product(g[i], tree[i * 2]);
        }
        vector<T> res(num_points);
        for(int i = 0; i < num_points; ++i) res[i] = g[siz + i][0];
        return res;
    }
};

template<typename T>
vector<T> multipoint_eval(FPS<T> &f, vector<T> &x){
    if(x.empty()) return {};
    SubProductTree<T> st(x);
    return st.eval(f);
}

template<typename T>
vector<T> multipoint_eval_geometric(const FPS<T> &f, const T &a, const T &r, int M){
    auto calc = [&](const T &r, int m) -> FPS<T> {
        FPS<T> res(m, T(1));
        T po = 1;
        for(int i = 0; i < m - 1; ++i) res[i + 1] = res[i] * po, po *= r;
        return res;
    };
    int N = (int)f.size();
    if(M == 0) return vector<T>();
    if(r == T(0)){
        vector<T> res(M);
        for(int i = 1; i < M; ++i) res[i] = f[0];
        res[0] = f(a);
        return res;
    }
    if(min(N, M) < 60){
        vector<T> res(M);
        T b = a;
        for(int i = 0; i < M; ++i) res[i] = f(b), b *= r;
        return res;
    }
    FPS<T> res = f;
    T po = 1;
    for(int i = 0; i < N; ++i) res[i] *= po, po *= a;
    FPS<T> A = calc(r, N + M - 1), B = calc(r.inv(), max(N, M));
    for(int i = 0; i < N; ++i) res[i] *= B[i];
    res = middle_product(A, res);
    for(int i = 0; i < M; ++i) res[i] *= B[i];
    return res;
}
