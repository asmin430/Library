// depend on Modint.cpp
template<typename Mint>
struct Convolution_base{
    static vector<Mint> roots, iroots, rate3, irate3;
    static int max_base;

    Convolution_base() = default;

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
vector<Mint> Convolution_base<Mint>::roots = vector<Mint>();
template <typename Mint>
vector<Mint> Convolution_base<Mint>::iroots = vector<Mint>();
template <typename Mint>
vector<Mint> Convolution_base<Mint>::rate3 = vector<Mint>();
template <typename Mint>
vector<Mint> Convolution_base<Mint>::irate3 = vector<Mint>();
template <typename Mint>
int Convolution_base<Mint>::max_base = 0;

template<typename T>
vector<T> ConvolutionNaive(const vector<T> &a, const vector<T> &b){
    const int n = a.size(), m = b.size();
    vector<T> c(n + m - 1);
    if(n < m){
        for(int j = 0; j < m; ++j) for(int i = 0; i < n; ++i) c[i + j] += a[i] * b[j];
    }else{
        for(int i = 0; i < n; ++i) for(int j = 0; j < m; ++j) c[i + j] += a[i] * b[j];
    }
    return c;
} 

template<typename T>
vector<T> Convolution(const vector<T> &a, const vector<T> &b){
    if(a.size() == 0 || b.size() == 0) return {};
    if(min(a.size(), b.size()) <= 120) return ConvolutionNaive(a, b);
    Convolution_base<T> conv;
    auto c = conv.multiply(a, b);
    return c;
}


constexpr pair<long long, long long> inv_gcd(long long a, long long b){
    a %= b;
    if(a < 0) a += b;
    if(a == 0) return {b, 0};
    long long s = b, t = a;
    long long m0 = 0, m1 = 1;
    while(t){
        long long u = s / t;
        s -= t * u;
        m0 -= m1 * u;
        long long tmp = s;
        s = t;
        t = tmp;
        tmp = m0;
        m0 = m1;
        m1 = tmp;
    }
    m0 %= (b / s);
    if(m0 < 0) m0 += b / s;
    return {s, m0};
}

template<typename T>
vector<T> ArbitraryConvolution(const vector<T> &a, const vector<T> &b){
    int p = T::mod() - 1;
    int cnt = 0;
    while(!(p & 1)){
        p >>= 1;
        ++cnt;
    }
    if(cnt >= 23){
        return Convolution(a, b);
    }
    if(a.size() == 0 || b.size() == 0) return {};
    if(min(a.size(), b.size()) <= 120) return ConvolutionNaive(a, b);
    int n = a.size(), m = b.size();
    static constexpr long long MOD1 = 754974721;
    static constexpr long long MOD2 = 167772161;
    static constexpr long long MOD3 = 469762049;
    static constexpr long long M1M2 = MOD1 * MOD2;
    static constexpr long long INV_M1_MOD2 = inv_gcd(MOD1, MOD2).second;
    static constexpr long long INV_M1M2_MOD3 = inv_gcd(M1M2, MOD3).second;
    using mint1 = MontgomeryModInt<MOD1>;
    using mint2 = MontgomeryModInt<MOD2>;
    using mint3 = MontgomeryModInt<MOD3>;
    vector<mint1> a21(n), b21(m);
    vector<mint2> a22(n), b22(m);
    vector<mint3> a23(n), b23(m);
    for(int i = 0; i < n; ++i){
        a21[i] = a[i].val();
        a22[i] = a[i].val();
        a23[i] = a[i].val();
    }
    for(int i = 0; i < m; ++i){
        b21[i] = b[i].val();
        b22[i] = b[i].val();
        b23[i] = b[i].val();
    }


    auto c1 = Convolution(a21, b21);
    auto c2 = Convolution(a22, b22);
    auto c3 = Convolution(a23, b23);
    const long long m1m2 = T(M1M2).val();
    vector<T> c(n + m - 1);
    for(int i = 0; i < n + m - 1; ++i){
        long long x1 = c1[i].val();
        long long x2 = (mint2(c2[i].val() - x1) * INV_M1_MOD2).val();
        long long x3 = (mint3(c3[i].val() - x1 - x2 * MOD1) * INV_M1M2_MOD3).val();
        c[i] = x1 + x2 * MOD1 + x3 * m1m2;
    }
    return c;   
}

vector<unsigned long long> Convolution_2_64(const vector<unsigned long long> &a, const vector<unsigned long long> &b){
    if(a.size() == 0 || b.size() == 0) return {};
    if(min(a.size(), b.size()) <= 120) return ConvolutionNaive(a, b);
    int n = a.size(), m = b.size();
    static constexpr long long MOD1 = 754974721;
    static constexpr long long MOD2 = 167772161;
    static constexpr long long MOD3 = 469762049;
    static constexpr long long MOD4 = 377487361;
    static constexpr long long MOD5 = 595591169;
    static constexpr long long MOD6 = 645922817;
    static constexpr long long M1_MOD2 = MOD1 % MOD2;
    static constexpr long long M1_MOD3 = MOD1 % MOD3;
    static constexpr long long M1_MOD4 = MOD1 % MOD4;
    static constexpr long long M1_MOD5 = MOD1 % MOD5;
    static constexpr long long M1_MOD6 = MOD1 % MOD6;
    static constexpr long long M1M2_MOD3 = M1_MOD3 * MOD2 % MOD3;
    static constexpr long long M1M2_MOD4 = M1_MOD4 * MOD2 % MOD4;
    static constexpr long long M1M2_MOD5 = M1_MOD5 * MOD2 % MOD5;
    static constexpr long long M1M2_MOD6 = M1_MOD6 * MOD2 % MOD6;
    static constexpr long long M1M2M3_MOD4 = M1M2_MOD4 * MOD3 % MOD4;
    static constexpr long long M1M2M3_MOD5 = M1M2_MOD5 * MOD3 % MOD5;
    static constexpr long long M1M2M3_MOD6 = M1M2_MOD6 * MOD3 % MOD6;
    static constexpr long long M1M2M3M4_MOD5 = M1M2M3_MOD5 * MOD4 % MOD5;
    static constexpr long long M1M2M3M4_MOD6 = M1M2M3_MOD6 * MOD4 % MOD6;
    static constexpr long long M1M2M3M4M5_MOD6 = M1M2M3M4_MOD6 * MOD5 % MOD6;
    static constexpr long long INV_M1_MOD2 = inv_gcd(M1_MOD2, MOD2).second;
    static constexpr long long INV_M1M2_MOD3 = inv_gcd(M1M2_MOD3, MOD3).second;
    static constexpr long long INV_M1M2M3_MOD4 = inv_gcd(M1M2M3_MOD4, MOD4).second;
    static constexpr long long INV_M1M2M3M4_MOD5 = inv_gcd(M1M2M3M4_MOD5, MOD5).second;
    static constexpr long long INV_M1M2M3M4M5_MOD6 = inv_gcd(M1M2M3M4M5_MOD6, MOD6).second;
    static constexpr unsigned long long M1 = MOD1;
    static constexpr unsigned long long M1M2 = M1 * MOD2;
    static constexpr unsigned long long M1M2M3 = M1M2 * MOD3;
    static constexpr unsigned long long M1M2M3M4 = M1M2M3 * MOD4;
    static constexpr unsigned long long M1M2M3M4M5 = M1M2M3M4 * MOD5;
    using mint1 = MontgomeryModInt<MOD1>;
    using mint2 = MontgomeryModInt<MOD2>;
    using mint3 = MontgomeryModInt<MOD3>;
    using mint4 = MontgomeryModInt<MOD4>;
    using mint5 = MontgomeryModInt<MOD5>;
    using mint6 = MontgomeryModInt<MOD6>;
    vector<mint1> a1(n), b1(m);
    vector<mint2> a2(n), b2(m);
    vector<mint3> a3(n), b3(m);
    vector<mint4> a4(n), b4(m);
    vector<mint5> a5(n), b5(m);
    vector<mint6> a6(n), b6(m);
    for(int i = 0; i < n; ++i){
        a1[i] = a[i] % MOD1;
        a2[i] = a[i] % MOD2;
        a3[i] = a[i] % MOD3;
        a4[i] = a[i] % MOD4;
        a5[i] = a[i] % MOD5;
        a6[i] = a[i] % MOD6;
    }
    for(int i = 0; i < m; ++i){
        b1[i] = b[i] % MOD1;
        b2[i] = b[i] % MOD2;
        b3[i] = b[i] % MOD3;
        b4[i] = b[i] % MOD4;
        b5[i] = b[i] % MOD5;
        b6[i] = b[i] % MOD6;
    }
    auto c1 = Convolution(a1, b1);
    auto c2 = Convolution(a2, b2);
    auto c3 = Convolution(a3, b3);
    auto c4 = Convolution(a4, b4);
    auto c5 = Convolution(a5, b5);
    auto c6 = Convolution(a6, b6);
    vector<unsigned long long> c(n + m - 1);
    for(int i = 0; i < n + m - 1; ++i){
        long long x1 = c1[i].val();
        long long x2 = (mint2((long long)c2[i].val() - x1) * INV_M1_MOD2).val();
        long long x3 = (mint3((long long)c3[i].val() - x1 - x2 * M1_MOD3) * INV_M1M2_MOD3).val();
        long long x4 = (mint4((long long)c4[i].val() - x1 - x2 * M1_MOD4 - x3 * M1M2_MOD4) * INV_M1M2M3_MOD4).val();
        long long x5 = (mint5((long long)c5[i].val() - x1 - x2 * M1_MOD5 - x3 * M1M2_MOD5 - x4 * M1M2M3_MOD5) * INV_M1M2M3M4_MOD5).val();
        long long x6 = (mint6((long long)c6[i].val() - x1 - x2 * M1_MOD6 - x3 * M1M2_MOD6 - x4 * M1M2M3_MOD6 - x5 * M1M2M3M4_MOD6) * INV_M1M2M3M4M5_MOD6).val();
        c[i] = x1 + x2 * M1 + x3 * M1M2 + x4 * M1M2M3 + x5 * M1M2M3M4 + x6 * M1M2M3M4M5;      
    }
    return c;
} 
