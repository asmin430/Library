class Mersenne61modint{
    static const long long md = (1LL << 61) - 1;
    long long _v;

    inline unsigned hi() const noexcept {return _v >> 31;}
    inline unsigned lo() const noexcept {return _v & ((1LL << 31) - 1);}

public:
    static long long mod() {return md;}

    Mersenne61modint() : _v(0){}

    explicit Mersenne61modint(long long x) : _v(x >= md ? x - md : x){}

    long long val() const noexcept {return _v;}

    Mersenne61modint operator+(const Mersenne61modint &x) const {
        return Mersenne61modint(_v + x._v);
    }

    Mersenne61modint operator-(const Mersenne61modint &x) const {
        return Mersenne61modint(_v + md - x._v);
    }

    Mersenne61modint operator*(const Mersenne61modint &x) const {
        using ull = unsigned long long;

        ull uu = (ull)hi() * x.hi() * 2;
        ull ll = (ull)lo() * x.lo();
        ull lu = (ull)hi() * x.lo() + (ull)lo() * x.hi();

        ull sum = uu + ll + ((lu & ((1ULL << 30) - 1)) << 31) + (lu >> 30);
        ull reduced = (sum >> 61) + (sum & ull(md));
        return Mersenne61modint(reduced);
    }

    Mersenne61modint pow(long long n) const {
        assert(n >= 0);
        Mersenne61modint ans(1), tmp = *this;
        while(n){
            if(n & 1) ans *= tmp;
            tmp *= tmp;
            n >>= 1;
        }
        return ans;
    }

    Mersenne61modint inv() const {return pow(md - 2);}

    Mersenne61modint operator/(const Mersenne61modint &x) const {return *this * x.inv();}

    Mersenne61modint operator-() const {return Mersenne61modint(md - _v);}
    Mersenne61modint operator+=(const Mersenne61modint &x){return *this = *this + x;}
    Mersenne61modint operator-=(const Mersenne61modint &x){return *this = *this - x;}
    Mersenne61modint operator*=(const Mersenne61modint &x){return *this = *this * x;}
    Mersenne61modint operator/=(const Mersenne61modint &x){return *this = *this / x;}
    
    Mersenne61modint operator+(unsigned &x) const {return Mersenne61modint(this->_v + x);}

    bool operator==(const Mersenne61modint &x) const {return _v == x._v;}
    bool operator!=(const Mersenne61modint &x) const {return _v != x._v;}
    bool operator<(const Mersenne61modint &x) const {return _v < x._v;}
};