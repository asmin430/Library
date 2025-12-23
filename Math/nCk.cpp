template<class T>
struct nCk{
    vector<T> fac;
    vector<T> finv;
    vector<T> inv;

    nCk(int M): fac(M), finv(M), inv(M){
        const int MOD = T::mod();
        fac[0] = fac[1] = 1;
        finv[0] = finv[1] = 1;
        inv[1] = 1;
        for(int i = 2; i < M; ++i){
            fac[i] = fac[i - 1] * i;
            inv[i] = (T)MOD - inv[MOD % i] * (MOD / i);
            finv[i] = finv[i - 1] * inv[i];
        }
    }

    T C(int n, int k){
        if(n < k) return 0;
        if(n < 0 || k < 0) return 0;
        return fac[n] * finv[k] * finv[n - k];
    }
};