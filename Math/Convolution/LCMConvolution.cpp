template<typename T>
void devisor_zeta_transform(vector<T> &f){
    int n = f.size();
    vector<int> is_prime(n, 1);
    for(int p = 2; p < n; ++p){
        if(is_prime[p]){
            for(int q = p * 2; q < n; q += p) is_prime[q] = 0;
            for(int j = 1; j * p < n; ++j) f[j * p] += f[j];
        }
    }
}

template<typename T>
void devisor_moebius_transform(vector<T> &f){
    int n = f.size();
    vector<int> is_prime(n, 1);
    for(int p = 2; p < n; ++p){
        if(is_prime[p]){
            for(int q = p * 2; q < n; q += p) is_prime[q] = 0;
            for(int j = (n - 1) / p; j > 0; --j) f[j * p] -= f[j];
        }
    }
}

template<typename T>
vector<T> LCMConvolution(vector<T> f, vector<T> g){
    assert(f.size() == g.size());
    devisor_zeta_transform(f);
    devisor_zeta_transform(g);
    int n = f.size();
    vector<T> h(n, 0);
    for(int i = 1; i < n; ++i) h[i] = f[i] * g[i];
    devisor_moebius_transform(h);
    return h;
};
