template<typename T>
void multiple_zeta_transform(vector<T> &f){
    int n = f.size();
    for(int k = 1; k < n; ++k){
        for(int i = k * 2; i < n; i += k){
            f[k] += f[i];
        }
    }
}

template<typename T>
void multiple_moebius_transform(vector<T> &f){
    int n = f.size();
    for(int k = n - 1; k >= 1; --k){
        for(int i = k * 2; i < n; i += k){
            f[k] -= f[i];
        }
    }
}

template<typename T>
vector<T> GCDConvolution(vector<T> f, vector<T> g){
    assert(f.size() == g.size());
    multiple_zeta_transform(f);
    multiple_zeta_transform(g);
    int n = f.size();
    vector<T> h(n, 0);
    for(int i = 1; i < n; ++i) h[i] = f[i] * g[i];
    multiple_moebius_transform(h);
    return h;
};
