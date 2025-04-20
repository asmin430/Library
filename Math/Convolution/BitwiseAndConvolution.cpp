template<typename T>
void superset_zeta_transform(vector<T> &f){
    const int n = (int)f.size();
    assert((n & (n - 1)) == 0);
    for(int i = 1; i < n; i <<= 1){
        for(int j = 0; j < n; j += (i << 1)){
            for(int k = 0; k < i; ++k){
                f[j + k] += f[j + k + i];
            }
        }
    }
}
template<typename T>
void superset_moebius_transform(vector<T> &f){
    const int n = (int)f.size();
    assert((n & (n - 1)) == 0);
    for(int i = 1; i < n; i <<= 1){
        for(int j = 0; j < n; j += (i << 1)){
            for(int k = 0; k < i; ++k){
                f[j + k] -= f[j + k + i];
            }
        }
    }
}

template<typename T>
vector<T> BitwiseAndConvolution(vector<T> f, vector<T> g){
    const int n = (int)f.size();
    assert(f.size() == g.size());
    assert((n & (n - 1)) == 0);
    superset_zeta_transform(f);
    superset_zeta_transform(g);
    for(int i = 0; i < n; ++i) f[i] *= g[i];
    superset_moebius_transform(f);
    return f;
}
