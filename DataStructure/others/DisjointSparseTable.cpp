template<typename T, typename F>
struct DisjointSparseTable{
    const F f;
    vector<vector<T>> st;
    vector<int> lookup;
    DisjointSparseTable(const vector<T> &v, const F &f): f(f){
        int b = 0;
        while((1 << b) <= v.size()) ++b;
        st.resize(b, vector<T>(v.size(), T()));
        for(int i = 0; i < v.size(); ++i) st[0][i] = v[i];
        for(int i = 1; i < b; ++i){
            int shift = (1 << i);
            for(int j = 0; j < v.size(); j += (shift << 1)){
                int t = min(j + shift, (int)v.size());
                st[i][t - 1] = v[t - 1];
                for(int k = t - 2; k >= j; --k) st[i][k] = f(v[k], st[i][k + 1]);
                if(v.size() <= t) break;
                st[i][t] = v[t];
                int r = min(t + shift, (int)v.size());
                for(int k = t + 1; k < r; ++k) st[i][k] = f(st[i][k - 1], v[k]);
            }
        }
        lookup.resize(1 << b);
        for(int i = 2; i < lookup.size(); ++i){
            lookup[i] = lookup[i >> 1] + 1;
        }
    }

    T query(int l, int r){
        if(l >= --r) return st[0][l];
        int p = lookup[l ^ r];
        return f(st[p][l], st[p][r]);
    }
};


