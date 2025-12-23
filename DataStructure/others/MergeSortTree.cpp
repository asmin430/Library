template<typename T>
struct MergeSortTree{

    int N;
    vector<vector<T>> x;
    vector<vector<T>> sum;
    vector<T> basis;

    MergeSortTree() = default;

    MergeSortTree(const vector<T> &vec) : N(vec.size()){
        x.resize(N * 2);
        sum.resize(N * 2);
        basis.resize(N + 1);
        basis[0] = 0;
        for(int i = 0; i < N; ++i) basis[i + 1] = basis[i] + vec[i];
        for(int i = 0; i < N; ++i) x[N + i] = {vec[i]};
        for(int i = 0; i < N; ++i) sum[N + i] = {0, vec[i]};
        for(int i = N - 1; i; --i){
            merge(all(x[2 * i]), all(x[2 * i + 1]), back_inserter(x[i]))
            sort(all(x[i]));
            sum[i].resize(x[i].size() + 1);
            sum[i][0] = 0;
            for(int j = 0; j < x[i].size(); ++j){
                sum[i][j + 1] = sum[i][j] + x[i][j];
            }
        }
    }

    // [l, r) の query 未満の要素の個数を出力．
    int countless(int l, int r, T query) const {
        l += N, r += N;
        int ret = 0;
        while(l < r){
            if(l & 1){
                ret += lower_bound(all(x[l]), query) - x[l].begin();
                ++l;
            }
            if(r & 1){
                --r;
                ret += lower_bound(all(x[r]), query) - x[r].begin();
            }
            l >>= 1;
            r >>= 1;
        }
        return ret;
    }

    // [l, r) の query より大きい要素の個数を出力．
    int countmore(int l, int r, T query) const {
        int tot = max(0, min(r, N) - max(l, 0));
        return tot - countless(l, r, query + 1);
    }

    // [l, r) の query 未満の要素の総和を出力．
    T sumless(int l, int r, T query) const {
        l += N, r += N;
        T ret = 0;
        while(l < r){
            if(l & 1){
                int s = lower_bound(all(x[l]), query) - x[l].begin();
                ret += sum[l][s];
                ++l;
            }
            if(r & 1){
                --r;
                int s = lower_bound(all(x[r]), query) - x[r].begin();
                ret += sum[r][s];
            }
            l >>= 1;
            r >>= 1;
        }
        return ret;
    }

    // [l, r) の query より大きい要素の総和を出力．
    T summore(int l, int r, T query) const {
        T tot = basis[r] - basis[l - 1];
        return tot - sumless(l, r, query + 1);
    }
};

