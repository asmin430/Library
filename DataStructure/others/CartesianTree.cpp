template<typename T, bool IS_MIN>
struct CartesianTree{
    int n;
    vector<T> &A;
    vector<pair<int, int>> range;
    vector<int> lch, rch, par;
    int root;

    // コンストラクタ
    CartesianTree(vector<T> &A): n(A.size()), A(A){
        range.assign(n, make_pair(-1, -1));
        lch.assign(n, -1);
        rch.assign(n, -1);
        par.assign(n, -1);
        if(n == 1){
            range[0] = make_pair(0, 1);
            root = 0;
            return;
        }
        auto is_sm = [&](int i, int j) -> bool {
            if(IS_MIN) return (A[i] < A[j]) || (A[i] == A[j] && i < j);
            return (A[i] > A[j]) || (A[i] == A[j] && i < j);
        };
        vector<int> st;
        for(int i = 0; i < n; ++i){
            while(!st.empty() && is_sm(i, st.back())){
                lch[i] = st.back();
                st.pop_back();
            }
            range[i].first = (st.empty() ? 0 : st.back() + 1);
            st.push_back(i);
        }
        st.clear();
        for(int i = n - 1; i >= 0; --i){
            while(!st.empty() && is_sm(i, st.back())){
                rch[i] = st.back();
                st.pop_back();
            }
            range[i].second = (st.empty() ? n : st.back());
            st.push_back(i);
        }
        for(int i = 0; i < n; ++i) if(lch[i] != -1) par[lch[i]] = i;
        for(int i = 0; i < n; ++i) if(rch[i] != -1) par[rch[i]] = i;
        for(int i = 0; i < n; ++i) if(par[i] == -1) root = i;
    }
    // i 番目の要素を値としたときの子孫の幅
    tuple<int, int, T> maximum_rectangle(int i){
        auto [l, r] = range[i];
        return make_tuple(l, r, A[i]);
    }
    // 子孫の幅と自身の値の積の最大値(実質的に最大の長方形)
    T maximum_rectangle_area(){
        assert(IS_MIN);
        T res = 0;
        for(int i = 0; i < n; ++i){
            auto [l, r, h] = maximum_rectangle(i);
            res = max(res, (r - l) * h);
        }
        return res;
    }

    ll count_subrestangle(bool baseline){
        assert(IS_MIN);
        ll res = 0;
        for(int i = 0; i < n; ++i){
            auto [l, r, h] = maximum_rectangle(i);
            ll x = (baseline ? h : h * (h + 1) / 2);
            res += x * (i - l + 1) * (r - i);
        }
        return res;
    }
};
