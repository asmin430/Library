template <typename T, auto op, auto e>
struct SegmentTree{
    vector<T> dat;
    int n, lg, sz;
    
    SegmentTree(){}
    SegmentTree(int n){build(vector<T> (n, e())); }
    SegmentTree(const vector<T> &v){build(v); }
    void build(const vector<T> &v){
        n = (int)v.size();
        lg = 1;
        while((1 << lg) < n) ++lg;
        sz = (1 << lg);
        dat.resize(sz << 1);
        for(int i = 0; i < n; ++i) dat[sz + i] = v[i];
        for(int i = sz - 1; i > 0; --i) update(i);
    }
    
    // i 番目の要素
    T operator[](int &i){return dat[sz + i]; }
    void update(int i){dat[i] = op(dat[2 * i], dat[2 * i + 1]); }
    // i 番目の要素を x に変更
    void set(int i, const T &x){
        assert(i < n && i >= 0);
        dat[i += sz] = x;
        while(i >>= 1) update(i);
    }
    // [l, r) のクエリ 
    T query(int l, int r){
        assert(0 <= l && l <= r && r <= n);
        T vl = e(), vr = e();
        l += sz;
        r += sz;
        while(l < r){
            if(l & 1) vl = op(vl, dat[l++]);
            if(r & 1) vr = op(dat[--r], vr);
            l >>= 1, r >>= 1;
        }
        return op(vl, vr);
    }
    // 全体クエリ
    T all_query(){return dat[1]; }
    // [l, x) が check を満たす最大の x
    template<class F>
    int max_right(F check, int l){
        assert(0 <= l && l <= n && check(e()));
        if(l == n) return n;
        l += sz;
        T sm = e();
        do{
            while(!(l & 1)) l >>= 1;
            if(!check(op(sm, dat[l]))){
                while(l < sz){
                    l = 2 * l;
                    if(check(op, dat[l])) sm = op(sm, dat[l++]);
                }
                return l - sz;
            }
            sm = op(sm, dat[l++]);
        }while((l & -l) != l);
        return n;
    }
    // [x, r) が check を満たす最小の x
    template<class F>
    int min_left(F check, int r){
        assert(0<= r && r <= n && check(e()));
        if(r == 0) return 0;
        r += sz;
        T sm = e();
        do{
            --r;
            while(r > 1 && (r & 1)) r >>= 1;
            if(!check(op(dat[r], sm))){
                while(r < sz){
                    r = 2 * r + 1;
                    if(check(op(dat[r], sm))) sm = op(dat[r], sm);
                }
                return r + 1 - sz;
            }
        }while((r & -r) != r);
        return 0;
    }
    //可換のときだけ A[i xor xor_val] (l <= i < r) のクエリを出力
    T xor_query(int l, int r, int xor_val){
        T x = e();
        for(int k = 0; k < lg + 1; ++k){
            if(l >= r) break;
            if(l & 1) x = op(x, dat[(sz >> k) + ((l++) ^ xor_val)]);
            if(r & 1) x = op(x, dat[(sz >> k) + ((--r) ^ xor_val)]);
            l >>= 1, r >>= 1, xor_val >>= 1;
        }
        return x;
    }
};
