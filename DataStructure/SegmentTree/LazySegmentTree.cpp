template <typename T, auto op, auto e, typename F, auto mapping, auto composition, auto id>
struct LazySegmentTree{
    int n, lg, sz;
    vector<T> dat;
    vector<F> laz;
    LazySegmentTree(){}
    LazySegmentTree(int n){build(vector<T>(n, e())); }
    LazySegmentTree(const vector<T> &v){build(v); }
    void build(const vector<T> &v){
        n = (int)v.size(), lg = 1;
        while((1 << lg) < n) ++lg;
        sz = (1 << lg);
        dat.assign(sz << 1, e());
        laz.assign(sz << 1, id());
        for(int i = 0; i < n; ++i) dat[sz + i] = v[i];
        for(int i = sz - 1; i >= 1; --i) update(i);
    }
    void update(int i){dat[i] = op(dat[2 * i], dat[2 * i + 1]); }
    // [p] に x をセット
    void set(int p, T x){
        assert(0 <= p && p < n);
        p += sz;
        for(int i = lg; i >= 1; --i) push(p >> i);
        dat[p] = x;
        for(int i = 1; i <= lg; ++i) update(p >> i);
    }
    T operator[](int p){
        assert(0 <= p && p < n);
        p += sz;
        for(int i = lg; i >= 1; --i) push(p >> i);
        return dat[p];
    }
    // [l, r) のクエリ
    T query(int l, int r){
        assert(0 <= l && l <= r && r <= n);
        if(l == r) return e();
        l += sz, r += sz;
        for(int i = lg; i >= 1; --i){
            if(((l >> i) << i) != l) push(l >> i);
            if(((r >> i) << i) != r) push((r - 1) >> i);
        }
        T xl = e(), xr = e();
        while(l < r){
            if(l & 1) xl = op(xl, dat[l++]);
            if(r & 1) xr = op(dat[--r], xr);
            l >>= 1, r >>= 1;
        }
        return op(xl, xr);
    }
    T all_query(){return dat[1]; }
    // [l, r) に f を作用
    void apply(int l, int r, F f){
        assert(0 <= l && l <= r && r <= n);
        if(l == r) return;
        l += sz, r += sz;
        for(int i = lg; i >= 1; --i){
            if(((l >> i) << i) != l) push(l >> i);
            if(((r >> i) << i) != r) push((r - 1) >> i);
        }
        int l2 = l, r2 = r;
        while(l < r){
            if(l & 1) apply_sub(l++, f);
            if(r & 1) apply_sub(--r, f);
            l >>= 1, r >>= 1;
        }
        l = l2, r = r2;
        for(int i = 1; i <= lg; ++i){
            if(((l >> i) << i) != l) update(l >> i);
            if(((r >> i) << i) != r) update((r - 1) >> i);
        }
    }
    template <typename C>
    int max_right(const C check, int l){
        assert(0 <= l && l <= n);
        assert(check(e()));
        if(l == n) return n;
        l += sz;
        for(int i = lg; i >= 1; --i) push(l >> i);
        T sm = e();
        do{
            while(!(l & 1)) l >>= 1;
            if(!check(op(sm, dat[l]))){
                while(l < sz){
                    push(l);
                    l <<= 1;
                    if(check(op(sm, dat[l]))) sm = op(sm, dat[l++]);
                }
                return l - sz;
            }
            sm = op(sm, dat[l++]);
        }while((l & -l) != l);
        return n;
    }
    template <typename C>
    int min_left(const C check, int r){
        assert(0 <= r && r <= n);
        assert(check(e()));
        if(r == 0) return 0;
        r += sz;
        for(int i = lg; i >= 1; --i) push((r - 1) >> i);
        T sm = e();
        do{
            r--;
            while(r > 1 && (r & 1)) r >>= 1;
            if(!check(op(dat[r], sm))){
                while(r < sz){
                    push(r);
                    r = (2 * r + 1);
                    if(check(op(dat[r], sm))) sm = op(dat[r--], sm);
                }
                return r + 1 - sz;
            }
            sm = op(dat[r], sm);
        }while((r & -r) != r);
        return 0;
    }
private:
    void apply_sub(int k, F f){
        dat[k] = mapping(dat[k], f);
        if(k < sz) laz[k] = composition(laz[k], f);
    }
    void push(int k){
        if(laz[k] == id()) return;
        apply_sub(2 * k, laz[k]);
        apply_sub(2 * k + 1, laz[k]);
        laz[k] = id();
    }
};






