template <typename T, auto composition, auto id>
struct DualSegmentTree{
    int n, lg, sz;
    vector<T> laz;
    DualSegmentTree(){}
    DualSegmentTree(int n){build(vector<T>(n, id())); }
    DualSegmentTree(const vector<T> &v){build(v); }
    void build(const vector<T> &v){
        n = (int)v.size();
        lg = 1;
        while((1 << lg) < n) ++lg;
        sz = (1 << lg);
        laz.assign(sz << 1, id());
        for(int i = 0; i < n; ++i) laz[sz + i] = v[i];
    }
    T operator[](int p){
        assert(0 <= p && p < n);
        p += sz;
        for(int i = lg; i >= 1; --i) push(p >> i);
        return laz[p];
    }
    void set(int p, T x){
        assert(0 <= p && p < n);
        p += sz;
        for(int i = lg; i >= 1; --i) push(p >> i);
        laz[p] = x;
    }
    void apply(int l, int r, const T &x){
        assert(0 <= l && l <= r && r <= n);
        if(l == r) return;
        l += sz, r += sz;
        for(int i = lg; i >= 1; --i){
            if(((l >> i) << i) != l) push(l >> i);
            if(((r >> i) << i) != r) push((r - 1) >> i);
        }
        while(l < r){
            if(l & 1) apply_sub(l++, x);
            if(r & 1) apply_sub(--r, x);
            l >>= 1, r >>= 1;
        }
    }
private:
    void push(int k){
        if(laz[k] == id()) return;
        apply_sub(2 * k, laz[k]);
        apply_sub(2 * k + 1, laz[k]);
        laz[k] = id();
    }
    void apply_sub(int k, T x){laz[k] = composition(laz[k], x); }
};
