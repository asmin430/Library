template<typename T, auto op, auto e>
struct SegmentTree2D{
    int H, W;
    vector<T> dat;
    SegmentTree2D(){}
    SegmentTree2D(int H, int W): H(H), W(W), dat(4 * H * W, e()){}
    SegmentTree2D(const vector<vector<T>> &v){
        H = v.size(), W = (H == 0 ? 0 : v[0].size());
        dat.assign(4 * H * W, e());
        for(int i = 0; i < H; ++i) for(int j = 0; j < W; ++j) dat[idx(H + i, W + j)] = v[i][j];
        for(int j = W; j < 2 * W; ++j) for(int i = H - 1; i >= 0; --i) dat[idx(i, j)] = op(dat[idx(2 * i, j)], dat[idx(2 * i + 1, j)]);
        for(int i = 0; i < 2 * H; ++i) for(int j = W - 1; j >= 0; --j) dat[idx(i, j)] = op(dat[idx(i, 2 * j)], dat[idx(i, 2 * j + 1)]);
    }
    //(i, j) に x をセット
    void set(int i, int j, T x){
        i += H, j += W;
        dat[idx(i, j)] = x;
        int k = i;
        while(k >>= 1) dat[idx(k, j)] = op(dat[idx(2 * k, j)], dat[idx(2 * k + 1, j)]);
        k = i;
        while(k){
            int l = j;
            while(l >>= 1) dat[idx(k, l)] = op(dat[idx(k, 2 * l)], dat[idx(k, 2 * l + 1)]);
            k >>= 1;
        }
    }
    // [il, ir), [jl, jr) のクエリ 
    T query(int il, int ir, int jl, int jr){
        assert(0 <= il && il <= ir && ir <= H);
        assert(0 <= jl && jl <= jr && jr <= W);
        T res = e();
        il += H, ir += H;
        while(il < ir){
            if(il & 1) res = op(res, query_sub(il++, jl, jr));
            if(ir & 1) res = op(res, query_sub(--ir, jl, jr));
            il >>= 1, ir >>= 1;
        }
        return res;
    }

private:
    inline int idx(int x, int y){ return x * 2 * W + y; }
    T query_sub(int i, int jl, int jr){
        T res = e();
        jl += W, jr += W;
        while(jl < jr){
            if(jl & 1) res = op(res, dat[idx(i, jl++)]);
            if(jr & 1) res = op(res, dat[idx(i, --jr)]);
            jl >>= 1, jr >>= 1;
        }
        return res;
    }
};
