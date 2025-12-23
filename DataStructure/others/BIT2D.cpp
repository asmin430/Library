template <typename T>
struct BIT2D { // 0-indexed
    int H, W;
    vector<vector<T>> bit;
    BIT2D(int H_, int W_) { init(H_, W_); }
    void init(int H_, int W_) {
        H = H_ + 1;
        W = W_ + 1;
        bit.assign(H, vector<T>(W, 0));
    }

    //(h, w) に x 加算
    void add(int h, int w, T x) {
        for (int i = h + 1; i < H; i += (i & -i)) {
            for (int j = w + 1; j < W; j += (j & -j)) {
                bit[i][j] += x;
            }
        }
    }

    T sum_sub(int h, int w) {
        T s(0);
        for (int i = h; i > 0; i -= (i & -i)) {
            for (int j = w; j > 0; j -= (j & -j)) {
                s += bit[i][j];
            }
        }
        return s;
    }

    // h1 ≦ i < h2 かつ w1 ≦ j < w2 の区間和
    T sum(int h1, int w1, int h2, int w2) {
        return sum_sub(h2, w2) - sum_sub(h2, w1) - sum_sub(h1, w2) + sum_sub(h1, w1);
    }
};