template<typename T>
struct BinaryIndexedTree{
private:
    int n;
    vector<T> dat;
public
    BinaryIndexedTree() = default;
    explicit BinaryIndexedTree(int n): n(n){dat.assign(n + 1, T());}
    explicit BinaryIndexedTree(const vector<T> &v): BinaryIndexedTree((int)v.size()){
        build(v);
    }
    void build(const vector<T> &v){
        assert(n == (int)v.size());
        for(int i = 1; i <= n; ++i) dat[i] = v[i - 1];
        for(int i = 1; i <= n; ++i){
            int j = i + (i & -i);
            if(j <= n) dat[j] += dat[i];
        }
    }
    void add(int k, const T &x){
        for(++k; k <= n; k += k & -k) dat[k] += x;
    }
    T sum(int r) const {
        T ret = T();
        for(; r > 0; r -= r & -r) ret += dat[r];
        return ret;
    }
    T sum(int l, int r) const {return sum(r) - sum(l);}
    int lower_bound(T x) const {
        int i = 0;
        for(int k = 1 << (__lg(n) + 1); k > 0; k >>= 1){
            if(i + k <= n && dat[i + k] < x){
                x -= dat[i + k];
                i += k;
            }
        }
        return i;
    }
    int upper_bound(T x) const {
        int i = 0;
        for(int k = 1 << (__lg(n) + 1); k > 0; k >>= 1){
            if(i + k <= n && dat[i + k] <= x){
                x -= dat[i + k];
                i += k;
            }
        }
        return i;
    }
};
