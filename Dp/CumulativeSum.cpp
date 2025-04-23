template<class T>
struct CumulativeSum{
    vector<T> dat;
    CumulativeSum() = default;
    explicit CumulativeSum(size_t sz): dat(sz + 1, 0){}
    explicit CumulativeSum(vector<T> &v): dat(v.size() + 1, 0){
        for(int i = 0; i < v.size(); ++i){
            dat[i + 1] = v[i];
        }
        build();
    }
    void add(int k, const T &x){ dat[k + 1] += x; }
    void build(){
        for(int i = 1; i < dat.size(); ++i){
            dat[i] += dat[i - 1];
        }
    }
    T sum(int r) const {
        if(r < 0) return 0;
        return dat[min(r, (int)dat.size() - 1)];
    }

    //[l, r) の総和
    T sum(int l, int r) const {
        return sum(r) - sum(l);
    }
};

