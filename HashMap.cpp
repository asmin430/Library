template<typename Val>
struct HashMap{
    HashMap(unsigned int n = 0){build(n); }
    void build(unsigned int n){
        unsigned int k = 8;
        while(k < n * 2){
            k *= 2;
        }
        cap = k / 2, msk = k - 1;
        key.resize(k), val.resize(k), used.assign(k, 0);
    }
    void clear(){used.assign(used.size(), 0); }
    int size(){return used.size() - cap; }
    int index(const unsigned long long& k){
        int i = 0;
        for(i = hash(k); used[i] && key[i] != k; i = (i + 1) & msk){
        }
        return i;
    }
    Val& operator[](const unsigned long long& k){
        if(cap == 0) extend();
        int i = index(k);
        if(!used[i]){used[i] = 1, key[i] = k, val[i] = Val{}, --cap; }
        return val[i];
    }
    Val get(const unsigned long long& k, Val default_value){
        int i = index(k);
        return (used[i] ? val[i]: default_value);
    }
    bool count(const unsigned long long& k){
        int i = index(k);
        return used[i] && key[i] == k;
    }
    template<typename F>
    void enumerate_all(F f){
        for(int i = 0; i < used.size(); ++i){
            if(used[i]){
                f(key[i], val[i]);
            }
        }
    }
private:
    unsigned int cap, msk;
    vector<unsigned long long> key;
    vector<Val> val;
    vector<bool> used;
    unsigned long long hash(unsigned long long x){
        static const unsigned long long FIXED_RANDOM = std::chrono::steady_clock::now().time_since_epoch().count();
        x += FIXED_RANDOM;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        return (x ^ (x >> 31)) & msk;
    }
    void extend(){
        vector<pair<unsigned long long, Val>> dat;
        dat.reserve(used.size() - cap);
        for(int i = 0; i < used.size(); ++i){
            if(used[i]) dat.emplace_back(key[i], val[i]);
        }
        build(2 * dat.size());
        for(auto& [a, b]: dat){
            (*this)[a] = b;
        }
    }
};