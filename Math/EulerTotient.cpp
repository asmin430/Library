vector<int> euler_totielnt_table(int N){
    vector<int> ret(N + 1);
    iota(begin(ret), end(ret), 0);
    for(int p = 2; p <= N; ++p){
        if(ret[p] == p){
            ret[p] = p - 1;
            for(int i = p * 2; i <= N; i += p){
                ret[i] = ret[i] / p * (p - 1);
            }
        }
    }
    return ret;
}

template<typename T>
T euler_totielnt_function(T n){
    T ret = n;
    if(n % 2 == 0){
        ret /= 2;
        do{n /= 2; }while(n % 2 == 0);
    }
    if(n % 3 == 0){
        ret = ret / 3 * 2;
        do{n /= 3; }while(n % 3 == 0);
    }
    for(T i = 5; i * i <= n; i += 4){
        if(n % i == 0){
            ret = ret / i * (i - 1);
            do{n /= i; }while(n % i == 0);
        }
    }
    if(n != 1) ret = ret / n * (n - 1);
    return ret;
}