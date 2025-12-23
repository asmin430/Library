vector<int> Z_algorithm(string S){
    int l = 0, r = 0;
    int N = S.size();
    vector<int> res(N, 0);
    res[0] = N;
    for(int i = 1; i < N; ++i){
        if(res[i - l] < r - i){
            res[i] = res[i - l];
        }else{
            r = max(r, i);
            while(r < N && S[r] == S[r - i]){
                ++r;
            }
            res[i] = r - i;
            l = i;
        }
    }
    return res;
}