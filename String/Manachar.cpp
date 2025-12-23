vector<int> manachar(string S){
    vector<int> res(S.size(), 0);
    int i = 0, j = 0;
    while(i < S.size()){
        while(i - j >= 0 && i + j < S.size() && S[i - j] == S[i + j]) ++j;
        res[i] = j;
        int k = 1;
        while(i - k >= 0 && k + res[i - k] < j) res[i + k] = res[i - k], ++k;
        i += k;
        j -= k;
    }
    return res;
}