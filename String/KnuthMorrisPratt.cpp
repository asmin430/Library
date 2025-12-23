int KMP(string S, string T){
    vector<int> table(T.size(), 0);
    table[0] = -1;
    int j = -1;
    for(int i = 0; i < T.size() - 1; ++i){
        while(j >= 0 && T[i] != T[j]){
            j = table[j];
        }
        table[i + 1] = j + 1;
        ++j;
    }
    int i = 0;
    j = 0;
    while(i + j < S.size()){
        if(S[i + j] == T[j]){
            ++j;
            if(j == T.size()){
                return i;
            }
        }else{
            i = i + j - table[j];
            if(j > 0){
                j = table[j];
            }
        }
    }
    return -1;
}