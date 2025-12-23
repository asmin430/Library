// depend on NTTFriendlyFormalPowerSeries.cpp
template<typename mint> mint BostanMori(const FPS<mint> &P, const FPS<mint>& Q, long long N){
    assert(!P.empty() && !Q.empty());
    if(N == 0 || Q.size() == 1){
        return P[0] / Q[0];
    }
    int qdeg = (int)Q.size();
    FPS<mint> P2{P}, minusQ{Q};
    P2.resize(qdeg - 1);
    for(int i = 1; i < (int)Q.size(); i += 2) minusQ[i] = -minusQ[i];
    P2 *= minusQ;
    FPS<mint> Q2 = Q * minusQ;
    FPS<mint> S(qdeg - 1), T(qdeg);
    for(int i = 0; i < (int)S.size(); ++i){
        if(!(N & 1)) S[i] = P2[i * 2];
        else S[i] = P2[i * 2] + P2[i * 2 + 1];
    }
    for(int i = 0; i < (int)T.size(); ++i){
        T[i] = Q2[i * 2];
    }
    return BostanMori(S, T, N >> 1);
}