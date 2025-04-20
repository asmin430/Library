struct LowLink{
    private:
    int n;
    vector<vector<int>> G;
    vector<bool> re;
    vector<int> ord;
    vector<int> low;
    
    public:
    vector<int> articulation; // 関節点
    vector<pair<int, int>> bridge; // 橋

    LowLink(const vector<vector<int>> &G): G(G), n(G.size()){
        re.assign(n, 0);
        ord.assign(n, 0);
        low.assign(n, 0);
    }

    void build(){
        int cnt = 0;
        for(int i = 0; i < n; ++i){
            if(!re[i]){
                dfs(i, -1, cnt);
            }
        }
    }

    private:
    void dfs(int cur, int par, int &cnt){
        re[cur] = true;
        ord[cur] = cnt;
        low[cur] = ord[cur];
        ++cnt;
        bool is_art = false;
        int son_cnt = 0;
        for(auto nxt: G[cur]){
            if(!re[nxt]){
                ++son_cnt;
                dfs(nxt, cur, cnt);
                if(nxt != par) low[cur] = min(low[cur], low[nxt]);
                if(par != -1 && ord[cur] <= low[nxt]) is_art = true;
                if(ord[cur] < low[nxt]) bridge.push_back(make_pair(min(cur, nxt), max(cur, nxt)));
            }else{
                if(nxt != par) low[cur] = min(low[cur], ord[nxt]);
            }
        }
        if(par == -1 && son_cnt >= 2) is_art = true;
        if(is_art) articulation.push_back(cur);
    }
};
