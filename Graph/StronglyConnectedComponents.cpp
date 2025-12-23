class StronglyConnectedComponents{
    private:
    int n;
    vector<vector<int>> G, rG;
    vector<int> order;
    vector<int> seen;

    public:
    vector<int> component;

    private:
    void dfs(int v){
        seen[v] = 1;
        for(auto nv: G[v]){
            if(seen[nv] == 0){
                dfs(nv);
            }
        }
        order.push_back(v);
    }

    void rdfs(int v, int k){
        component[v] = k;
        for(auto nv: rG[v]){
            if(component[nv] < 0){
                rdfs(nv, k);
            }
        }
    }

    public:
    StronglyConnectedComponents(int _n, vector<vector<int>> &_G, vector<vector<int>> &_rG){
        n = _n;
        G = _G;
        rG = _rG;
        component.assign(n, -1);
        seen.assign(n, 0);

        for(int v = 0; v < n; ++v){
            if(seen[v] == 0){
                dfs(v);
            }
        }

        int k = 0;
        reverse(begin(order), end(order));
        for(auto v: order){
            if(component[v] == -1){
                rdfs(v, k);
                ++k;
            }
        }
    }

    bool same(int u, int v){
        return component[u] == component[v];
    }

    int representative(int v){
        return component[v];
    }
};