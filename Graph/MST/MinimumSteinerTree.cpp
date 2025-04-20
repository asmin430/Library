template<class T>
struct MinimumSteinerTree{
private:
    struct edge{
        int to;
        T cost;

        edge(){}
        edge(int _to, T _cost): to(_to), cost(_cost){}
    };
    int V;
    vector<vector<edge>> G;
    vector<vector<T>> dp;
    const T infty = numeric_limits<T>::max() / 10;

public:
    // node_size = 頂点数
    MinimumSteinerTree(int node_size): V(node_size), G(V){}

    // u, v にコスト cost の辺を追加
    void add_edge(int u, int v, T cost){
        G[u].push_back(edge(v, cost));
        G[v].push_back(edge(u, cost));
    }


    vector<vector<T>> solve(vector<int> &terminal){
        int t = (int)terminal.size();
        assert(t > 0);
        dp.clear();
        dp.resize((1 << t), vector<T>(V, infty));
        for(int i = 0; i < t; ++i){
            assert(0 <= terminal[i] && terminal[i] < V);
            dp[(1 << i)][terminal[i]] = 0;
        }
        priority_queue<pair<T, int>, vector<pair<T, int>>, greater<pair<T, int>>> que;
        for(int i = 1; i < (1 << t); ++i){
            for(int j = 0; j < V; ++j){
                for(int k = i; k > 0; k = (k - 1) & i){
                    if(dp[k][j] != infty && dp[i ^ k][j] != infty){
                        T nxt = dp[k][j] + dp[i ^ k][j];
                        if(nxt < dp[i][j]){
                            dp[i][j] = nxt;
                        }
                    }
                }
            }
            for(int j = 0; j < V; ++j){
                if(dp[i][j] != infty) que.push(make_pair(dp[i][j], j));
            }
            while(!que.empty()){
                auto [c, v] = que.top();
                que.pop();
                if(dp[i][v] < c) continue;
                for(int j = 0; j < (int)G[v].size(); ++j){
                    const auto &e = G[v][j];
                    if(dp[i][e.to] > c + e.cost){
                        dp[i][e.to] = c + e.cost;
                        que.push(make_pair(dp[i][e.to], e.to));
                    }
                }
            }
        }
        return dp;
    }
};
