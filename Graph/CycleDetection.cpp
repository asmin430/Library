template<class T>
struct Edge{
    int from, to;
    T val;
    Edge(int f = -1, int t = -1, T v = -1): from(f), to(t), val(v){}
};

template<class T> using Graph = vector<vector<Edge<T>>>;

template<class T>
struct CycleDetection{
    Graph<T> G;
    vector<bool> seen, finished;
    vector<Edge<T>> history;
    CycleDetection(){}
    CycleDetection(const Graph<T> &g){ init(g); }
    void init(const Graph<T> &g){
        G = g;
        seen.assign(G.size(), false);
        finished.assign(G.size(), false);
    }

    int dfs(int v, Edge<T> e, bool t = true){
        seen[v] = true;
        history.push_back(e);
        for(const Edge<T> &e2: G[v]){
            if(t && e2.to == e.from) continue;
            if(finished[e2.to]) continue;
            if(seen[e2.to] && !finished[e2.to]){
                history.push_back(e2);
                return e2.to;
            }

            int pos = dfs(e2.to, e2, t);
            if(pos != -1) return pos;
        }
        finished[v] = true;
        history.pop_back();
        return -1;
    }

    vector<Edge<T>> reconstruct(int pos){
        vector<Edge<T>> cycle;
        while(!history.empty()){
            const Edge<T> &e = history.back();
            cycle.push_back(e);
            history.pop_back();
            if(e.from == pos) break;
        }

        reverse(begin(cycle), end(cycle));
        return cycle;
    }

    vector<Edge<T>> detect(bool t = true){
        int pos = -1;
        for(int v = 0; v < (int)G.size() && pos == -1; ++v){
            if(seen[v]) continue;
            history.clear();
            pos = dfs(v, Edge<T>(), t);
            if(pos != -1) return reconstruct(pos);
        }
        return vector<Edge<T>>();
    }
};