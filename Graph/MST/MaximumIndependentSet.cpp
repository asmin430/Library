struct MaximumIndependentSet{
    std::vector<long long> conn;
    int V;         // # of vertices
    int nret;      // Largest possible size of independent set
    long long ret; // Result is saved here: use (1), don't use (0)
    long long _avail;
    long long _tmp_state;

    void mis_dfs() {
        bool retry = true;
        std::stack<int> st;
        while (retry) {
            retry = false;
            for (int i = 0; i < V; i++)
                if (_avail & (1LL << i)) {
                    int nb = __builtin_popcountll(_avail & conn[i]);
                    if (nb <= 1) {
                        st.emplace(i), _avail -= 1LL << i, _tmp_state |= 1LL << i;
                        retry = true;
                        if (nb == 1) {
                            int j = __builtin_ctzll(_avail & conn[i]);
                            st.emplace(j), _avail &= ~(1LL << j);
                        }
                    }
                }
        }

        int t = __builtin_popcountll(_tmp_state);
        if (t > nret) nret = t, ret = _tmp_state;

        int d = -1, n = -1;
        for (int i = 0; i < V; i++)
            if (_avail & (1LL << i)) {
                int c = __builtin_popcountll(_avail & conn[i]);
                if (c > d) d = c, n = i;
            }

        if (d > 0) {
            long long nxt = _avail & conn[n];
            _avail -= 1LL << n;
            mis_dfs();
            _tmp_state |= 1LL << n;
            _avail &= ~nxt;
            mis_dfs();
            _avail |= nxt;
            _avail |= 1LL << n;
            _tmp_state &= ~(1LL << n);
        }

        while (st.size()) {
            int i = st.top();
            _avail |= 1LL << i;
            _tmp_state &= ~(1LL << i);
            st.pop();
        }
    }
    MaximumIndependentSet(const std::vector<std::vector<int>> &e)
        : conn(e.size()), V(e.size()), nret(-1), _avail((1LL << V) - 1), _tmp_state(0) {
        assert(V <= 63);
        for (int i = 0; i < V; i++)
            for (auto &j : e[i])
                if (i != j) conn[i] |= 1LL << j, conn[j] |= 1LL << i;
        mis_dfs();
    }
};