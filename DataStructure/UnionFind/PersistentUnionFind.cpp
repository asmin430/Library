template <typename T, int LOG>
struct PersistentArray{
    struct Node{
        T dat;
        Node *child[1 << LOG] = {};

        Node() {}

        Node(const T &dat) : dat(dat) {}
    };

    Node *root;

    PersistentArray() : root(nullptr){}

    T get(Node *t, int k){
        if(k == 0) return t->dat;
        return get(t->child[k & ((1 << LOG) - 1)], k >> LOG);
    }

    T get(const int &k){ return get(root, k); }

    pair<Node *, T *> mutable_get(Node *t, int k){
        t = t ? new Node(*t) : new Node();
        if(k == 0) return {t, &t->dat};
        auto p = mutable_get(t->child[k & ((1 << LOG) - 1)], k >> LOG);
        t->child[k & ((1 << LOG) - 1)] = p.first;
        return {t, p.second};
    }

    T *mutable_get(const int &k){
        auto ret = mutable_get(root, k);
        root = ret.first;
        return ret.second;
    }

    Node *build(Node *t, const T &dat, int k){
        if(!t) t = new Node();
        if(k == 0){
            t->dat = dat;
            return t;
        }
        auto p = build(t->child[k & ((1 << LOG) - 1)], dat, k >> LOG);
        t->child[k & ((1 << LOG) - 1)] = p;
        return t;
    }
    void build(const vector<T> &v){
        root = nullptr;
        for(int i = 0; i < v.size(); i++){
            root = build(root, v[i], i);
        }
    }
};

struct PersistentUnionFind {
    PersistentArray<int, 3> dat;
    PersistentUnionFind(){}
    PersistentUnionFind(int sz){ dat.build(vector<int>(sz, -1)); }
    int find(int k){
        int p = dat.get(k);
        return p >= 0 ? find(p) : k;
    }
    int size(int k){ return (-dat.get(find(k))); }
    bool unite(int x, int y){
        x = find(x);
        y = find(y);
        if(x == y) return false;
        auto u = dat.get(x);
        auto v = dat.get(y);
        if(u < v){
            auto a = dat.mutable_get(x);
            *a += v;
            auto b = dat.mutable_get(y);
            *b = x;
        }else{
            auto a = dat.mutable_get(y);
            *a += u;
            auto b = dat.mutable_get(x);
            *b = y;
        }
        return true;
    }
};
