template<typename TreeDPInfo>
struct LinkCutTree{
    using Path = typename TreeDPInfo::Path;
    using Info = typename TreeDPInfo::Info;
    private:
    struct Node{
        int id;
        Node *l, *r, *p;
        Info info;
        Path sum, mus;
        bool rev;
        bool is_root() const { return (!p || (p->l != this && p->r != this)); }
        Node(const Info &info, int id): id(id), info(info), l(nullptr), r(nullptr), p(nullptr), rev(false){}
    };
    public:
    using NP = Node *;
    private:
    vector<NP> nodes;
    void toggle(NP t){
        swap(t->l, t->r);
        swap(t->sum, t->mus);
        t->rev ^= true;
    }
    void rotr(NP t){
        NP x = t->p, y = x->p;
        push(x), push(t);
        if((x->l = t->r)) t->r->p = x;
        t->r = x, x->p = t;
        update(x), update(t);
        if((t->p = y)){
            if(y->l == x) y->l = t;
            if(y->r == x) y->r = t;
        }
    }
    void rotl(NP t){
        NP x = t->p, y = x->p;
        push(x), push(t);
        if((x->r = t->l)) t->l->p = x;
        t->l = x, x->p = t;
        update(x), update(t);
        if((t->p = y)){
            if(y->l == x) y->l = t;
            if(y->r == x) y->r = t;
        }
    }
    private:
    void push(NP t){
        if(t->rev){
            if(t->l) toggle(t->l);
            if(t->r) toggle(t->r);
            t->rev = false;
        }
    }
    void push_rev(NP t){
        if(t->rev){
            if(t->l) toggle(t->l);
            if(t->r) toggle(t->r);
            t->rev = false;
        }
    }
    void update(NP t){
        Path key = TreeDPInfo::vertex(t->info);
        t->sum = key;
        t->mus = key;
        if(t->l){
            t->sum = TreeDPInfo::compress(t->l->sum, t->sum);
            t->mus = TreeDPInfo::compress(t->mus, t->l->mus);
        }
        if(t->r){
            t->sum = TreeDPInfo::compress(t->sum, t->r->sum);
            t->mus = TreeDPInfo::compress(t->r->mus, t->mus);
        }
    }
    void splay(NP t){
        push(t);
        while(not t->is_root()){
            NP q = t->p;
            if(q->is_root()){
                push_rev(q), push_rev(t);
                if(q->l == t) rotr(t);
                else rotl(t);
            }else{
                NP r = q->p;
                push_rev(r), push_rev(q), push_rev(t);
                if(r->l == q){
                    if(q->l == t) rotr(q), rotr(t);
                    else rotl(t), rotr(t);
                }else{
                    if(q->r == t) rotl(q), rotl(t);
                    else rotr(t), rotl(t);
                }
            }
        }
    }
    NP expose(NP t){
        NP rp = nullptr;
        for(NP cur = t; cur; cur = cur->p){
            splay(cur);
            cur->r = rp;
            update(cur);
            rp = cur;
        }
        splay(t);
        return rp;
    }
    void link(NP child, NP parent){
        if(is_connected(child, parent)){
            throw runtime_error("child and parent must be different connected components");
        }
        if(child->l){
            throw runtime_error("child must be root");
        }
        child->p = parent;
        parent->r = child;
        update(parent);
    }
    void cut(NP child){
        expose(child);
        NP parent = child->l;
        if(not parent){
            throw runtime_error("child must not be root");
        }
        child->l = nullptr;
        parent->p = nullptr;
        update(child);
    }
    void evert(NP t){
        expose(t);
        toggle(t);
        push(t);
    }
    bool is_connected(NP u, NP v){
        expose(u), expose(v);
        return u == v or u->p;
    }
    vector<NP> build(vector<Info> &vs){
        vector<NP> nodes(vs.size());
        for(int i = 0; i < (int)vs.size(); i++){
            nodes[i] = alloc(vs[i]);
        }
        return nodes;
    }
    NP lca(NP u, NP v){
        if(not is_connected(u, v)) return nullptr;
        expose(u);
        return expose(v);
    }
    void set_key(NP t, const Info &v){
        expose(t);
        t->info = std::move(v);
        update(t);
    }
    const Path &query_path(NP u){
        expose(u);
        return u->sum;
    }
    const Path &query_path(NP u, NP v){
        evert(u);
        return query_path(v);
    }
    template <typename C>
    pair<NP, Path> find_first(NP u, const C &check){
        expose(u);
        Path sum = TreeDPInfo::vertex(u->info);
        if(check(sum)) return {u, sum};
        u = u->l;
        while(u){
            push(u);
            if(u->r){
                Path nxt = TreeDPInfo::compress(u->r->sum, sum);
                if(check(nxt)){
                    u = u->r;
                    continue;
                }
                sum = nxt;
            }
            Path nxt = TreeDPInfo::compress(TreeDPInfo::vertex(u->info), sum);
            if(check(nxt)){
                splay(u);
                return {u, nxt};
            }
            sum = nxt;
            u = u->l;
        }
        return {nullptr, sum};
    }
    public:
    LinkCutTree() = default;
    void exopse(int u){
        assert(u >= 0 && u < (int)nodes.size());
        expose(nodes[u]);
    }
    void link(int u, int v){
        assert(u >= 0 && u < (int)nodes.size());
        assert(v >= 0 && v < (int)nodes.size());
        evert(nodes[u]);
        link(nodes[u], nodes[v]);
    }
    void cut(int u){
        assert(u >= 0 && u < (int)nodes.size());
        cut(nodes[u]);
    }
    void evert(int u){
        assert(u >= 0 && u < (int)nodes.size());
        evert(nodes[u]);
    }
    void alloc(const Info &info){
        NP t = new Node(info, nodes.size());
        update(t);
        nodes.push_back(t);
    }
    bool is_connected(int u, int v){
        assert(u >= 0 && u < (int)nodes.size());
        assert(v >= 0 && v < (int)nodes.size());
        return is_connected(nodes[u], nodes[v]);
    }
    int lca(int u, int v){
        assert(u >= 0 && u < (int)nodes.size());
        assert(v >= 0 && v < (int)nodes.size());
        NP lca = lca(nodes[u], nodes[v]);
        if(lca == nullptr) return -1;
        return lca->id;
    }
    void set_key(int u, const Info &v){
        assert(u >= 0 && u < (int)nodes.size());
        set_key(nodes[u], v);
    }
    const Path &query_path(int u){
        assert(u >= 0 && u < (int)nodes.size());
        return query_path(nodes[u]);
    }
    const Path &query_path(int u, int v){
        assert(u >= 0 && u < (int)nodes.size());
        assert(v >= 0 && v < (int)nodes.size());
        return query_path(nodes[u], nodes[v]);
    }
    template<typename C>
    pair<int, Path> find_first(int u, const C &check){
        assert(u >= 0 && u < (int)nodes.size());
        auto [t, path] = find_first(nodes[u], check);
        if(t == nullptr) return {-1, path};
        return {t->id, path};
    }
};


struct TreeDPInfo {
    struct Path{ mint a, b; };
    struct Info{ mint a, b; };
    static Path vertex(const Info& u){ return {u.a, u.b};}
    static Path compress(const Path& p, const Path& c){ return {c.a * p.a, c.a * p.b + c.b}; }
};
