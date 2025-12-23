template<typename TreeDPInfo>
struct LinkCutTreeSubtree{
    using Point = typename TreeDPInfo::Point;
    using Path = typename TreeDPInfo::Path;
    using Info = typename TreeDPInfo::Info;
    private:
    struct Node{
        int id;
        Node *l, *r, *p;
        Info info;
        Path sum, mus;
        Point point;
        bool rev;
        bool is_root() const { return (!p || (p->l != this && p->r != this)); }
        Node(const Info &info, int id): info(info), id(id), l(nullptr), r(nullptr), p(nullptr), rev(false), point(Point::id()){}
    };
    public:
    using NP = Node *;
    vector<NP> nodes;
    private:
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
    void push(NP t){
        if(t->rev) {
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
        Path key = TreeDPInfo::add_vertex(t->point, t->info);
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
            if(cur->r){
                cur->point = TreeDPInfo::rake(cur->point, TreeDPInfo::add_edge(cur->r->sum));
            }
            cur->r = rp;
            if(cur->r){
                cur->point = TreeDPInfo::rake(cur->point, TreeDPInfo::add_edge(cur->r->sum).inv());
            }
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
    const Path &query(NP u){
        evert(u);
        return u->sum;
    }
    const Path &query_path(NP u){
        expose(u);
        return u->sum;
    }
    const Path &query_path(NP u, NP v){
        evert(u);
        return query_path(v);
    }
    Path query_subtree(NP u){
        expose(u);
        NP l = u->l;
        u->l = nullptr;
        update(u);
        auto ret = u->sum;
        u->l = l;
        update(u);
        return ret;
    }
    Path query_subtree(NP r, NP u){
        evert(r);
        return query_subtree(u);
    }
    public:
    LinkCutTreeSubtree() = default;
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
    const Path query_path(int u){
        assert(u >= 0 && u < (int)nodes.size());
        return query_path(nodes[u]);
    }
    const Path query_path(int u, int v){
        assert(u >= 0 && u < (int)nodes.size());
        assert(v >= 0 && v < (int)nodes.size());
        return query_path(nodes[u], nodes[v]);
    }
    const Path query_subtree(int u){
        assert(u >= 0 && u < (int)nodes.size());
        return query_subtree(nodes[u]);
    }
    const Path query_subtree(int r, int u){
        assert(r >= 0 && r < (int)nodes.size());
        assert(u >= 0 && u < (int)nodes.size());
        return query_subtree(nodes[r], nodes[u]);
    }

};

struct TreeDPInfo {
    struct Point {
        ll sum;
        static constexpr Point id(){return {0LL};}
        Point inv() const { return {-sum}; }
    };
    struct Path{ll sum;};
    struct Info{ll sum;};
    static Path add_vertex(const Point& d, const Info& u){return {d.sum + u.sum}; }
    static Point add_edge(const Path& d){return {d.sum};}
    static Point rake(const Point& s, const Point& p){ return {s.sum + p.sum}; }
    static Path compress(const Path& p, const Path& c){ return {p.sum + c.sum}; }
};