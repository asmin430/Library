template<class Type>
class PotentialUnionFind{
public:
    PotentialUnionFind() = default;

    explicit PotentialUnionFind(size_t n) : m_parentsOrSize(n, -1), m_diffWeights(n){}

    //i の代表元
    int root(int i){
        if(m_parentsOrSize[i] < 0) return i;
        const int r = root(m_parentsOrSize[i]);
        m_diffWeights[i] += m_diffWeights[m_parentsOrSize[i]];
        return (m_parentsOrSize[i] = r);
    }

    //a, b を重み w で結ぶ
    void merge(int a, int b, Type w){
        w += weight(a);
        w -= weight(b);
        a = root(a);
        b = root(b);
        if(a != b){
            if(-m_parentsOrSize[a] < -m_parentsOrSize[b]){
                swap(a, b);
                w = -w;
            }
            m_parentsOrSize[a] += m_parentsOrSize[b];
            m_parentsOrSize[b] = a;
            m_diffWeights[b] = w;
        }
    }

    //(b の重み) - (a の重み)
    Type diff(int a, int b){
        return (weight(b) - weight(a));
    }

    //a, b が連結か
    bool same(int a, int b){
        return (root(a) == root(b));
    }

    //i を含む連結成分の要素数
    int size(int i){
        return -m_parentsOrSize[root(i)];
    }

private:
    vector<int> m_parentsOrSize;
    vector<Type> m_diffWeights;
    Type weight(int i){
        root(i);
        return m_diffWeights[i];
    }
};
