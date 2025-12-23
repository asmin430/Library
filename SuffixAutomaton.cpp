struct SuffixAutomaton{
    struct Node{
        unordered_map<char, int> nxt;
        int link, len;
    };
    vector<Node> nodes;
    int last;
    SuffixAutomaton(){
        nodes.push_back({{}, -1, 0});
        last = 0;
    }
    void push(char c){
        int new_node = (int)nodes.size();
        nodes.push_back({{}, -1, nodes[last].len + 1});
        int p = last;
        while(p != -1 && nodes[p].nxt.find(c) == nodes[p].nxt.end()){
            nodes[p].nxt[c] = new_node;
            p = nodes[p].link;
        }
        int q = (p == -1 ? 0 : nodes[p].nxt[c]);
        if(p == -1 || nodes[p].len + 1 == nodes[q].len){
            nodes[new_node].link = q;
        }else{
            int new_q = (int)nodes.size();
            nodes.push_back({nodes[q].nxt, nodes[q].link, nodes[p].len + 1});
            nodes[q].link = new_q;
            nodes[new_node].link = new_q;
            while(p != -1 && nodes[p].nxt[c] == q){
                nodes[p].nxt[c] = new_q;
                p = nodes[p].link;
            }
        }
        last = new_node;
    }
};