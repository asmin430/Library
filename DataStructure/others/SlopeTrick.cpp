template<typename T>
struct SlopeTrick{
    private:
    T min_f;
    priority_queue<T> L;
    priority_queue<T, vector<T>, greater<T>> R;
    const T INFT = numeric_limits<T>::max() / 3;
    T add_l, add_r;

    public:
    SlopeTrick(): min_f(T(0)), add_l(T(0)), add_r(T(0)){
        L.push(-INFT);
        R.push(INFT);
    }

    tuple<T, T, T> min_get(){
        return make_tuple(min_f, L.top() + add_l, R.top() + add_r);
    }
    void add_a(T a){
        min_f += a;
    }
    void add_x_a(T a){
        min_f += max(T(0), L.top() + add_l - a);
        L.push(a - add_l);
        R.push(L.top() + add_l - add_r);
        L.pop();
    }
    void add_a_x(T a){
        min_f += max(T(0), a - (R.top() + add_r));
        R.push(a - add_r);
        L.push(R.top() + add_r - add_l);
        R.pop();
    }
    void add_abs(T a){
        add_a_x(a);
        add_x_a(a);
    }
    void right_clear(){
        R = priority_queue<T, vector<T>, greater<T>>{};
        R.push(INFT);
    }
    void left_clear(){
        L = priority_queue<T>{};
        L.push(-INFT);
    }
    void shift(T a, T b){
        assert(a <= b);
        add_l += a;
        add_r += b;
    }
    void shift(T a){
        shift(a, a);
    }
};
