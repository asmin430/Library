namespace geometry{
    using Real = long double;
    using Point = complex<Real>;
    const Real EPS = 1e-12;
    const Real PI = acos(-1);
    // a == b かどうか
    inline bool equal(const Real& a, const Real& b){
        return fabs(a - b) < EPS;
    }
    // r の符号
    inline int sign(const Real &r) { return r <= -EPS ? -1 : r >= EPS ? 1 : 0; }
    // 実部の比較
    bool compare_x(const Point &a, const Point &b) {
        return equal(real(a), real(b)) ? imag(a) < imag(b) : real(a) < real(b);
    }
    // 虚部の比較
    bool compare_y(const Point &a, const Point &b) {
        return equal(imag(a), imag(b)) ? real(a) < real(b) : imag(a) < imag(b);
    }
    // a の単位ベクトル
    Point unitvector(const Point& a){
        return a / abs(a);
    }
    // a の法線ベクトル
    Point normalvector(const Point& a){
        return a * Point(0, 1);
    }
    // a,b の内積
    Real dot(const Point &a, const Point &b){
        return (a.real() * b.real() + a.imag() * b.imag());
    }
    // a,b の外積
    Real cross(const Point &a, const Point &b){
        return (a.real() * b.imag() - a.imag() * b.real());
    }
    // p を rad ((単位)rad) 回転
    Point rotate(const Point &p, const Real &rad){
        return Point(p.real() * cos(rad) - p.imag() * sin(rad), p.real() * sin(rad) + p.imag() * cos(rad));
    }
    // 角度を rad -> ° に変換
    Real radtodeg(const Real &rad){
        return rad * 180 / PI;
    }
    // 角度を °　-> rad に変換
    Real degtorad(const Real &deg){
        return deg * PI / 180;
    }
    // 直電の構造体
    struct Line{
        Point a, b;
        Line() = default;
        Line(const Point &a, const Point &b): a(a), b(b){}
        Line(Real A, Real B, Real C){
            if(equal(A, 0)){
                a = Point(0, C / B);
                b = Point(1, C / B);
            }else if(equal(B, 0)){
                a = Point(C / A, 0);
                b = Point(C / A, 1);
            }else{
                a = Point(0, C / B);
                b = Point(C / A, 0);
            }
        }
    };
    // 線分の構造体
    struct Segment: Line{
        Segment() = default;
        Segment(Point a, Point b): Line(a, b){}
    };
    // 円の構造体
    struct Circle{
        Point p;
        Real r;
        Circle() = default;
        Circle(const Point &p, const Real &r): p(p), r(r){}
    };
    // 直線 l に点 p からおろした垂線の足
    Point projection(const Line &l, const Point &p){
        Real t = dot(p - l.a, l.a - l.b);
        t /= norm(l.a - l.b);
        return l.a + (l.a - l.b) * t;
    }
    // 直線 l について点 p と対称な点
    Point reflection(const Line &l, const Point &p){
        return p + (projection(l, p) - p) * (Real)2;
    }
    // a, b, c の位置関係
    int ccw(const Point &a, const Point &b, const Point &c){
        Point d = b - a;
        Point e = c - a;
        if(cross(d, e) > EPS) return 1; // counter clockwise
        if(cross(d, e) < -EPS) return -1; // clockwise
        if(dot(d, e) < 0) return 2; // c-b-a
        if(norm(d) < norm(e)) return -2; // a-b-c
        return 0; // a-c-b
    }
    // 直線 l1, l2 が垂直か？
    bool is_orthogonal(const Line &l1, const Line &l2){
        return equal(dot(l1.b - l1.a, l2.b - l2.a), 0);
    }
    // 直線 l1, l2 が平行か？
    bool is_parallel(const Line &l1, const Line &l2){
        return equal(cross(l1.b - l1.a, l2.b - l2.a), 0);
    }
    // 線分 l1, l2 か交差するか？
    bool is_intersect(const Segment &s1, const Segment &s2){
        return (ccw(s1.a, s1.b, s2.a) * ccw(s1.a, s1.b, s2.b) <= 0 && ccw(s2.a, s2.b, s1.a) * ccw(s2.a, s2.b, s1.b) <= 0);
    }
    // 直線 l1, l2 の交点
    Point crosspoint(const Line &l1, const Line &l2){
        Real A = cross(l1.b - l1.a, l2.b - l2.a);
        Real B = cross(l1.b - l1.a, l1.b - l2.a);
        if(equal(abs(A), 0) && equal(abs(B), 0)){
            return l1.a;
        }
        return (l2.b - l2.a) * (B / A) + l2.a;
    }
    // 線分 l1, l2 の交点
    Point crosspoint(const Segment &s1, const Segment &s2){
        return crosspoint(Line(s1), Line(s2));
    }
    // 直線 l1 と点 p の距離
    Real distanceLineAndPoint(const Line &l, const Point &p){
        return abs(cross(l.b - l.a, p - l.a)) / abs(l.b - l.a);
    }
    // 線分 l1 と点 p の距離
    Real distanceSegmentAndPoint(const Segment &l, const Point &p){
        if(dot(l.b - l.a, p - l.a) < EPS){
            return abs(p - l.a);
        }
        if(dot(l.a - l.b, p - l.b) < EPS){
            return abs(p - l.b);
        }
        return abs(cross(l.b - l.a, p - l.a)) / abs(l.b - l.a);
    }
    // 線分 l1, l2 の距離
    Real distanceSegmentAndSegment(const Segment &s1, const Segment &s2){
        if(is_intersect(s1, s2)){
            return (Real)(0);
        }
        return min({distanceSegmentAndPoint(s1, s2.a), distanceSegmentAndPoint(s1, s2.b), distanceSegmentAndPoint(s2, s1.a), distanceSegmentAndPoint(s2, s1.b)});
    }
    // 多角形の面積
    Real PolygonArea(const vector<Point> &p){
        Real res = 0;
        int n = p.size();
        for(int i = 0; i < n - 1; ++i){
            res += cross(p[i], p[i + 1]);
        }
        res += cross(p[n - 1], p[0]);
        return res * 0.5;
    }
    // p は凸多角形であるか？
    bool isConvex(const vector<Point> &p){
        int n = p.size();
        int now, pre, nxt;
        for(int i = 0; i < n; ++i){
            pre = (i - 1 + n) % n;
            nxt = (i + 1) % n;
            now = i;
            if(ccw(p[pre], p[now], p[nxt]) == -1){
                return false;
            }
        }
        return true;
    }
    // (多角形 p に) 含まれる:2, 辺上にある:1, 含まれない:0
    int isContained(const vector<Point> &g, const Point &p){
        bool in = false;
        int n = (int)g.size();
        for(int i = 0; i < n; ++i){
            Point a = g[i] - p;
            Point b = g[(i + 1) % n] - p;
            if(imag(a) > imag(b)) swap(a, b);
            if(imag(a) <= EPS && EPS < imag(b) && cross(a, b) < -EPS) in = !in;
            if(cross(a, b) == 0 && dot(a, b) <= 0) return 1;
        }
        return (in ? 2 : 0);
    }
    // 凸包を出力
    vector<Point> ConvexHull(vector<Point> &p, bool strict = true){
        long double si = -1;
        if(strict) si = 1;
        int n = (int)p.size(), k = 0;
        sort(begin(p), end(p), [](const Point &a, const Point &b){
            return (real(a) != real(b) ? real(a) < real(b) : imag(a) < imag(b));
        });
        vector<Point> ch(2 * n);

        for(int i = 0; i < n; ch[k++] = p[i++]){
            while(k >= 2 && cross(ch[k - 1] - ch[k - 2], p[i] - ch[k - 1]) < EPS * si) --k;
        }
        for(int i = n - 2, t = k + 1; i >= 0; ch[k++] = p[i--]){
            while(k >= t && cross(ch[k - 1] - ch[k - 2], p[i] - ch[k - 1]) < EPS * si) --k;
        }
        ch.resize(k - 1);
        return ch;
    }
    // 凸多角形の直径
    pair<int, int> ConvexDiameter(const vector<Point> &p){
        int n = (int)p.size();
        int is = 0, js = 0;
        for(int i = 1; i < n; ++i){
            if(imag(p[i]) > imag(p[is])) is = i;
            if(imag(p[i]) < imag(p[js])) js = i;
        }
        Real maxdis = norm(p[is] - p[js]);
        int maxi, maxj, i, j;
        i = maxi = is;
        j = maxj = js;
        do{
            if(cross(p[(i + 1) % n] - p[i], p[(j + 1) % n] - p[j]) >= 0){
                j = (j + 1) % n;
            }else{
                i = (i + 1) % n;
            }
            if(norm(p[i] - p[j]) > maxdis){
                maxdis = norm(p[i] - p[j]);
                maxi = i;
                maxj = j;
            }
        }while(i != is || j != js);
        return make_pair(maxi, maxj);
    }
    //　(凸多角形 p を直線 l で分割したときの)左側を出力
    vector<Point> ConvexCut(const vector<Point> &p, const Line &l){
        vector<Point> ret;
        for(int i = 0; i < p.size(); ++i){
            const Point &now = p[i];
            const Point &nxt = p[(i + 1) % p.size()];
            auto cf = cross(l.a - now, l.b - now);
            auto cs = cross(l.a - nxt, l.b - nxt);
            if(sign(cf) >= 0){
                ret.push_back(now);
            }
            if(sign(cf) * sign(cs) < 0){
                ret.push_back(crosspoint(Line(now, nxt), l));
            }
        }
        return ret;
    }
    // 最も近い点のペアの距離
    Real ClosestPair(vector<Point> p){
        sort(begin(p), end(p), [](const Point &a, const Point &b){
            return real(a) < real(b);
        });
        int n = p.size();
        if(n <= 1){
            return 1e18;
        }
        int m = n / 2;
        Real x = real(p[m]);
        Real res = min(ClosestPair(vector<Point>(begin(p), begin(p) + m)), ClosestPair(vector<Point>(begin(p) + m, end(p))));
        sort(begin(p), end(p), [](const Point &a, const Point &b){
            return imag(a) < imag(b);
        });
        deque<int> deq;
        for(int i = 0; i < n; ++i){
            if(res - abs(real(p[i]) - x) < EPS) continue;
            while(!deq.empty() && imag(p[deq.front()]) - imag(p[i]) + res < EPS) deq.pop_front();
            for(auto e: deq){
                if(res > sqrt(norm(p[i] - p[e]))){
                    res = sqrt(norm(p[i] - p[e]));
                }
            }
            deq.push_back(i);
        }
        return res;
    }
    // 点集合 p を完全に含む最小の円
    Circle MinimumEnclosingCircle(vector<Point> p, int seed = random_device()()){
        int n = p.size();
        assert(n >= 1);
        if(n == 1){
            return Circle(p[0], 0.0);
        }
        std::mt19937 mt(seed);
        std::shuffle(begin(p), end(p), mt);
        auto ps = begin(p);
        auto make_circle_3 = [](const Point &a, const Point &b, const Point &c) -> Circle {
            Real A = norm(b - c), B = norm(c - a), C = norm(a - b);
            Real S = cross(b - a, c - a);
            Point p = (A * (B + C - A) * a + B * (C + A - B) * b + C * (A + B - C) * c) / (4 * S * S); 
            Real r2 = norm(p - a);
            return Circle(p, r2);
        };
        auto make_circle_2 = [](const Point &a, const Point &b) -> Circle {
            Point c = (a + b) / (Real)2;
            Real r2 = norm(a - c);
            return Circle(c, r2);
        };
        auto in_circle = [](const Point &a, const Circle &c) -> bool {
            return norm(a - c.p) <= c.r + EPS;
        };
        Circle c = make_circle_2(ps[0], ps[1]);
        for(int i = 2; i < n; ++i){
            if(!in_circle(ps[i], c)){
                c = make_circle_2(ps[0], ps[i]);
                for(int j = 1; j < i; ++j){
                    if(!in_circle(ps[j], c)){
                        c = make_circle_2(ps[i], ps[j]);
                        for(int k = 0; k < j; ++k){
                            if(!in_circle(ps[k], c)){
                                c = make_circle_3(ps[i], ps[j], ps[k]);
                            }
                        }
                    }
                }
            }
        }
        c.r = sqrt(c.r);
        return c;
    }
    
    int CircleIntersect(const Circle &c1, const Circle &c2){
        Real dist = (real(c1.p) - real(c2.p)) * (real(c1.p) - real(c2.p)) + (imag(c1.p) - imag(c2.p)) * (imag(c1.p) - imag(c2.p));
        if(equal(dist, (c1.r - c2.r) * (c1.r - c2.r))) return 1;
        if(equal(dist, (c1.r + c2.r) * (c1.r + c2.r))) return 3;
        if(dist > (c1.r + c2.r) * (c1.r + c2.r)) return 4;
        if(dist < (c1.r - c2.r) * (c1.r - c2.r)) return 0;
        return 2;
    }
    // 三角形 abc の内接円
    Circle InCircle(const Point &a, const Point &b, const Point &c){
        Real A = abs(b - c), B = abs(c - a), C = abs(a - b);
        Point p(A * real(a) + B * real(b) + C * real(c), A * imag(a) + B * imag(b) + C * imag(c));
        p /= (A + B + C);
        Real r = distanceLineAndPoint(Line(a, b), p);
        return Circle(p, r);
    }
    // 三角形 abc の外接円
    Circle CircumCircle(const Point &a, const Point &b, const Point &c){
        Real A = norm(b - c), B = norm(c - a), C = norm(a - b);
        Real S = cross(b - a, c - a);
        Point p = (A * (B + C - A) * a + B * (C + A - B) * b + C * (A + B - C) * c) / (4 * S * S); 
        Real r = abs(p - a);
        return Circle(p, r);
    }
    // 円 c と直線 l の交点
    vector<Point> CrossPointLineAndCircle(const Circle &c, const Line &l){
        vector<Point> res;
        Real d = distanceLineAndPoint(l, c.p);
        if(d > c.r + EPS) return res;
        Point h = projection(l, c.p);
        if(equal(d, c.r)){
            res.push_back(h);
            return res;
        }
        Point e = unitvector(l.b - l.a);
        Real ph = sqrt(c.r * c.r - d * d);
        res.push_back(h - e * ph);
        res.push_back(h + e * ph);
        return res;
    }
    // 円 c と線分 l の交点
    vector<Point> CrossPointSegmentAndCircle(const Circle &c, const Segment &l){
        auto res1 = CrossPointLineAndCircle(c, (Line)l);
        vector<Point> res;
        for(auto p: res1){
            if(abs(ccw(l.a, p, l.b)) == 2) res.push_back(p);
        }
        return res;
    }
    // 円 c1, c2 の交点
    vector<Point> CrossPointCircleAndCircle(const Circle &c1, const Circle &c2){
        vector<Point> res;
        int mode = CircleIntersect(c1, c2);
        Real d = abs(c1.p - c2.p);
        if(mode == 4){
            return res;
        }
        if(mode == 0){
            return res;
        }
        if(mode == 3){
            Real t = c1.r / (c1.r + c2.r);
            res.push_back(c1.p + (c2.p - c1.p) * t);
            return res;
        }
        if(mode == 1){
            if(c2.r < c1.r - EPS){
                res.push_back(c1.p + (c2.p - c1.p) * (c1.r / d));
            }else{
                res.push_back(c2.p + (c1.p - c2.p) * (c2.r / d));
            }
            return res;
        }
        Real rc1 = (c1.r * c1.r + d * d - c2.r * c2.r) / (2 * d);
        Real rs1 = sqrt(c1.r * c1.r - rc1 * rc1);
        if(c1.r - abs(rc1) < EPS){
            rs1 = 0;
        }
        Point e12 = (c2.p - c1.p) / abs(c2.p - c1.p);
        res.push_back(c1.p + rc1 * e12 + rs1 * e12 * Point(0, 1));
        res.push_back(c1.p + rc1 * e12 + rs1 * e12 * Point(0, -1));
        return res;
    }
    // 点 p を通る円 c の接線
    vector<Point> tangentToCircle(const Point &p, const Circle &c){
        return CrossPointCircleAndCircle(c, Circle(p, sqrt(norm(c.p - p) - c.r * c.r)));
    }
    // 円 a, b の共通接線
    vector<Line> commonTangent(const Circle &a, const Circle &b){
        vector<Line> ret;
        double g = abs(a.p - b.p);
        if(equal(g, 0)){
            return ret;
        }
        Point u = unitvector(b.p - a.p);
        Point v = rotate(u, PI / 2);
        for(int s : {-1, 1}){
            Real h = (a.r + b.r * s) / g;
            if(equal(h * h, 1)){
                ret.emplace_back(a.p + (h > 0 ? u : -u) * a.r, a.p + (h > 0 ? u : -u) * a.r + v);
            }else if(1 - h * h > 0){
                Point U = u * h, V = v * sqrt(1 - h * h);
                ret.emplace_back(a.p + (U + V) * a.r, b.p - (U + V) * (b.r * s));
                ret.emplace_back(a.p + (U - V) * a.r, b.p - (U - V) * (b.r * s));
            }
        }
        return ret;
    }
    // 円 c と多角形 ps の共通部分の面積
    Real CommonAreaCircleAndPolygon(const Circle &c, const vector<Point> &ps){
        const int n = ps.size();
        vector<Point> nps;
        for(int i = 0; i < n; ++i){
            Point p1 = ps[i], p2 = ps[(i + 1) % n];
            nps.push_back(p1 - c.p);
            auto res1 = CrossPointLineAndCircle(c, Line(p1, p2));
            if(res1.size() == 2){
                Point dp1 = res1[0], dp2 = res1[1];
                Real ratio1 = dot(dp1 - p1, p2 - p1) / norm(p2 - p1);
                Real ratio2 = dot(dp2 - p1, p2 - p1) / norm(p2 - p1);
                if(0 < ratio1 && ratio1 < 1){
                    if(0 < ratio2 && ratio2 < 1){
                        if(ratio1 > ratio2) swap(dp1, dp2);
                        nps.push_back(dp1 - c.p);
                        nps.push_back(dp2 - c.p);
                    }else{
                        nps.push_back(dp1 - c.p);
                    }
                }else if(0 < ratio2 && ratio2 < 1){
                    nps.push_back(dp2 - c.p);
                }
            }
        }
        int m = nps.size();
        Real res = 0;
        for(int i = 0; i < m; ++i){
            Point p1 = nps[i], p2 = nps[(i + 1) % m];
            Point mid = (p1 + p2) / (Real)2.0;
            if(norm(mid) > c.r * c.r - EPS){
                Real th = atan2(cross(p1, p2), dot(p1, p2));
                res += c.r * c.r * th / 2.0;
            }else{
                res += cross(p1, p2) / 2.0;
            }
        }
        return res;
    }
    // 円 c1, c2 の共通部分の面積
    Real CommonAreaCircleAndCircle(const Circle &c1, const Circle &c2){
        Real d = norm(c1.p - c2.p);
        Real res = 0.0;
        if(d <= (c1.r - c2.r) * (c1.r - c2.r)) res = min(c1.r, c2.r) * min(c1.r, c2.r) * acos(-1);
        else if(d < (c1.r + c2.r) * (c1.r + c2.r)){
            Real t1 = acosl((Real)(c1.r * c1.r - c2.r * c2.r + d) / 2.0 / c1.r / sqrtl(d));
            Real t2 = acosl((Real)(c2.r * c2.r - c1.r * c1.r + d) / 2.0 / c2.r / sqrtl(d));
            res = c1.r * c1.r * t1 + c2.r * c2.r * t2 - c1.r * c1.r * sinl(2 * t1) / 2.0 - c2.r * c2.r * sinl(2 * t2) / 2.0;
        }
        return res;
    }
    

}

using namespace geometry;
