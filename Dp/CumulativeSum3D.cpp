template<class T>
struct CumulativeSum3D{
    vector<vector<vector<T>>> dat;
    CumulativeSum2D(int a, int b, int c): dat(a + 1, vector<vector<T>>(b + 1, vector<T>(c + 1, 0))){}
    CumulativeSum2D(vector<vector<vector<T>>> &v): dat(v.size() + 1, vector<vector<T>>(v[0].size() + 1, vector<T>(v[0][0].size() + 1, 0))){
        for(int i = 0; i < v.size(); ++i){
            for(int j = 0; j < v[i].size(); ++j){
                for(int k = 0; k < v[i][j].size(); ++k){
                    dat[i + 1][j + 1][k + 1] += v[i][j][k];
                }
            }
        }
        build();
    }
    void add(int x, int y, int z, T w){
        ++x, ++y, ++z;
        if(x >= dat.size() || y >= dat[0].size(), z >= dat[0][0].size()) return;
        dat[x][y][z] += w;
    }
    void build(){
        for(int i = 1; i < dat.size(); ++i){
            for(int j = 1; j < dat[i].size(); ++j){
                for(int k = 1; k < dat[i][j].size(); ++k){
                    dat[i][j][k] += dat[i][j][k - 1] + dat[i][j - 1][k] + dat[i - 1][j][k] - dat[i][j - 1][k - 1] - dat[i - 1][j - 1][k] - dat[i - 1][j][k - 1] + dat[i - 1][j - 1][k - 1];
                }
            }
        }
    }

    //[sx, gx),[sy, gy),[sz, gz) の総和
    T sum(int sx, int sy, int sz, int gx, int gy, int gz) const {
        return (dat[gx][gy][gz] - dat[sx][gy][gz] - dat[gx][sy][gz] - dat[gx][gy][sz] + dat[gx][sy][sz] + dat[sx][gy][sz] + dat[sx][sy][gz] - dat[sx][sy][sz]);
    }
};
