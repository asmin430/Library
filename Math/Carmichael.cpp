long long carmichael(long long x){
    long long r = 1;
    int b = 0;
    while((x & 1) == 0){
        ++b;
        x >>= 1;
        r <<= 1;
    }
    if(b > 1){
        if(b == 2) r >>= 1;
        else r >>= 2;
    }
    long long y = 3;
    while(y * y <= x){
        int c = 0;
        long long s = 1;
        while(x % y == 0){
            x /= y;
            ++c;
            s *= y;
        }
        if(c > 0){
            r = lcm(r, s / y * (y - 1));
        }
        ++y;
    }
    if(x > 1){
        r = lcm(r, x - 1);
    }
    return r;
}