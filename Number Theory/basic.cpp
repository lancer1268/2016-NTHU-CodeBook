typedef long long int LL;
template<typename T>
void gcd(const T &a,const T &b,T &d,T &x,T &y){
    if(!b) d=a,x=1,y=0;
    else gcd(b,a%b,d,y,x), y-=x*(a/b);
}

const int MAXPRIME = 1000000;
int isp[MAXPRIME], prime[MAXPRIME], primecnt;
void sieve() {
    memset(isp,0,sizeof(isp));
    primecnt = 0;
    for(int i=2;i<MAXPRIME;++i) {
        if(!isp[i]) prime[primecnt++] = i;
        for(int j=0;j<primecnt;++j) {
            if(i*prime[j]>=MAXPRIME) break;
            isp[i*prime[j]] = prime[j];
            if(i%prime[j]==0) break;
        }
    }
}

bool g_test(const LL &g, const LL &p, const vector<LL> &v) {
    for(int i=0;i<v.size();++i)
        if(modexp(g,(p-1)/v[i],p)==1)
            return false;
    return true;
}
LL primitive_root(const LL &p) {
    if(p==2) return 1;
    vector<LL> v;
    Factor(p-1,v);
    v.erase(unique(v.begin(), v.end()), v.end());
    for(LL g=2;g<p;++g)
        if(g_test(g,p,v))
            return g;
    puts("primitive_root NOT FOUND");
    return -1;
}

int Legendre(const LL &a, const LL &p) { return modexp(a%p,(p-1)/2,p); }

LL inv(const LL &a, const LL &n) {
    LL d,x,y;
    gcd(a,n,d,x,y);
    return d==1 ? (x+n)%n : -1;
}

LL china(const int &n, const int a[], const int m[]) {
    // x = a[i] (mod m[i])
    LL M=1, d, y, x=0;
    for(int i=0;i<n;++i) M *= m[i];
    for(int i=0;i<n;++i) {
        LL w = M/m[i];
        gcd(m[i], w, d, d, y);
        x = (x + y*w*a[i]) % M;
    }
    return (x+M)%M;
}

LL log_mod(const LL &a, const LL &b, const LL &p) {
    // a ^ x = b ( mod p )
    int m=sqrt(p+.5), e=1;
    LL v=inv(modexp(a,m,p), p);
    map<LL,int> x;
    x[1]=0;
    for(int i=1;i<m;++i) {
        e = LLmul(e,a,p);
        if(!x.count(e)) x[e] = i;
    }
    for(int i=0;i<m;++i) {
        if(x.count(b)) return i*m + x[b];
        b = LLmul(b,v,p);
    }
    return -1;
}

LL Tonelli_Shanks(const LL &n, const LL &p) {
    // x^2 = n ( mod p )
    if(n==0) return 0;
    if(Legendre(n,p)!=1) while(1) { puts("SQRT ROOT does not exist"); }
    int S = 0;
    LL Q = p-1;
    while( !(Q&1) ) { Q>>=1; ++S; }
    if(S==1) return modexp(n%p,(p+1)/4,p);
    LL z = 2;
    for(;Legendre(z,p)!=-1;++z)
    LL c = modexp(z,Q,p);
    LL R = modexp(n%p,(Q+1)/2,p), t = modexp(n%p,Q,p);
    int M = S;
    while(1) {
        if(t==1) return R;
        LL b = modexp(c,1L<<(M-i-1),p);
        R = LLmul(R,b,p);
        t = LLmul( LLmul(b,b,p), t, p);
        c = LLmul(b,b,p);
        M = i;
    }
    return -1;
}
