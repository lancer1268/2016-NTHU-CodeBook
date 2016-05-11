

LL LLmul(LL a,LL b,const LL &mod)
{
    LL ans = 0;
    while(b) {
        if(b&1) {
            ans += a;
            if(ans>=mod) ans -= mod;
        }
        b >>= 1;
        a <<= 1;
        if(a>=mod) a-=mod;
    }
    return ans;
}

LL modexp(LL x,LL e,const LL &mod)
{
    LL ans = 1;
    while(e) {
        if(e&1) ans = LLmul(ans,x,mod);
        e >>= 1;
        x = LLmul(x,x,mod);
    }
    return ans%mod;
}

bool Miller(const int &base,const LL &n)
{
    if(base>=n) return true;
    LL d = n-1, s = 0;
    while( !(d&1) ) { s++; d>>=1; }
    LL x = modexp(base,d,n);
    if(x==1) return true;
    for(int r=0;r<s;r++,x=LLmul(x,x,n))
        if(x==n-1)
            return true;
    return false;
}

bool Isprime(const LL &n)
{
    if(n<MAXPRIME) {
        return !iscom[(int)n];
    }
    for(int i=0;i<12;++i) {
        if(n==prime[i]) return true;
        if(n%prime[i]==0) return false;
    }
    for(int i=0;i<12;++i)
        if(!Miller(prime[i],n))
            return false;
    return true;
}

LL func(const LL n,const LL mod,const int c)
{
    return (LLmul(n,n,mod) + c +mod)%mod;
}

LL pollorrho(const LL n,const int c)
{
    LL a = 1, b = 1;
    a = func(a,n,c)%n;
    b = func(b,n,c)%n; b = func(b,n,c)%n;
    while(gcd(abs(a-b),n)==1) {
        a = func(a,n,c)%n;
        b = func(b,n,c)%n; b = func(b,n,c)%n;
    }
    return gcd(abs(a-b),n);
}

void prefactor(LL &n, vector<LL> &v)
{
    for(int i=0;i<12;++i) {
        while(n%prime[i]==0) {
            v.push_back(prime[i]);
            n /= prime[i];
        }
    }
}

void smallfactor(LL n, vector<LL> &v)
{
    if(n<MAXPRIME) {
        while(iscom[(int)n]) {
            v.push_back(iscom[(int)n]);
            n /= iscom[(int)n];
        }
        v.push_back(n);
    } else {
        for(int i=0;i<primecnt&&prime[i]*prime[i]<=n;++i) {
            while(n%prime[i]==0) {
                v.push_back(prime[i]);
                n /= prime[i];
            }
        }
        if(n!=1) v.push_back(n);
    }
}

void comfactor(const LL &n, vector<LL> &v)
{
    if(n<1e9) {
        smallfactor(n,v);
        return;
    }
    if(Isprime(n)) {
        v.push_back(n);
        return;
    }
    LL d;
    for(int c=3;;++c) {
        d = pollorrho(n,c);
        if(d!=n) break;
    }
    comfactor(d,v);
    comfactor(n/d,v);
}

void Factor(const LL &x, vector<LL> &v)
{
    LL n = x;
    if(n==1) { puts("Factor 1"); return; }
    prefactor(n,v);
    if(n==1) return;
    comfactor(n,v);
    sort(v.begin(),v.end());
}

void AllFactor(const LL &n,vector<LL> &v)
{
    vector<LL> tmp;
    Factor(n,tmp);
    v.clear();
    v.push_back(1);
    int len;
    LL now = 1;
    for(int i=0;i<(int)tmp.size();++i) {
        if(i==0 || tmp[i]!=tmp[i-1]) {
            len = v.size();
            now = 1;
        }
        now *= tmp[i];
        for(int j=0;j<len;++j)
            v.push_back(v[j]*now);
    }
}
