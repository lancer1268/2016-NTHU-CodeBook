Codebook
=======

Code:.\Computational Geometry\Geometry.cpp
================

```cpp
#include<cmath>
#include<vector>
#include<algorithm>
template<typename T>
struct point{
	T x,y;
	point(){}
	point(const T&dx,const T&dy):x(dx),y(dy){}
	inline const point operator+(const point &b)const{
		return point(x+b.x,y+b.y);
	}
	inline const point operator-(const point &b)const{
		return point(x-b.x,y-b.y);
	}
	inline const point operator*(const point &b)const{
		return point(x*b.x,y*b.y);
	}
	inline const point operator/(const point &b)const{
		return point(x/b.x,y/b.y);
	}
	inline const point operator+(const T &b)const{
		return point(x+b,y+b);
	}
	inline const point operator-(const T &b)const{
		return point(x-b,y-b);
	}
	inline const point operator*(const T &b)const{
		return point(x*b,y*b);
	}
	inline const point operator/(const T &b)const{
		return point(x/b,y/b);
	}
	inline bool operator==(const point &b)const{
		return x==b.x&&y==b.y;
	}
	inline const T dot(const point &b)const{
		return x*b.x+y*b.y;
	}
	inline const T cross(const point &b)const{
		return x*b.y-y*b.x;
	}
	inline point normal()const{/*求法向量*/
		return point(-y,x);
	}
	inline const T abs2()const{/*向量長度的平方*/
		return dot(*this);
	}
};
template<typename T>
struct line{
	line(){}
	point<T> p1,p2;
	T a,b,c;/*ax+by+c=0*/
	line(const point<T>&x,const point<T>&y):p1(x),p2(y){}
	inline void pton(){/*轉成一般式*/ 
		a=p1.y-p2.y;
		b=p2.x-p1.x;
		c=-a*p1.x-b*p1.y;
	}
	inline T cross(const point<T> &p)const{/*點和有向直線的關係，>0左邊、=0在線上<0右邊*/
		return (p2-p1).cross(p-p1);
	}
	inline bool point_on_segment(const point<T>&p)const{/*點是否線段上*/ 
		return cross(p)==0&&(p1-p).dot(p2-p)<=0; 
	}
	inline T dis2(const point<T> &p,bool is_segment=0)const{/*點跟直線/線段的距離平方*/ 
		point<T> v=p2-p1,v1=p-p1;
		if(is_segment){
			point<T> v2=p-p2;
			if(v.dot(v1)<=0)return v1.abs2();
			if(v.dot(v2)>=0)return v2.abs2();
		}
		T tmp=v.cross(v1);
		return tmp*tmp/v.abs2();
	}
	inline point<T> mirror(const point<T> &p)const{/*點對直線的鏡射*/ 
		/*要先呼叫pton轉成一般式*/
		point<T> ans;
		T d=a*a+b*b;
		ans.x=(b*b*p.x-a*a*p.x-2*a*b*p.y-2*a*c)/d;
		ans.y=(a*a*p.y-b*b*p.y-2*a*b*p.x-2*b*c)/d;
		return ans;
	}
	inline bool equal(const line &l)const{/*直線相等*/ 
		return cross(l.p1)==0&&cross(l.p2)==0;
	}
	inline bool parallel(const line &l)const{/*直線平行*/ 
		return (p1-p2).cross(l.p1-l.p2)==0;
	}
	inline bool cross_seg(const line &l)const{/*直線是否交線段*/
		return (p2-p1).cross(l.p1)*(p2-p1).cross(l.p2)<=0;
	}
	inline char line_intersect(const line &l)const{/*直線相交情況，-1無限多點、1交於一點、0不相交*/ 
		return parallel(l)?(cross(l.p1)==0?-1:0):1;
	}
	inline char seg_intersect(const line &l)const{/*線段相交情況，-1無限多點、1交於一點、0不相交*/ 
		T c1=(p2-p1).cross(l.p1-p1);
		T c2=(p2-p1).cross(l.p2-p1);
		T c3=(l.p2-l.p1).cross(p1-l.p1);
		T c4=(l.p2-l.p1).cross(p2-l.p1);
		if(c1*c2<0&&c3*c4<0)return 1;
		if(c1==0&&c2==0){
			if(p1==l.p1&&(p2-p1).dot(l.p2)<=0)return 1;
			if(p1==l.p2&&(p2-p1).dot(l.p1)<=0)return 1;
			if(p2==l.p1&&(p1-p2).dot(l.p2)<=0)return 1;
			if(p2==l.p2&&(p1-p2).dot(l.p1)<=0)return 1;
			return -1;
		}
		if(c1==0&&point_on_segment(l.p1))return 1;
		if(c2==0&&point_on_segment(l.p2))return 1;
		if(c3==0&&l.point_on_segment(p1))return 1;
		if(c4==0&&l.point_on_segment(p2))return 1;
		return 0;
	}
	inline point<T> line_intersection(const line &l)const{/*直線交點*/ 
		point<T> a=p2-p1,b=l.p2-l.p1,s=l.p1-p1;
		//if(a.cross(b)==0)return INF;
		return p1+a*s.cross(b)/a.cross(b);
	}
	inline point<T> seg_intersection(const line &l)const{/*線段交點*/ 
		T c1=(p2-p1).cross(l.p1-p1);
		T c2=(p2-p1).cross(l.p2-p1);
		T c3=(l.p2-l.p1).cross(p1-l.p1);
		T c4=(l.p2-l.p1).cross(p2-l.p1);
		if(c1*c2<0&&c3*c4<0)return line_intersection(l);
		if(c1==0&&c2==0){
			if(p1==l.p1&&(p2-p1).dot(l.p2)<=0)return p1;
			if(p1==l.p2&&(p2-p1).dot(l.p1)<=0)return p1;
			if(p2==l.p1&&(p1-p2).dot(l.p2)<=0)return p2;
			if(p2==l.p2&&(p1-p2).dot(l.p1)<=0)return p2;
		}else{
			if(c1==0&&point_on_segment(l.p1))return l.p1;
			if(c2==0&&point_on_segment(l.p2))return l.p2;
			if(c3==0&&l.point_on_segment(p1))return p1;
			if(c4==0&&l.point_on_segment(p2))return p2;
		}
		//return INF;
	}
};
template<typename T>
struct polygon{
	polygon(){}
	std::vector<point<T> > p;
	inline const point<T>& operator[](int id)const{
		return p[id];
	}
	inline T area()const{/*多邊形面積*/ 
		T ans=0;
		for(int i=p.size()-1,j=0;j<(int)p.size();i=j++)
			ans+=p[i].cross(p[j]);
		return ans/2;
	}
	inline point<T> gravity()const{/*多邊形重心*/
		T cx=0,cy=0,w=0;
		for(int i=p.size()-1,j=0;j<(int)p.size();i=j++){
			T a=p[i].cross(p[j]);
			cx+=(p[i].x+p[j].x)*a;
			cy+=(p[i].y+p[j].y)*a;
			w+=a;
		}
		return point<T>(cx/3/w,cy/3/w);
	}
	inline char chas(const point<T> &x)const{/*點是否在凸多邊形內，是的話回傳1、在邊上回傳-1、否則回傳0*/
		T tp=0,np;
		for(int i=p.size()-1,j=0;j<(int)p.size();i=j++){
			if(!(np=(p[j]-x).cross(p[i]-x)))return -1;
			if(tp*np<0)return 0;
			tp=(np!=0)?np:tp;
		}
		return 1;
	}
	inline char ahas(const point<T>& t)const{/*點是否在簡單多邊形內，是的話回傳1、在邊上回傳-1、否則回傳0*/ 
		bool c=0;
		for(int i=0,j=p.size()-1;i<p.size();j=i++)
			if(line<T>(p[i],p[j]).point_on_segment(t))return -1;
			else if((p[i].y>t.y)!=(p[j].y>t.y)&&
			t.x<(p[j].x-p[i].x)*(t.y-p[i].y)/(p[j].y-p[i].y)+p[i].x)
				c=!c;
		return c;
	}
	inline polygon cut(const line<T> &l)const{/*凸包對直線切割，得到直線l左側的凸包*/
		polygon ans;
		for(int n=p.size(),i=n-1,j=0;j<n;i=j++){
			if(l.cross(p[i])>=0){
				ans.p.push_back(p[i]); 
				if(l.cross(p[j])<0)
					ans.p.push_back(l.line_intersection(line<T>(p[i],p[j])));
			}else if(l.cross(p[j])>0)
				ans.p.push_back(l.line_intersection(line<T>(p[i],p[j])));
		}
		return ans;
	}
	inline static bool graham_cmp(const point<T>& a,const point<T>& b){
		return (a.x<b.x)||(a.x==b.x&&a.y<b.y);/*凸包排序函數*/
	}
	inline void graham(std::vector<point<T> > &s){/*凸包*/
		sort(s.begin(),s.end(),graham_cmp);
		p.resize(s.size()+1);
		int m=0;
		for(int i=0;i<(int)s.size();++i){
			while(m>=2&&(p[m-1]-p[m-2]).cross(s[i]-p[m-2])<=0)--m;
			p[m++]=s[i];
		}
		for(int i=s.size()-2,t=m+1;i>=0;--i){
			while(m>=t&&(p[m-1]-p[m-2]).cross(s[i]-p[m-2])<=0)--m;
			p[m++]=s[i];
		}
		if(s.size()>1)--m; 
		p.resize(m);
	}
	inline static char sign(const T&x){
		return x>=0?1:-1;
	}
	inline static bool angle_cmp(const line<T>& A,const line<T>& B){
		point<T>a=A.p2-A.p1,b=B.p2-B.p1;
		//return atan2(a.y,a.x)<atan2(b.y,b.x); 
		char ay=sign(a.y),by=sign(b.y),ax=sign(a.x),bx=sign(b.x);
		return ay>by||(ay==by&&(ax*ay>bx*by||(ax*ay==bx*by&&a.cross(b)>0)));
	}
	inline int halfplane_intersection(std::vector<line<T> > &s){
		sort(s.begin(),s.end(),angle_cmp);
		int L,R,n=s.size();
		std::vector<point<T> > px(n);
		std::vector<line<T> > q(n);
		q[L=R=0]=s[0];
		for(int i=1;i<n;++i){
			while(L<R&&s[i].cross(px[R-1])<=0)--R;
			while(L<R&&s[i].cross(px[L])<=0)++L;
			q[++R]=s[i];
			if(q[R].parallel(q[R-1])){
				--R;
				if(q[R].cross(s[i].p1)>0)q[R]=s[i];
			}
			if(L<R)px[R-1]=q[R-1].line_intersection(q[R]);
		}
		while(L<R&&q[L].cross(px[R-1])<=0)--R;
		p.clear();
		if(R-L<=1)return 0;
		px[R]=q[R].line_intersection(q[L]);
		for(int i=L;i<=R;++i)p.push_back(px[i]);
		return R-L+1;
	}
};
template<typename T>
struct triangle{
	point<T> a,b,c;
	triangle(){}
	triangle(const point<T> &a,const point<T> &b,const point<T> &c):a(a),b(b),c(c){}
	inline T area()const{
		T t=(b-a).cross(c-a)/2;
		return t>0?t:-t;
	}
	inline point<T> barycenter()const{/*重心*/
		return (a+b+c)/3;
	}
	inline point<T> circumcenter()const{/*外心*/
		static line<T> u,v;
		u.p1.x=(a.x+b.x)/2;
		u.p1.y=(a.y+b.y)/2;
		u.p2.x=u.p1.x-a.y+b.y;
		u.p2.y=u.p1.y+a.x-b.x;
		v.p1.x=(a.x+c.x)/2;
		v.p1.y=(a.y+c.y)/2;
		v.p2.x=v.p1.x-a.y+c.y;
		v.p2.y=v.p1.y+a.x-c.x;
		return u.line_intersection(v);
	}
	inline point<T> incenter()const{/*內心，用到根號*/
		T A=sqrt((b-c).abs2()),B=sqrt((a-c).abs2()),C=sqrt((a-b).abs2());
		return point<T>(A*a.x+B*b.x+C*c.x,A*a.y+B*b.y+C*c.y)/(A+B+C);
	}
	inline point<T> perpencenter()const{/*垂心*/
		return barycenter()*3-circumcenter()*2;
	}
};

```

Code:.\Computational Geometry\SmallestCircle.cpp
================

```cpp
#include"Geometry.cpp"
#include<vector>
struct Circle{
    typedef point<double> p;
    typedef const point<double> cp;
    p x;
    double r2;
    bool incircle(cp &c)const{
        return (x-c).abs2() <= r2;
    }
};

Circle TwoPointCircle(Circle::cp &a,Circle::cp &b)
{
    Circle::p m=(a+b)/2;
    return (Circle){m,(a-m).abs2()};
}

Circle outcircle(Circle::p a,Circle::p b,Circle::p c)
{
    if(TwoPointCircle(a,b).incircle(c)) return TwoPointCircle(a,b);
    if(TwoPointCircle(b,c).incircle(a)) return TwoPointCircle(b,c);
    if(TwoPointCircle(c,a).incircle(b)) return TwoPointCircle(c,a);

    Circle::p ret;
    double a1=b.x-a.x, b1=b.y-a.y, c1=(a1*a1+b1*b1)/2;
    double a2=c.x-a.x, b2=c.y-a.y, c2=(a2*a2+b2*b2)/2;
    double d = a1*b2 - a2*b1;
    ret.x = a.x + (c1*b2-c2*b1)/d;
    ret.y = a.y + (a1*c2-a2*c1)/d;
    return (Circle){ret,(ret-a).abs2()};
}
//rand required
Circle SmallestCircle(std::vector<Circle::p> &p) // save Points in p
{
    int n = p.size();
    if(n==1) return (Circle){p[0],0.0};
    if(n==2) return TwoPointCircle(p[0],p[1]);

    random_shuffle(p.begin(),p.end());
    Circle c = {p[0],0.0};
    for(int i=0;i<n;++i){
        if(c.incircle(p[i])) continue;
        c = Circle{p[i],0.0};
        for(int j=0;j<i;++j){
            if(c.incircle(p[j])) continue;
            c = TwoPointCircle(p[i],p[j]);
            for(int k=0;k<j;++k){
                if(c.incircle(p[k])) continue;
                c = outcircle(p[i],p[j],p[k]);
            }
        }
    }
    return c;
}

```

Code:.\Data_Structure\space.cpp
================

```cpp

double norm(void) {
    return hypot(x,y);
}

double dis(Point a,Point b)
{
    return (a-b).norm();
}

```

Code:.\Data_Structure\treap.cpp
================

```cpp
#ifndef TREAP
#define TREAP
template<typename T>
class treap{
	private:
		struct node{
			T data;
			unsigned fix;
			int s;
			node *ch[2];
			node(const T&d):data(d),s(1){}
			node():s(0){ch[0]=ch[1]=this;}
		}*nil,*root;
		unsigned x;
		inline unsigned ran(){
			return x=x*0xdefaced+1;
		}
		inline void rotate(node *&a,bool d){
			node *b=a;
			a=a->ch[!d];
			a->s=b->s;
			b->ch[!d]=a->ch[d];
			a->ch[d]=b;
			b->s=b->ch[0]->s+b->ch[1]->s+1;
		}
		void insert(node *&o,const T &data){
			if(!o->s){
				o=new node(data),o->fix=ran();
				o->ch[0]=o->ch[1]=nil;
			}else{
				o->s++;
				bool d=o->data<data;
				insert(o->ch[d],data);
				if(o->ch[d]->fix>o->fix)rotate(o,!d);
			}
		}
		node *merge(node *a,node *b){
			if(!a->s||!b->s)return a->s?a:b;
			if(a->fix>b->fix){
				a->ch[1]=merge(a->ch[1],b);
				a->s=a->ch[0]->s+a->ch[1]->s+1;
				return a;
			}else{
				b->ch[0]=merge(a,b->ch[0]);
				b->s=b->ch[0]->s+b->ch[1]->s+1;
				return b;
			}
		}
		bool erase(node *&o,const T &data){
			if(!o->s)return 0;
			if(o->data==data){
				node *t=o;
				o=merge(o->ch[0],o->ch[1]);
				delete t;
				return 1;
			}
			if(erase(o->ch[o->data<data],data)){
				o->s--;return 1;
			}else return 0;
		}
		void clear(node *&o){
			if(o->s)clear(o->ch[0]),clear(o->ch[1]),delete o;
		}
	public:
		treap(unsigned s=20150119):nil(new node),root(nil),x(s){}
		~treap(){clear(root),delete nil;}
		inline void clear(){clear(root),root=nil;}
		inline void insert(const T &data){
			insert(root,data);
		}
		inline bool erase(const T &data){
			return erase(root,data);
		}
		inline bool find(const T&data){
			for(node *o=root;o->s;)
			if(o->data==data)return 1;
			else o=o->ch[o->data<data];
			return 0;
		}
		inline int rank(const T&data){
			int cnt=0;
			for(node *o=root;o->s;)
			if(o->data<data)cnt+=o->ch[0]->s+1,o=o->ch[1];
			else o=o->ch[0];
			return cnt;
		}
		inline const T&kth(int k){
			for(node *o=root;;)
			if(k<=o->ch[0]->s)o=o->ch[0];
			else if(k==o->ch[0]->s+1)return o->data;
			else k-=o->ch[0]->s+1,o=o->ch[1];
		}
		inline const T&operator[](int k){
			return kth(k);
		}
		inline const T&preorder(const T&data){
			node *x=root,*y=0;
			while(x->s)
			if(x->data<data)y=x,x=x->ch[1];
			else x=x->ch[0];
			if(y)return y->data;
			return data;
		}
		inline const T&successor(const T&data){
			node *x=root,*y=0;
			while(x->s)
			if(data<x->data)y=x,x=x->ch[0];
			else x=x->ch[1];
			if(y)return y->data;
			return data;
		}
		inline int size(){return root->s;}
};
#endif

```

Code:.\Graph\Arborescence_EV.cpp
================

```cpp
#include <bits/stdc++.h>
using namespace std;

struct node {
    int from, to, cost;
    node(int from=0,int to=0,int cost=0):from(from),to(to),cost(cost){};
} edge[M];

int m, n, m, c;
int far[N], In[N], ID[N], vis[N];

bool MST(int cost,int n,int root)
{
    long long int ans=0;
    while(true)
    {
        for(int i=0;i<n;++i) IN[i].first = INF;
        for(int i=0;i<m;++i)
            if(edge[i].from!=edge[i].to)
                IN[edge[i].to] = min(IN[edge[i].to],make_pair(edge[i].cost,edge[i].from));
        for(int i=0;i<n;++i)
            if(i!=root && IN[i].first==INF)
                return false; // NO Arborescence

        int cntnode = 0;
        memset(ID,-1,sizeof(ID));
        memset(vis,-1,sizeof(vis));
        In[root] = 0;
        for(int i=0;i<n;++i) ans += IN[i].first;
        for(int i=0;i<n;++i) {
            int x;
            for(x=i;vis[x]!=i&&ID[x]==-1&&x!=root;x=IN[x].second)
                vis[x] = i;
            if(ID[x]==-1 && x!=root) {
                for(int i=IN[x].second;u!=x;u=IN[u].second)
                    ID[u] = cntnode;
                ++cntnode;
            }
        }
        if(cntnode==0)  break; // END

        for(int i=0;i<n;++i)
            if(ID[i]==-1)
                ID[i] = cntnode++;

        for(int i=0;i<m;++i) {
            int v = edge[i].to;
            edge[i].from = ID[edge[i].from];
            edge[i].to = ID[edge[i].to];
            if(edge[i].from!=edge[i].to)
                edge[i].cost -= IN[edge[i].to].first;
        }
        n=cntnode;
        root=ID[root];
    }
    return ans<=cost;
}

```

Code:.\Number_Theory\basic.cpp
================

```cpp
typedef long long int LL;

LL gcd(const LL &a, const LL &b) { return b==0 ? a : gcd(b,a%b); }
void gcd(const LL &a, const LL &b, LL &d, LL &x, LL &y) {
    if(!b) d=a,x=1,y=0;
    else gcd(b,a%b,d,y,x), y-=x*(a/b);
}

const int MAXPRIME = 1000000;
int iscom[MAXPRIME], prime[MAXPRIME], primecnt;
void sieve() {
    memset(iscom,0,sizeof(iscom));
    primecnt = 0;
    for(int i=2;i<MAXPRIME;++i) {
        if(!iscom[i]) prime[primecnt++] = i;
        for(int j=0;j<primecnt;++j) {
            if(i*prime[j]>=MAXPRIME) break;
            iscom[i*prime[j]] = prime[j];
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

// Harmonic series H_n = ln(n) + gamma + 1/(2n) - 1/(12nn) + 1/(120nnnn)
// gamma = 0.57721566490153286060651209008240243104215933593992

```

Code:.\Number_Theory\enumerate.cpp
================

```cpp

void all_divdown(const LL &n) { // all n/x
    for(LL a=1;a<=n;a=n/(n/(a+1))) {
        dosomething;
    }
}

LL Catalan(const LL &n) { return (2n)! / n! / n! / (n+1); }

```

Code:.\Number_Theory\eulerphi.cpp
================

```cpp
int eulerPhi(int n){
	int m = sqrt(n+0.5);
	int res=n;
	for(int i=2; i<=m; i++){
		if(n%i==0){
			res = res*(i-1)/i;
			while(n%i==0)n/=i;
		}
	}
	if(n>1) res = res*(n-1)/n;
	return res;
}

vector<int> phiTable(int n){
	vector<int>phi(n+1, 0);
	phi[1] = 1;
	for(int i=2; i<=n; i++) if(!phi[i])
		for(int j=i; j<=n; j+=i){
			if(!phi[j])phi[j] = j;
			phi[j] = phi[j]*(i-1)/i;
		}
	return phi;
}
```

Code:.\Number_Theory\Factor.cpp
================

```cpp
LL LLmul(LL a, LL b, const LL &mod) {
    LL ans=0;
    while(b) {
        if(b&1) {
            ans+=a;
            if(ans>=mod) ans-=mod;
        }
        a<<=1, b>>=1;
        if(a>=mod) a-=mod;
    }
    return ans;
}

LL modexp(LL x, LL e, const LL &mod) {
    LL ans=1;
    while(e) {
        if(e&1) ans=LLmul(ans,x,mod);
        e>>=1;
        x=LLmul(x,x,mod);
    }
    return ans%mod;
}

bool Miller(const int &base, const LL &n) {
    if(base>=n) return true;
    LL d=n-1, s=0;
    while(!(d&1)) ++s, d>>=1;
    LL x=modexp(base,d,n);
    if(x==1) return true;
    for(int r=0;r<s;r++,x=LLmul(x,x,n))
        if(x==n-1)
            return true;
    return false;
}

bool Isprime(const LL &n) {
    if(n<MAXPRIME) return !iscom[(int)n];
    for(int i=0;i<12;++i) {
        if(n==prime[i]) return true;
        if(n%prime[i]==0) return false;
    }
    for(int i=0;i<12;++i)
        if(!Miller(prime[i],n))
            return false;
    return true;
}

LL func(const LL n,const LL mod,const int c) {
    return (LLmul(n,n,mod)+c+mod)%mod;
}

LL pollorrho(const LL n, const int c) {
    LL a=1, b=1;
    a=func(a,n,c)%n;
    b=func(b,n,c)%n; b=func(b,n,c)%n;
    while(gcd(abs(a-b),n)==1) {
        a=func(a,n,c)%n;
        b=func(b,n,c)%n; b=func(b,n,c)%n;
    }
    return gcd(abs(a-b),n);
}

void prefactor(LL &n, vector<LL> &v) {
    for(int i=0;i<12;++i) {
        while(n%prime[i]==0) {
            v.push_back(prime[i]);
            n/=prime[i];
        }
    }
}

void smallfactor(LL n, vector<LL> &v) {
    if(n<MAXPRIME) {
        while(iscom[(int)n]) {
            v.push_back(iscom[(int)n]);
            n/=iscom[(int)n];
        }
        v.push_back(n);
    } else {
        for(int i=0;i<primecnt&&prime[i]*prime[i]<=n;++i) {
            while(n%prime[i]==0) {
                v.push_back(prime[i]);
                n/=prime[i];
            }
        }
        if(n!=1) v.push_back(n);
    }
}

void comfactor(const LL &n, vector<LL> &v) {
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

void Factor(const LL &x, vector<LL> &v) {
    LL n = x;
    if(n==1) { puts("Factor 1"); return; }
    prefactor(n,v);
    if(n==1) return;
    comfactor(n,v);
    sort(v.begin(),v.end());
}

void AllFactor(const LL &n,vector<LL> &v) {
    vector<LL> tmp;
    Factor(n,tmp);
    v.clear();
    v.push_back(1);
    int len;
    LL now=1;
    for(int i=0;i<tmp.size();++i) {
        if(i==0 || tmp[i]!=tmp[i-1]) {
            len = v.size();
            now = 1;
        }
        now*=tmp[i];
        for(int j=0;j<len;++j)
            v.push_back(v[j]*now);
    }
}

```

Code:.\Number_Theory\find_real_root.cpp
================

```cpp
// an*x^n + ... + a1x + a0 = 0;
int sign(double x){
    return x < -eps ? -1 : x > eps;
}

double get(const vector<double>&coef, double x){
    double e = 1, s = 0;
    for(auto i : coef) s += i*e, e *= x;
    return s;
}

double find(const vector<double>&coef, int n, double lo, double hi){
    double sign_lo, sign_hi;
    if( !(sign_lo = sign(get(coef,lo))) ) return lo;
    if( !(sign_hi = sign(get(coef,hi))) ) return hi;
    if(sign_lo * sign_hi > 0) return INF;
    for(int stp = 0; stp < 100 && hi - lo > eps; ++stp){
        double m = (lo+hi)/2.0;
        int sign_mid = sign(get(coef,m));
        if(!sign_mid) return m;
        if(sign_lo*sign_mid < 0) hi = m;
        else lo = m;
    }
    return (lo+hi)/2.0;
}

vector<double> cal(vector<double>coef, int n){
    vector<double>res;
    if(n == 1){
        if(sign(coef[1])) res.pb(-coef[0]/coef[1]);
        return res;
    }
    vector<double>dcoef(n);
    for(int i = 0; i < n; ++i) dcoef[i] = coef[i+1]*(i+1);
    vector<double>droot = cal(dcoef, n-1);
    droot.insert(droot.begin(), -INF);
    droot.pb(INF);
    for(int i = 0; i+1 < droot.size(); ++i){
        double tmp = find(coef, n, droot[i], droot[i+1]);
        if(tmp < INF) res.pb(tmp);
    }
    return res;
}

int main () {
    vector<double>ve;
    vector<double>ans = cal(ve, n);
    // 視情況把答案 +eps，避免 -0
}

```

Code:.\Tree_problem\link cut tree.cpp
================

```cpp
#include<vector>
struct splay_tree{
	int ch[2],pa;/*子節點跟父母*/ 
	bool rev;/*反轉的懶惰標記*/ 
	splay_tree():pa(0),rev(0){ch[0]=ch[1]=0;}
};
std::vector<splay_tree> node;
/*
有的時候用vector會TLE，要注意 
這邊以node[0]作為null節點 
*/
inline bool isroot(int x){/*判斷是否為這棵splay tree的根*/ 
	return node[node[x].pa].ch[0]!=x&&node[node[x].pa].ch[1]!=x;
}
inline void down(int x){/*懶惰標記下推*/ 
	if(node[x].rev){
		if(node[x].ch[0])node[node[x].ch[0]].rev^=1;
		if(node[x].ch[1])node[node[x].ch[1]].rev^=1;
		std::swap(node[x].ch[0],node[x].ch[1]);
		node[x].rev^=1;
	}
}
void push_down(int x){/*將所有祖先的懶惰標記下推*/ 
	if(!isroot(x))push_down(node[x].pa);
	down(x);
}
inline void up(int x){}/*將子節點的資訊向上更新*/ 
inline void rotate(int x){/*旋轉，會自行判斷轉的方向*/ 
	int y=node[x].pa,z=node[y].pa,d=(node[y].ch[1]==x);
	node[x].pa=z;
	if(!isroot(y))node[z].ch[node[z].ch[1]==y]=x;
	node[y].ch[d]=node[x].ch[d^1];
	node[node[y].ch[d]].pa=y;
	node[y].pa=x,node[x].ch[d^1]=y;
	up(y);
	up(x);
}
inline void splay(int x){/*將節點x伸展到所在splay tree的根*/ 
	push_down(x);
	while(!isroot(x)){
		int y=node[x].pa;
		if(!isroot(y)){
			int z=node[y].pa;
			if((node[z].ch[0]==y)^(node[y].ch[0]==x))rotate(y);
			else rotate(x);
		}
		rotate(x);
	}
}
inline int access(int x){
	int last=0;
	while(x){
		splay(x);
		node[x].ch[1]=last;
		up(x);
		last=x;
		x=node[x].pa;
	}
	return last;/*回傳access後splay tree的根*/
}
inline void access(int x,bool is=0){/*is=0就是一般的access*/
	int last=0;
	while(x){
		splay(x);
		if(is&&!node[x].pa){
			//printf("%d\n",max(node[last].ma,node[node[x].ch[1]].ma));
		}
		node[x].ch[1]=last;
		up(x);
		last=x;
		x=node[x].pa;
	}
}
inline void query_edge(int u,int v){
	access(u);
	access(v,1);
}
inline void make_root(int x){
	access(x),splay(x);
	node[x].rev^=1;
}
inline void make_root(int x){
	node[access(x)].rev^=1;
	splay(x);
}
inline void cut(int x,int y){
	make_root(x);
	access(y);
	splay(y);
	node[y].ch[0]=0;
	node[x].pa=0;
}
inline void cut_parents(int x){
	access(x);
	splay(x);
	node[node[x].ch[0]].pa=0;
	node[x].ch[0]=0;
}
inline void link(int x,int y){
	make_root(x); 
	node[x].pa=y;
}
inline int find_root(int x){
	x=access(x);
	while(node[x].ch[0])x=node[x].ch[0];
	splay(x);
	return x;
}
inline int query(int u,int v){
/*
傳回uv路徑splay tree的根結點
這種寫法無法求LCA
*/
	make_root(u);
	return access(v);
}
inline int query_lca(int u,int v){
/*假設求鏈上點權的總和，sum是子樹的權重和，data是節點的權重*/
	access(u);
	int lca=access(v);
	splay(u);
	if(u==lca){
		//return node[lca].data+node[node[lca].ch[1]].sum;
	}else{
		//return node[lca].data+node[node[lca].ch[1]].sum+node[u].sum;
	}
}
struct EDGE{
	int a,b,w;
}e[10005]; 
int n;
std::vector<std::pair<int ,int > >G[10005];
/*first表示子節點，second表示邊的編號*/
int pa[10005],edge_node[10005];
/*pa是父母節點，暫存用的，edge_node是每個編被存在哪個點裡面的陣列*/ 
inline void bfs(int root){
/*在建構的時候把每個點都設成一個splay tree，不會壞掉*/ 
	std::queue<int > q;
	for(int i=1;i<=n;++i)pa[i]=0;
	q.push(root);
	while(q.size()){
		int u=q.front();
		q.pop();
		for(int i=0;i<(int)G[u].size();++i){
			int v=G[u][i].first;
			if(v!=pa[u]){
				pa[v]=u;
				node[v].pa=u;
				node[v].data=e[G[u][i].second].w;
				edge_node[G[u][i].second]=v;
				up(v);
				q.push(v);
			}
		}
	}
}
inline void change(int x,int b){
	splay(x);
	//node[x].data=b;
	up(x);
}

```

Code:.\������.cpp
================

```cpp
#include<vector>
#define MAXN 100005
typedef std::vector<int >::iterator VIT;
int siz[MAXN],max_son[MAXN],pa[MAXN],dep[MAXN];
/*節點大小、大小最大的孩子、父母節點、深度*/
int link_top[MAXN],link[MAXN],cnt;
/*每個點所在鍊的鏈頭、樹鏈剖分的DFS序、時間戳*/ 
std::vector<int >G[MAXN];/*用vector存樹*/ 
void find_max_son(int x){
	siz[x]=1;
	max_son[x]=-1;
	for(VIT i=G[x].begin();i!=G[x].end();++i){
		if(*i==pa[x])continue;
		pa[*i]=x;
		dep[*i]=dep[x]+1;
		find_max_son(*i);
		if(max_son[x]==-1||siz[*i]>siz[max_son[x]])max_son[x]=*i;
		siz[x]+=siz[*i];
	}
}
void build_link(int x,int top){
	link[x]=++cnt;/*記錄x點的時間戳*/
	link_top[x]=top;
	if(max_son[x]==-1)return;
	build_link(max_son[x],top);/*優先走訪最大孩子*/ 
	for(VIT i=G[x].begin();i!=G[x].end();++i){
		if(*i==max_son[x]||*i==pa[x])continue;
		build_link(*i,*i);
	}
}
inline int find_lca(int a,int b){
	/*求LCA，可以在過程中對區間進行處理*/ 
	int ta=link_top[a],tb=link_top[b];
	while(ta!=tb){
		if(dep[ta]<dep[tb]){ 
			std::swap(ta,tb);
			std::swap(a,b);
		}
		//這裡可以對a所在的鏈做區間處理  
		//區間為(link[ta],link[a]) 
		ta=link_top[a=pa[ta]];
	}
	/*最後a,b會在同一條鏈，若a!=b還要在進行一次區間處理*/ 
	return dep[a]<dep[b]?a:b;
}

```

