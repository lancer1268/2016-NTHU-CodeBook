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
	point operator+(const point &b)const{return point(x+b.x,y+b.y);}
	point operator-(const point &b)const{return point(x-b.x,y-b.y);}
	point operator*(const point &b)const{return point(x*b.x,y*b.y);}
	point operator/(const point &b)const{return point(x/b.x,y/b.y);}
	point operator+(const T &b)const{return point(x+b,y+b);}
	point operator-(const T &b)const{return point(x-b,y-b);}
	point operator*(const T &b)const{return point(x*b,y*b);}
	point operator/(const T &b)const{return point(x/b,y/b);}
	bool operator==(const point &b)const{return x==b.x&&y==b.y;}
	T dot(const point &b)const{return x*b.x+y*b.y;}
	T cross(const point &b)const{return x*b.y-y*b.x;}
	// 法向量
	point normal()const{return point(-y,x);}
	// 向量長度的平方
	T abs2()const{return dot(*this);}
};
template<typename T>
struct line{
	line(){}
	point<T> p1,p2;
	T a,b,c;/*ax+by+c=0*/
	line(const point<T>&x,const point<T>&y):p1(x),p2(y){}
	void pton(){/*轉成一般式*/ 
		a=p1.y-p2.y;
		b=p2.x-p1.x;
		c=-a*p1.x-b*p1.y;
	}
	/*點和有向直線的關係，>0左邊、=0在線上<0右邊*/
	T cross(const point<T> &p)const{return (p2-p1).cross(p-p1);}
	bool point_on_segment(const point<T>&p)const{/*點是否線段上*/ 
		return cross(p)==0&&(p1-p).dot(p2-p)<=0; 
	}
	T dis2(const point<T> &p,bool is_segment=0)const{/*點跟直線/線段的距離平方*/ 
		point<T> v=p2-p1,v1=p-p1;
		if(is_segment){
			point<T> v2=p-p2;
			if(v.dot(v1)<=0)return v1.abs2();
			if(v.dot(v2)>=0)return v2.abs2();
		}
		T tmp=v.cross(v1);
		return tmp*tmp/v.abs2();
	}
	point<T> mirror(const point<T> &p)const{/*點對直線的鏡射*/ 
		/*要先呼叫pton轉成一般式*/
		point<T> ans;
		T d=a*a+b*b;
		ans.x=(b*b*p.x-a*a*p.x-2*a*b*p.y-2*a*c)/d;
		ans.y=(a*a*p.y-b*b*p.y-2*a*b*p.x-2*b*c)/d;
		return ans;
	}
	/*直線相等*/ 
	bool equal(const line &l)const{return cross(l.p1)==0&&cross(l.p2)==0;}
	/*直線平行*/
	bool parallel(const line &l)const{return (p1-p2).cross(l.p1-l.p2)==0;}
	bool cross_seg(const line &l)const{/*直線是否交線段*/
		return (p2-p1).cross(l.p1)*(p2-p1).cross(l.p2)<=0;
	}
	/*直線相交情況，-1無限多點、1交於一點、0不相交*/
	int line_intersect(const line &l)const{return parallel(l)?(cross(l.p1)==0?-1:0):1;}
	/*線段相交情況，-1無限多點、1交於一點、0不相交*/
	int seg_intersect(const line &l)const{
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
	point<T> line_intersection(const line &l)const{/*直線交點*/ 
		point<T> a=p2-p1,b=l.p2-l.p1,s=l.p1-p1;
		//if(a.cross(b)==0)return INF;
		return p1+a*s.cross(b)/a.cross(b);
	}
	point<T> seg_intersection(const line &l)const{/*線段交點*/ 
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
	const point<T>& operator[](int id)const{return p[id];}
	T area()const{/*多邊形面積*/ 
		T ans=0;
		size_t psize=p.size();
		for(size_t i= psize-1,j=0;j<psize;i=j++)
			ans+=p[i].cross(p[j]);
		return ans/2;
	}
	point<T> center_of_mass()const{/*多邊形重心*/
		T cx=0,cy=0,w=0;
		size_t psize=p.size();
		for(size_t i=psize-1,j=0;j<psize;i=j++){
			T a=p[i].cross(p[j]);
			cx+=(p[i].x+p[j].x)*a;
			cy+=(p[i].y+p[j].y)*a;
			w+=a;
		}
		return point<T>(cx/3/w,cy/3/w);
	}
	/*點是否在凸多邊形內，是的話回傳1、在邊上回傳-1、否則回傳0*/
	int chas(const point<T> &x)const{
		T tp=0,np;
		size_t psize = p.size();
		for(size_t i=psize-1,j=0;j<psize;i=j++){
			if(!(np=(p[j]-x).cross(p[i]-x)))return -1;
			if(tp*np<0)return 0;
			tp=(np!=0)?np:tp;
		}
		return 1;
	}
	/*點是否在簡單多邊形內，是的話回傳1、在邊上回傳-1、否則回傳0*/
	int ahas(const point<T>& t)const{
		bool c=0;
		size_t psize = p.size();
		for(size_t i=0,j=psize-1;i<psize;j=i++)
			if(line<T>(p[i],p[j]).point_on_segment(t))return -1;
			else if((p[i].y>t.y)!=(p[j].y>t.y)&&
			t.x<(p[j].x-p[i].x)*(t.y-p[i].y)/(p[j].y-p[i].y)+p[i].x)
				c=!c;
		return c;
	}
	polygon cut(const line<T> &l)const{/*凸包對直線切割，得到直線l左側的凸包*/
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
	static bool graham_cmp(const point<T>& a,const point<T>& b){
		return (a.x<b.x)||(a.x==b.x&&a.y<b.y);/*凸包排序函數*/
	}
	void graham(std::vector<point<T> > &s){/*凸包*/
		sort(s.begin(),s.end(),graham_cmp);
		p.resize(s.size()+1);
		int m=0;
		for(size_t i=0;i<s.size();++i){
			while(m>=2&&(p[m-1]-p[m-2]).cross(s[i]-p[m-2])<=0)--m;
			p[m++]=s[i];
		}
		for(int i=(int)s.size()-2,t=m+1;i>=0;--i){
			while(m>=t&&(p[m-1]-p[m-2]).cross(s[i]-p[m-2])<=0)--m;
			p[m++]=s[i];
		}
		if(s.size()>1)--m; 
		p.resize(m);
	}
	static int sign(const T&x){return x>=0?1:-1;}
	static bool angle_cmp(const line<T>& A,const line<T>& B){
		point<T>a=A.p2-A.p1,b=B.p2-B.p1;
		//return atan2(a.y,a.x)<atan2(b.y,b.x); 
		int ay=sign(a.y),by=sign(b.y),ax=sign(a.x),bx=sign(b.x);
		return ay>by||(ay==by&&(ax*ay>bx*by||(ax*ay==bx*by&&a.cross(b)>0)));
	}
	int halfplane_intersection(std::vector<line<T> > &s){
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
	T area()const{
		T t=(b-a).cross(c-a)/2;
		return t>0?t:-t;
	}
	/*重心*/
	point<T> center_of_mass()const{return (a+b+c)/3;}
	point<T> circumcenter()const{/*外心*/
		static line<T> u,v;
		u.p1=(a+b)/2;
		u.p2.x=u.p1.x-a.y+b.y,u.p2.y=u.p1.y+a.x-b.x;
		v.p1=(a+c)/2;
		v.p2.x=v.p1.x-a.y+c.y,v.p2.y=v.p1.y+a.x-c.x;
		return u.line_intersection(v);
	}
	point<T> incenter()const{/*內心，用到根號*/
		T A=sqrt((b-c).abs2()),B=sqrt((a-c).abs2()),C=sqrt((a-b).abs2());
		return point<T>(A*a.x+B*b.x+C*c.x,A*a.y+B*b.y+C*c.y)/(A+B+C);
	}
	/*垂心*/
	point<T> perpencenter()const{return center_of_mass()*3-circumcenter()*2;}
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
    bool incircle(cp &c)const{return (x-c).abs2()<=r2;}
};

Circle TwoPointCircle(Circle::cp &a, Circle::cp &b) {
    Circle::p m=(a+b)/2;
    return (Circle){m,(a-m).abs2()};
}

Circle outcircle(Circle::p a, Circle::p b, Circle::p c) {
    if(TwoPointCircle(a,b).incircle(c)) return TwoPointCircle(a,b);
    if(TwoPointCircle(b,c).incircle(a)) return TwoPointCircle(b,c);
    if(TwoPointCircle(c,a).incircle(b)) return TwoPointCircle(c,a);
    Circle::p ret;
    double a1=b.x-a.x, b1=b.y-a.y, c1=(a1*a1+b1*b1)/2;
    double a2=c.x-a.x, b2=c.y-a.y, c2=(a2*a2+b2*b2)/2;
    double d = a1*b2 - a2*b1;
    ret.x=a.x+(c1*b2-c2*b1)/d;
    ret.y=a.y+(a1*c2-a2*c1)/d;
    return (Circle){ret,(ret-a).abs2()};
}
//rand required
Circle SmallestCircle(std::vector<Circle::p> &p){
    int n=p.size();
    if(n==1) return (Circle){p[0],0.0};
    if(n==2) return TwoPointCircle(p[0],p[1]);
    random_shuffle(p.begin(),p.end());
    Circle c = {p[0],0.0};
    for(int i=0;i<n;++i){
        if(c.incircle(p[i])) continue;
        c=Circle{p[i],0.0};
        for(int j=0;j<i;++j){
            if(c.incircle(p[j])) continue;
            c=TwoPointCircle(p[i],p[j]);
            for(int k=0;k<j;++k){
                if(c.incircle(p[k])) continue;
                c=outcircle(p[i],p[j],p[k]);
            }
        }
    }
    return c;
}

```

Code:.\Computational Geometry\最近點對.cpp
================

```cpp
#define INF LLONG_MAX/*預設是long long最大值*/
template<typename T>
T closest_pair(vector<point<T> >&v,vector<point<T> >&t,int l,int r){
	T dis=INF,tmd;
	if(l>=r)return dis;
	int mid=(l+r)/2;
	if((tmd=closest_pair(v,t,l,mid))<dis)dis=tmd;
	if((tmd=closest_pair(v,t,mid+1,r))<dis)dis=tmd;
	t.clear();
	for(int i=l;i<=r;++i)
		if((v[i].x-v[mid].x)*(v[i].x-v[mid].x)<dis)t.push_back(v[i]);
	sort(t.begin(),t.end(),point<T>::y_cmp);/*如果用merge_sort的方式可以O(n)*/
	for(int i=0;i<(int)t.size();++i)
		for(int j=1;j<=3&&i+j<(int)t.size();++j)
			if((tmd=(t[i]-t[i+j]).abs2())<dis)dis=tmd;
	return dis;
}
template<typename T>
inline T closest_pair(vector<point<T> > &v){
	vector<point<T> >t;
	sort(v.begin(),v.end(),point<T>::x_cmp);
	return closest_pair(v,t,0,v.size()-1);/*最近點對距離*/
}

```

Code:.\Data Structure\Dynemic KD tree.cpp
================

```cpp
#ifndef SUNMOON_DYNEMIC_KD_TREE
#define SUNMOON_DYNEMIC_KD_TREE
#include<algorithm>
#include<vector>
#include<queue>
#include<cmath>
template<typename T,size_t kd>//kd���ܦ��X�Ӻ��� 
class kd_tree{
	public:
		struct point{
			T d[kd];
			inline T dist(const point &x)const{
				T ret=0;
				for(size_t i=0;i<kd;++i)ret+=std::abs(d[i]-x.d[i]);
				return ret;
			}
			inline bool operator==(const point &p){
				for(size_t i=0;i<kd;++i){
					if(d[i]!=p.d[i])return 0;
				}
				return 1;
			}
			inline bool operator<(const point &b)const{
				return d[0]<b.d[0];
			}
		};
	private:
		struct node{
			node *l,*r;
			point pid;
			int s;
			node(const point &p):l(0),r(0),pid(p),s(1){}
			inline void up(){
				s=(l?l->s:0)+1+(r?r->s:0);
			}
		}*root;
		const double alpha,loga;
		const T INF;//�O�o�n��INF�A���ܷ��j�� 
		int maxn;
		struct __cmp{
			int sort_id;
			inline bool operator()(const node*x,const node*y)const{
				return operator()(x->pid,y->pid);
			}
			inline bool operator()(const point &x,const point &y)const{
				if(x.d[sort_id]!=y.d[sort_id])
					return x.d[sort_id]<y.d[sort_id];
				for(size_t i=0;i<kd;++i){
					if(x.d[i]!=y.d[i])return x.d[i]<y.d[i];
				}
				return 0;
			}
		}cmp;
		void clear(node *o){
			if(!o)return;
			clear(o->l);
			clear(o->r);
			delete o;
		}
		inline int size(node *o){
			return o?o->s:0;
		}
		std::vector<node*> A;
		node* build(int k,int l,int r){
			if(l>r)return 0;
			if(k==kd)k=0;
			int mid=(l+r)/2;
			cmp.sort_id=k;
			std::nth_element(A.begin()+l,A.begin()+mid,A.begin()+r+1,cmp);
			node *ret=A[mid];
			ret->l=build(k+1,l,mid-1);
			ret->r=build(k+1,mid+1,r);
			ret->up();
			return ret;
		}
		inline bool isbad(node*o){
			return size(o->l)>alpha*o->s||size(o->r)>alpha*o->s;
		}
		void flatten(node *u,typename std::vector<node*>::iterator &it){
			if(!u)return;
			flatten(u->l,it);
			*it=u;
			flatten(u->r,++it);
		}
		inline void rebuild(node*&u,int k){
			if((int)A.size()<u->s)A.resize(u->s);
			typename std::vector<node*>::iterator it=A.begin();
			flatten(u,it);
			u=build(k,0,u->s-1);
		}
		bool insert(node*&u,int k,const point &x,int dep){
			if(!u){
				u=new node(x);
				return dep<=0;
			}
			++u->s;
			cmp.sort_id=k;
			if(insert(cmp(x,u->pid)?u->l:u->r,(k+1)%kd,x,dep-1)){
				if(!isbad(u))return 1;
				rebuild(u,k);
			}
			return 0;
		}
		node *findmin(node*o,int k){
			if(!o)return 0; 
			if(cmp.sort_id==k)return o->l?findmin(o->l,(k+1)%kd):o;
			node *l=findmin(o->l,(k+1)%kd);
			node *r=findmin(o->r,(k+1)%kd);
			if(l&&!r)return cmp(l->pid,o->pid)?l:o;
			if(!l&&r)return cmp(r->pid,o->pid)?r:o;
			if(!l&&!r)return o;
			if(cmp(l->pid,r->pid))return cmp(l->pid,o->pid)?l:o;
			return cmp(r->pid,o->pid)?r:o;
		}
		bool erase(node *&u,int k,const point &x){
			if(!u)return 0;
			if(u->pid==x){
				if(u->r);
				else if(u->l){
					u->r=u->l;
					u->l=0;
				}else{
					delete u;
					u=0;
					return 1;
				}
				--u->s;
				cmp.sort_id=k;
				u->pid=findmin(u->r,(k+1)%kd)->pid;
				return erase(u->r,(k+1)%kd,u->pid);
			}
			cmp.sort_id=k;
			if(erase(cmp(x,u->pid)?u->l:u->r,(k+1)%kd,x)){
				--u->s;return 1;
			}else return 0;
		}
		inline T heuristic(const T h[])const{
			T ret=0;
			for(size_t i=0;i<kd;++i)ret+=h[i];
			return ret;
		}
		int qM;
		std::priority_queue<std::pair<T,point > >pQ;
		void nearest(node *u,int k,const point &x,T *h,T &mndist){
			if(u==0||heuristic(h)>=mndist)return;
			T dist=u->pid.dist(x),old=h[k];
			/*mndist=std::min(mndist,dist);*/
			if(dist<mndist){
				pQ.push(std::make_pair(dist,u->pid));
				if((int)pQ.size()==qM+1){
					mndist=pQ.top().first,pQ.pop();
				}
			}
			if(x.d[k]<u->pid.d[k]){
				nearest(u->l,(k+1)%kd,x,h,mndist);
				h[k]=abs(x.d[k]-u->pid.d[k]);
				nearest(u->r,(k+1)%kd,x,h,mndist);
			}else{
				nearest(u->r,(k+1)%kd,x,h,mndist);
				h[k]=abs(x.d[k]-u->pid.d[k]);
				nearest(u->l,(k+1)%kd,x,h,mndist);
			}
			h[k]=old;
		}
		std::vector<point>in_range;
		void range(node *u,int k,const point&mi,const point&ma){
			if(!u)return;
			bool is=1;
			for(int i=0;i<kd;++i)
				if(u->pid.d[i]<mi.d[i]||ma.d[i]<u->pid.d[i]){
					is=0;break;
				}
			if(is)in_range.push_back(u->pid);
			if(mi.d[k]<=u->pid.d[k])range(u->l,(k+1)%kd,mi,ma);
			if(ma.d[k]>=u->pid.d[k])range(u->r,(k+1)%kd,mi,ma);
		}
	public:
		kd_tree(const T &INF,double a=0.75):alpha(a),loga(log2(1.0/a)),INF(INF),maxn(1){}
		inline void clear(){
			clear(root),root=0,maxn=1;
		}
		inline void build(int n,const point *p){
			clear(root),A.resize(maxn=n);
			for(int i=0;i<n;++i)A[i]=new node(p[i]);
			root=build(0,0,n-1);
		}
		inline void insert(const point &x){
			insert(root,0,x,std::__lg(size(root))/loga);
			if(root->s>maxn)maxn=root->s;
		}
		inline bool erase(const point &p){
			bool d=erase(root,0,p);
			if(root&&root->s<alpha*maxn)rebuild();
			return d;
		}
		inline void rebuild(){
			if(root)rebuild(root,0);
			maxn=root->s;
		}
		inline T nearest(const point &x,int k){
			qM=k;
			T mndist=INF,h[kd]={};
			nearest(root,0,x,h,mndist);
			mndist=pQ.top().first;
			pQ=std::priority_queue<std::pair<T,point > >();
			return mndist;/*�^����x��k���I���Z��*/ 
		}
		inline const std::vector<point> &range(const point&mi,const point&ma){
			in_range.clear();
			range(root,0,mi,ma);
			return in_range;/*�^�Ǥ���mi��ma�������Ivector*/ 
		}
		inline int size(){return root?root->s:0;}
};
#endif

```

Code:.\Data Structure\kd_tree_replace_segment_tree.cpp
================

```cpp
/*kd��N�������u�q��*/
struct node{
	node *l,*r;
	point pid,mi,ma;
	int s;
	int data;
	node(const point &p,int d):l(0),r(0),pid(p),mi(p),ma(p),s(1),data(d),dmin(d),dmax(d){}
	inline void up(){
		mi=ma=pid;
		s=1;
		if(l){
			for(int i=0;i<kd;++i){
				mi.d[i]=min(mi.d[i],l->mi.d[i]);
				ma.d[i]=max(ma.d[i],l->ma.d[i]);
			}
			s+=l->s;
		}
		if(r){
			for(int i=0;i<kd;++i){
				mi.d[i]=min(mi.d[i],r->mi.d[i]);
				ma.d[i]=max(ma.d[i],r->ma.d[i]);
			}
			s+=r->s;
		}
	}
	inline void up2(){
		//��L�i�k�аO�V�W��s 
	}
	inline void down(){
		//��L�i�k�аO�U�� 
	}
}*root;

/*�ˬd�϶��]�t�Ϊ����*/
inline bool range_include(node *o,const point &L,const point &R){
	for(int i=0;i<kd;++i){ 
		if(L.d[i]>o->ma.d[i]||R.d[i]<o->mi.d[i])return 0;
	}//�u�n(L,R)�϶����Mo���϶����涰�N�^��true
	return 1;
}
inline bool range_in_range(node *o,const point &L,const point &R){
	for(int i=0;i<kd;++i){
		if(L.d[i]>o->mi.d[i]||o->ma.d[i]>R.d[i])return 0;
	}//�p�G(L,R)�϶������]�to���϶��N�^��true
	return 1;
}
inline bool point_in_range(node *o,const point &L,const point &R){
	for(int i=0;i<kd;++i){
		if(L.d[i]>o->pid.d[i]||R.d[i]<o->pid.d[i])return 0;
	}//�p�G(L,R)�϶������]�to->pid�o���I�N�^��true
	return 1;
}

/*���I�ק�A�H���I��Ȭ���*/
void update(node *u,const point &x,int data,int k=0){
	if(!u)return;
	u->down();
	if(u->pid==x){
		u->data=data;
		u->up2();
		return;
	}
	cmp.sort_id=k;
	update(cmp(x,u->pid)?u->l:u->r,x,data,(k+1)%kd);
	u->up2();
}

/*�϶��ק�*/ 
void update(node *o,const point &L,const point &R,int data){
	if(!o)return;
	o->down();
	if(range_in_range(o,L,R)){
		//�϶��i�k�аO�ק� 
		o->down();
		return;
	}
	if(point_in_range(o,L,R)){
		//�o���I�b(L,R)�϶��A���O�L�����k�l�𤣤@�w�b�϶���
		//���I�i�k�аO�ק� 
	}
	if(o->l&&range_include(o->l,L,R))update(o->l,L,R,data);
	if(o->r&&range_include(o->r,L,R))update(o->r,L,R,data);
	o->up2();
}

/*�϶��d�ߡA�H�`�M����*/ 
int query(node *o,const point &L,const point &R){
	if(!o)return 0;
	o->down();
	if(range_in_range(o,L,R))return o->sum;
	int ans=0;
	if(point_in_range(o,L,R))ans+=o->data;
	if(o->l&&range_include(o->l,L,R))ans+=query(o->l,L,R);
	if(o->r&&range_include(o->r,L,R))ans+=query(o->r,L,R);
	return ans;
}

```

Code:.\Data Structure\skew_heap.cpp
================

```cpp
#ifndef SKEW_HEAP
#define SKEW_HEAP
template<typename T,typename _Compare=std::less<T> >
class skew_heap{
	private:
		struct node{
			T data;
			node *l,*r;
			node(const T&d):data(d),l(0),r(0){}
		}*root;
		int _size;
		_Compare cmp;
		node *merge(node *a,node *b){
			if(!a||!b)return a?a:b;
			if(cmp(a->data,b->data))return merge(b,a);
			node *t=a->r;
			a->r=a->l;
			a->l=merge(b,t);
			return a;
		}
		void _clear(node *&o){
			if(o)_clear(o->l),_clear(o->r),delete o;
		}
	public:
		skew_heap():root(0),_size(0){}
		~skew_heap(){_clear(root);}
		inline void clear(){
			_clear(root);root=0;_size=0;
		}
		inline void join(skew_heap &o){
			root=merge(root,o.root);
			o.root=0;
			_size+=o._size;
			o._size=0;
		}
		inline void swap(skew_heap &o){
			node *t=root;
			root=o.root;
			o.root=t;
			int st=_size;
			_size=o._size;
			o._size=st;
		}
		inline void push(const T&data){
			_size++;
			root=merge(root,new node(data));
		}
		inline void pop(){
			if(_size)_size--;
			node *tmd=merge(root->l,root->r);
			delete root;
			root=tmd;
		}
		inline const T& top(){return root->data;}
		inline int size(){return _size;}
		inline bool empty(){return !_size;}
};
#endif

```

Code:.\Data Structure\split_merge.cpp
================

```cpp
void split(node *o,node *&a,node *&b,int k){
	if(!o)a=b=0;
	else{
		o->down();
		if(k<=size(o->l)){
			b=o;
			split(o->l,a,b->l,k);
		}else{
			a=o;
			split(o->r,a->r,b,k-size(o->l)-1);
		}
		o->up();
	}
}
node *merge(node *a,node *b){
	if(!a||!b)return a?a:b;
	static int x;
	if(x++%(a->s+b->s)<a->s){
		a->down();
		a->r=merge(a->r,b);
		a->up();
		return a;
	}else{
		b->down();
		b->l=merge(a,b->l);
		b->up();
		return b;
	}
}

```

Code:.\Data Structure\treap.cpp
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

Code:.\Flow\dinic.cpp
================

```cpp
#define MAXN 105
#define INF INT_MAX
int n;/*�I��*/
int level[MAXN],cur[MAXN];/*�h���B���e���u��*/
struct edge{
	int v,pre;
	long long cap,flow,r;
	edge(int v,int pre,long long cap):v(v),pre(pre),cap(cap),flow(0),r(cap){}
};
int g[MAXN];
std::vector<edge> e;
inline void init(){
	memset(g,-1,sizeof(int)*(n+1));
	e.clear();
}
inline void add_edge(int u,int v,long long cap,bool directed=false){
	e.push_back(edge(v,g[u],cap));
	g[u]=e.size()-1;
	e.push_back(edge(u,g[v],directed?0:cap));
	g[v]=e.size()-1;
}
inline int bfs(int s,int t){
	memset(level,0,sizeof(int)*(n+1));
	memcpy(cur,g,sizeof(int)*(n+1));
	std::queue<int >q;
	q.push(s);
	level[s]=1;
	while(q.size()){
		int u=q.front();q.pop();
		for(int i=g[u];~i;i=e[i].pre){
			if(!level[e[i].v]&&e[i].r){
				level[e[i].v]=level[u]+1;
				q.push(e[i].v);
				if(e[i].v==t)return 1;
			}
		}
	}
	return 0;
}
long long dfs(int u,int t,long long cur_flow=INF){
	if(u==t||!cur_flow)return cur_flow;
	long long df,tf=0;
	for(int &i=cur[u];~i;i=e[i].pre){
		if(level[e[i].v]==level[u]+1&&e[i].r){
			if(df=dfs(e[i].v,t,std::min(cur_flow,e[i].r))){
				e[i].flow+=df;
				e[i^1].flow-=df;
				e[i].r-=df;
				e[i^1].r+=df;
				tf+=df;
				if(!(cur_flow-=df))break;
			}
		}
	}
	if(!df)level[u]=0;
	return tf;
}
inline long long dinic(int s,int t,bool clean=true){
	if(clean){
		for(size_t i=0;i<e.size();++i){
			e[i].flow=0;
			e[i].r=e[i].cap;
		}
	}
	long long ans=0;
	while(bfs(s,t))ans+=dfs(s,t);
	return ans;
}

```

Code:.\Flow\ISAP.cpp
================

```cpp
#define MAXN 105
#define INF INT_MAX
int n;/*�I��*/
int d[MAXN],gap[MAXN],cur[MAXN];
/*�h���Bgap[i]=�h����i���I���ӼơB���e���u��*/ 
struct edge{
	int v,pre;
	long long cap,flow,r;
	edge(int v,int pre,long long cap):v(v),pre(pre),cap(cap),flow(0),r(cap){}
};
int g[MAXN];
std::vector<edge> e;
inline void init(){
	memset(g,-1,sizeof(int)*(n+1));
	e.clear();
}
inline void add_edge(int u,int v,long long cap,bool directed=false){
	e.push_back(edge(v,g[u],cap));
	g[u]=e.size()-1;
	e.push_back(edge(u,g[v],directed?0:cap));
	g[v]=e.size()-1;
}
long long dfs(int u,int s,int t,long long cur_flow=INF){
	if(u==t)return cur_flow;
	long long tf=cur_flow,df;
	for(int &i=cur[u];~i;i=e[i].pre){
		if(e[i].r&&d[u]==d[e[i].v]+1){
			df=dfs(e[i].v,s,t,std::min(tf,e[i].r));
			e[i].flow+=df;
			e[i^1].flow-=df;
			e[i].r-=df;
			e[i^1].r+=df;
			if(!(tf-=df)||d[s]==n)return cur_flow-tf;
		}
	}
	int minh=n;
	for(int i=cur[u]=g[u];~i;i=e[i].pre){
		if(e[i].r&&d[e[i].v]<minh)minh=d[e[i].v];
	}
	if(!--gap[d[u]])d[s]=n;
	else ++gap[d[u]=++minh];
	return cur_flow-tf;
}
inline long long isap(int s,int t,bool clean=true){
	memset(d,0,sizeof(int)*(n+1));
	memset(gap,0,sizeof(int)*(n+1));
	memcpy(cur,g,sizeof(int)*(n+1));
	if(clean){
		for(size_t i=0;i<e.size();++i){
			e[i].flow=0;
			e[i].r=e[i].cap;
		}
	}
	long long max_flow=0;
	for(gap[0]=n;d[s]<n;)max_flow+=dfs(s,s,t);
	return max_flow;
}

```

Code:.\Flow\MinCostMaxFlow.cpp
================

```cpp
#define MAXN 440
#define INF 999999999
struct edge{
	int v,pre;
	int cap,cost;
	edge(int v,int pre,int cap,int cost):v(v),pre(pre),cap(cap),cost(cost){}
};
int n,S,T;
int dis[MAXN],piS,ans;
bool vis[MAXN];
std::vector<edge> e;
int g[MAXN];
inline void init(){
	memset(g,-1,sizeof(int)*n);
	e.clear();
}
inline void add_edge(int u,int v,int cost,int cap,bool directed=false){
	e.push_back(edge(v,g[u],cap,cost));
	g[u]=e.size()-1;
	e.push_back(edge(u,g[v],directed?0:cap,-cost));
	g[v]=e.size()-1;
}
int augment(int u,int cur_flow){
	if(u==T||!cur_flow)return ans+=piS*cur_flow,cur_flow;
	vis[u]=1;
	int r=cur_flow,d;
	for(int i=g[u];~i;i=e[i].pre){
		if(e[i].cap&&!e[i].cost&&!vis[e[i].v]){
			d=augment(e[i].v,std::min(r,e[i].cap));
			e[i].cap-=d;
			e[i^1].cap+=d;
			if(!(r-=d))break;
		}
	}
	return cur_flow-r;
}
inline bool modlabel(){
	for(int i=0;i<n;++i)dis[i]=INF;
	dis[T]=0;
	static std::deque<int>q;
	q.push_back(T);
	while(q.size()){
		int u=q.front();
		q.pop_front();
		int dt;
		for(int i=g[u];~i;i=e[i].pre){
			if(e[i^1].cap&&(dt=dis[u]-e[i].cost)<dis[e[i].v]){
				if((dis[e[i].v]=dt)<=dis[q.size()?q.front():S]){
					q.push_front(e[i].v);
				}else q.push_back(e[i].v);
			}
		}
	}
	for(int u=0;u<n;++u){
		for(int i=g[u];~i;i=e[i].pre){
			e[i].cost+=dis[e[i].v]-dis[u];
		}
	}
	piS+=dis[S];
	return dis[S]<INF;
}
inline int mincost(){
	piS=ans=0;
	while(modlabel()){
		do memset(vis,0,sizeof(bool)*n);
		while(augment(S,INF));
	}
	return ans;
}

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

Code:.\Graph\blossom matching.cpp
================

```cpp
#define MAXN 505
vector<int>g[MAXN];
int pa[MAXN],match[MAXN],st[MAXN],S[MAXN],v[MAXN];
int t,n;
inline int lca(int x,int y){
	for(++t;;swap(x,y)){
		if(x==0)continue;
		if(v[x]==t)return x;
		v[x]=t;
		x=st[pa[match[x]]];
	}
}
#define qpush(x) q.push(x),S[x]=0
inline void flower(int x,int y,int l,queue<int> &q){
	while(st[x]!=l){
		pa[x]=y;
		if(S[y=match[x]]==1)qpush(y);
		st[x]=st[y]=l,x=pa[y];
	}
}
inline bool bfs(int x){
	for(int i=1;i<=n;++i)st[i]=i;
	memset(S+1,-1,sizeof(int)*n);
	queue<int>q;qpush(x);
	while(q.size()){
		x=q.front(),q.pop();
		for(size_t i=0;i<g[x].size();++i){
			int y=g[x][i];
			if(S[y]==-1){
				pa[y]=x,S[y]=1;
				if(!match[y]){
					for(int lst;x;y=lst,x=pa[y])
						lst=match[x],match[x]=y,match[y]=x;
					return 1;
				}
				qpush(match[y]);
			}else if(!S[y]&&st[y]!=st[x]){
				int l=lca(y,x);
				flower(y,x,l,q),flower(x,y,l,q);
			}
		}
	}
	return 0;
}
inline int blossom(){
	int ans=0;
	for(int i=1;i<=n;++i)
		if(!match[i]&&bfs(i))++ans;
	return ans;
}

```

Code:.\Graph\KM.cpp
================

```cpp
#define MAXN 100
int n;
int g[MAXN][MAXN],lx[MAXN],ly[MAXN],slack_y[MAXN];
int match_y[MAXN];
bool vx[MAXN],vy[MAXN];//要保證g是完全二分圖
bool dfs(int x,bool adjust=1){//DFS找增廣路，is=1表示要交換邊 
	if(vx[x])return 0;
	vx[x]=1;
	for(int y=0;y<n;++y){
		if(vy[y])continue;
		int t=lx[x]+ly[y]-g[x][y];
		if(t==0){
			vy[y]=1;
			if(match_y[y]==-1||dfs(match_y[y],adjust)){
				if(adjust)match_y[y]=x;
				return 1;
			}
		}else if(slack_y[y]>t)slack_y[y]=t;
	}
	return 0;
}
inline int km(){
	memset(ly,0,sizeof(int)*n);
	memset(match_y,-1,sizeof(int)*n);
	for(int x=0;x<n;++x){
		lx[x]=0;
		for(int y=0;y<n;++y){
			lx[x]=max(lx[x],g[x][y]);
		}
	}
	for(int x=0;x<n;++x){
		for(int y=0;y<n;++y)slack_y[y]=INT_MAX;
		memset(vx,0,sizeof(bool)*n);
		memset(vy,0,sizeof(bool)*n);
		if(dfs(x))continue;
		bool flag=1;
		while(flag){
			int cut=INT_MAX;
			for(int y=0;y<n;++y){
				if(!vy[y]&&cut>slack_y[y])cut=slack_y[y];
			}
			for(int j=0;j<n;++j){
				if(vx[j])lx[j]-=cut;
				if(vy[j])ly[j]+=cut;
				else slack_y[j]-=cut;
			}
			for(int y=0;y<n;++y){
				if(!vy[y]&&slack_y[y]==0){
					vy[y]=1;
					if(match_y[y]==-1||dfs(match_y[y],0)){
						flag=0;//測試成功，有增廣路 
						break;
					}
				}
			}
		}
		memset(vx,0,sizeof(bool)*n);
		memset(vy,0,sizeof(bool)*n);
		dfs(x);//最後要記得將邊翻反轉 
	}
	int ans=0;
	for(int y=0;y<n;++y)ans+=g[match_y[y]][y];
	return ans;
}

```

Code:.\Graph\Rectilinear_Steiner_tree.cpp
================

```cpp
//�����ҫ��y�̤p�ͦ���c�y��(�h���D���n��) 
#include<vector>
#include<algorithm>
#define T int
#define INF 0x3f3f3f3f
struct point{
	T x,y;
	int id;//�C���I���s�����n���@�ˡA�q0�}�l�s�� 
	point(){}
	T dist(const point &p)const{
		return std::abs(x-p.x)+std::abs(y-p.y);
	}
};
inline bool cmpx(const point &a,const point &b){
	return a.x<b.x||(a.x==b.x&&a.y<b.y);
}
struct edge{
	int u,v;
	T cost;
	edge(int u,int v,const T&c):u(u),v(v),cost(c){}
	inline bool operator<(const edge&e)const{
		return cost<e.cost;
	}
};
struct bit_node{
	T mi;
	int id;
	bit_node(const T&mi=INF,int id=-1):mi(mi),id(id){}
};
std::vector<bit_node> bit;
inline void bit_update(int i,const T&data,int id){
	for(;i;i-=i&(-i)){
		if(data<bit[i].mi)bit[i]=bit_node(data,id);
	}
}
inline int bit_find(int i,int m){
	bit_node x;
	for(;i<=m;i+=i&(-i)){
		if(bit[i].mi<x.mi)x=bit[i];
	}
	return x.id;
}
inline std::vector<edge> build_graph(int n,point p[]){
	std::vector<edge> e;//�^�Ǫ���N�i�H�ΨӨD�̤p�ͦ��� 
	for(int dir=0;dir<4;++dir){//4�خy���ܴ� 
		if(dir%2){
			for(int i=0;i<n;++i)std::swap(p[i].x,p[i].y);
		}else if(dir==2){
			for(int i=0;i<n;++i)p[i].x=-p[i].x;
		}
		std::sort(p,p+n,cmpx);
		std::vector<T>ga(n),gb;
		for(int i=0;i<n;++i)ga[i]=p[i].y-p[i].x;
		gb=ga;
		std::sort(gb.begin(),gb.end());
		gb.resize(std::unique(gb.begin(),gb.end())-gb.begin());
		int m=gb.size();
		bit=std::vector<bit_node>(m+1);
		for(int i=n-1;i>=0;--i){
			int pos=std::lower_bound(gb.begin(),gb.end(),ga[i])-gb.begin()+1;
			int ans=bit_find(pos,m);
			if(~ans){
				e.push_back(edge(p[i].id,p[ans].id,p[i].dist(p[ans])));
			}
			bit_update(pos,p[i].x+p[i].y,i);
		}
	}
	return e;
}
```

Code:.\Graph\一般圖最大權匹配.cpp
================

```cpp
#include<bits/stdc++.h>
using namespace std;
#define INF INT_MAX
#define MAXN 400
struct edge{
	int u,v,w;
	edge(){}
	edge(int u,int v,int w):u(u),v(v),w(w){}
};
int n,n_x;
edge g[MAXN*2+1][MAXN*2+1];
int lab[MAXN*2+1];
int match[MAXN*2+1],slack[MAXN*2+1],st[MAXN*2+1],pa[MAXN*2+1];
int flower_from[MAXN*2+1][MAXN+1],S[MAXN*2+1],vis[MAXN*2+1];
vector<int> flower[MAXN*2+1];
queue<int> q;
inline int e_delta(const edge &e){ // does not work inside blossoms
	return lab[e.u]+lab[e.v]-g[e.u][e.v].w*2;
}
inline void update_slack(int u,int x){
	if(!slack[x]||e_delta(g[u][x])<e_delta(g[slack[x]][x]))slack[x]=u;
}
inline void set_slack(int x){
	slack[x]=0;
	for(int u=1;u<=n;++u)
		if(g[u][x].w>0&&st[u]!=x&&S[st[u]]==0)update_slack(u,x);
}
void q_push(int x){
	if(x<=n)q.push(x);
	else for(size_t i=0;i<flower[x].size();i++)q_push(flower[x][i]);
}
inline void set_st(int x,int b){
	st[x]=b;
	if(x>n)for(size_t i=0;i<flower[x].size();++i)
			set_st(flower[x][i],b);
}
inline int get_pr(int b,int xr){
	int pr=find(flower[b].begin(),flower[b].end(),xr)-flower[b].begin();
	if(pr%2==1){//檢查他在前一層圖是奇點還是偶點
		reverse(flower[b].begin()+1,flower[b].end());
		return (int)flower[b].size()-pr;
	}else return pr;
}
inline void set_match(int u,int v){
	match[u]=g[u][v].v;
	if(u>n){
		edge e=g[u][v];
		int xr=flower_from[u][e.u],pr=get_pr(u,xr);
		for(int i=0;i<pr;++i)set_match(flower[u][i],flower[u][i^1]);
		set_match(xr,v);
		rotate(flower[u].begin(),flower[u].begin()+pr,flower[u].end());
	}
}
inline void augment(int u,int v){
	for(;;){
		int xnv=st[match[u]];
		set_match(u,v);
		if(!xnv)return;
		set_match(xnv,st[pa[xnv]]);
		u=st[pa[xnv]],v=xnv;
	}
}
inline int get_lca(int u,int v){
	static int t=0;
	for(++t;u||v;swap(u,v)){
		if(u==0)continue;
		if(vis[u]==t)return u;
		vis[u]=t;//這種方法可以不用清空v陣列 
		u=st[match[u]];
		if(u)u=st[pa[u]];
	}
	return 0;
}
inline void add_blossom(int u,int lca,int v){
	int b=n+1;
	while(b<=n_x&&st[b])++b;
	if(b>n_x)++n_x;
	lab[b]=0,S[b]=0;
	match[b]=match[lca];
	flower[b].clear();
	flower[b].push_back(lca);
	for(int x=u,y;x!=lca;x=st[pa[y]])
		flower[b].push_back(x),flower[b].push_back(y=st[match[x]]),q_push(y);
	reverse(flower[b].begin()+1,flower[b].end());
	for(int x=v,y;x!=lca;x=st[pa[y]])
		flower[b].push_back(x),flower[b].push_back(y=st[match[x]]),q_push(y);
	set_st(b,b);
	for(int x=1;x<=n_x;++x)g[b][x].w=g[x][b].w=0;
	for(int x=1;x<=n;++x)flower_from[b][x]=0;
	for(size_t i=0;i<flower[b].size();++i){
		int xs=flower[b][i];
		for(int x=1;x<=n_x;++x)
			if(g[b][x].w==0||e_delta(g[xs][x])<e_delta(g[b][x]))
				g[b][x]=g[xs][x],g[x][b]=g[x][xs];
		for(int x=1;x<=n;++x)
			if(flower_from[xs][x])flower_from[b][x]=xs;
	}
	set_slack(b);
}
inline void expand_blossom(int b){ // S[b] == 1
	for(size_t i=0;i<flower[b].size();++i)
		set_st(flower[b][i],flower[b][i]);
	int xr=flower_from[b][g[b][pa[b]].u],pr=get_pr(b,xr);
	for(int i=0;i<pr;i+=2){
		int xs=flower[b][i],xns=flower[b][i+1];
		pa[xs]=g[xns][xs].u;
		S[xs]=1,S[xns]=0;
		slack[xs]=0,set_slack(xns);
		q_push(xns);
	}
	S[xr]=1,pa[xr]=pa[b];
	for(size_t i=pr+1;i<flower[b].size();++i){
		int xs=flower[b][i];
		S[xs]=-1,set_slack(xs);
	}
	st[b]=0;
}
inline bool on_found_edge(const edge &e){
	int u=st[e.u],v=st[e.v];
	if(S[v]==-1){
		pa[v]=e.u,S[v]=1;
		int nu=st[match[v]];
		slack[v]=slack[nu]=0;
		S[nu]=0,q_push(nu);
	}else if(S[v]==0){
		int lca=get_lca(u,v);
		if(!lca){
			augment(u,v),augment(v,u);
			return true;
		}else add_blossom(u,lca,v);
	}
	return false;
}
inline bool matching(){
	memset(S+1,-1,sizeof(int)*n_x);
	memset(slack+1,0,sizeof(int)*n_x);
	q=queue<int>();
	for(int x=1;x<=n_x;++x)
		if(st[x]==x&&!match[x])pa[x]=0,S[x]=0,q_push(x);
	if(q.empty())return false;
	for(;;){
		while(q.size()){
			int u=q.front();q.pop();
			if(S[st[u]]==1)continue;
			for(int v=1;v<=n;++v)
				if(g[u][v].w>0&&st[u]!=st[v]){
					if(e_delta(g[u][v])==0){
						if(on_found_edge(g[u][v]))return true;
					}else update_slack(u,st[v]);
				}
		}
		int d=INF;
		for(int b=n+1;b<=n_x;++b)
			if(st[b]==b&&S[b]==1)d=min(d,lab[b]/2);
		for(int x=1;x<=n_x;++x)
			if(st[x]==x&&slack[x]){
				if(S[x]==-1)d=min(d,e_delta(g[slack[x]][x]));
				else if(S[x]==0)d=min(d,e_delta(g[slack[x]][x])/2);
			}
		for(int u=1;u<=n;++u){
			if(S[st[u]]==0){
				if(lab[u]<=d)return 0;
				lab[u]-=d;
			}else if(S[st[u]]==1)lab[u]+=d;
		}
		for(int b=n+1;b<=n_x;++b)
			if(st[b]==b){
				if(S[st[b]]==0)lab[b]+=d*2;
				else if(S[st[b]]==1)lab[b]-=d*2;
			}
		q=queue<int>();
		for(int x=1;x<=n_x;++x)
			if(st[x]==x&&slack[x]&&st[slack[x]]!=x&&e_delta(g[slack[x]][x])==0)
				if(on_found_edge(g[slack[x]][x]))return true;
		for(int b=n+1;b<=n_x;++b)
			if(st[b]==b&&S[b]==1&&lab[b]==0)expand_blossom(b);
	}
	return false;
}
inline pair<long long,int> weight_blossom(){
	memset(match+1,0,sizeof(int)*n);
	n_x=n;
	int n_matches=0;
	long long tot_weight=0;
	for(int u=0;u<=n;++u)st[u]=u,flower[u].clear();
	int w_max=0;
	for(int u=1;u<=n;++u)
		for(int v=1;v<=n;++v){
			flower_from[u][v]=(u==v?u:0);
			w_max=max(w_max,g[u][v].w);
		}
	for(int u=1;u<=n;++u)lab[u]=w_max;
	while(matching())++n_matches;
	for(int u=1;u<=n;++u)
		if(match[u]&&match[u]<u)
			tot_weight+=g[u][match[u]].w;
	return make_pair(tot_weight,n_matches);
}
inline void init_weight_graph(){
	for(int u=1;u<=n;++u)
		for(int v=1;v<=n;++v)
			g[u][v]=edge(u,v,0);
}
int main(){
	int m;
	scanf("%d%d",&n,&m);
	init_weight_graph();
	for(int i=0;i<m;++i){
		int u,v,w;
		scanf("%d%d%d",&u,&v,&w);
		g[u][v].w=g[v][u].w=w;
	}
	printf("%lld\n",weight_blossom().first);
	for(int u=1;u<=n;++u)printf("%d ",match[u]);puts("");
	return 0;
}

```

Code:.\Number Theory\basic.cpp
================

```cpp
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

```

Code:.\Number Theory\bit_set.cpp
================

```cpp
void sub_set(int S){
	int sub=S;
	do{
		//對某集合的子集合的處理 
		sub=(sub-1)&S;
	}while(sub!=S);
}
void k_sub_set(int k,int n){
	int comb=(1<<k)-1,S=1<<n;
	while(comb<S){
		//對大小為k的子集合的處理 
		int x=comb&-comb,y=comb+x;
		comb=((comb&~y)/x>>1)|y;
	}
}

```

Code:.\Number Theory\enumerate.cpp
================

```cpp
void all_divdown(const LL &n) { // all n/x
    for(LL a=1;a<=n;a=n/(n/(a+1))) {
        // dosomething;
    }
}

```

Code:.\Number Theory\eulerphi.cpp
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

Code:.\Number Theory\Factor.cpp
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
inline long long mod_mul(long long a,long long b,long long m){
	a%=m,b%=m;
	long long y=(long long)((double)a*b/m+0.5);/* fast for m < 2^58 */
	long long r=(a*b-y*m)%m;
	return r<0?r+m:r;
}
template<typename T>
inline T pow(T a,T b,T mod){//a^b%mod
	T ans=1;
	for(;b;a=mod_mul(a,a,mod),b>>=1)
		if(b&1)ans=mod_mul(ans,a,mod);
	return ans;
}
int sprp[3]={2,7,61};//int範圍可解
int llsprp[7]={2,325,9375,28178,450775,9780504,1795265022};//至少unsigned long long範圍
template<typename T>
inline bool isprime(T n,int *sprp,int num){
	if(n==2)return 1;
	if(n<2||n%2==0)return 0;
	int t=0;
	T u=n-1;
	for(;u%2==0;++t)u>>=1;
	for(int i=0;i<num;++i){
		T a=sprp[i]%n;
		if(a==0||a==1||a==n-1)continue;
		T x=pow(a,u,n);
		if(x==1||x==n-1)continue;
		for(int j=0;j<t;++j){
			x=mod_mul(x,x,n);
			if(x==1)return 0;
			if(x==n-1)break;
		}
		if(x==n-1)continue;
		return 0;
	}
	return 1;
}

LL func(const LL n,const LL mod,const int c) {
    return (LLmul(n,n,mod)+c+mod)%mod;
}

LL pollorrho(const LL n, const int c) {//循環節長度 
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
        while(isp[(int)n]) {
            v.push_back(isp[(int)n]);
            n/=isp[(int)n];
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

Code:.\Number Theory\FFT.cpp
================

```cpp
#ifndef SUNMOON_FFT
#define SUNMOON_FFT
#include<vector>
#include<complex>
#include<algorithm>
template<typename T,typename VT=std::vector<std::complex<T> > >
struct FFT{
	const T pi;
	FFT(const T pi=acos((T)-1)):pi(pi){}
	inline unsigned int bit_reverse(unsigned int a,int len){
		a=((a&0x55555555U)<<1)|((a&0xAAAAAAAAU)>>1);
		a=((a&0x33333333U)<<2)|((a&0xCCCCCCCCU)>>2);
		a=((a&0x0F0F0F0FU)<<4)|((a&0xF0F0F0F0U)>>4);
		a=((a&0x00FF00FFU)<<8)|((a&0xFF00FF00U)>>8);
		a=((a&0x0000FFFFU)<<16)|((a&0xFFFF0000U)>>16);
		return a>>(32-len);
	}
	inline void fft(bool is_inv,VT &in,VT &out,int N){
		int bitlen=std::__lg(N),num=is_inv?-1:1;
		for(int i=0;i<N;++i)out[bit_reverse(i,bitlen)]=in[i];
		for(int step=2;step<=N;step<<=1){
			const int mh=step>>1;
			for(int i=0;i<mh;++i){
				std::complex<T> wi=exp(std::complex<T>(0,i*num*pi/mh));
				for(int j=i;j<N;j+=step){
					int k=j+mh;
					std::complex<T> u=out[j],t=wi*out[k];
					out[j]=u+t;
					out[k]=u-t;
				}
			}
		}
		if(is_inv)for(int i=0;i<N;++i)out[i]/=N;
	}
};
#endif

```

Code:.\Number Theory\find_real_root.cpp
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

Code:.\Number Theory\formula
================

```cpp
Sigma_{d|n} phi(n) = n
Sigma_{d|n} mu(n) = (n==1)
g(n) = Sigma_{d|n} f(d) => f(n) = Sigma_{d|n} mu(d)*g(n/d)
Catalan number: (2n)!/n!/n!/(n+1)
Harmonic series H_n = ln(n) + gamma + 1/(2n) - 1/(12nn) + 1/(120nnnn)
gamma = 0.57721566490153286060651209008240243104215933593992
i-th gray code:i^(i>>1)

```

Code:.\Number Theory\Gauss_Elimination.cpp
================

```cpp
const int MAX = 300;
const double EPS = 1e-8;

double mat[MAX][MAX];
void Gauss(int n) {
	for(int i=0; i<n; i++) {
		bool ok = 0;
		for(int j=i; j<n; j++) {
			if(fabs(mat[j][i]) > EPS) {
				swap(mat[j], mat[i]);
				ok = 1;
				break;
			}
		}
		if(!ok) continue;

		double fs = mat[i][i];
		for(int j=i+1; j<n; j++) {
			double r = mat[j][i] / fs;
			for(int k=i; k<n; k++) {
				mat[j][k] -= mat[i][k] * r;
			}
		}
	}
}

```

Code:.\Number Theory\NTT.cpp
================

```cpp
2615053605667*(2^18)+1,3
15*(2^27)+1,31
479*(2^21)+1,3
7*17*(2^23)+1,3
3*3*211*(2^19)+1,5
25*(2^22)+1,3
#ifndef SUNMOON_NTT
#define SUNMOON_NTT
#include<vector>
#include<algorithm>
template<typename T,typename VT=std::vector<T> >
struct NTT{
	const T P,G;
	NTT(T p=(1<<23)*7*17+1,T g=3):P(p),G(g){}
	inline unsigned int bit_reverse(unsigned int a,int len){
		a=((a&0x55555555U)<<1)|((a&0xAAAAAAAAU)>>1);
		a=((a&0x33333333U)<<2)|((a&0xCCCCCCCCU)>>2);
		a=((a&0x0F0F0F0FU)<<4)|((a&0xF0F0F0F0U)>>4);
		a=((a&0x00FF00FFU)<<8)|((a&0xFF00FF00U)>>8);
		a=((a&0x0000FFFFU)<<16)|((a&0xFFFF0000U)>>16);
		return a>>(32-len);
	}
	inline T pow_mod(T n,T k,T m){
		T ans=1;
		for(n=(n>=m?n%m:n);k;k>>=1){
			if(k&1)ans=ans*n%m;
			n=n*n%m;
		}
		return ans;
	}
	inline void ntt(bool is_inv,VT &in,VT &out,int N){
		int bitlen=std::__lg(N);
		for(int i=0;i<N;++i)out[bit_reverse(i,bitlen)]=in[i];
		for(int step=2,id=1;step<=N;step<<=1,++id){
			T wn=pow_mod(G,(P-1)>>id,P),wi=1,u,t;
			const int mh=step>>1;
			for(int i=0;i<mh;++i){
				for(int j=i;j<N;j+=step){
					u=out[j],t=wi*out[j+mh]%P;
					out[j]=u+t;
					out[j+mh]=u-t;
					if(out[j]>=P)out[j]-=P;
					if(out[j+mh]<0)out[j+mh]+=P;
				}
				wi=wi*wn%P;
			}
		}
		if(is_inv){
			for(int i=1;i<N/2;++i)std::swap(out[i],out[N-i]);
			T invn=pow_mod(N,P-2,P);
			for(int i=0;i<N;++i)out[i]=out[i]*invn%P;
		}
	}
};
#endif

```

Code:.\Number Theory\random.cpp
================

```cpp
inline int random_int(){
	static int seed=20160424;
	return seed+=(seed<<16)+0x1db3d743;
}
inline long long random_long_long(){
	static long long seed=20160424;
	return seed+=(seed<<32)+0xdb3d742c265539d;
}

```

Code:.\Number Theory\外星模運算.cpp
================

```cpp
//a[0]^(a[1]^a[2]^...)
#include<bits/stdc++.h>
using namespace std;
#define maxn 1000000
int euler[maxn+5];
bool is_prime[maxn+5];
inline void init_euler(){
	is_prime[1]=1;//一不是質數
	for(int i=1;i<=maxn;i++)euler[i]=i;
	for(int i=2;i<=maxn;i++){
		if(!is_prime[i]){//是質數
			euler[i]--;
			for(int j=i<<1;j<=maxn;j+=i){
				is_prime[j]=1;
				euler[j]=euler[j]/i*(i-1);
			}
		}
	}
}
inline long long pow(long long a,long long b,long long mod){//a^b%mod
	long long ans=1;
	for(;b;a=a*a%mod,b>>=1)
		if(b&1)ans=ans*a%mod;
	return ans;
}
bool isless(long long *a,int n,int k){
	if(*a==1)return k>1;
	if(--n==0)return *a<k;
	int next=0;
	for(long long b=1;b<k;++next)
		b*=*a;
	return isless(a+1,n,next);
}
long long high_pow(long long *a,int n,long long mod){
	if(*a==1||--n==0)return *a%mod;
	int k=0,r=euler[mod];
	for(long long tma=1;tma!=pow(*a,k+r,mod);++k)
		tma=tma*(*a)%mod;
	if(isless(a+1,n,k))return pow(*a,high_pow(a+1,n,k),mod);
	int tmd=high_pow(a+1,n,r);
	int t=(tmd-k+r)%r;
	return pow(*a,k+t,mod);
}
long long a[1000005];
int t,mod;
int main(){
	init_euler();
	scanf("%d",&t);
	#define n 4 
	while(t--){
		for(int i=0;i<n;++i)scanf("%lld",&a[i]);
		scanf("%d",&mod);
		printf("%lld\n",high_pow(a,n,mod));
	}
	return 0;
}

```

Code:.\STL IO\input.cpp
================

```cpp
﻿inline int read()
{
    int x=0; bool f=0; char c=getchar();
    while(ch<'0'||'9'<ch)f|=ch=='-',ch=getchar();
    while('0'<=ch&&ch<='9')x=x*10-'0'+ch,ch=getchar();
    return f?-x:x;
}
inline int read()
{//輸入不可以包含 : ; < > = ? 
    int x=0; bool f=0; char c=getchar();
    while((c>>4)&3!=3)f|=ch=='-',c=getchar();
    while((c>>4)&3==3)x=x*10-'0'+c,c=getchar();
    return f?-x:x;
}



```

Code:.\String\AC自動機.cpp
================

```cpp
#ifndef SUNMOON_AHO_CORASICK_AUTOMATON
#define SUNMOON_AHO_CORASICK_AUTOMATON
#include<queue>
#include<vector>
template<char L='a',char R='z'>
class ac_automaton{
	private:
		struct joe{
			int next[R-L+1],fail,efl,ed,cnt_dp,vis;
			joe():ed(0),cnt_dp(0),vis(0){
				for(int i=0;i<=R-L;++i)next[i]=0;
			}
		};
	public:
		std::vector<joe> S;
		std::vector<int> q;
		int qs,qe,vt;
		ac_automaton():S(1),qs(0),qe(0),vt(0){}
		inline void clear(){
			q.clear();
			S.resize(1);
			for(int i=0;i<=R-L;++i)S[0].next[i]=0;
			S[0].cnt_dp=S[0].vis=qs=qe=vt=0;
		}
		inline void insert(const char *s){
			int o=0;
			for(int i=0,id;s[i];++i){
				id=s[i]-L;
				if(!S[o].next[id]){
					S.push_back(joe());
					S[o].next[id]=S.size()-1;
				}
				o=S[o].next[id];
			}
			++S[o].ed;
		}
		inline void build_fail(){
			S[0].fail=S[0].efl=-1;
			q.clear();
			q.push_back(0);
			++qe;
			while(qs!=qe){
				int pa=q[qs++],id,t;
				for(int i=0;i<=R-L;++i){
					t=S[pa].next[i];
					if(!t)continue;
					id=S[pa].fail;
					while(~id&&!S[id].next[i])id=S[id].fail;
					S[t].fail=~id?S[id].next[i]:0;
					S[t].efl=S[S[t].fail].ed?S[t].fail:S[S[t].fail].efl;
					q.push_back(t);
					++qe;
				}
			}
		}
		/*DP出每個前綴在字串s出現的次數並傳回所有字串被s匹配成功的次數O(N+M)*/
		inline int match_0(const char *s){
			int ans=0,id,p=0,i;
			for(i=0;s[i];++i){
				id=s[i]-L;
				while(!S[p].next[id]&&p)p=S[p].fail;
				if(!S[p].next[id])continue;
				p=S[p].next[id];
				++S[p].cnt_dp;/*匹配成功則它所有後綴都可以被匹配(DP計算)*/
			}
			for(i=qe-1;i>=0;--i){
				ans+=S[q[i]].cnt_dp*S[q[i]].ed;
				if(~S[q[i]].fail)S[S[q[i]].fail].cnt_dp+=S[q[i]].cnt_dp;
			}
			return ans;
		}
		/*多串匹配走efl邊並傳回所有字串被s匹配成功的次數O(N*M^1.5)*/ 
		inline int match_1(const char *s)const{
			int ans=0,id,p=0,t;
			for(int i=0;s[i];++i){
				id=s[i]-L;
				while(!S[p].next[id]&&p)p=S[p].fail;
				if(!S[p].next[id])continue;
				p=S[p].next[id];
				if(S[p].ed)ans+=S[p].ed;
				for(t=S[p].efl;~t;t=S[t].efl){
					ans+=S[t].ed;/*因為都走efl邊所以保證匹配成功*/
				}
			}
			return ans;
		}
		/*枚舉(s的子字串∩A)的所有相異字串各恰一次並傳回次數O(N*M^(1/3))*/
		inline int match_2(const char *s){
			int ans=0,id,p=0,t;
			++vt;
			/*把戳記vt+=1，只要vt沒溢位，所有S[p].vis==vt就會變成false
			這種利用vt的方法可以O(1)歸零vis陣列*/ 
			for(int i=0;s[i];++i){
				id=s[i]-L;
				while(!S[p].next[id]&&p)p=S[p].fail;
				if(!S[p].next[id])continue;
				p=S[p].next[id];
				if(S[p].ed&&S[p].vis!=vt){
					S[p].vis=vt;
					ans+=S[p].ed;
				}
				for(t=S[p].efl;~t&&S[t].vis!=vt;t=S[t].efl){
					S[t].vis=vt;
					ans+=S[t].ed;/*因為都走efl邊所以保證匹配成功*/
				}
			}
			return ans;
		}
		/*把AC自動機變成真的自動機*/
		inline void evolution(){
			for(qs=1;qs!=qe;){
				int p=q[qs++];
				for(int i=0;i<=R-L;++i)
					if(S[p].next[i]==0)S[p].next[i]=S[S[p].fail].next[i];
			}
		}
};
#endif

```

Code:.\String\hash.cpp
================

```cpp
#define MAXN 1000000
#define prime_mod 1073676287
/*prime_mod 必須要是質數*/
typedef long long T;
char s[MAXN+5];
T h[MAXN+5];/*hash陣列*/ 
T h_base[MAXN+5];/*h_base[n]=(prime^n)%prime_mod*/ 
inline void hash_init(int len,T prime=0xdefaced){
	h_base[0]=1;
	for(int i=1;i<=len;++i){
		h[i]=(h[i-1]*prime+s[i-1])%prime_mod;
		h_base[i]=(h_base[i-1]*prime)%prime_mod;
	}
}
inline T get_hash(int l,int r){/*閉區間寫法，設編號為0 ~ len-1*/
	return (h[r+1]-(h[l]*h_base[r-l+1])%prime_mod+prime_mod)%prime_mod;
}

```

Code:.\String\KMP.cpp
================

```cpp
/*產生fail function*/ 
inline void kmp_fail(char *s,int len,int *fail){
	int id=-1;
	fail[0]=-1;
	for(int i=1;i<len;++i){
		while(~id&&s[id+1]!=s[i])id=fail[id];
		if(s[id+1]==s[i])++id;
		fail[i]=id;
	}
}
/*以字串B匹配字串A，傳回匹配成功的數量(用B的fail)*/
inline int kmp_match(char *A,int lenA,char *B,int lenB,int *fail){
	int id=-1,ans=0;
	for(int i=0;i<lenA;++i){
		while(~id&&B[id+1]!=A[i])id=fail[id];
		if(B[id+1]==A[i])++id;
		if(id==lenB-1){/*匹配成功*/
			++ans;
			id=fail[id];
		}
	}
	return ans;
}

```

Code:.\String\manacher.cpp
================

```cpp
//原字串: asdsasdsa 
//先把字串變成這樣: @a#s#d#s#a#s#d#s#a#
inline void manacher(char *s,int len,int *z){
	int l=0,r=0;
	for(int i=1;i<len;++i){
		z[i]=r>i?min(z[2*l-i],r-i):1;
		while(s[i+z[i]]==s[i-z[i]])++z[i];
		if(z[i]+i>r)r=z[i]+i,l=i;
	}
}

```

Code:.\String\suffix array lcp.cpp
================

```cpp
#define radix_sort(x,y){\
	for(i=0;i<A;++i)c[i]=0;\
	for(i=0;i<len;++i)c[x[y[i]]]++;\
	for(i=1;i<A;++i)c[i]+=c[i-1];\
	for(i=len-1;i>=0;--i)sa[--c[x[y[i]]]]=y[i];\
}
inline void suffix_array(const char *s,int len,int *sa,int *rank,int *tmp,int *c){
	int A='z'+1,i,k,id,*t;
	for(i=0;i<len;++i){
		tmp[i]=i;
		rank[i]=s[i];
	}
	radix_sort(rank,tmp);
	for(k=1;id<len-1;k<<=1){
		id=0;
		for(i=len-k;i<len;++i)tmp[id++]=i;
		for(i=0;i<len;++i){
			if(sa[i]>=k)tmp[id++]=sa[i]-k;
		}
		radix_sort(rank,tmp);
		t=rank;rank=tmp;tmp=t;
		id=0;
		rank[sa[0]]=0;
		for(i=1;i<len;++i){
			if(tmp[sa[i-1]]!=tmp[sa[i]]||sa[i-1]+k>=len||tmp[sa[i-1]+k]!=tmp[sa[i]+k])++id;
			rank[sa[i]]=id;
		}
		A=id+1;
	}
}
#undef radix_sort
//h:高度數組 sa:後綴數組 rank:排名 
inline void suffix_array_lcp(const char *s,int len,int *h,int *sa,int *rank){
	for(int i=0;i<len;++i)rank[sa[i]]=i;
	for(int i=0,k=0;i<len;++i){
		if(rank[i]==0)continue;
		if(k)--k;
		while(s[i+k]==s[sa[rank[i]-1]+k])++k;
		h[rank[i]]=k;
	}
	h[0]=0;
}

```

Code:.\String\Z.cpp
================

```cpp
inline void z_alg(char *s,int len,int *z){
	int l=0,r=0;
	z[0]=len;
	for(int i=1;i<len;++i){
		z[i]=i>r?0:(i-l+z[i-l]<z[l]?z[i-l]:r-i+1);
		while(i+z[i]<len&&s[i+z[i]]==s[z[i]])++z[i];
		if(i+z[i]-1>r)r=i+z[i]-1,l=i;
	}
}

```

Code:.\Tarjan\橋連通分量.cpp
================

```cpp
#include<vector>
#include<algorithm>
#define N 50005
std::vector<int> s[N];
int low[N],v[N]={0},Time=0,ans=0,cnt=0;
int st[N],top=0,contract[N];/*BCC用*/ 
void dfs(int x,int p){/*x當前點，p父親*/
	low[x]=v[x]=++Time;
	st[top++]=x;/*BCC用*/ 
	for(int i=0,r;i<(int)s[x].size();++i){
		if(!v[r=s[x][i]])dfs(r,x);
		if(r!=p)low[x]=std::min(low[x],low[r]);
		if(v[x]<low[r])++ans;/*這條邊是橋*/
	}/*傳回橋的數量*/ 
	if(v[x]==low[x]){/*處理BCC*/ 
		int u;
		do{
			contract[u=st[--top]]=cnt;/*每個點所在的BCC*/
		}while(x!=u);
		++cnt;
	}
}

```

Code:.\Tarjan\雙連通分量&割點.cpp
================

```cpp
#include<vector>
#include<stack>
#include<algorithm>
#define N 500005 
struct Edge{
	int a,b;
};
std::vector<int> s[N];
std::vector<int> bcc[N];/*存每塊雙連通分量的點*/
int low[N],v[N]={0},Time=0;
int bcc_id[N],bcc_cnt;/*割點的bcc_id沒意義*/
bool is_cut[N];/*是否為割點*/ 
std::stack<Edge,std::vector<Edge> > st;
void dfs(int x,int p){/*x當前點，p父親*/
	int i,r,is=0,child=0;
	low[x]=v[x]=++Time;
	for(i=0;i<(int)s[x].size();++i){
		Edge e=(Edge){x,r=s[x][i]};
		if(!v[r]){
			st.push(e);
			dfs(r,x),++child;
			low[x]=std::min(low[x],low[r]);
			if(v[x]<=low[r]){
				is_cut[x]=1;
				++bcc_cnt;
				bcc[bcc_cnt].clear();
				for(;;){
					Edge u=st.top();
					st.pop();
					if(bcc_id[u.a]!=bcc_cnt){
						bcc[bcc_cnt].push_back(u.a);
						bcc_id[u.a]=bcc_cnt;
					}
					if(bcc_id[u.b]!=bcc_cnt){
						bcc[bcc_cnt].push_back(u.b);
						bcc_id[u.b]=bcc_cnt;
					}
					if(u.a==x&&u.b==r)break; 
				}
			}
		}else if(v[r]<v[x]&&r!=p){/*反向邊*/ 
			low[x]=std::min(low[x],v[r]);
		}
	}
	if(x==p&&child<2)is_cut[x]=0;/*x是dfs樹的根要特判*/ 
}
inline void bcc_init(int n){/*使用前如果有東西要清掉*/
	Time=bcc_cnt=0;
	for(int i=1;i<=n;++i){
		v[i]=0;
		is_cut[i]=0;
		bcc_id[i]=0;
	}
}

```

Code:.\Tree problem\HeaveLight.cpp
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

Code:.\Tree problem\link cut tree.cpp
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

