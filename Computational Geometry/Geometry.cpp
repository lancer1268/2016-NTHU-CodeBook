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
