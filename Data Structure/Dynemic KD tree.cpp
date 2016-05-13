#ifndef SUNMOON_DYNEMIC_KD_TREE
#define SUNMOON_DYNEMIC_KD_TREE
#include<algorithm>
#include<queue>
#include<cmath>
#define T int
#define kd 2 /*有kd維*/
#define MAXN 100005 /*最大可以插入的數量*/
#define INF 999999999
const double alpha=0.75,loga=log2(1.0/alpha);
int sort_id;
struct point{
	T d[kd];
	int id;
	inline T dist(const point &x)const{
		T ret=0;
		for(int i=0;i<kd;++i)ret+=std::abs(d[i]-x.d[i]);
		return ret;
	}
	inline bool operator<(const point &b)const{
		return d[sort_id]<b.d[sort_id];
	}
};
class kd_tree{
	private:
		struct node{
			node *l,*r;
			point pid;
			int s;
			node(const point &p):l(0),r(0),pid(p),s(1){}
			inline void up(){
				s=(l?l->s:0)+1+(r?r->s:0);
			}
		}*root,*A[MAXN];
		int qM;
		std::priority_queue<std::pair<T,point > >pQ;
		struct cmp{
			inline bool operator()(const node*x,const node*y)const{
				return x->pid<y->pid;
			}
		};
		void clear(node *o){
			if(!o)return;
			clear(o->l);
			clear(o->r);
			delete o;
		}
		inline int size(node *o){
			return o?o->s:0;
		}
		node* build(int k,int l,int r){
			if(l>r)return 0;
			if(k==kd)k=0;
			int mid=(l+r)/2;
			sort_id=k;
			std::nth_element(A+l,A+mid,A+r+1,cmp());
			node *ret=A[mid];
			ret->l=build(k+1,l,mid-1);
			ret->r=build(k+1,mid+1,r);
			ret->up();
			return ret;
		}
		inline bool isbad(node*o){
			return size(o->l)>alpha*o->s||size(o->r)>alpha*o->s;
		}
		void flatten(node *u,node **&buf){
			if(!u)return;
			flatten(u->l,buf);
			*buf=u,++buf;
			flatten(u->r,buf);
		}
		bool insert(node*&u,int k,const point &x,int dep){
			if(!u){
				u=new node(x);
				return dep<=0;
			}
			++u->s;
			if(insert(x.d[k]<u->pid.d[k]?u->l:u->r,(k+1)%kd,x,dep-1)){
				if(!isbad(u))return 1;
				node **ptr=&A[0];
				flatten(u,ptr);
				u=build(k,0,u->s-1);
			}
			return 0;
		}
		inline int heuristic(const int h[])const{
			int ret=0;
			for(int i=0;i<kd;++i)ret+=h[i];
			return ret;
		}
		void nearest(node *u,int k,const point &x,T *h,T &mndist){
			if(u==0||heuristic(h)>=mndist)return;
			point now=u->pid;
			int dist=u->pid.dist(x),old=h[k];
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
			if(mi.d[k]>=u->pid.d[k])range(u->r,(k+1)%kd,mi,ma);
		}
	public:
		inline void clear(){
			clear(root),root=0;
		}
		inline void build(int n,const point *p){
			for(int i=0;i<n;++i)A[i]=new node(p[i]);
			root=build(0,0,n-1);
		}
		inline void insert(const point &x){
			insert(root,0,x,std::__lg(size(root))/loga);
		}
		inline T nearest(const point &x,int k){
			qM=k;
			T mndist=INF,h[kd]={};
			nearest(root,0,x,h,mndist);
			mndist=pQ.top().first;
			pQ=std::priority_queue<std::pair<T,point > >();
			return mndist;/*回傳離x第k近的點的距離*/ 
		}
		inline const std::vector<point> &range(const point&mi,const point&ma){
			in_range.clear();
			range(root,0,mi,ma);
			return in_range;/*回傳介於mi到ma之間的點vector*/ 
		}
};
#undef T
#undef kd
#undef MAXN
#undef INF
#endif
