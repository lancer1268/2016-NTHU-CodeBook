#ifndef SUNMOON_KD_TREE
#define SUNMOON_KD_TREE
#include<algorithm>
#include<queue>
#define T int
#define kd 2 /*有kd維*/
#define MAXN 200005 /*最大可以插入的數量*/
#define INF 999999999
int sort_id;
struct point{
	T d[kd];
	int id;
	inline T dist(const point &x)const{
		T ret=0;
		for(int i=0;i<kd;++i)ret+=std::abs(d[i]-x.d[i]);
		return ret;
	}
	inline bool operator==(const point &b)const{
		for(int i=0;i<kd;++i)if(d[i]!=b.d[i])return 0;
		return 1;
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
			node(const point &p):l(0),r(0),pid(p){}
		}*root,*A[MAXN];
		int s;
		struct cmp{
			inline bool operator()(const node*x,const node*y)const{
				return x->pid<y->pid;
			}
		};
		node* build(int k,int l,int r){
			if(l>r)return 0;
			if(k==kd)k=0;
			int mid=(l+r)/2;
			sort_id=k;
			std::nth_element(A+l,A+mid,A+r+1,cmp());
			node *ret=A[mid];
			ret->l=build(k+1,l,mid-1);
			ret->r=build(k+1,mid+1,r);
			return ret;
		}
		node **mnp;
		int mnk;
		void findmin(node*&o,int d,int k){
			if(!o)return;
			if(!mnp||o->pid.d[d]<(*mnp)->pid.d[d]){
				mnp=&o;
				mnk=k;
			}
			findmin(o->l,d,(k+1)%kd);
			if(d==k)return;
			findmin(o->r,d,(k+1)%kd);
		}
		inline int heuristic(const int h[])const{
			int ret=0;
			for(int i=0;i<kd;++i)ret+=h[i];
			return ret;
		}
		void nearest(node *&u,int k,const point &x,T *h,T &mndist){
			if(u==0||heuristic(h)>=mndist)return;
			point now=u->pid;
			int dist=u->pid.dist(x),old=h[k];
			if(dist<mndist){
				mnp=&u;
				mnk=k;
				if(!(mndist=dist))return;
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
	public:
		inline void clear(){
			for(int i=0;i<s;++i)delete A[i];
			root=0;
			s=0;
		}
		inline void build(int n,const point *p){
			s=n;
			for(int i=0;i<n;++i)A[i]=new node(p[i]);
			root=build(0,0,n-1);
		}
		inline bool erase(point p){
			T mndist=1,h[kd]={};
			nearest(root,0,p,h,mndist);
			if(mndist)return 0;
			for(node **o=mnp;;){
				if((*o)->r);
				else if((*o)->l){
					(*o)->r=(*o)->l;
					(*o)->l=0;
				}else{
					(*o)=0;
					return 1;
				}
				mnp=0;
				findmin((*o)->r,mnk,(mnk+1)%kd);
				(*o)->pid=(*mnp)->pid;
				o=mnp;
			}
		}
		inline T nearest(const point &x){
			T mndist=INF,h[kd]={};
			nearest(root,0,x,h,mndist);
			return mndist;/*回傳離x第k近的點的距離*/ 
		}
};
#undef T
#undef kd
#undef MAXN
#undef INF
#endif
