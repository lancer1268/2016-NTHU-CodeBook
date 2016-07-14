#define MAXN 805
#define INF 0x3f3f3f3f3f3f3f3fll
struct edge{
	int u,v;
	long long cap,flow,cost;
	edge(int u,int v,long long cap,long long cost):u(u),v(v),cap(cap),flow(0),cost(cost){}
};
int n,S,T;
long long piS,ans;
long long dis[MAXN];
bool vis[MAXN];
vector<edge> e;
vector<unsigned> g[MAXN];
inline void init(){
	for(int i=0;i<n;++i)g[i].clear();
	e.clear();
}
inline void add_edge(int u,int v,long long cost,long long cap){
	e.push_back(edge(u,v,cap,cost));
	e.push_back(edge(v,u,0,-cost));
	g[u].push_back(e.size()-2);
	g[v].push_back(e.size()-1);
}
long long augment(int u,long long cap){
	if(u==T)return ans+=piS*cap,cap;
	vis[u]=1;
	long long r=cap;
	for(size_t i=0;i<g[u].size();++i){
		int id=g[u][i];
		if(e[id].cap&&!e[id].cost&&!vis[e[id].v]){
			long long d=augment(e[id].v,r<e[id].cap?r:e[id].cap);
			e[id].cap-=d;
			e[id].flow+=d;
			e[id^1].cap+=d;
			e[id^1].flow-=d;
			if(!(r-=d))return cap;
		}
	}
	return cap-r;
}
inline bool modlabel(){
	for(int i=0;i<=n;++i)dis[i]=INF;
	dis[T]=0;
	deque<int>q;
	q.push_back(T);
	while(q.size()){
		long long dt;
		int u=q.front();
		q.pop_front();
		for(size_t i=0;i<g[u].size();++i){
			int id=g[u][i];
			if(e[id^1].cap&&(dt=dis[u]-e[id].cost)<dis[e[id].v]){
				if((dis[e[id].v]=dt)<=dis[q.size()?q.front():0]){
					q.push_front(e[id].v);
				}else q.push_back(e[id].v);
			}
		}
	}
	for(int u=0;u<n;++u){
		for(size_t i=0;i<g[u].size();++i){
			e[g[u][i]].cost+=dis[e[g[u][i]].v]-dis[u];
		}
	}
	piS+=dis[S];
	return dis[S]<INF;
}
inline long long mincost(){
	piS=ans=0;
	while(modlabel())
		do memset(vis,0,sizeof(bool)*n);
		while(augment(S,INF));
	return ans;
}
