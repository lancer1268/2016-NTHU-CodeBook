#include<string.h>
#include<queue>
#include<algorithm>
#define MAXN 100
int n,m;/*點數、邊數*/
int g[MAXN+5][MAXN+5],f[MAXN+5][MAXN+5],r[MAXN+5][MAXN+5];
/*容量上限、流量、剩餘容量*/ 
int level[MAXN+5];/*層次*/ 
inline int bfs(int s,int t){
	memset(level,0,sizeof(level));
	std::queue<int >q;
	q.push(s);
	level[s]=1;
	while(q.size()){
		int x=q.front();q.pop();
		for(int i=1;i<=n;i++){
			if(!level[i]&&r[x][i]){
				level[i]=level[x]+1;
				q.push(i);
			}
		}
	}
	return level[t]!=0;
}
int dfs(int x,int flow,int t){
	if(x==t)return flow;
	int tmd=flow,tmp;
	for(int i=1;i<=n&&tmd;++i){
		if(level[i]==level[x]+1&&r[x][i]){
			tmp=dfs(i,std::min(tmd,r[x][i]),t);
			f[x][i]+=tmp;
			f[i][x]-=tmp;
			r[x][i]=g[x][i]-f[x][i];
			r[i][x]=g[i][x]-f[i][x];
			if(!(tmd-=tmp))return flow;
		}
	}
	return flow-tmd;
}
inline int dinic(int s,int t){
	memset(f,0,sizeof(f));
	memcpy(r,g,sizeof(g));
	int ans=0,tmd=0;
	while(bfs(s,t)){
		while((tmd=dfs(s,INT_MAX,t)))ans+=tmd;
	}
	return ans;
}
