#include<string.h>
#include<queue>
#include<algorithm>
#define MAXN 100+5
int n,m;/*點數、邊數*/
/*容量上限、流量、剩餘容量*/ 
int g[MAXN][MAXN],f[MAXN][MAXN],r[MAXN][MAXN];
int level[MAXN];/*層次*/ 
int bfs(int s,int t){
	memset(level,0,sizeof(level));
	std::queue<int>q;
	q.push(s);
	level[s]=1;
	while(!q.empty()){
		int x=q.front();q.pop();
		for(int i=1;i<=n;++i){
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
	int rf=flow,tmp;
	for(int i=1;i<=n&&rf;++i){
		if(level[i]==level[x]+1&&r[x][i]){
			tmp=dfs(i,std::min(rf,r[x][i]),t);
			f[x][i]+=tmp;
			f[i][x]-=tmp;
			r[x][i]=g[x][i]-f[x][i];
			r[i][x]=g[i][x]-f[i][x];
			if(!(rf-=tmp))return flow;
		}
	}
	return flow-rf;
}
int dinic(int s,int t){
	memset(f,0,sizeof(f));
	memcpy(r,g,sizeof(g));
	int ret=0,f=0;
	while(bfs(s,t)){
		while((f=dfs(s,INT_MAX,t)))ret+=f;
	}
	return ret;
}
