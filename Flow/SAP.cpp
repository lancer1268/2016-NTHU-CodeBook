#include<string.h>
#include<algorithm>
#define MAXN 100
int n,m;/*點數、邊數*/
int g[MAXN+5][MAXN+5],f[MAXN+5][MAXN+5],r[MAXN+5][MAXN+5];
/*容量上限、流量、剩餘容量*/ 
int d[MAXN+5],gap[MAXN+5];
/*層次、gap[i]=層次為i的點之個數*/ 
int dfs(int x,int flow,int s,int t){
	if(x==t)return flow;
	int tmd=flow,tmp;
	for(int i=1;i<=n;++i){
		if(r[x][i]&&d[x]==d[i]+1){
			tmp=dfs(i,std::min(tmd,r[x][i]),s,t);
			f[x][i]+=tmp;
			f[i][x]-=tmp;
			r[x][i]=g[x][i]-f[x][i];
			r[i][x]=g[i][x]-f[i][x];
			if(!(tmd-=tmp))return flow;
		}
	}
	if(!--gap[d[x]])d[s]=n;
	++gap[++d[x]];
	return flow-tmd;
}
inline int sap(int s,int t){
	memset(f,0,sizeof(f));
	memset(d,0,sizeof(d));
	memset(gap,0,sizeof(gap));
	memcpy(r,g,sizeof(g));
	int ans=0;
	for(gap[0]=n;d[s]<n;)ans+=dfs(s,INT_MAX,s,t);
	return ans;
}
