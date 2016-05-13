#include<string.h>
#include<algorithm>
#define MAXN 100+5
int n,m;/*點數、邊數*/
/*容量上限、流量、剩餘容量*/ 
int g[MAXN][MAXN],f[MAXN][MAXN],r[MAXN][MAXN];
/*層次、gap[i]=層次為i的點之個數*/ 
int d[MAXN],gap[MAXN];
int dfs(int x,int flow,int s,int t){
	if(x==t)return flow;
	int rf=flow,tmp;
	for(int i=1;i<=n;++i){
		if(r[x][i]&&d[x]==d[i]+1){
			tmp=dfs(i,std::min(rf,r[x][i]),s,t);
			f[x][i]+=tmp;
			f[i][x]-=tmp;
			r[x][i]=g[x][i]-f[x][i];
			r[i][x]=g[i][x]-f[i][x];
			if(!(rf-=tmp))return flow;
		}
	}
	if(!--gap[d[x]])d[s]=n;
	++gap[++d[x]];
	return flow-rf;
}
int sap(int s,int t){
	memset(f,0,sizeof(f));
	memset(d,0,sizeof(d));
	memset(gap,0,sizeof(gap));
	memcpy(r,g,sizeof(g));
	int ans=0;
	for(gap[0]=n;d[s]<n;)ans+=dfs(s,INT_MAX,s,t);
	return ans;
}
