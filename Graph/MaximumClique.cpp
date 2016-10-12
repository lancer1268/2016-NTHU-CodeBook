const int MAXN=105;
int N;
bool G[MAXN][MAXN];
int Set[MAXN],DP[MAXN],Ans;
inline bool is_clique(const int end,const int point){
	for(int i=1;i<end;++i){
		if(!G[Set[i]][point])return false;
	}
	return true;
}
void dfs(int depth,int now){
	if(depth+N-now+1<=Ans||depth+DP[now]<=Ans)return;
	for(int i=now;i<=N;++i){
		if(is_clique(depth+1,i)){
			Set[depth+1]=i;
			dfs(depth+1,i+1);
		}
	}
	if(depth>Ans)Ans=depth;
}
inline int max_clique(){
	memset(DP,0,sizeof(DP));
	Ans=0;
	DP[N]=1;
	for(int i=N-1;i>=1;--i){
		Set[1]=i;
		dfs(1,i+1);
		DP[i]=Ans;
	}
	return DP[1];
}
