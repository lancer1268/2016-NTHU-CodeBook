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
