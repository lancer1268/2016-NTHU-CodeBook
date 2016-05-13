#include<vector>
#include<stack>
#include<algorithm>
#define N 500005 
struct Edge{
	int a,b;
};
std::vector<int> s[N];
std::vector<int> bcc[N];/*存每塊雙連通分量的點*/
int low[N],v[N]={0},Time=0;
int bcc_id[N],bcc_cnt;/*割點的bcc_id沒意義*/
bool is_cut[N];/*是否為割點*/ 
std::stack<Edge,std::vector<Edge> > st;
void dfs(int x,int p){/*x當前點，p父親*/
	int i,r,is=0,child=0;
	low[x]=v[x]=++Time;
	for(i=0;i<(int)s[x].size();++i){
		Edge e=(Edge){x,r=s[x][i]};
		if(!v[r]){
			st.push(e);
			dfs(r,x),++child;
			low[x]=std::min(low[x],low[r]);
			if(v[x]<=low[r]){
				is_cut[x]=1;
				++bcc_cnt;
				bcc[bcc_cnt].clear();
				for(;;){
					Edge u=st.top();
					st.pop();
					if(bcc_id[u.a]!=bcc_cnt){
						bcc[bcc_cnt].push_back(u.a);
						bcc_id[u.a]=bcc_cnt;
					}
					if(bcc_id[u.b]!=bcc_cnt){
						bcc[bcc_cnt].push_back(u.b);
						bcc_id[u.b]=bcc_cnt;
					}
					if(u.a==x&&u.b==r)break; 
				}
			}
		}else if(v[r]<v[x]&&r!=p){/*反向邊*/ 
			low[x]=std::min(low[x],v[r]);
		}
	}
	if(x==p&&child<2)is_cut[x]=0;/*x是dfs樹的根要特判*/ 
}
inline void bcc_init(int n){/*使用前如果有東西要清掉*/
	Time=bcc_cnt=0;
	for(int i=1;i<=n;++i){
		v[i]=0;
		is_cut[i]=0;
		bcc_id[i]=0;
	}
}
