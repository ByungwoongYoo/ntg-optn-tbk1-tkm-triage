#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
inline int cd(int a,int b){int d=abs(a-b);return min(d,100-d);} 
struct Search{
 vector<int>Nv,O; long long nodes[15]{}; long long remsets=0; bool found=false,timed=false; vector<int>sol; int K,first; double limit; chrono::steady_clock::time_point t0;
 bool timeout(){if(chrono::duration<double>(chrono::steady_clock::now()-t0).count()>limit){timed=true;return true;}return false;}
 bool adddfs(const vector<int>&R,const unsigned char cv[100][51],vector<int>&ad,int st,int need,int cnt[51]){
  if(timeout())return false;if(!need){sol=R;sol.insert(sol.end(),ad.begin(),ad.end());sort(sol.begin(),sol.end());found=true;return true;}if((int)O.size()-st<need)return false;
  for(int ii=st;ii<=(int)O.size()-need&&!found;ii++){int x=O[ii],nc[51];memcpy(nc,cnt,sizeof(nc));bool ok=true;for(int d=1;d<=50;d++){nc[d]+=cv[x][d];if(nc[d]>(d==50?1:2)){ok=false;break;}}if(!ok)continue;for(int y:ad){int d=cd(x,y);if(++nc[d]>(d==50?1:2)){ok=false;break;}}if(!ok)continue;nodes[ad.size()+1]++;ad.push_back(x);if(adddfs(R,cv,ad,ii+1,need-1,nc))return true;ad.pop_back();if(timeout())return false;}return false;
 }
 void process(const vector<int>&remidx){remsets++;bool rem[14]{};for(int i:remidx)rem[i]=true;vector<int>R;for(int i=0;i<14;i++)if(!rem[i])R.push_back(Nv[i]);int cnt[51]{};bool ok=true;for(int i=0;i<(int)R.size()&&ok;i++)for(int j=i+1;j<(int)R.size();j++){int d=cd(R[i],R[j]);if(++cnt[d]>(d==50?1:2)){ok=false;break;}}if(!ok)return;static unsigned char cv[100][51];memset(cv,0,sizeof(cv));for(int x:O)for(int y:R)cv[x][cd(x,y)]++;vector<int>ad;adddfs(R,cv,ad,0,K,cnt);}
 void remdfs(int pos,int left,vector<int>&v){if(found||timed)return;if(left==0){process(v);return;}for(int i=pos;i<=14-left;i++){v.push_back(i);remdfs(i+1,left-1,v);v.pop_back();if(found||timed)return;}}
 void run(){bool in[100]{};for(int x:Nv)in[x]=true;for(int x=0;x<100;x++)if(!in[x])O.push_back(x);t0=chrono::steady_clock::now();vector<int>v;if(first>=0){v.push_back(first);remdfs(first+1,K-1,v);}else remdfs(0,K,v);}
};
int main(int ac,char**av){if(ac<6){cerr<<"usage SET(1|2) K FIRST(-1|0..7) LIMIT OUT\n";return 2;}int setno=atoi(av[1]),K=atoi(av[2]),first=atoi(av[3]);double limit=atof(av[4]);string out=av[5];vector<int>N1={0,5,7,31,43,58,61,62,63,72,80,84,91,97};vector<int>N2={0,3,9,11,33,46,52,62,63,64,67,77,84,91};Search s{setno==1?N1:N2,{}, {},0,false,false,{},K,first,limit};s.run();double e=chrono::duration<double>(chrono::steady_clock::now()-s.t0).count();ofstream f(out);f<<"{\n\"near_miss\":\"N"<<setno<<"\",\n\"k\":"<<K<<",\n\"first_removed_index\":"<<first<<",\n\"completed_exhaustively\":"<<(s.timed?"false":"true")<<",\n\"timed_out\":"<<(s.timed?"true":"false")<<",\n\"witness_found\":"<<(s.found?"true":"false")<<",\n\"elapsed_seconds\":"<<e<<",\n\"removal_subsets_processed\":"<<s.remsets<<",\n\"addition_nodes\":[";for(int i=1;i<=K;i++){if(i>1)f<<",";f<<s.nodes[i];}f<<"],\n\"witness\":[";for(size_t i=0;i<s.sol.size();i++){if(i)f<<",";f<<s.sol[i];}f<<"]\n}\n";cout<<"N"<<setno<<" k="<<K<<" first="<<first<<" timed="<<s.timed<<" found="<<s.found<<" remsets="<<s.remsets<<" elapsed="<<e<<"\n";return 0;}
