#include <array>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
static constexpr int MODS[7]={2,4,5,10,20,25,50};
struct Search{
 int a2; double limit; chrono::steady_clock::time_point start; array<int,14>a{}; array<unsigned char,51>dc{}; array<array<unsigned char,50>,7>qd{}; array<array<unsigned char,50>,7>rc{}; uint64_t nodes=0,prunes=0,qprunes=0,lbprunes=0; bool found=false,timed=false; array<int,14>w{};
 int cap(int mi,int d){int m=MODS[mi];return d==0?2*(100/m-1):2*(100/m);} bool timeout(){if((nodes&((1u<<20)-1))!=0)return false;if(chrono::duration<double>(chrono::steady_clock::now()-start).count()>limit){timed=true;return true;}return false;}
 bool lower(int depth){int rem=14-depth;for(int mi=0;mi<7;mi++){int m=MODS[mi];vector<int>c(m);for(int r=0;r<m;r++)c[r]=rc[mi][r];for(int k=0;k<rem;k++){auto it=min_element(c.begin(),c.end());(*it)++;}int s=0;for(int n:c)s+=n*(n-1);if(s>cap(mi,0)){lbprunes++;return false;}}return true;}
 bool add(int x,int depth,vector<int>&dt,vector<array<int,3>>&qt){for(int i=0;i<depth;i++){int raw=x-a[i],d=min(raw,100-raw),c=(d==50?1:2);if(++dc[d]>c){--dc[d];for(int z:dt)--dc[z];for(auto z:qt){--qd[z[0]][z[1]];--qd[z[0]][z[2]];}return false;}dt.push_back(d);for(int mi=0;mi<7;mi++){int m=MODS[mi],d1=raw%m,d2=(m-d1)%m;++qd[mi][d1];++qd[mi][d2];qt.push_back({mi,d1,d2});if(qd[mi][d1]>cap(mi,d1)||qd[mi][d2]>cap(mi,d2)){for(int z:dt)--dc[z];for(auto z:qt){--qd[z[0]][z[1]];--qd[z[0]][z[2]];}qprunes++;return false;}}}for(int mi=0;mi<7;mi++)rc[mi][x%MODS[mi]]++;if(!lower(depth+1)){for(int mi=0;mi<7;mi++)rc[mi][x%MODS[mi]]--;for(int z:dt)--dc[z];for(auto z:qt){--qd[z[0]][z[1]];--qd[z[0]][z[2]];}return false;}return true;}
 void undo(int x,const vector<int>&dt,const vector<array<int,3>>&qt){for(int mi=0;mi<7;mi++)rc[mi][x%MODS[mi]]--;for(int z:dt)--dc[z];for(auto z:qt){--qd[z[0]][z[1]];--qd[z[0]][z[2]];}}
 void dfs(int depth,int last,int evens){if(found||timed)return;++nodes;if(timeout())return;int rem=14-depth;if(evens>9||evens+rem<5){prunes++;return;}if(depth==14){int wrap=100-a[13];if(wrap<1||a[2]-1>wrap){prunes++;return;}found=true;w=a;return;}int minx=last+1;int maxx=100-rem; // leave rem-1 future positive gaps and wrap >=1
 // reflection constraint requires final wrap >= a2-1, so final a13 <= 101-a2
 maxx=min(maxx,101-a2-(rem-1));
 for(int x=minx;x<=maxx&&!found&&!timed;x++){vector<int>dt;vector<array<int,3>>qt;if(!add(x,depth,dt,qt)){prunes++;continue;}a[depth]=x;dfs(depth+1,x,evens+(x%2==0));undo(x,dt,qt);}}
 void run(){start=chrono::steady_clock::now();a[0]=0;for(int mi=0;mi<7;mi++)rc[mi][0]++;vector<int>d1;vector<array<int,3>>q1;if(!add(1,1,d1,q1))return;a[1]=1;vector<int>d2;vector<array<int,3>>q2;if(add(a2,2,d2,q2)){a[2]=a2;dfs(3,a2,1+0+(a2%2==0));undo(a2,d2,q2);}undo(1,d1,q1);for(int mi=0;mi<7;mi++)rc[mi][0]--;}
};
static bool verify(const array<int,14>&a){array<int,100>r{};for(int i=0;i<14;i++)for(int j=0;j<14;j++)if(i!=j)r[(a[i]-a[j]+100)%100]++;for(int d=1;d<100;d++)if(r[d]>2)return false;return true;}
int main(int argc,char**argv){if(argc<4){cerr<<"usage A2 LIMIT OUT\n";return 2;}int a2=atoi(argv[1]);double lim=atof(argv[2]);string out=argv[3];Search s{a2,lim};s.run();double wall=chrono::duration<double>(chrono::steady_clock::now()-s.start).count();bool ok=s.found&&verify(s.w);ofstream f(out);f<<"{\n\"normalized_unit_pair\":[0,1],\n\"fixed_third_element\":"<<a2<<",\n\"completed_exhaustively\":"<<(s.timed?"false":"true")<<",\n\"timed_out\":"<<(s.timed?"true":"false")<<",\n\"witness_found\":"<<(s.found?"true":"false")<<",\n\"witness_verified\":"<<(ok?"true":"false")<<",\n\"nodes\":"<<s.nodes<<",\n\"prunes\":"<<s.prunes<<",\n\"quotient_prunes\":"<<s.qprunes<<",\n\"lower_bound_prunes\":"<<s.lbprunes<<",\n\"wall_seconds\":"<<wall<<",\n\"witness\":[";if(s.found){for(int i=0;i<14;i++){if(i)f<<",";f<<s.w[i];}}f<<"]\n}\n";cout<<"a2="<<a2<<" nodes="<<s.nodes<<" timed="<<s.timed<<" found="<<s.found<<" ok="<<ok<<" wall="<<wall<<"\n";return 0;}
