#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
using namespace std;

static constexpr int MODS[7]={2,4,5,10,20,25,50};
static atomic<bool> GLOBAL_FOUND(false),GLOBAL_TIMED(false);
static chrono::steady_clock::time_point GLOBAL_START;
static double GLOBAL_LIMIT=3200.0;

struct BranchResult{
 int gap=0,third=0,fourth=0; bool completed=false,timed=false,found=false,verified=false;
 uint64_t nodes=0,prunes=0,qprunes=0,lbprunes=0; double wall=0.0; array<int,14>w{};
};

struct Search{
 int gap,third,fourth; array<int,14>a{}; array<unsigned char,51>dc{};
 array<array<unsigned char,50>,7>qd{}; array<array<unsigned char,50>,7>rc{};
 uint64_t nodes=0,prunes=0,qprunes=0,lbprunes=0; bool found=false,timed=false; array<int,14>w{};
 int cap(int mi,int d) const {int m=MODS[mi];return d==0?2*(100/m-1):2*(100/m);}
 bool timeout(){
   if((nodes&((1u<<20)-1))!=0)return false;
   if(GLOBAL_FOUND.load(memory_order_relaxed))return true;
   if(chrono::duration<double>(chrono::steady_clock::now()-GLOBAL_START).count()>GLOBAL_LIMIT){timed=true;GLOBAL_TIMED.store(true);return true;}
   return false;
 }
 bool lower(int depth){
   int rem=14-depth;
   for(int mi=0;mi<7;mi++){
     int m=MODS[mi]; vector<int>c(m);
     for(int r=0;r<m;r++)c[r]=rc[mi][r];
     for(int k=0;k<rem;k++){auto it=min_element(c.begin(),c.end());(*it)++;}
     int s=0;for(int n:c)s+=n*(n-1);
     if(s>cap(mi,0)){lbprunes++;return false;}
   }
   return true;
 }
 bool add(int x,int depth,vector<int>&dt,vector<array<int,3>>&qt){
   for(int i=0;i<depth;i++){
     int raw=x-a[i],d=min(raw,100-raw),c=(d==50?1:2);
     if(++dc[d]>c){--dc[d];for(int z:dt)--dc[z];for(auto z:qt){--qd[z[0]][z[1]];--qd[z[0]][z[2]];}return false;}
     dt.push_back(d);
     for(int mi=0;mi<7;mi++){
       int m=MODS[mi],d1=raw%m,d2=(m-d1)%m;
       ++qd[mi][d1];++qd[mi][d2];qt.push_back({mi,d1,d2});
       if(qd[mi][d1]>cap(mi,d1)||qd[mi][d2]>cap(mi,d2)){
         for(int z:dt)--dc[z];for(auto z:qt){--qd[z[0]][z[1]];--qd[z[0]][z[2]];}qprunes++;return false;
       }
     }
   }
   for(int mi=0;mi<7;mi++)rc[mi][x%MODS[mi]]++;
   if(!lower(depth+1)){
     for(int mi=0;mi<7;mi++)rc[mi][x%MODS[mi]]--;
     for(int z:dt)--dc[z];for(auto z:qt){--qd[z[0]][z[1]];--qd[z[0]][z[2]];}return false;
   }
   return true;
 }
 void undo(int x,const vector<int>&dt,const vector<array<int,3>>&qt){
   for(int mi=0;mi<7;mi++)rc[mi][x%MODS[mi]]--;
   for(int z:dt)--dc[z];for(auto z:qt){--qd[z[0]][z[1]];--qd[z[0]][z[2]];}
 }
 void dfs(int depth,int last,int evens){
   if(found||timed||GLOBAL_FOUND.load(memory_order_relaxed))return;
   ++nodes;if(timeout())return;
   int rem=14-depth;
   if(evens>9||evens+rem<5){prunes++;return;}
   if(depth==14){
     int wrap=100-a[13],h=third-gap;
     if(wrap<gap||h>wrap){prunes++;return;}
     found=true;w=a;GLOBAL_FOUND.store(true);return;
   }
   int h=third-gap;
   int minx=last+gap;
   // After choosing x there remain rem-1 interior gaps of size >=gap and one wrap gap >=h.
   int maxx=100-h-gap*(rem-1);
   for(int x=minx;x<=maxx&&!found&&!timed&&!GLOBAL_FOUND.load(memory_order_relaxed);x++){
     vector<int>dt;vector<array<int,3>>qt;
     if(!add(x,depth,dt,qt)){prunes++;continue;}
     a[depth]=x;dfs(depth+1,x,evens+(x%2==0));undo(x,dt,qt);
   }
 }
 static bool verify(const array<int,14>&a,int g,int t){
   if(a[0]!=0||a[1]!=g||a[2]!=t)return false;
   for(int i=1;i<14;i++)if(a[i]-a[i-1]<g)return false;
   if(100-a[13]<g||t-g>100-a[13])return false;
   array<int,100>r{};
   for(int i=0;i<14;i++)for(int j=0;j<14;j++)if(i!=j)r[(a[i]-a[j]+100)%100]++;
   for(int d=1;d<100;d++)if(r[d]>2)return false;
   return true;
 }
 BranchResult run(){
   auto st=chrono::steady_clock::now();
   a[0]=0;for(int mi=0;mi<7;mi++)rc[mi][0]++;
   vector<int>d1;vector<array<int,3>>q1;
   if(add(gap,1,d1,q1)){
     a[1]=gap;vector<int>d2;vector<array<int,3>>q2;
     if(add(third,2,d2,q2)){
       a[2]=third;vector<int>d3;vector<array<int,3>>q3;
       if(add(fourth,3,d3,q3)){
         a[3]=fourth;dfs(4,fourth,1+(gap%2==0)+(third%2==0)+(fourth%2==0));undo(fourth,d3,q3);
       }
       undo(third,d2,q2);
     }
     undo(gap,d1,q1);
   }
   for(int mi=0;mi<7;mi++)rc[mi][0]--;
   BranchResult b;b.gap=gap;b.third=third;b.fourth=fourth;b.completed=!timed;b.timed=timed;b.found=found;b.verified=found&&verify(w,gap,third);
   b.nodes=nodes;b.prunes=prunes;b.qprunes=qprunes;b.lbprunes=lbprunes;b.wall=chrono::duration<double>(chrono::steady_clock::now()-st).count();b.w=w;return b;
 }
};

static void write_json(const string&path,int gap,int third,const vector<BranchResult>&R,double wall){
 bool all=true,any=false,ok=false;uint64_t nodes=1,pr=0,qp=0,lb=0;array<int,14>w{};
 for(auto&r:R){all&=r.completed;any|=r.found;ok|=r.verified;nodes+=r.nodes;pr+=r.prunes;qp+=r.qprunes;lb+=r.lbprunes;if(r.verified)w=r.w;}
 int umin=third+gap,umax=100-(third-gap)-10*gap;
 ofstream f(path);f<<"{\n";
 f<<"  \"problem\": \"14-element B_2[2] subset of Z_100\",\n";
 f<<"  \"normalization\": \"minimum cyclic gap g mapped to {0,g}; third selected element fixed; reflection enforces third-g <= wrap\",\n";
 f<<"  \"gap\": "<<gap<<",\n  \"third\": "<<third<<",\n";
 f<<"  \"fourth_range\": ["<<umin<<", "<<umax<<"],\n";
 f<<"  \"expected_fourth_branches\": "<<max(0,umax-umin+1)<<",\n";
 f<<"  \"observed_fourth_branches\": "<<R.size()<<",\n";
 f<<"  \"completed_exhaustively\": "<<(all?"true":"false")<<",\n  \"timed_out\": "<<((!all)?"true":"false")<<",\n";
 f<<"  \"witness_found\": "<<(any?"true":"false")<<",\n  \"witness_verified\": "<<(ok?"true":"false")<<",\n";
 f<<"  \"aggregate_nodes_including_third_root\": "<<nodes<<",\n  \"aggregate_prunes\": "<<pr<<",\n";
 f<<"  \"aggregate_quotient_prunes\": "<<qp<<",\n  \"aggregate_lower_bound_prunes\": "<<lb<<",\n";
 f<<"  \"wall_seconds\": "<<setprecision(12)<<wall<<",\n";
 f<<"  \"coverage\": \"For normalized minimum gap g and fixed third t, the fourth selected element ranges from t+g through 100-(t-g)-10g. Every later consecutive gap is at least g and the final wrap is at least t-g by reflection. The disjoint fourth-element branches therefore cover the full normalized t case.\",\n";
 f<<"  \"witness\": [";if(ok){for(int i=0;i<14;i++){if(i)f<<",";f<<w[i];}}f<<"],\n";
 f<<"  \"branches\": [\n";
 for(size_t k=0;k<R.size();k++){auto&r=R[k];f<<"    {\"fourth\":"<<r.fourth<<",\"completed\":"<<(r.completed?"true":"false")<<",\"timed\":"<<(r.timed?"true":"false")<<",\"found\":"<<(r.found?"true":"false")<<",\"verified\":"<<(r.verified?"true":"false")<<",\"nodes\":"<<r.nodes<<",\"prunes\":"<<r.prunes<<",\"qprunes\":"<<r.qprunes<<",\"lbprunes\":"<<r.lbprunes<<",\"wall_seconds\":"<<setprecision(12)<<r.wall<<"}"<<(k+1<R.size()?",":"")<<"\n";}
 f<<"  ]\n}\n";
}

int main(int argc,char**argv){
 if(argc<7){cerr<<"usage: GAP THIRD LIMIT_SECONDS THREADS RESULT_JSON BRANCH_TSV\n";return 2;}
 int gap=atoi(argv[1]),third=atoi(argv[2]);GLOBAL_LIMIT=atof(argv[3]);int nth=max(1,atoi(argv[4]));string jout=argv[5],tout=argv[6];
 int tmin=2*gap,tmax=50-5*gap;
 if(gap<1||gap>7||third<tmin||third>tmax){cerr<<"invalid gap/third; expected third in ["<<tmin<<","<<tmax<<"]\n";return 2;}
 int umin=third+gap,umax=100-(third-gap)-10*gap;
 vector<int>us;for(int u=umin;u<=umax;u++)us.push_back(u);
 vector<BranchResult>R(us.size());atomic<size_t>next(0);GLOBAL_START=chrono::steady_clock::now();
 vector<thread>ts;
 for(int q=0;q<nth;q++)ts.emplace_back([&](){while(!GLOBAL_FOUND.load()&&!GLOBAL_TIMED.load()){size_t k=next.fetch_add(1);if(k>=us.size())break;Search s;s.gap=gap;s.third=third;s.fourth=us[k];R[k]=s.run();}});
 for(auto&t:ts)t.join();double wall=chrono::duration<double>(chrono::steady_clock::now()-GLOBAL_START).count();
 for(size_t k=0;k<R.size();k++)if(R[k].fourth==0){R[k].gap=gap;R[k].third=third;R[k].fourth=us[k];R[k].completed=false;R[k].timed=GLOBAL_TIMED.load();}
 write_json(jout,gap,third,R,wall);
 ofstream o(tout);o<<"gap\tthird\tfourth\tcompleted\ttimed\tfound\tverified\tnodes\tprunes\tqprunes\tlbprunes\twall_seconds\n";
 for(auto&r:R)o<<r.gap<<'\t'<<r.third<<'\t'<<r.fourth<<'\t'<<r.completed<<'\t'<<r.timed<<'\t'<<r.found<<'\t'<<r.verified<<'\t'<<r.nodes<<'\t'<<r.prunes<<'\t'<<r.qprunes<<'\t'<<r.lbprunes<<'\t'<<setprecision(12)<<r.wall<<'\n';
 bool all=true,ok=false;for(auto&r:R){all&=r.completed;ok|=r.verified;}
 cout<<"gap="<<gap<<" third="<<third<<" branches="<<R.size()<<" all="<<all<<" witness="<<ok<<" wall="<<wall<<"\n";
 return ok?10:(all?0:124);
}
