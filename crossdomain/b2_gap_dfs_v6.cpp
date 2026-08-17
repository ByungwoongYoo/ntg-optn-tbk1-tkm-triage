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
struct Search {
    int g; double limit_sec; chrono::steady_clock::time_point start;
    array<int,14> a{}; array<unsigned char,51> distcnt{};
    array<array<unsigned char,50>,7> qdiff{}; // oriented difference counts mod m
    array<array<unsigned char,50>,7> rcount{};
    uint64_t nodes=0, prunes=0, qprunes=0, lbprunes=0; bool found=false,timed=false; array<int,14>witness{};

    bool timeout(){ if((nodes&((1u<<20)-1))!=0)return false; if(chrono::duration<double>(chrono::steady_clock::now()-start).count()>limit_sec){timed=true;return true;}return false; }
    int cap(int mi,int d) const {int m=MODS[mi]; return d==0 ? 2*(100/m-1) : 2*(100/m);}

    bool min_zero_class_feasible(int depth){
        int rem=14-depth;
        for(int mi=0;mi<7;mi++){
            int m=MODS[mi]; vector<int> c(m); for(int r=0;r<m;r++)c[r]=rcount[mi][r];
            // optimistic minimum final same-residue ordered pairs by placing remaining points in currently least-loaded classes.
            for(int k=0;k<rem;k++){auto it=min_element(c.begin(),c.end());(*it)++;}
            int s=0;for(int n:c)s+=n*(n-1);
            if(s>cap(mi,0)){lbprunes++;return false;}
        }
        return true;
    }
    bool add_point(int x,int depth,vector<int>&dt,vector<array<int,3>>&qt){
        for(int i=0;i<depth;i++){
            int raw=x-a[i]; int cd=raw>50?100-raw:raw; int dcap=cd==50?1:2;
            if(++distcnt[cd]>dcap){--distcnt[cd]; for(int d:dt)--distcnt[d]; for(auto z:qt){--qdiff[z[0]][z[1]];--qdiff[z[0]][z[2]];} return false;}
            dt.push_back(cd);
            for(int mi=0;mi<7;mi++){
                int m=MODS[mi], d1=raw%m, d2=(m-d1)%m;
                ++qdiff[mi][d1]; ++qdiff[mi][d2]; qt.push_back({mi,d1,d2});
                if(qdiff[mi][d1]>cap(mi,d1) || qdiff[mi][d2]>cap(mi,d2)){
                    for(int d:dt)--distcnt[d]; for(auto z:qt){--qdiff[z[0]][z[1]];--qdiff[z[0]][z[2]];} qprunes++; return false;
                }
            }
        }
        for(int mi=0;mi<7;mi++)rcount[mi][x%MODS[mi]]++;
        if(!min_zero_class_feasible(depth+1)){
            for(int mi=0;mi<7;mi++)rcount[mi][x%MODS[mi]]--;
            for(int d:dt)--distcnt[d]; for(auto z:qt){--qdiff[z[0]][z[1]];--qdiff[z[0]][z[2]];} return false;
        }
        return true;
    }
    void undo_point(int x,const vector<int>&dt,const vector<array<int,3>>&qt){
        for(int mi=0;mi<7;mi++)rcount[mi][x%MODS[mi]]--;
        for(int d:dt)--distcnt[d]; for(auto z:qt){--qdiff[z[0]][z[1]];--qdiff[z[0]][z[2]];}
    }
    void dfs(int depth,int last,int evens){
        if(found||timed)return; ++nodes;if(timeout())return;
        int rem=14-depth; if(evens>9||evens+rem<5){prunes++;return;}
        if(depth==14){int wrap=100-a[13];if(wrap<g){prunes++;return;}int second=a[2]-a[1];if(second>wrap){prunes++;return;}found=true;witness=a;return;}
        int maxx=100-g*rem, minx=last+g;
        for(int x=minx;x<=maxx&&!found&&!timed;x++){
            vector<int>dt;dt.reserve(depth);vector<array<int,3>>qt;qt.reserve(depth*7);
            if(!add_point(x,depth,dt,qt)){prunes++;continue;}a[depth]=x;dfs(depth+1,x,evens+(x%2==0));undo_point(x,dt,qt);
        }
    }
    void run(){start=chrono::steady_clock::now();a[0]=0;for(int mi=0;mi<7;mi++)rcount[mi][0]++;vector<int>dt;vector<array<int,3>>qt;if(!add_point(g,1,dt,qt))return;a[1]=g;dfs(2,g,1+(g%2==0));undo_point(g,dt,qt);for(int mi=0;mi<7;mi++)rcount[mi][0]--;}
};
static bool verify(const array<int,14>&a){array<int,100>r{};for(int i=0;i<14;i++)for(int j=0;j<14;j++)if(i!=j)r[(a[i]-a[j]+100)%100]++;for(int d=1;d<100;d++)if(r[d]>2)return false;return true;}
int main(int argc,char**argv){if(argc<4){cerr<<"usage GAP LIMIT OUT\n";return 2;}int g=atoi(argv[1]);double lim=atof(argv[2]);string out=argv[3];Search s{g,lim};s.run();double wall=chrono::duration<double>(chrono::steady_clock::now()-s.start).count();bool ok=s.found&&verify(s.witness);ofstream f(out);f<<"{\n\"gap\":"<<g<<",\n\"algorithm\":\"exact DFS + quotient-difference capacity pruning mod 2,4,5,10,20,25,50\",\n\"completed_exhaustively\":"<<(s.timed?"false":"true")<<",\n\"timed_out\":"<<(s.timed?"true":"false")<<",\n\"witness_found\":"<<(s.found?"true":"false")<<",\n\"witness_verified\":"<<(ok?"true":"false")<<",\n\"nodes\":"<<s.nodes<<",\n\"prunes\":"<<s.prunes<<",\n\"quotient_prunes\":"<<s.qprunes<<",\n\"lower_bound_prunes\":"<<s.lbprunes<<",\n\"wall_seconds\":"<<wall<<",\n\"witness\":[";if(s.found){for(int i=0;i<14;i++){if(i)f<<",";f<<s.witness[i];}}f<<"]\n}\n";cout<<"gap="<<g<<" nodes="<<s.nodes<<" prunes="<<s.prunes<<" q="<<s.qprunes<<" lb="<<s.lbprunes<<" timed="<<s.timed<<" found="<<s.found<<" ok="<<ok<<" wall="<<wall<<"\n";return 0;}
