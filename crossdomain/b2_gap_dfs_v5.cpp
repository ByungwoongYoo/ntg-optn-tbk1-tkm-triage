#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

struct Search {
    int g;
    double limit_sec;
    chrono::steady_clock::time_point start;
    array<int,14> a{};
    array<unsigned char,51> cnt{};
    uint64_t nodes=0, prunes=0;
    bool found=false, timed=false;
    array<int,14> witness{};

    bool timeout(){
        if ((nodes & ((1u<<20)-1))!=0) return false;
        double s=chrono::duration<double>(chrono::steady_clock::now()-start).count();
        if(s>limit_sec){timed=true; return true;} return false;
    }
    bool adddist(int x,int depth, vector<int>& touched){
        for(int i=0;i<depth;i++){
            int d=x-a[i]; if(d>50)d=100-d;
            int cap=(d==50?1:2);
            if(++cnt[d]>cap){
                for(int t:touched) --cnt[t];
                --cnt[d];
                return false;
            }
            touched.push_back(d);
        }
        return true;
    }
    void undo(const vector<int>& touched){ for(int d:touched)--cnt[d]; }
    void dfs(int depth,int last,int evens){
        if(found||timed) return;
        ++nodes; if(timeout())return;
        int remaining=14-depth;
        if(evens>9 || evens+remaining<5){++prunes;return;}
        if(depth==14){
            int wrap=100-a[13]; if(wrap<g){++prunes;return;}
            // Reflection breaker around fixed first minimum gap: second cyclic gap <= last cyclic gap.
            int secondgap=a[2]-a[1]; if(secondgap>wrap){++prunes;return;}
            found=true; witness=a; return;
        }
        // Need room for the remaining points and a final wrap gap >=g.
        int maxx=100-g*(remaining); // after choosing x, remaining-1 future gaps + wrap
        int minx=last+g;
        for(int x=minx;x<=maxx && !found && !timed;x++){
            vector<int> touched; touched.reserve(depth);
            if(!adddist(x,depth,touched)){++prunes;continue;}
            a[depth]=x;
            dfs(depth+1,x,evens+(x%2==0));
            undo(touched);
        }
    }
    void run(){
        start=chrono::steady_clock::now();
        a[0]=0; a[1]=g;
        vector<int> t;
        if(!adddist(g,1,t)){cerr<<"bad seed\n";return;}
        dfs(2,g,1+(g%2==0));
        undo(t);
    }
};

static bool verify(const array<int,14>& a){
    array<int,100> r{};
    for(int i=0;i<14;i++) for(int j=0;j<14;j++) if(i!=j){int d=(a[i]-a[j]+100)%100; r[d]++;}
    for(int d=1;d<100;d++) if(r[d]>2)return false;
    return true;
}

int main(int argc,char**argv){
    if(argc<4){cerr<<"usage: prog GAP LIMIT_SECONDS OUT_JSON\n";return 2;}
    int g=atoi(argv[1]); double lim=atof(argv[2]); string out=argv[3];
    Search s{g,lim}; s.run();
    double wall=chrono::duration<double>(chrono::steady_clock::now()-s.start).count();
    bool ok=s.found && verify(s.witness);
    ofstream f(out);
    f<<"{\n  \"gap\": "<<g<<",\n  \"algorithm\": \"C++ incremental exact DFS on sorted cyclic ruler with distance multiplicity pruning\",\n";
    f<<"  \"completed_exhaustively\": "<<(s.timed?"false":"true")<<",\n  \"timed_out\": "<<(s.timed?"true":"false")<<",\n";
    f<<"  \"witness_found\": "<<(s.found?"true":"false")<<",\n  \"witness_verified\": "<<(ok?"true":"false")<<",\n";
    f<<"  \"nodes\": "<<s.nodes<<",\n  \"prunes\": "<<s.prunes<<",\n  \"wall_seconds\": "<<wall<<",\n  \"witness\": [";
    if(s.found){for(int i=0;i<14;i++){if(i)f<<",";f<<s.witness[i];}}
    f<<"]\n}\n"; f.close();
    cout<<"gap="<<g<<" nodes="<<s.nodes<<" prunes="<<s.prunes<<" timed="<<s.timed<<" found="<<s.found<<" verified="<<ok<<" wall="<<wall<<"\n";
    if(s.found){for(int x:s.witness)cout<<x<<" ";cout<<"\n";}
    return 0;
}
