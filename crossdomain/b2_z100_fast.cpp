#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <mutex>
#include <numeric>
#include <random>
#include <thread>
#include <vector>
using namespace std;
constexpr int N=100,K=14;

struct Eval{int cost,maxr;};
Eval evalset(const array<int,K>& a){
    int cnt[N]{};
    for(int i=0;i<K;i++)for(int j=0;j<K;j++)if(i!=j){int d=(a[i]-a[j]+N)%N;if(d)cnt[d]++;}
    int c=0,m=0;for(int d=1;d<N;d++){m=max(m,cnt[d]);int e=max(0,cnt[d]-2);c+=e*e;}
    return {c,m};
}
bool verify(const array<int,K>&a){
    array<bool,N> seen{};for(int x:a){if(x<0||x>=N||seen[x])return false;seen[x]=true;}
    return evalset(a).cost==0;
}
int main(int argc,char**argv){
    int seconds=argc>1?stoi(argv[1]):1500;int nt=argc>2?stoi(argv[2]):8;
    atomic<bool> found(false);atomic<int> globalbest(999999);atomic<unsigned long long> moves(0);
    mutex mu;array<int,K> best{};auto st=chrono::steady_clock::now();
    auto worker=[&](int tid){
        mt19937_64 rng(0x20260817ULL+tid*0x9e3779b97f4a7c15ULL);uniform_real_distribution<double> U(0,1);
        while(!found && chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-st).count()<seconds){
            array<int,K>a{};a[0]=0;vector<int>pool(99);iota(pool.begin(),pool.end(),1);shuffle(pool.begin(),pool.end(),rng);for(int i=1;i<K;i++)a[i]=pool[i-1];sort(a.begin()+1,a.end());
            auto ev=evalset(a);double T=2.2;int stale=0;
            for(int step=0;step<500000 && !found;step++){
                if((step&4095)==0 && chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-st).count()>=seconds)break;
                bool in[N]{};for(int x:a)in[x]=true;
                int pos=1+(rng()%(K-1));int nv=1+(rng()%(N-1));if(in[nv])continue;
                int old=a[pos];a[pos]=nv;sort(a.begin()+1,a.end());auto ne=evalset(a);moves++;
                int dc=ne.cost-ev.cost;
                if(dc<=0 || U(rng)<exp(-dc/max(.02,T))){ev=ne;stale=dc<0?0:stale+1;}else{
                    // restore: locate nv and replace with old
                    auto it=find(a.begin()+1,a.end(),nv);*it=old;sort(a.begin()+1,a.end());stale++;
                }
                T=max(.025,T*.99996);
                if(ev.cost<globalbest.load()){
                    lock_guard<mutex>lk(mu);if(ev.cost<globalbest.load()){globalbest=ev.cost;best=a;cerr<<"best="<<ev.cost<<" maxr="<<ev.maxr<<" moves="<<moves.load()<<"\n";}
                }
                if(ev.cost==0 && verify(a)){lock_guard<mutex>lk(mu);best=a;found=true;break;}
                if(stale>50000){
                    // targeted restart/perturbation of three positions
                    for(int q=0;q<3;q++){bool used[N]{};for(int x:a)used[x]=true;int pos2=1+rng()%(K-1),v;do v=1+rng()%(N-1);while(used[v]);a[pos2]=v;sort(a.begin()+1,a.end());}
                    ev=evalset(a);T=1.2;stale=0;
                }
            }
        }
    };
    vector<thread>ts;for(int t=0;t<nt;t++)ts.emplace_back(worker,t);for(auto&th:ts)th.join();
    cout<<"status="<<(found?"WITNESS":"NO_WITNESS_WITHIN_BUDGET")<<"\n";
    cout<<"best_cost="<<globalbest.load()<<"\n";cout<<"moves="<<moves.load()<<"\n";cout<<"set=";for(int i=0;i<K;i++){if(i)cout<<",";cout<<best[i];}cout<<"\n";cout<<"verified="<<(verify(best)?"true":"false")<<"\n";
}
