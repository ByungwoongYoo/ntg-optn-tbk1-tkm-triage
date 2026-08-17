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
constexpr int N=32;

int violation(const array<int,N>& p){
    // displacement keys: dx 1..31, dy -(31)..31; for fixed dx, duplicates are forbidden.
    int c[N][2*N-1]{};
    int bad=0;
    for(int i=0;i<N;i++) for(int j=i+1;j<N;j++){
        int dx=j-i, dy=p[j]-p[i]+N-1;
        int &z=c[dx][dy]; if(z>=1) bad++; z++;
    }
    return bad;
}
bool verify(const array<int,N>&p){
    array<int,N> q=p; sort(q.begin(),q.end()); for(int i=0;i<N;i++) if(q[i]!=i)return false;
    return violation(p)==0;
}

int main(int argc,char**argv){
    int seconds=argc>1?stoi(argv[1]):1200;
    int threads=argc>2?stoi(argv[2]):max(1u,thread::hardware_concurrency());
    atomic<bool> found(false); atomic<long long> evals(0); atomic<int> globalBest(1000000);
    mutex mtx; array<int,N> bestPerm{}; iota(bestPerm.begin(),bestPerm.end(),0);
    auto start=chrono::steady_clock::now();
    auto worker=[&](int tid){
        mt19937_64 rng(20260817ULL + 0x9e3779b97f4a7c15ULL*(tid+1));
        uniform_real_distribution<double> U(0,1);
        array<int,N> p;
        while(!found){
            if(chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-start).count()>=seconds) break;
            iota(p.begin(),p.end(),0); shuffle(p.begin(),p.end(),rng);
            int c=violation(p); double T=2.0;
            int stagnant=0;
            for(long long step=0; step<300000 && !found; ++step){
                if((step&4095)==0 && chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-start).count()>=seconds) break;
                int i=rng()%N,j=rng()%N; if(i==j)continue;
                swap(p[i],p[j]); int c2=violation(p); evals++;
                int delta=c2-c;
                if(delta<=0 || U(rng)<exp(-delta/max(0.02,T))){c=c2; stagnant=(delta<0?0:stagnant+1);} else {swap(p[i],p[j]); stagnant++;}
                T=max(0.03,T*0.99996);
                if(c<globalBest.load()){
                    lock_guard<mutex> lk(mtx);
                    if(c<globalBest.load()){globalBest=c;bestPerm=p; cerr<<"best "<<c<<" thread "<<tid<<" evals "<<evals.load()<<"\n";}
                }
                if(c==0 && verify(p)){
                    lock_guard<mutex> lk(mtx); bestPerm=p; globalBest=0; found=true; break;
                }
                if(stagnant>30000){ // perturb without changing permutation validity
                    for(int k=0;k<5;k++){int a=rng()%N,b=rng()%N;swap(p[a],p[b]);}
                    c=violation(p); T=1.0; stagnant=0;
                }
            }
        }
    };
    vector<thread> ts; for(int t=0;t<threads;t++)ts.emplace_back(worker,t); for(auto&th:ts)th.join();
    cout<<"status="<<(found?"WITNESS":"NO_WITNESS_WITHIN_BUDGET")<<"\n";
    cout<<"best_violations="<<globalBest.load()<<"\n";
    cout<<"evaluations="<<evals.load()<<"\n";
    cout<<"permutation="; for(int i=0;i<N;i++){if(i)cout<<",";cout<<bestPerm[i];} cout<<"\n";
    cout<<"verified="<<(verify(bestPerm)?"true":"false")<<"\n";
    return 0;
}
