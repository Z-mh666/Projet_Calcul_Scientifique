#include <cstdlib>
#include <iostream>
#include <ostream>
#include <vector>
#include <list>
#include <fstream>
#include <cmath>
#include <time.h>
#include <ctime>
#include <math.h>
#include "projet.hpp"

using namespace std;


Point operator+(const Point & P, const Point &Q)
{
Point R(P ) ;
return R+=Q;
}
Point operator-(const Point & P, const Point &Q)
{
Point R(P ) ;
return R-=Q;
}
Point operator *( const Point & P, double a )
{
Point R(P ) ;
return R*=a ;
}
Point operator* (double a , const Point & P)
{
Point R(P ) ;
return R*=a ;
}
Point operator/ ( const Point & P, double a )
{
Point R(P ) ;
return R/=a ;
}

ostream & operator <<(ostream & os , const Point &P)
{
os<<"("<<P.x<<","<<P.y<<")" ;
return os ;
}


vector<double> Subdiv(double a,int N)
{
    vector<double> xi = vector<double>(N+1);
    for (int i=0;i<N+1;i++){
        xi[i] = -a +i*2*a/N;
    }
    return xi;
}


void print_v(vector<double> A)    //affichage de vecteur
{
    cout<<"[";
    for( int i =0; i<A.size() ; i++){
        cout<<A[i]<<"," ;
    }
    cout<<A[A.size()-1]<<"]";
}


int numgb(int N, const Point &P)
{
    return (N+1)*P.y+P.x;
}


Point Point::invnumgb(int N, int s){
    x = s%(N+1);// i: le reste de la division s par (N+1)
    y = s/(N+1);// j: le quotient de la division s par (N+1)
    return *this;  
}


int numint(int N, const Point &P)
{
    return (P.y-1)*(N-1)+P.x-1; 
}

Point Point::invnumint(int N, int k){
    x = k%(N+1)+1; // i-1 est le reste de la division k par (N+1)
    y = k/(N+1)+1; // j-1 est le quotient de la division k par (N+1)
    return *this;
}

int num_int_gb(int N, int k){
    Point P;
    P = P.invnumint(N,k);
    int s=numgb(N,P);
    return s;
}

int num_gb_int(int N, int s){
    Point P;
    P = P.invnumgb(N,s);
    int k= numint(N,P);
    return k;
}



ostream & operator <<(ostream &os, const Triangle &T){
    os<<"("<<T.a<<","<<T.b<<","<<T.c<<")" ;
    return os ;
}


vector<Triangle>  maillageTR(int N,int M)
{
    vector<Triangle> TRG = vector<Triangle>(2*N*M);
    Point P;
    int i=0;
    for(int s=0;s<(N+1)*(M)-1;s++){
        P.invnumgb(N,s);
        if (P.x == N){
            continue;
        }
        else{
            TRG[2*i].a = s;
            TRG[2*i].b = s+N+1;
            TRG[2*i].c = s+N+2;
            TRG[2*i+1].a = s;
            TRG[2*i+1].b = s+1;
            TRG[2*i+1].c = s+N+2;
            i++;
        }
    }
    return TRG;
}





int main(){

    Point P ;
    P.invnumgb(6,0);
    
    vector<Triangle> TRG =  maillageTR(4,2);
    for (int i=0;i<16;i++){
    cout<<TRG[i]<<endl;
    }












    return 0;
}







