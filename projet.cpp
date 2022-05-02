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

void print_m(vector<vector<double>> A)    //affichage de matrice
{
    for( int i =0; i<A.size() ; i++){
        cout<<"ligne"<<i<<": " ;
        for ( int j =0; j<A.size() ; j++){
            cout<<A[i][j]<<" " ;
        }
    cout<<endl ; 
    }
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

Point point_exact(double a,double b,int N,int M,int s){
    vector<double> xi = Subdiv(a,N);
    vector<double> yi = Subdiv(b,M);
    Point R;
    Point P;
    R.invnumgb(N, s);
    P.x = xi[R.x];
    P.y = yi[R.y];
    return P;
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

vector<Point> sommetTR(double a,double b,int N,int M,const Triangle &T){
    vector<Point> St = vector<Point>(3);
    St[0] = point_exact(a,b,N,M,T.a);
    St[1] = point_exact(a,b,N,M,T.b);
    St[2] = point_exact(a,b,N,M,T.c);
    return St;
}


vector<vector<double>> CalcMatBT(double a,double b,int M,const Triangle& t)
{
    vector<vector<double>> BT(2,vector<double>(2));
    int N = t.c - t.a -2;
    vector<Point> St = sommetTR(a,b,N,M,t);
    Point p0 = St[0];
    Point p1 = St[1];
    Point p2 = St[2];
    BT[0][0] = p1.x-p0.x;
    BT[0][1] = p2.x-p0.x;
    BT[1][0] = p1.y-p0.y;
    BT[1][1] = p2.y-p0.y;
    return BT;
}



double DetBT(double a,double b,int M,const Triangle&t)    // Derminant de BT 
{
    vector<vector<double>> BT = CalcMatBT(a,b,M,t);
    return (BT[0][0]*BT[1][1]-BT[0][1]*BT[1][0]);
}


vector<vector<double>> GradGrad(double a,double b,int M,const Triangle &t,double (*eta)(double x,double y)){
    int N = t.c - t.a -2;
    vector<Point> St = sommetTR(a,b,N,M,t);
    double D = abs(DetBT(a,b,M,t));
    vector<vector<double>> GrdGrd(3,vector<double>(3));
    double val = (eta(St[0].x,St[0].y)+eta(St[1].x,St[1].y)+eta(St[2].x,St[2].y));
    GrdGrd[0][0] = 1./6/D*(pow(St[1].y-St[2].y,2)+pow(St[2].x-St[1].x,2))*val;
    GrdGrd[0][1] = 1./6/D*((St[1].y-St[2].y)*(St[2].y-St[0].y)+(St[2].x-St[1].x)*(St[0].x-St[2].x))*val;
    GrdGrd[0][2] = 1./6/D*((St[1].y-St[2].y)*(St[0].y-St[1].y)+(St[2].x-St[1].x)*(St[1].x-St[0].x))*val;
    GrdGrd[1][0] = 1*GrdGrd[0][1];
    GrdGrd[1][1] = 1./6/D*(pow(St[2].y-St[0].y,2)+pow(St[0].x-St[2].x,2))*val;
    GrdGrd[1][2] = 1./6/D*((St[2].y-St[0].y)*(St[0].y-St[1].y)+(St[0].x-St[2].x)*(St[1].x-St[0].x))*val;
    GrdGrd[2][0] = 1*GrdGrd[0][2];
    GrdGrd[2][1] = 1*GrdGrd[1][2];
    GrdGrd[2][2] = 1./6/D*(pow(St[0].y-St[1].y,2)+pow(St[1].x-St[0].x,2))*val;
    return GrdGrd;
}

double eta(double x,double y){
    return 2-sin(x+2*y);
}

int DansTrg(double a,double b,int N,int M,const Triangle &t,double x,double y)
{
    int info = 0;
    vector<Point> St = sommetTR(a,b,N,M,t);
    if ((x<=St[2].x)and(x>=St[0].x)and(y<=St[2].y)and(y>=St[0].y)){
        if (St[0].x==St[1].x){
            if((y-St[0].y)/(x-St[0].x) >= 1.){
                info = 1;
            }
        }
        else{
            if((y-St[0].y)/(x-St[0].x) <= 1.){
                info = 1;
            }
        }
    }
    return info;
}






int main(){
    int a = 1;
    int b = 1;
    int N = 5;
    int M = 5;
    vector<Triangle> TRG =  maillageTR(5,5);
    cout<<TRG[2]<<","<<TRG[3]<<endl;
    vector<vector<double>> BT = CalcMatBT(a,b,M,TRG[2]);
    Point Q;
    Q = point_exact(1,1,5,4,8);
    cout<<Q<<endl;
    vector<Point> St = sommetTR(a,b,N,M,TRG[2]);
    cout<<St[0]<<St[1]<<St[2]<<endl;
    BT = CalcMatBT(a,b,M,TRG[2]);
    print_m(BT);
    double Det = DetBT(a,b,M,TRG[2]);
    cout<<Det<<endl;
    vector<vector<double>>GrdGrd =  GradGrad(a,b,M,TRG[3],&eta);
    print_m(GrdGrd);
    int info = DansTrg(a,b,N,M,TRG[2],-0.5,-0.5);
    cout<<info<<endl;





    return 0;
}







