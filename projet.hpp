#ifndef PROJET_H
#define PROJET_H

#include <iostream>
#include <ostream>
#include <vector>
#include <list>

using namespace std;

class Point
{
private:
    
public:
    double x;
    double y;
    Point(double x0=0,double y0 =0):x(x0),y(y0){}
    Point& operator+=( const Point & P) {x+=P.x ; y+=P.y ; return *this;}
    Point& operator-=( const Point & P) {x-=P.x ; y-=P.y ; return *this;}
    Point& operator *= (double a ){ x*=a ; y*=a ; return *this;}
    Point& operator/= (double a ){ x/=a ; y/=a ; return *this;}
    Point invnumgb(int N, int s);
    Point invnumint(int N, int k);
};
Point operator + ( const Point& , const Point &);
Point operator - ( const Point& , const Point &);
Point operator * ( const Point& ,double a);
Point operator * (double a , const Point &);
Point operator / ( const Point& ,double a);
//bool operator == ( const P oin t& , const P oin t& ) ;
//bool operator != ( const P oin t& , const P oin t& ) ;
ostream & operator <<(ostream &, const Point &);
int numgb(int N, const Point &);
int numint(int N, const Point&);
int num_int_gb(int N, int k);
int num_gb_int(int N, int s);


vector<double> Subdiv(double a,int N);
Point point_exact(double a,double b,int N,int M,int s);
void print_v(vector<double> A);
void print_m(vector<vector<double>> A);


class Triangle
{
    private:

    public:
    double a;
    double b;
    double c;
    Triangle(double a0 = 0,double b0 = 0,double c0 = 0):a(a0),b(b0),c(c0){} 
};
ostream & operator <<(ostream &, const Triangle &);
vector<Triangle>  maillageTR(int N,int M);
vector<Point> sommetTR(double a,double b,int N,int M,const Triangle &);
vector<vector<double>> CalcMatBT(double a,double b,int M,const Triangle&);
double DetBT(double a,double b,int M,const Triangle&);
vector<vector<double>> GradGrad(double a,double b,int M,const Triangle &t,double (*eta)(double x,double y));
int DansTrg(double a,double b,int N,int M,const Triangle &t,double x,double y);






#endif