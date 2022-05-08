#ifndef PROJET_H
#define PROJET_H

#include <iostream>
#include <ostream>
#include <vector>
#include <list>
#include <variant>

using namespace std;

class Point
{
private:
    
public:
    double x;
    double y;
    Point(double x0=0,double y0 =0):x(x0),y(y0){}
    Point invnumgb(int N, int s);
    Point invnumint(int N, int k);
};

vector<double> operator + ( vector<double> u , vector<double> v); // addition vectorielle
vector<double> operator - ( vector<double> u , vector<double> v); //soustraction vectorielle
vector<double> operator * (double a , vector<double> u); // multiplication reel vecteur
ostream & operator <<(ostream &, const Point &);


int numgb(int N, const Point &);
int numint(int N, const Point&);
int num_int_gb(int N, int k);
int num_gb_int(int N, int s);


vector<double> Subdiv(double a,int N);   //subdivision de l'intervalle
Point point_exact(double a,double b,int N,int M,int s);
void print_v(vector<double> A);   //affichage du vecteur
void print_m(vector<vector<double>> A);     //affichage de la matrice
double ps(vector<double> u,vector<double> v);     //produit scalaire
vector<double> pmv(vector<vector<double>> A, vector<double> v);   // produit matrice vecteur
double norme(vector<double> u);   // norme 2
double norm_inf(vector<double> u); //norme infinie



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
void print_t(vector<Triangle> A);   //affichage du vecteur de traingles
vector<Triangle>  maillageTR(int N,int M);
vector<Point> sommetTR(double a,double b,int N,int M,const Triangle &);
vector<vector<double>> CalcMatBT(double a,double b,int M,const Triangle&);
double DetBT(double a,double b,int M,const Triangle&);
vector<vector<double>> GradGrad(double a,double b,int M,const Triangle &t,double (*eta)(double x,double y));
int DansTrg(double a,double b,int N,int M,const Triangle &t,double x,double y);
vector<double> extendVec(vector<double> V,int N,int M);
vector<double> IntVec(vector<double> VV,int N,int M);
vector<double> matvec(vector<double> V, double a,double b,int N,int M);
vector<vector<double>> mat_elem(double a,double b,int N,int M, const Triangle &t);
vector<vector<double>> mat_assem(double a,double b,int N,int M);
vector<double> scdmembre(double a,double b,int N,int M,double (*rhsf)(double x,double y));
double normL2Grad(vector<double> V,double a,double b,int N,int M);
double normL2(vector<double> V,double a,double b,int N,int M);
vector<double> Grad_conj(vector<double> B, double a,double b,int N,int M);
double rhsf(double x,double y);
double u(double x,double y);
vector<double> solExa(double a,double b,int N,int M);
Triangle erreurs(double a,double b,int N,int M, double(*solExa)(double a,double b,int N,int M));

#endif
