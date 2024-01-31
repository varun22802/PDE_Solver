#pragma once
#include <bits/stdc++.h>
using namespace std;

template <class T>
class NLsolver{
    private:
        int n,choice;
        Matrix<T> X;
        Matrix<T> M;
        Matrix<T> J;
        Discretizer<T> D1;
    public:
        NLsolver();
        NLsolver(Matrix<T>,Matrix<T>,Discretizer<T>,int,int);
        Matrix<T> solve();
};

template <class T>
NLsolver<T> :: NLsolver(){
    Matrix<T> N(1,1);
    N.setval(0,0,0);
    M = N;
}

template <class T>
NLsolver<T> :: NLsolver(Matrix<T> N, Matrix<T> O,Discretizer<T> P, int num, int choice_n){
    M = N;
    J = O;
    D1 = P;
    n = num;
    choice = choice_n;
    Matrix<T> X1(n,1);
    for(int i=0;i<n;i++){
        X1.setval(i,0,1);
    }
    X = X1;
}

template <class T>
Matrix<T> NLsolver<T> :: solve(){
    Matrix<T> F;
    F = M;
    vector<T> F1(n);
    for(int i=0;i<n;i++){
        F1[i] = F.getVal(i,0);
    }
    if(choice==1){
        GaussElimination<double> G;
        Matrix<T> Y;
        Y = G.solve(J*(-1),F1);
        double norm;
        norm = diff_norm(X,X+Y,n);
        X = X+Y;
        D1.discretize(X);
        M = D1.return_F();
        J = D1.return_J();
        while(norm>=0.000001){
            F = M;
            for(int i=0;i<n;i++){
                F1[i] = F.getVal(i,0);
            }
            Y = G.solve(J*(-1),F1);
            norm = diff_norm(X,X+Y,n);
            X=X+Y;
            D1.discretize(X);
            M = D1.return_F();
            J = D1.return_J();
        }
        return X;
    }    
    else if(choice==2){
        LU_Decomposition<double> G;
        Matrix<T> Y;
        Y = G.solve(J*(-1),F1);
        double norm;
        norm = diff_norm(X,X+Y,n);
        X = X+Y;
        D1.discretize(X);
        M = D1.return_F();
        J = D1.return_J();
        while(norm>=0.000001){
            F = M;
            for(int i=0;i<n;i++){
                F1[i] = F.getVal(i,0);
            }
            Y = G.solve(J*(-1),F1);
            norm = diff_norm(X,X+Y,n);
            X=X+Y;
            D1.discretize(X);
            M = D1.return_F();
            J = D1.return_J();
        }
        return X;
    }
    else if(choice==3){
        TriDiagonal<double> G;
        Matrix<T> Y;
        Y = G.solve(J*(-1),F1);
        double norm;
        norm = diff_norm(X,X+Y,n);
        X = X+Y;
        D1.discretize(X);
        M = D1.return_F();
        J = D1.return_J();
        while(norm>=0.000001){
            F = M;
            for(int i=0;i<n;i++){
                F1[i] = F.getVal(i,0);
            }
            Y = G.solve(J*(-1),F1);
            norm = diff_norm(X,X+Y,n);
            X=X+Y;
            D1.discretize(X);
            M = D1.return_F();
            J = D1.return_J();
        }
        return X;
    }    
    else if(choice==4){
        Gauss_Jacobi<double> G;
        Matrix<T> Y;
        Y = G.solve(J*(-1),F1);
        double norm;
        norm = diff_norm(X,X+Y,n);
        X = X+Y;
        D1.discretize(X);
        M = D1.return_F();
        J = D1.return_J();
        while(norm>=0.000001){
            F = M;
            for(int i=0;i<n;i++){
                F1[i] = F.getVal(i,0);
            }
            Y = G.solve(J*(-1),F1);
            norm = diff_norm(X,X+Y,n);
            X=X+Y;
            D1.discretize(X);
            M = D1.return_F();
            J = D1.return_J();
        }
        return X;
    }    
    else if(choice==5){
        Gauss_Seidal<double> G;
        Matrix<T> Y;
        Y = G.solve(J*(-1),F1);
        double norm;
        norm = diff_norm(X,X+Y,n);
        X = X+Y;
        D1.discretize(X);
        M = D1.return_F();
        J = D1.return_J();
        while(norm>=0.000001){
            F = M;
            for(int i=0;i<n;i++){
                F1[i] = F.getVal(i,0);
            }
            Y = G.solve(J*(-1),F1);
            norm = diff_norm(X,X+Y,n);
            X=X+Y;
            D1.discretize(X);
            M = D1.return_F();
            J = D1.return_J();
        }
        return X;
    }    
    else if(choice==6){
        SOR<double> G;
        Matrix<T> Y;
        Y = G.solve(J*(-1),F1);
        double norm;
        norm = diff_norm(X,X+Y,n);
        X = X+Y;
        D1.discretize(X);
        M = D1.return_F();
        J = D1.return_J();
        while(norm>=0.000001){
            F = M;
            for(int i=0;i<n;i++){
                F1[i] = F.getVal(i,0);
            }
            Y = G.solve(J*(-1),F1);
            norm = diff_norm(X,X+Y,n);
            X=X+Y;
            D1.discretize(X);
            M = D1.return_F();
            J = D1.return_J();
        }
        return X;
    }    
    else{
        cout<<"Did not chose a valid method"<<endl;
        exit(0);
    }
}