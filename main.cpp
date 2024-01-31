#include<bits/stdc++.h>
#include "AD.h"
#include "Discretizer.h"
#include "Matrix.h"
#include "LinearSolver.h"
#include "NonLInearSolver.h"
using namespace std;

template<class T>
T diff_norm(Matrix<T>,Matrix<T>,int);

int main(){
  
// ******************** Start driver code from here *******************
    string A,B,C,D,E;
    string f1,f2,g1,g2;
    int nx,ny;
    double ax,ay;
    double bx,by;
    
    cout<<"Enter the boundary points a,b in the order ax,ay,bx,by separated by spaces:"<<endl;
    cin>>ax;
    cin>>ay;
    cin>>bx;
    cin>>by;
    
    cout<<"Enter nx and ny, the number of partitions you want to divide x and y in:"<<endl;
    cin>>nx;
    cin>>ny;
    
    cout<<"Enter A,B,C,D,E in the form of equations separated by spaces:"<<endl;
    cin>>A;
    if(A[0]=='-')A = "0"+A;
    cin>>B;
    if(B[0]=='-')B = "0"+B;
    cin>>C;
    if(C[0]=='-')C = "0"+C;
    cin>>D;
    if(D[0]=='-')D = "0"+D;
    cin>>E;
    if(E[0]=='-')E = "0"+E;
    
    cout<<"Enter f1,f2,g1,g2 the border functions (f1 is function at x=ax, f2 at x=bx, g1 at y = ay and g2 at y=by) separated by spaces:"<<endl;
    cin>>f1;
    cin>>f2;
    cin>>g1;
    cin>>g2;
    
    int choice;
    cout<<"Enter 1 if you want to use Gauss Elimination, 2 for LU decomp, 3 for TriDiagonal, 4 for Gauss Jacobi, 5 for Gauss Seidal, 6 for SOR:"<<endl;
    cin>>choice;
    
    int n=(nx-1)*(ny-1);
    Discretizer<double> D1(A,B,C,D,E,ax,ay,bx,by);
    D1.set_boundary(f1,f2,g1,g2);
    D1.set_partitions(nx,ny);
    Matrix<double> M,J;
    Matrix<double> u(n,1);
    for(int i=0;i<n;i++){
        u.setval(i,0,1);
    }
    D1.discretize(u);
    M = D1.return_F();
    //cout<<M;
    J = D1.return_J();
    //cout<<J;
    NLsolver<double> Nl(M,J,D1,n,choice);
    Matrix<double> X = Nl.solve();
    vector<double> ui(n);
    cout<<"The final u in grid form with top left as (0,0) is:"<<endl;
    int count = 0;
    for(int i=0;i<nx-1;i++){
        for(int j=0;j<ny-1;j++){
            cout<<fixed<<setprecision(6);
            ui[count] = X.getVal(count,0);
            cout<<ui[count]<<"\t";
            count++;
        }
        cout<<endl;
    }
    
    
return 0;
}

template<class T>
T diff_norm(Matrix<T> X,Matrix<T> Y,int n){
    T norm = 0;
    for(int i=0;i<n;i++){
        norm += pow((X.getVal(i,0)-Y.getVal(i,0)),2);
    }
    norm = sqrt(norm);
    return norm;
}
