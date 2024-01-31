#pragma once
#include "Matrix.h"
using namespace std;

template<class T>
class GaussElimination{
    private:
        int n;
    public:
        GaussElimination();
        void set_n(int);
        Matrix<T> solve(Matrix<T>, vector<T>);
};

template<class T>
class LU_Decomposition{
    private:
        int n;
    public:
        LU_Decomposition();
        void set_n(int);
        Matrix<T> solve(Matrix<T>, vector<T>);
};

template<class T>
class TriDiagonal{
    private:
        int n;
    public:
        TriDiagonal();
        void set_n(int);
        Matrix<T> solve(Matrix<T>, vector<T>);
};

template<class T>
class Gauss_Jacobi{
    private:
        int n;
    public:
        Gauss_Jacobi();
        void set_n(int);
        Matrix<T> solve(Matrix<T>, vector<T>);
};

template<class T>
class Gauss_Seidal{
    private:
        int n;
    public:
        Gauss_Seidal();
        void set_n(int);
        Matrix<T> solve(Matrix<T>, vector<T>);
};

template<class T>
class SOR{
    private:
        int n;
    public:
        SOR();
        void set_n(int);
        Matrix<T> solve(Matrix<T>, vector<T>);
};


//GaussElimination Method
template<class T>
GaussElimination<T> :: GaussElimination(){
    n = 1;
}

template<class T>
void GaussElimination<T> :: set_n(int num){
    n = num;
}

template<class T>
Matrix<T> GaussElimination<T>::solve(Matrix<T> M, vector<T> u)
{
    vector<T> x(u.size());
    T temp;
    int i,j,k;
    for(i=0;i<u.size();i++)
    {
        temp = fabs(M.getVal(i,i));
        k = i;
        for (j=i+1;j<u.size();j++)
            if(temp<fabs(M.getVal(j,i)))
            {
                temp = fabs(M.getVal(j,i));
                k = j;
            }
        if (fabs(M.getVal(k,i))<0.00001)
        {
            cout << "The matrix is singular: The system has either no solution or infinitely many solution";
            exit(0);
        }
        if(k!=i)
        {
            for(j=0;j<u.size();j++)
            {
                temp = M.getVal(k,j);
                M.setval(k,j,M.getVal(i,j));
                M.setval(i,j,temp);
            }
            temp = u[k];
            u[k] = u[i];
            u[i] = temp;
        }
        for(j=i+1;j<u.size();j++)
        {
            temp = M.getVal(j,i)/M.getVal(i,i);
            for(k=0;k<u.size();k++)
                M.setval(j,k, M.getVal(j,k) - temp*M.getVal(i,k));
            u[j] = u[j] - temp*u[i];
        }
    }
    x[u.size()-1] = u[u.size()-1] / M.getVal(u.size()-1, u.size()-1);
    for(i=u.size()-2;i>=0;i--)
    {
        x[i] = u[i];
        for(j=i+1;j<u.size();j++)
            x[i] = x[i] - M.getVal(i,j)*x[j];
        x[i] = x[i] / M.getVal(i,i);
    }
    
    Matrix<T> Y(u.size(),1);
    for(int i=0;i<u.size();i++){
        Y.setval(i,0,x[i]);
    }
    return Y;
}


//LU_Decomposition Method
template<class T>
LU_Decomposition<T> :: LU_Decomposition(){
    n = 1;
}

template<class T>
void LU_Decomposition<T> :: set_n(int num){
    n = num;
}

template<class T>
Matrix<T> LU_Decomposition<T>::solve(Matrix<T> M, vector<T> u){
    vector<T> x(u.size()), z(u.size());
    Matrix<T> l(u.size(),u.size()), u1(u.size(),u.size());
    int i,j,k;
    
    for(i=0;i<u.size();i++)
        l.setval(i,0, M.getVal(i,0));
    for(j=1;j<u.size();j++)
        u1.setval(0,j, M.getVal(0,j)/l.getVal(0,0));
    for(i=0;i<u.size();i++)
        u1.setval(i,i,1);
    for(i=1;i<u.size();i++)
        for(j=1;j<u.size();j++)
            if(i>=j)
            {
                l.setval(i,j,M.getVal(i,j));
                for(k=0;k<=j-1;k++){
                    
                    l.setval(i,j,l.getVal(i,j)-l.getVal(i,k)*u1.getVal(k,j));
                }
            }
            else
            {
                u1.setval(i,j,M.getVal(i,j));
                for(k=0;k<=i-1;k++)
                    u1.setval(i, j, u1.getVal(i,j)-l.getVal(i,k)*u1.getVal(k,j));
                u1.setval(i,j, u1.getVal(i,j)/l.getVal(i,i));
            }
    /*        
    cout<<"The lower triangular matrix L:"<<endl;
    for(i=0;i<u.size();i++)
    {
        for(j=0;j<=i;j++)
            cout<<"\t"<<l(i,j);
        cout<<endl;
    }
    cout<<"\nThe upper triangular matrix U:"<<endl;
    for(i=0;i<u.size();i++)
    {
        for(j=0;j<i;j++)
            cout<<"\t";
        for(j=i;j<u.size();j++)
            cout<<"\t"<<u1(i,j);
        cout<<endl;
    }
    */
    z[0]=u[0]/l.getVal(0,0);
    for(i=1;i<u.size();i++)
    {
        z[i]=u[i];
        for(j=0;j<=i-1;j++)
            z[i]-=l.getVal(i,j)*z[j];
        z[i]/=l.getVal(i,i);
    }
    x[u.size()-1]=z[u.size()-1];
    for(i=u.size()-2;i>=0;i--)
    {
        x[i]=z[i];
        for(j=i+1;j<u.size();j++)
            x[i]-=u1.getVal(i,j)*x[j];
    }
    
    Matrix<T> Y(u.size(),1);
    for(int i=0;i<u.size();i++){
        Y.setval(i,0,x[i]);
    }
    return Y;
}

//TriDiagonal Method
template<class T>
TriDiagonal<T> :: TriDiagonal(){
    n = 1;
}

template<class T>
void TriDiagonal<T> :: set_n(int num){
    n = num;
}

template<class T>
Matrix<T> TriDiagonal<T>::solve(Matrix<T> M, vector<T> u)
{
    vector<T> x(u.size());
    T A[u.size()-1],B[u.size()],C[u.size()-1],D[u.size()];
    T C_star[u.size()-1],D_star[u.size()];
    int i,j;
    
    for(i=0;i<u.size()-2;i++)
        for(j=i;j<u.size()-2;j++)
    {
        if(M.getVal(i,j+2)!=0)
        {
            cout<<"Method can't be applied";
            exit(0);
        }

    }
    for(i=u.size()-1;i>1;i--)
        for(j=i;j>1;j--)
    {
        if(M.getVal(i,j-2)!=0)
        {
            cout<<"Method can't be applied";
            exit(0);
        }
    }
    for(i=0;i<u.size();i++)
    {
        B[i]=M.getVal(i,i);
        D[i]=u[i];
    }
    for(i=1;i<u.size();i++)
    {
       A[i]=M.getVal(i,i-1);
    }
    for(i=0;i<u.size()-1;i++)
    {
       C[i]=M.getVal(i,i+1);
    }

    C_star[0]=C[0]/B[0];
    D_star[0]=D[0]/B[0];
    for(i=1;i<u.size();i++)
    {
        C_star[i]=C[i]/(B[i]-A[i]*C_star[i-1]);
        D_star[i]=(D[i]-A[i]*D_star[i-1])/(B[i]-A[i]*C_star[i-1]);
    }
    x[u.size()-1]=D_star[u.size()-1];
    for(i=u.size()-2;i>=0;i--)
    {
        x[i]=D_star[i]-C_star[i]*x[i+1];
    }
    Matrix<T> Y(u.size(),1);
    for(int i=0;i<u.size();i++){
        Y.setval(i,0,x[i]);
    }
    return Y;
}

// Gauss Jacobi Method
template<class T>
Gauss_Jacobi<T> :: Gauss_Jacobi(){
    n = 1;
}

template<class T>
void Gauss_Jacobi<T> :: set_n(int num){
    n = num;
}

template<class T>
Matrix<T> Gauss_Jacobi<T>::solve(Matrix<T> M, vector<T> u)
{
    vector<T> x,xn;
    int i,j,flag;
    T sum,eps;
    
    x=vector<T>(u.size());
    xn=vector<T>(u.size());
    //cout<<"Enter the accuracy u want\n";
    //cin>>eps;
    eps = 0.0001;
    for(i=0;i<u.size();i++)
        x[i]=0;
    flag=0;
    for(i=0;i<u.size();i++)
    {
        sum=0;
        for(j=0;j<u.size();j++)
            if(i!=j)
                sum+=fabs(M.getVal(i,j));
        if(sum>fabs(M.getVal(i,i)))
            flag=1;
    }
    if(flag==1)
    {
        flag=0;
        for(j=0;j<u.size();j++)
        {
            sum=0;
            for(i=0;i<u.size();i++)
                if(i!=j)
                    sum+=fabs(M.getVal(i,j));
            if(sum>fabs(M.getVal(j,j)))
                flag=1;
        }
    }
    if(flag==1)
    {
        cout<<"The co-efficient matrix is not diagonally dominant\n";
        cout<<"The Gauss-Jacobi method doesn't converge surely\n";
        exit(0);
    }
    do
    {
        for(i=0;i<u.size();i++)
        {
            sum=u[i];
            for(j=0;j<u.size();j++)
                if(j!=i)
                    sum-=M.getVal(i,j)*x[j];
            xn[i]=sum/M.getVal(i,i);
        }
        flag=0;
        for(i=0;i<u.size();i++)
            if(fabs(x[i]-xn[i])>eps)
                flag=1;
        if(flag==1)
            for(i=0;i<u.size();i++)
                x[i]=xn[i];
    }while(flag==1);
    
    Matrix<T> Y(u.size(),1);
    for(int i=0;i<u.size();i++){
        Y.setval(i,0,xn[i]);
    }
    return Y;
}

//Gauss Seidal method
template<class T>
Gauss_Seidal<T> :: Gauss_Seidal(){
    n = 1;
}

template<class T>
void Gauss_Seidal<T> :: set_n(int num){
    n = num;
}


template<class T>
Matrix<T> Gauss_Seidal<T>::solve(Matrix<T> M, vector<T> u)
{
    vector<T> x,xn;
    int test,i,j;
    T g,sum,eps=0.0001;
    /*
    for(i=0;i<u.size();i++){
       for(j=0;j<u.size();j++){
         cout<<M(i,j)<<" ";
       }
       cout<<"\n";
    }
    */  

    x=vector<T>(u.size());
    xn=vector<T>(u.size());
    for(i=0;i<u.size();i++)
    {
        x[i] = 0;
    }
    
    for(i=0;i<u.size();i++)
        x[i]=0;

        
    test=0;
    for(i=0;i<u.size();i++)
    {
        sum=0;
        for(j=0;j<u.size();j++)
            if(i!=j)
                sum+=fabs(M.getVal(i,j));
        if(sum>fabs(M.getVal(i,i)))
            test=1;
    }
 // Diagonally dominance verification by columns.
    if(test==1)
    {
        test=0;
        for(j=0;j<u.size();j++)
        {
            sum=0;
            for(i=0;i<u.size();i++)
                if(i!=j)
                    sum+=fabs(M.getVal(i,j));
            if(sum>fabs(M.getVal(j,j)))
                test=1;
        }
    }

    if(test==1)
    {
        cout<<"The co-efficient matrix is not diagonally dominant\n";
        cout<<"The Gauss-Seidel method doesn't converge surely\n";
        exit(0);
    }
    do
    {
        for(i=0;i<u.size();i++)
        {
            sum=u[i];
            for(j=0;j<u.size();j++)
            {
                if(j<i)
                    sum-=M.getVal(i,j)*xn[j];
                else if(j>i)
                    sum-=M.getVal(i,j)*x[j];
            }
            xn[i]=sum/M.getVal(i,i);
        }
        test=0;
        for(i=0;i<u.size();i++)
            if(fabs(x[i]-xn[i])>eps)
                test=1;
        if(test==1)
            for(i=0;i<u.size();i++)
                x[i] = xn[i];
    }while(test==1);

    Matrix<T> Y(u.size(),1);
    for(int i=0;i<u.size();i++){
        Y.setval(i,0,xn[i]);
    }
    return Y;
}

// Successive Over Relaxation method
template<class T>
SOR<T> :: SOR(){
    n = 1;
}

template<class T>
void SOR<T> :: set_n(int num){
    n = num;
}


template<class T>
Matrix<T> SOR<T>::solve(Matrix<T> M, vector<T> u){
    vector<T> x,xn;
    int i,j,flag;
    T sum,eps=0.001,w=1;
    
    x=vector<T>(u.size());
    xn=vector<T>(u.size());
    
    
    for(i=0;i<u.size();i++)
        x[i]=0;
    do
    {
        for(i=0;i<u.size();i++)
        {
            sum=u[i]*w+M.getVal(i,i)*x[i];
            for(j=0;j<u.size();j++)
            {
                if(j<i)
                    sum-=M.getVal(i,j)*xn[j]*w;
                else if(j>=i)
                    sum-=M.getVal(i,j)*x[j]*w;
                xn[i]=sum/M.getVal(i,i);
            }
        }
        flag=0;
        for(i=0;i<u.size();i++)
            if(fabs(x[i]-xn[i])>eps)
                flag=1;
        if(flag==1)
            for(i=0;i<u.size();i++)
                x[i]=xn[i];
    }while(flag==1);

    Matrix<T> Y(u.size(),1);
    for(int i=0;i<u.size();i++){
        Y.setval(i,0,xn[i]);
    }
    return Y;
}