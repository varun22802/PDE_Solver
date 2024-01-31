#pragma once
#include <bits/stdc++.h>
using namespace std;

template<class T>
class Matrix{
    private:
        T **M;
        int num_rows, num_cols;
    public:
        Matrix();
        Matrix(int, int);
        void set_nrc(int,int);
        void input();
        void display();
        int get_nrow();
        int get_ncol();
        Matrix<T> num_mat_sub(T);
        T getVal(int,int);
        void setval(int,int,T);
        Matrix<T> operator+(Matrix<T>);
        Matrix<T> operator-(Matrix<T>);
        Matrix<T> operator*(Matrix<T>);
        Matrix<T> operator+(T);
        Matrix<T> operator-(T);
        Matrix<T> operator*(T);
        
        template<class U>
        friend ostream& operator<<(ostream&, Matrix<T>&);
        
        template<class U>
        friend istream& operator>>(istream&, Matrix<T>&);
        // Define other operators as per requirement.
};

template<class T>
Matrix<T> :: Matrix(){
    num_rows = 0;
    num_cols = 0;
}

template<class T>
Matrix<T> :: Matrix(int m, int n)
{
    num_rows = m;
    num_cols = n;
    M = new T*[num_rows];
    for(int i=0;i<num_rows;i++)
        M[i] = new T[num_cols];
}

template<class T>
void Matrix<T> :: set_nrc(int r,int c){
    num_rows = r;
    num_cols = c;
}

template<class T>
int Matrix<T> :: get_nrow(){
    return num_rows;
}

template<class T>
int Matrix<T> :: get_ncol(){
    return num_cols;
}

template<class T>
void Matrix<T> :: input(){
    for(int j=0;j<num_rows;j++){
        for(int k=0;k<num_cols;k++){
            T n;
            cin>>n;
            M[j][k] = n;
        }
    }
}

template<class T>
void Matrix<T> :: display(){
    for(int j=0;j<num_rows;j++){
        for(int k=0;k<num_cols;k++){
            T n;
            n = getVal(j,k);
            cout<<n<<"\t";
        }
        cout<<endl;
    }
}

template<class T>
ostream& operator<<(ostream& out, Matrix<T>& Mat){
   Mat.display();
}

template<class T>
istream& operator>>(istream& in, Matrix<T>& Mat){
    Mat.input();
}

template<class T>
Matrix<T> Matrix<T> :: operator+(Matrix<T> A){
    Matrix<T> P;
    if(num_rows != A.num_rows || num_cols != A.num_cols){
        cout << "Invalid\n";
        exit(0);
    }
    else{
        P = Matrix(num_rows, num_cols);
        for(int i=0;i<num_rows;i++)
            for(int j=0;j<num_cols;j++)
                P.M[i][j] = M[i][j] + A.M[i][j];
        return P;
    }
}

template<class T>
Matrix<T> Matrix<T> :: operator-(Matrix<T> A){
    Matrix<T> P;
    if(num_rows != A.num_rows || num_cols != A.num_cols){
        cout << "Invalid\n";
        exit(0);
    }
    else{
        P = Matrix(num_rows, num_cols);
        for(int i=0;i<num_rows;i++)
            for(int j=0;j<num_cols;j++)
                P.M[i][j] = M[i][j] - A.M[i][j];
        return P;
    }
}

template<class T>
Matrix<T> Matrix<T> :: operator*(Matrix<T> A){
    Matrix<T> P;
    if (num_cols != A.num_rows){
        cout << "Invalid\n";
        exit(0);
    }
    else{
        P = Matrix(num_rows, A.num_cols);
        for(int i=0;i<num_rows;i++)
            for(int j=0;j<A.num_cols;j++){
                double zero = 0;
                P.M[i][j]=zero;
                for(int k=0;k<num_cols;k++)
                    P.M[i][j] = P.M[i][j] + M[i][k]*A.M[k][j];
            }
        return P;
    }
}

template<class T>
Matrix<T> Matrix<T> :: operator+(T num){
    Matrix<T> P;
    P = Matrix(num_rows, num_cols);
    for(int i=0 ; i<num_rows ; i++)
        for(int j=0 ; j<num_cols ; j++)
            P.M[i][j] = M[i][j] + num;
    return P;
    
}

template<class T>
Matrix<T> Matrix<T> :: operator-(T num){
    Matrix<T> P;
    P = Matrix(num_rows, num_cols);
    for(int i=0 ; i<num_rows ; i++)
        for(int j=0 ; j<num_cols ; j++)
            P.M[i][j] = M[i][j] - num;
    return P;
    
}

template<class T>
Matrix<T> Matrix<T> :: operator*(T num){
    Matrix<T> P;
    P = Matrix(num_rows, num_cols);
    for(int i=0 ; i<num_rows ; i++)
        for(int j=0 ; j<num_cols ; j++)
            P.M[i][j] = M[i][j]*num;
    return P;
}

template<class T>
Matrix<T> Matrix<T> :: num_mat_sub(T num){
    Matrix<T> P;
    P = Matrix(num_rows, num_cols);
    for(int i=0 ; i<num_rows ; i++)
        for(int j=0 ; j<num_cols ; j++)
            P.M[i][j] = num - M[i][j];
    return P;
}

template<class T>
T Matrix<T> :: getVal(int i,int j){
    return M[i][j];
}

template<class T>
void Matrix<T> :: setval(int i,int j,T num){
    M[i][j] = num;
}
