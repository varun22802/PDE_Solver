#pragma once
#include "Matrix.h"
#include "AD.h"
#include <stack>
using namespace std;

bool IsOperator(char);  
bool IsOperand(char);  
bool eqlOrhigher(char, char);
string convert(string);
string infix_to_postfix(string);


template<class T>
class Discretizer{
    private:
        int nx,ny,n;
        double ax,ay,bx,by;
        string A,B,C,D,E;
        string f1,f2,g1,g2;
        Matrix<T> M,J;
        AD *F;
    public:
        Discretizer();
        Discretizer(string,string,string,string,string,double,double,double,double);
        void set_boundary(string,string,string,string);
        void set_partitions(int,int);
        void discretize(Matrix<T> u);
        Matrix<T> return_F();
        Matrix<T> return_J();
};

template<class T>
Discretizer<T> :: Discretizer(){
    nx=1;
    ny=1;
}

template<class T>
Discretizer<T> :: Discretizer(string n1,string n2,string n3,string n4,string n5, double n6, double n7, double n8, double n9){
    A = n1;
    B = n2;
    C = n3;
    D = n4;
    E = n5;
    ax = n6;
    ay = n7;
    bx = n8;
    by = n9;
}

template<class T>
void Discretizer<T> :: set_boundary(string a,string b,string c,string d){
    f1 = a;
    f2 = b;
    g1 = c;
    g2 = d;
}

template<class T>
void Discretizer<T> :: set_partitions(int x, int y){
    nx = x;
    ny = y;
    n = (nx-1)*(ny-1);
    F = new AD[n];
}

template<class T>
void Discretizer<T> :: discretize(Matrix<T> Mu1){
    double dx,dy,l,w;
    setNum_Var(n);
    double Mu[n];
    for(int i=0;i<n;i++){
        Mu[i] = Mu1.getVal(i,0);
    }
    l = abs(ax-bx);
    w = abs(ay-by);
    dx = l/nx;
    dy = w/ny;
    
    string temp[9];
    temp[0] = A;
    temp[1] = B;
    temp[2] = C;
    temp[3] = D;
    temp[4] = E;
    temp[5] = f1;
    temp[6] = f2;
    temp[7] = g1;
    temp[8] = g2;
    
    
    for(int ic=0;ic<n;ic++){
        AD *functions;
        functions = new AD[9];
        for(int jf=0;jf<9;jf++){
            //cout<<jf<<" count"<<ic<<endl;
            string infix_expression = convert(temp[jf]);
            int len = infix_expression.length();
            double inttrack[len];
            int j=0;
            double number = 0;
            string numb = "";
            int intcount = 0;
            while(j<len){
                char c = infix_expression[j];
                if((c>='0' && c<= '9') || c=='.'){
                    numb += c;
                    j++;
                }
                else{
                    if(numb != ""){
                        number = stod(numb);
                        inttrack[intcount] = number;
                        number = 0;
                        numb = "";
                        intcount++;
                    }    
                    j++;
                }
            }
            if(numb!=""){
                inttrack[intcount]=stod(numb);
                intcount++;
            }
            
            string postfix = infix_to_postfix(infix_expression);
            stack <AD> funcStack;
            int i = 0;
            int k = 0;
            double num = 0;
            numb = "";
            while(postfix[i] != '\0'){
                if(postfix[i] == 'u'){
                    AD y(Mu[ic],ic);
                    funcStack.push(y);
                }
                if(postfix[i] == 'x'){
                    double temp1 = (ic+1)%(nx-1);
                    if(temp1==0){
                        temp1 = ax + (nx-1)*dx;
                    }    
                    else{
                        temp1 = temp1*dx + ax; 
                    }
                    AD y(temp1,-1);
                    funcStack.push(y);
                }
                if(postfix[i] == 'y'){
                    int temp2 = (ic)/(nx-1);
                    double temp1 = (double)temp2;
                    temp1 = temp1 +1;
                    temp1 = ay + temp1*dy;
                    AD y(temp1,-1);
                    funcStack.push(y);
                }
                if((postfix[i] <= '9' && postfix[i] >='0') || postfix[i] == '.'){
                    numb += postfix[i];
                    if(stod(numb) == inttrack[k]){
                        num = stod(numb);
                        numb = "";
                        AD y(num,-1);
                        funcStack.push(y);
                        num = 0;
                        k++;
                    }
                }
                
                AD func;
                if(postfix[i]=='+'||postfix[i]=='-'||postfix[i]=='*'||postfix[i]=='/'||postfix[i]=='^'){
                    AD func1 = funcStack.top();
                    funcStack.pop();
                    AD func2 = funcStack.top();
                    funcStack.pop();
                    switch(postfix[i])
                    {
                        case '+' : funcStack.push(func2+func1);
                                    break;
                        case '-' : funcStack.push(func2-func1);
                                    break;
                        case '*' : funcStack.push(func2*func1);
                                    break;
                        case '/' : funcStack.push(func2/func1);
                                    break;
                        case '^' : funcStack.push(func2^func1);
        
                    }
                }
                
                if(postfix[i] >= 'A' && postfix[i] <= 'O'){
                    AD func = funcStack.top();
                    funcStack.pop();
                    switch(postfix[i])
                    {
                        case 'A' : funcStack.push(sin(func));
                                    break;
                        case 'B' : funcStack.push(cos(func));
                                    break;
                        case 'C' : funcStack.push(tan(func));
                                    break;
                        case 'D' : funcStack.push(cosec(func));
                                    break;
                        case 'E' : funcStack.push(sec(func));
                                    break;
                        case 'F' : funcStack.push(cot(func));
                                    break;
                        case 'G' : funcStack.push(sinh(func));
                                    break;
                        case 'H' : funcStack.push(cosh(func));
                                    break;
                        case 'I' : funcStack.push(tanh(func));
                                    break;
                        case 'J' : funcStack.push(arcsin(func));
                                    break;
                        case 'K' : funcStack.push(arccos(func));
                                    break;
                        case 'L' : funcStack.push(arctan(func));
                                    break;
                        case 'M' : funcStack.push(abs(func));
                                    break;
                        case 'N' : funcStack.push(log(func));
                                    break;
                        case 'O' : funcStack.push(exp(func));
                                    break;
                    }
                }
                i++;
            }   
            AD f = funcStack.top();
            functions[jf] = f;
            //cout<<f.get_f()<<endl;
        }
        AD Ad,Bd,Cd,Dd,Ed,f1d,f2d,g1d,g2d;
        Ad = functions[0]/(dx*dx);
        Bd = functions[1]/(dy*dy);
        Cd = functions[2]/(dx);
        Dd = functions[3]/(dy);
        Ed = functions[4];
        f1d = functions[5];
        f2d = functions[6];
        g1d = functions[7];
        g2d = functions[8];
        
        
        AD ui(Mu[ic],ic);
        double cond;
        cond = (ic+1);
        cond /= (nx-1);
        if(ic==0){
            AD ui1(Mu[ic+1],ic+1);
            AD ui2(Mu[ic+nx-1],ic+nx-1);
            F[ic] = ((0-2*Ad+0-2*Bd+0-Cd+0-Dd)*ui)+((Ad+Cd)*ui1)+((Bd+Dd)*ui2)+(Ad*f1d+Bd*g1d+Ed);
        }
        else if(ic==nx-2){
            AD ui1(Mu[ic-1],ic-1);
            AD ui2(Mu[ic+nx-1],ic+nx-1);
            F[ic] = (((0-(2*Ad))+(0-(2*Bd))+(0-(Cd))+(0-(Dd)))*ui)+((Ad)*ui1)+(((Bd)+(Dd))*ui2)+((Ad)*f2d+Bd*g1d+(Cd)*f2d+Ed);
        }
        else if(ic==n-nx+1){
            AD ui1(Mu[ic+1],ic+1);
            AD ui2(Mu[ic-nx+1],ic-nx+1);
            F[ic] = (((0-(2*Ad))+(0-(2*Bd))+(0-(Cd))+(0-(Dd)))*ui)+(((Ad)+(Cd))*ui1)+((Bd)*ui2)+((Ad)*f1d+Bd*g2d+(Dd)*g2d+Ed);
        }
        else if(ic==n-1){
            AD ui1(Mu[ic-1],ic-1);
            AD ui2(Mu[ic-nx+1],ic-nx+1);
            F[ic] = (((0-(2*Ad))+(0-(2*Bd))+(0-(Cd))+(0-(Dd)))*ui)+((Ad)*ui1)+((Bd)*ui2)+((Ad)*f2d+Bd*g2d+(Cd)*f2d+(Dd)*g2d+Ed);
        }
        
        else if(ic>0 && ic<nx-2){
            AD ui1(Mu[ic-1],ic-1);
            AD ui2(Mu[ic+1],ic+1);
            AD ui3(Mu[ic+nx-1],ic+nx-1);
            F[ic] = (((0-(2*Ad))+(0-(2*Bd))+(0-(Cd))+(0-(Dd)))*ui)+(((Ad)+(Cd))*ui2)+((Ad)*ui1)+(((Bd)+(Dd))*ui3)+(Bd*g1d+Ed);
        }
        else if(ic>n-nx+1 && ic<n-1){
            AD ui1(Mu[ic-1],ic-1);
            AD ui2(Mu[ic+1],ic+1);
            AD ui3(Mu[ic-nx+1],ic-nx+1);
            F[ic] = (((0-(2*Ad))+(0-(2*Bd))+(0-(Cd))+(0-(Dd)))*ui)+(((Ad)+(Cd))*ui2)+((Ad)*ui1)+((Bd)*ui3)+(Bd*g2d+(Dd)*g2d+Ed);
        }
        else if((ic+1)%(nx-1)==1 && cond>1 && cond<ny-2){
            AD ui1(Mu[ic+1],ic+1);
            AD ui2(Mu[ic+nx-1],ic+nx-1);
            AD ui3(Mu[ic-nx+1],ic-nx+1);
            F[ic] = (((0-(2*Ad))+(0-(2*Bd))+(0-(Cd))+(0-(Dd)))*ui)+(((Ad)+(Cd))*ui1)+(((Bd)+(Dd))*ui2)+((Bd)*ui3)+(Ad*f1d+Ed);
        }
        else if((ic+1)%(nx-1)==0 && (ic+1)/(nx-1)>1 && (ic+1)/(nx-1)<ny-1){
            AD ui1(Mu[ic-1],ic-1);
            AD ui2(Mu[ic+nx-1],ic+nx-1);
            AD ui3(Mu[ic-nx+1],ic-nx+1);
            F[ic] = (((0-(2*Ad))+(0-(2*Bd))+(0-(Cd))+(0-(Dd)))*ui)+((Ad)*ui1)+(((Bd)+(Dd))*ui2)+((Bd)*ui3)+(Ad*f2d+(Cd)*f2d+Ed);
        }
        else{
            AD ui1(Mu[ic-1],ic-1);
            AD ui2(Mu[ic+1],ic+1);
            AD ui3(Mu[ic-nx+1],ic-nx+1);
            AD ui4(Mu[ic+nx-1],ic+nx-1);
            F[ic] = (((0-(2*Ad))+(0-(2*Bd))+(0-(Cd))+(0-(Dd)))*ui)+((Ad)*ui1)+(((Ad)+(Cd))*ui2)+(((Bd)+(Dd))*ui4)+((Bd)*ui3)+Ed;
        }
    }
}


template <class T>
Matrix<T> Discretizer<T> :: return_F(){
    Matrix<T> M(n,1);
    for(int i=0;i<n;i++){
        double f = F[i].get_f();
        M.setval(i,0,f);
    }
    return M;
}

template <class T>
Matrix<T> Discretizer<T> :: return_J(){
    Matrix<T> J(n,n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            double df = F[i].get_Df(j);
            J.setval(i,j,df);
        }
    }
    return J;
}


bool IsOperator(char c){  
    if(c == '+' || c == '-' || c == '*' || c == '/' || c == '^' )  
        return true;
    if(c >= 'A' && c<='O')
        return true;
    return false;  
}  
  
bool IsOperand(char c){
    if( c >= 'a' )  
        return true;  
    if((c >= '0' && c <= '9') || c=='.') 
        return true;  
    return false;  
}

int precedence(char op){  
    if(op == '+' || op == '-')                     
        return 1;  
    if (op == '*' || op == '/')  
        return 2;  
    if(op == '^')                                  
        return 3;
    if(op >= 'A' && op <= 'O')
        return 4;
    return 0; 
} 
 
bool eqlOrhigher (char op1, char op2){  
    int p1 = precedence(op1);  
    int p2 = precedence(op2);  
    if (p1 == p2)  
    {  
    if (op1 == '^' )  
    return false;  
    return true;  
    }  
    return  (p1>p2 ? true : false);  
    }  

string convert(string input){
    string infix = "";
    int len;
    len = input.length();
    for(int i = 0;i<len;i++){
        if(input[i]=='s' && input[i+1] =='i' && input[i+2] =='n' && input[i+3]!='h'){
            infix += 'A';
            i=i+2;
        }    
        else if(input[i]=='c' && input[i+1]=='o' && input[i+3]!='e' && input[i+3]!='h'){
            infix += 'B';
            i=i+2;
        }    
        else if(input[i]=='t' && input[i+1] =='a' && input[i+3]!='h'){
            infix += 'C';
            i=i+2;
        }    
        else if(input[i]=='c' && input[i+1]=='o' && input[i+3]=='e'){
            infix += 'D';
            i=i+4;
        }    
        else if(input[i]=='s' && input[i+1]=='e' && input[i+2]=='c'){
            infix += 'E';
            i=i+2;
        }    
        else if(input[i]=='c' && input[i+1]=='o' && input[i+2]=='t'){
            infix += 'F';
            i=i+2;
        }
        else if(input[i]=='s' && input[i+1] =='i' && input[i+3]=='h'){
            infix += 'G';
            i=i+3;
        }    
        else if(input[i]=='c' && input[i+1]=='o' && input[i+3]=='h'){
            infix += 'H';
            i=i+3;
        }    
        else if(input[i]=='t' && input[i+1]=='a' && input[i+3]=='h'){
            infix += 'I';
            i=i+3;
        }    
        else if(input[i]=='a' && input[i+1]=='r' && input[i+3]=='s'){
            infix += 'J';
            i=i+5;
        }
        else if(input[i]=='a' && input[i+1]=='r' && input[i+3]=='c'){
            infix += 'K';
            i=i+5;
        }
        else if(input[i]=='a' && input[i+1]=='r' && input[i+3]=='t'){
            infix += 'L';
            i=i+5;
        }
        else if(input[i]=='a' && input[i+1]=='b' && input[i+2]=='s'){
            infix += 'M';
            i=i+2;
        }
        else if(input[i]=='l' && input[i+1]=='o' && input[i+2]=='g'){
            infix += 'N';
            i=i+2;
        }
        else if(input[i]=='e' && input[i+1]=='x' && input[i+2]=='p'){
            infix += 'O';
            i=i+2;
        }
        else{
            infix += input[i];
        }
    }
    return infix;
}
    
string infix_to_postfix(string infix){  
    stack <char> S;  
    string postfix ="";    
    char ch;  
      
    S.push( '(' );  
    infix += ')';  
    
      
      
    for(int i = 0; i<infix.length(); i++){  
        ch = infix[i];  
        if(ch == ' ')  
        continue;  
        else if(ch == '(')  
        S.push(ch);  
        else if(IsOperand(ch))  
        postfix += ch;  
        else if(IsOperator(ch)){  
        while(!S.empty() && eqlOrhigher(S.top(), ch)){  
        postfix += S.top();  
        S.pop();  
        }  
        S.push(ch);  
        }  
        else if(ch == ')'){  
        while(!S.empty() && S.top() != '('){  
        postfix += S.top();  
        S.pop();  
        }  
        S.pop();  
        }
    }  
    return postfix;  
}
