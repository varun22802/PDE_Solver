#pragma once
using namespace std;

int Num_Var = 0;
void setNum_Var(int n){
    Num_Var = n;
}

class AD{
private:
	double f;
	vector<double> df;
	int id;
 public:
	AD();
	AD(double,int);
    void print();
	double get_f();
	double get_Df(int);
	int get_id();
	void set_id(int);
	void setIndVar();
	void set_f(double);
	void set_df(int, double);
	vector<double> getGradient();
	friend double** getJacobian(AD*);
	AD operator+(AD); // f + g
	AD operator-(AD); // f - g
	AD operator*(AD); // f * g
	AD operator/(AD); // f / g
	AD operator^(double); // f^c
	// Do for the rest of all possible operations like,
	// c+f, c-f, c*f, c/f, c^f.

	
	friend AD operator+(double,AD);

	
	friend AD operator-(double,AD);

	
	friend AD operator*(double,AD);

	
	friend AD operator/(double,AD);

	
	friend AD operator^(double,AD);

	// f+c, f-c, f*c, f/c.

	AD operator+(double); // f + c
	AD operator-(double); // f - c
	AD operator*(double); // f*c
	AD operator/(double); // f/c
    AD operator^(AD); // f^g
    
	// sin(f), cos(f), tan(f), cosec(f), sec(f), cot(f).
	// arcsin(f), arccos(f), arctan(f), sinh(f), cosh(f), tanh(f).
	// log(f), exp(f), abs(f)

	friend AD sin(AD);
	
    friend AD cos(AD);
	
	friend AD tan(AD);


	friend AD cosec(AD);
	
	friend AD sec(AD);
	
	friend AD cot(AD);

	
	friend AD arcsin(AD);
	
	friend AD arccos(AD);
	
	friend AD arctan(AD);
	
	
	friend AD sinh(AD);
	
	friend AD cosh(AD);

	friend AD tanh(AD);


	friend AD log(AD);

	friend AD exp(AD);
	
	friend AD abs(AD);

};

AD :: AD(){
    this->f = 0;
    this->id = -1;
}

AD :: AD(double val,int n){
    this->f = val;
    this->id = n;
    vector<double> DF(Num_Var,0);
    if(id != -1){
        DF[id] = 1;
    }
    df = DF;
}

double AD :: get_f(){
    return f;
} 

int AD:: get_id(){
    return id;
}

void AD::set_id(int n){
    id = n;
}

void AD::setIndVar(){
    for(int i=0;i<this->id;i++)this->df[i] = 0;
    this->df[this->id] = 1;
    for(int i=this->id+1;i<Num_Var;i++)this->df[i] = 0;
}

double AD :: get_Df(int index){
    return df[index];
}


void AD :: set_f(double t){
    f = t;
} 


void AD :: set_df(int n, double t){
    df[n] = t;
}

AD AD :: operator+(AD g){
    AD h;
    h.f = this->f + g.f;
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = this->df[i]+g.df[i];
    h.df = DF;
    return h;
} 


AD AD :: operator-(AD g){
    AD h;
    h.f = this->f - g.f;
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = this->df[i]-g.df[i];
    h.df = DF;
    return h;
}


AD AD :: operator*(AD g){
    AD h;
    h.f = this->f * g.f;
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] =(this->f*g.df[i])+(this->df[i]*g.f);
    h.df = DF;
    return h;
} 


AD AD :: operator/(AD g){
    AD h;
    h.f = this->f / g.f;
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = (df[i]*g.f - g.df[i]*f)/(g.f*g.f);
    h.df = DF;
    return h;
}

AD AD :: operator+(double a){
    f += a;
    return *this;
}


AD AD :: operator-(double a){
    f -= a;
    return *this;
}


AD AD :: operator*(double a){
    f *= a;
    for(int i=0;i<Num_Var;i++)df[i] *= a; 
    return *this;
}


AD AD :: operator/(double a){
    f /= a;
    for(int i=0;i<Num_Var;i++)df[i] /= a;
    return *this;
}


AD AD :: operator^(double n){
    AD h;
    h.f = pow(f,n);
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = n*(df[i])*pow(f,n-1);
    h.df = DF;
    return h;
}


AD AD :: operator^(AD g){
    AD h;
    h.f = pow(f,g.f);
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = (h.f)*(log(pow(f,g.df[i])) + ((g.f*df[i])/f) );
    h.df = DF;
    return h;
}


AD operator+(double a, AD g){
    AD h;
    h = g+a;
    return h;
}


AD operator-(double a, AD g){
    AD h;
    h.set_f(a-g.f);
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = (-1)*g.df[i];
    h.df = DF;
    return h;
}


AD operator*(double a, AD g){
    AD h;
    h = g*a;
    return h;
}


AD operator/(double a, AD g){
    AD h;
    h = g^(-1);
    h = h*a;
    return h;
}


AD operator^(double a, AD g){
    AD h;
    h.f = pow(a,g.f);
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = g.df[i]*h.f*log(a);
    h.df = DF;
    return h;
}

AD sin(AD g){
    AD h;
    h.f = sin(g.f);
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = cos(g.f)*g.df[i];
    h.df = DF;
    return h;
}


AD cos(AD g){
    AD h;
    h.f = cos(g.f);
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = -1*sin(g.f)*g.df[i];
    h.df = DF;
    return h;
}


AD tan(AD g){
    AD h;
    AD i;
    h = sin(g);
    i = cos(g);
    return h/i;
}


AD cot(AD g){
    AD h;
    AD i;
    h = sin(g);
    i = cos(g);
    return i/h;
}


AD sec(AD g){
    AD h;
    h = cos(g);
    return 1.0/h;
}


AD cosec(AD g){
    AD h;
    h = sin(g);
    return 1.0/h;
}


AD log(AD g){
    AD h;
    h.f = log(g.get_f());
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = g.df[i]/g.f;
    h.df = DF;
    return h;
}


AD exp(AD g){
    AD h;
    h.f = exp(g.f);
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = g.df[i]*h.f;
    h.df = DF;
    return h;
}


AD arcsin(AD g){
    AD h;
    h.f = asin(g.f);
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = g.df[i]/(pow((1-pow(g.f,2)),.5));
    h.df = DF;
    return h;
}


AD arccos(AD g){
    AD h;
    h.f = acos(g.f);
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = (-1)*g.df[i]/(pow((1-pow(g.f,2)),.5));
    h.df = DF;
    return h;
}


AD arctan(AD g){
    AD h;
    h.f = atan(g.f);
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = g.df[i]/(1+(pow(g.f,2)));
    h.df = DF;
    return h;
}


AD sinh(AD g){
    AD h;
    h.f = sinh(g.f);
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = g.df[i]*cosh(g.f);
    h.df = DF;
    return h;
}



AD cosh(AD g){
    AD h;
    h.f = cosh(g.f);
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = g.df[i]*sinh(g.f);
    h.df = DF;
    return h;
}



AD tanh(AD g){
    AD h;
    h.f = tanh(g.f);
    vector<double> DF(Num_Var,0);
    for(int i=0;i<Num_Var;i++)DF[i] = g.df[i]*(1-(pow(h.f,2)));
    h.df = DF;
    return h;
}


AD abs(AD g){
    AD h;
    if(g.get_f() > 0 ){
        h.f = g.f;
        vector<double> DF(Num_Var,0);
        for(int i=0;i<Num_Var;i++)DF[i] = g.df[i];
        h.df = DF;  
    }
    else if(g.get_f() < 0){
        h.f = (-1)*g.f;
        vector<double> DF(Num_Var,0);
        for(int i=0;i<Num_Var;i++)DF[i] = (-1)*g.df[i];
        h.df = DF;
    }
    else{
        h.f = 0;
        vector<double> DF(Num_Var,0);
        for(int i=0;i<Num_Var;i++)DF[i] = log(-1);
        h.df = DF;
    }
    return h;
}


void AD :: print(){
    cout<<fixed<<setprecision(4);
    if(this->f == -0)f = 0;
    for(int i = 0;i<Num_Var;i++){
        if(df[i] == -0)df[i] == 0;
        cout<<df[i]<<" ";
    }
}