#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include "dssCPP.h"


using namespace std;

#define MYTIME 0.5
#define ERROR 1
#define SIZE 41

#include <boost/numeric/odeint.hpp>

double analytical_solution (double, double);
double analytical_solution_green (double , double);
ostream& operator<<(ostream& os, const array_type& a);
double u(double, double);
double v(double, double);


struct output_observer_append
{
    string filename_;
    size_t count_;
    output_observer_append( const string &filename ) : filename_( filename ) , count_( 0 ) { }

    void operator()( const array_type &x , double t )
    {
    	std::ofstream fout(filename_, std::ios_base::app);
    	fout<<t<<",";
        for( int i=0 ; i<SIZE-1 ; ++i ){
       		fout <<x[i] << ",";
        }
        fout<<x[SIZE-1];
        fout<<"\n";
        ++count_;
    }
};

struct output_observer_append_error
{
    string filename_;
    size_t count_;
    output_observer_append_error( const string &filename ) : filename_( filename ) , count_( 0 ) { }

    void operator()( const array_type &x , double t )
    {
    	std::ofstream fout(filename_, std::ios_base::app);
    	fout<<t<<",";
        for( int i=0 ; i<(SIZE-1)/2 ; ++i ){
        	double s_x = i*(1.0/40.0);
        	double error_u = u(s_x,t);
        	cout<<error_u- x[i]<<"\n";
        	double error_v = v(s_x,t) - x[i];
       		fout <<error_u << ",";
        }
        fout<<x[SIZE-1];
        fout<<"\n";
        ++count_;
    }
};

struct streaming_observer
{
    std::ostream& m_out;

    streaming_observer( std::ostream &out ) : m_out( out ) { }

    struct write_element
    {
        std::ostream &m_out;
        write_element( std::ostream &out ) : m_out( out ) { };

        template< class T >
        void operator()( const T &t ) const
        {
            m_out << "\t" << t;
        }
    };

    template< class State , class Time >
    void operator()( const State &x , const Time &t ) const
    {

		m_out << t << ' '<<"\n";
		for (auto i : x)
			m_out<<t<<" : "<<i<<endl;
    }
};
void system_of_equation_PDE_v2 ( const array_type &u , array_type &ut , const double t/* t */ )
{
	//                                                                  -4x
	// u  = ((v − 1)u )  + (16xt − 2t − 16(v − 1)(u − 1))(u − 1) + 10x e
	//  t            x x
    //                                         -4x
	// vt = v   + u  + 4u − 4 + x2 − 2t − 10t e
	//		 xx    x

	double e4=pow(exp(1.0),4);

	array_type U(SIZE), V(SIZE);
	array_type Ut(SIZE),Vt(SIZE);
	array_type Vx(SIZE),Ux(SIZE);
	array_type Vxx(SIZE), Uxx(SIZE);
	array_type X(SIZE);

	// Spatial increment
	double  dx=1.0/(SIZE-1);

	// Values of x along the spatial grid
	for (int k=0; k<SIZE; k++)
		X[k]=dx*k;

	double Xl=X[0];
	double Xu=X[SIZE-1];

	for (int i=0; i<SIZE; i++) {
	    U[i]=u[i];
	    V[i]=u[i+SIZE];
	}

	// Note that the derivatives in t, u (0) and v (0), are set to zero to maintain the
	//                                  t         t
	// constant values of u(0) and v(0) at the boundary x = 0.

	U[0]=1.0;
	Ut[0]=0.0;
	V[0]=1.0;
	Vt[0]=0.0;

	Ux=dss02(Xl,Xu,SIZE,U);
	Vx=dss02(Xl,Xu,SIZE,V);

	//Boundary conditions at x = 1
	Ux[SIZE-1]=3.0-3.0*U[SIZE-1];
	Vx[SIZE-1]=e4*(U[SIZE-1]-1.0)/5.0;

	Vxx=dss02(Xl,Xu,SIZE,Vx);

	for (int i=0; i<SIZE; i++)
		Vx[i]=(V[i] - 1.0)*Ux[i];

	Uxx=dss02(Xl,Xu,SIZE,Vx);

	for (int i=1; i<SIZE ; i++) {
		double ex = exp(-4.0*X[i]);
	    Ut[i]=Uxx[i]+(16.0*X[i]*t-2.0*t-16.0*(V[i]-1.0))*(U[i]-1.0)+10.*X[i]*ex;
	    Vt[i]=Vxx[i]+Ux[i]+4.0*U[i]-4.0+X[i]*X[i]-2.0*t-10.0*t*ex;
	}

	for (int i=0; i<SIZE; i++) {
	    ut[i]   =Ut[i];
	    ut[i+SIZE] =Vt[i];
	}
}


void system_of_equation_PDE_v4 ( const array_type &u , array_type &ut , const double t/* t */ )
{
	//                                                                  -4x
	// u  = ((v − 1)u )  + (16xt − 2t − 16(v − 1)(u − 1))(u − 1) + 10x e
	//  t            x x
    //                                         -4x
	// vt = v   + u  + 4u − 4 + x2 − 2t − 10t e
	//		 xx    x

	double e4=pow(exp(1.0),4);

	array_type U(SIZE), V(SIZE);
	array_type Ut(SIZE),Vt(SIZE);
	array_type Vx(SIZE),Ux(SIZE);
	array_type Vxx(SIZE), Uxx(SIZE);
	array_type X(SIZE);

	// Spatial increment
	double  dx=1.0/(SIZE-1);

	// Values of x along the spatial grid
	for (int k=0; k<SIZE; k++)
		X[k]=dx*k;

	double Xl=X[0];
	double Xu=X[SIZE-1];

	for (int i=0; i<SIZE; i++) {
	    U[i]=u[i];
	    V[i]=u[i+SIZE];
	}

	// Note that the derivatives in t, u (0) and v (0), are set to zero to maintain the
	//                                  t         t
	// constant values of u(0) and v(0) at the boundary x = 0.

	U[0]=1.0;
	Ut[0]=0.0;
	V[0]=1.0;
	Vt[0]=0.0;

	Ux=dss04(Xl,Xu,SIZE,U);
	Vx=dss04(Xl,Xu,SIZE,V);

	//Boundary conditions at x = 1
	Ux[SIZE-1]=3.0-3.0*U[SIZE-1];
	Vx[SIZE-1]=e4*(U[SIZE-1]-1.0)/5.0;

	Vxx=dss04(Xl,Xu,SIZE,Vx);

	for (int i=0; i<SIZE; i++)
		Vx[i]=(V[i] - 1.0)*Ux[i];

	Uxx=dss04(Xl,Xu,SIZE,Vx);

	for (int i=1; i<SIZE ; i++) {
		double ex = exp(-4.0*X[i]);
	    Ut[i]=Uxx[i]+(16.0*X[i]*t-2.0*t-16.0*(V[i]-1.0))*(U[i]-1.0)+10.*X[i]*ex;
	    Vt[i]=Vxx[i]+Ux[i]+4.0*U[i]-4.0+X[i]*X[i]-2.0*t-10.0*t*ex;
	}

	for (int i=0; i<SIZE; i++) {
	    ut[i]   =Ut[i];
	    ut[i+SIZE] =Vt[i];
	}

}

void system_of_equation_PDE_v6 ( const array_type &u , array_type &ut , const double t/* t */ )
{
	//                                                                  -4x
	// u  = ((v − 1)u )  + (16xt − 2t − 16(v − 1)(u − 1))(u − 1) + 10x e
	//  t            x x
    //                                         -4x
	// vt = v   + u  + 4u − 4 + x2 − 2t − 10t e
	//		 xx    x

	double e4=pow(exp(1.0),4);

	array_type U(SIZE), V(SIZE);
	array_type Ut(SIZE),Vt(SIZE);
	array_type Vx(SIZE),Ux(SIZE);
	array_type Vxx(SIZE), Uxx(SIZE);
	array_type X(SIZE);

	// Spatial increment
	double  dx=1.0/(SIZE-1);

	// Values of x along the spatial grid
	for (int k=0; k<SIZE; k++)
		X[k]=dx*k;

	double Xl=X[0];
	double Xu=X[SIZE-1];

	for (int i=0; i<SIZE; i++) {
	    U[i]=u[i];
	    V[i]=u[i+SIZE];
	}
	// Note that the derivatives in t, u (0) and v (0), are set to zero to maintain the
	//                                  t         t
	// constant values of u(0) and v(0) at the boundary x = 0.

	U[0]=1.0;
	Ut[0]=0.0;
	V[0]=1.0;
	Vt[0]=0.0;

	Ux=dss06(Xl,Xu,SIZE,U);
	Vx=dss06(Xl,Xu,SIZE,V);

	//Boundary conditions at x = 1
	Ux[SIZE-1]=3.0-3.0*U[SIZE-1];
	Vx[SIZE-1]=e4*(U[SIZE-1]-1.0)/5.0;

	Vxx=dss06(Xl,Xu,SIZE,Vx);

	for (int i=0; i<SIZE; i++)
		Vx[i]=(V[i] - 1.0)*Ux[i];

	Uxx=dss06(Xl,Xu,SIZE,Vx);

	for (int i=1; i<SIZE ; i++) {
		double ex = exp(-4.0*X[i]);
	    Ut[i]=Uxx[i]+(16.0*X[i]*t-2.0*t-16.0*(V[i]-1.0))*(U[i]-1.0)+10.*X[i]*ex;
	    Vt[i]=Vxx[i]+Ux[i]+4.0*U[i]-4.0+X[i]*X[i]-2.0*t-10.0*t*ex;
	}

	  for (int i=0; i<SIZE; i++) {
	    ut[i]   =Ut[i];
	    ut[i+SIZE] =Vt[i];
	  }
}

void system_of_equation_PDE_v10 ( const array_type &u , array_type &ut , const double t/* t */ )
{
	//                                                                  -4x
	// u  = ((v − 1)u )  + (16xt − 2t − 16(v − 1)(u − 1))(u − 1) + 10x e
	//  t            x x
    //                                         -4x
	// vt = v   + u  + 4u − 4 + x2 − 2t − 10t e
	//		 xx    x

	double e4=pow(exp(1.0),4);

	array_type U(SIZE), V(SIZE);
	array_type Ut(SIZE),Vt(SIZE);
	array_type Vx(SIZE),Ux(SIZE);
	array_type Vxx(SIZE), Uxx(SIZE);
	array_type X(SIZE);

	// Spatial increment
	double  dx=1.0/(SIZE-1);

	// Values of x along the spatial grid
	for (int k=0; k<SIZE; k++)
		X[k]=dx*k;

	double Xl=X[0];
	double Xu=X[SIZE-1];

	for (int i=0; i<SIZE; i++) {
	    U[i]=u[i];
	    V[i]=u[i+SIZE];
	}

	// Note that the derivatives in t, u (0) and v (0), are set to zero to maintain the
	//                                  t         t
	// constant values of u(0) and v(0) at the boundary x = 0.

	U[0]=1.0;
	Ut[0]=0.0;
	V[0]=1.0;
	Vt[0]=0.0;

	Ux=dss10(Xl,Xu,SIZE,U);
	Vx=dss10(Xl,Xu,SIZE,V);

	//Boundary conditions at x = 1
	Ux[SIZE-1]=3.0-3.0*U[SIZE-1];
	Vx[SIZE-1]=e4*(U[SIZE-1]-1.0)/5.0;

	Vxx=dss10(Xl,Xu,SIZE,Vx);

	for (int i=0; i<SIZE; i++)
		Vx[i]=(V[i] - 1.0)*Ux[i];

	Uxx=dss10(Xl,Xu,SIZE,Vx);

	for (int i=1; i<SIZE ; i++) {
		double ex = exp(-4.0*X[i]);
	    Ut[i]=Uxx[i]+(16.0*X[i]*t-2.0*t-16.0*(V[i]-1.0))*(U[i]-1.0)+10.*X[i]*ex;
	    Vt[i]=Vxx[i]+Ux[i]+4.0*U[i]-4.0+X[i]*X[i]-2.0*t-10.0*t*ex;
	}

	  for (int i=0; i<SIZE; i++) {
	    ut[i]   =Ut[i];
	    ut[i+SIZE] =Vt[i];
	  }



}


int main(int /* argc */ , char** /* argv */) {

    using namespace std;
    using namespace boost::numeric::odeint;

    double t0=0.0;
    /*
    Function inital_1 is called by the main program to define
    the initial conditions in the MOL solution of two nonlinear PDEs
    */
    array_type u(2*SIZE,1.0);

    string fn= "/home/julio/temp/2PDEs.csv";

    std::ofstream out;
    out.open(fn); // append instead of overwrite
    out <<"t,";
    for (int i=0;i<SIZE;i++)
    	out <<i*(1.0/(SIZE-1))<<",";

    for (int i=0;i<SIZE;i++)
    	out <<i*(1.0/(SIZE-1))<<",";

    out << "\n";
    out.close();

    int mf = 1;

    void (*FuncPtr) (const array_type&, array_type&, const double);
    if (mf == 2)
    	FuncPtr = system_of_equation_PDE_v2;
    if (mf == 4)
    	FuncPtr = system_of_equation_PDE_v4;
    if (mf == 6)
        FuncPtr = system_of_equation_PDE_v6;
    else
    	FuncPtr = system_of_equation_PDE_v10;


    //size_t steps = integrate( FuncPtr, u, 0.0 , 1.0 , 0.1, streaming_observer( cout ));
    //size_t steps = integrate( FuncPtr, u, 0.0 , 0.05 , 0.1, output_observer_append( fn ));
    size_t steps = integrate( FuncPtr, u, 0.0 , 1.0, 0.1, output_observer_append_error( fn ));

    runge_kutta4< state_type > stepper;

	//[ harm_iterator_const_step]
//	std::for_each( make_const_step_time_iterator_begin( stepper , FuncPtr, u , 0.01 , 0.5, 0.01 ) ,
//				   make_const_step_time_iterator_end( stepper , FuncPtr, u ) ,
//				   []( std::pair< const state_type & , const double & > x ) {
//						//if (x.second == MYTIME)
//							for (size_t k=0; k<x.first.SIZE(); k++) {
//								double x_v = -10.0+0.2*k;
//								double temp = analytical_solution_green(MYTIME,x_v);
//								double error = (temp-x.first[k]);
//								cout << x.second << " "<<x_v<<" "<< x.first[k]<<" "<<error<< "\n";
//							}
//	} );
//	//]



    cout << "!!!Hello World!!! with " <<steps<< endl; // prints !!!Hello World!!!


	return 0;
}

double analytical_solution (double x, double time) {

	return (exp(-(M_PI*M_PI/4.0)*time) *  sin(M_PI_2*x));
}

double analytical_solution_green (double time, double x) {

	return (exp(-x*x/(4.0*time))/(2.0*sqrt(M_PI*time)));

}

double u(double x, double t){
	return (1 + 10.0*x*t*exp(-4.0*x));
}

double v(double x, double t) {
	return (1 + x*x*t);
}


ostream& operator<<(ostream& os, const array_type& a)
{
	int count=0;
	os<<"Array type:";
	for (auto i : a) {
		cout<<i<<"\n";
	count++;
	}

		cout<<"\n";

    return os;
}

