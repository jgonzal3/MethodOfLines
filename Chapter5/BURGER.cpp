#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include "dssCPP.h"


using namespace std;

#define ERROR 1
#define CASE 10
#define BC 2

#include <boost/numeric/odeint.hpp>

double analytical_solution (double, double);
double analytical_solution_green (double , double);
ostream& operator<<(ostream& os, const array_type& a);
double phi(double x, double t);


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
        for( int i=0 ; i<SIZE-1 ; ++i ){
        	double s_x = i*(1.0/200.0);
        	double real_u = phi(s_x,t);
        	cout<<x[i]-real_u<<"\n";
       		fout <<x[i]-real_u << ",";
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

void burgers_solution ( const array_type &u , array_type &ut , const double t/* t */ )
{

	// u + uxu - µ u   = 0	, VIS = µ/ρ
	//  t     x  ρ  xx

	// u = VISu  - uxu
	//  t      xx     x

	array_type U = u;
	array_type Ux(SIZE);
	array_type Uxx(SIZE);

	double Xl=0.0;
	double Xu=1.0;

	// BC at x = 0
	U[0] = phi(Xl,t);

	// BC at x = 1
	U[SIZE-1] = phi(Xu,t);

	switch (CASE){
	case 2:
		Ux  = dss02(Xl,Xu,SIZE,U); 	Uxx = dss02(Xl,Xu,SIZE,Ux);
		break;
	case 4:
		Ux  = dss04(Xl,Xu,SIZE,U); 	Uxx = dss04(Xl,Xu,SIZE,Ux);
		break;
	case 6:
		Ux  = dss06(Xl,Xu,SIZE,U); 	Uxx = dss06(Xl,Xu,SIZE,Ux);
		break;
	case 8:
		Ux  = dss08(Xl,Xu,SIZE,U); 	Uxx = dss08(Xl,Xu,SIZE,Ux);
		break;
	case 10:
		Ux  = dss10(Xl,Xu,SIZE,U); 	Uxx = dss10(Xl,Xu,SIZE,Ux);
		break;
	default:
		exit(0);
	}

	for (int i=0; i<SIZE ; i++)
	    ut[i]=VIS*Uxx[i]-U[i]*Ux[i];
}

void burgers_solution_BC_Newman ( const array_type &u , array_type &ut , const double t/* t */ )
{

	// u + uxu - µ u   = 0	, VIS = µ/ρ
	//  t     x  ρ  xx

	// u = VISu  - uxu
	//  t      xx     x

	array_type U = u;
	array_type Ux(SIZE);
	array_type Uxx(SIZE);

	double Xl=0.0;
	double Xu=1.0;

	switch (CASE){
	case 2:
		Ux  = dss02(Xl,Xu,SIZE,U);
		//BC at x = 0
		Ux[0]=0.0;
		//BC at x = 1
		Ux[SIZE-1]=0.0;
		Uxx = dss02(Xl,Xu,SIZE,Ux);
		break;
	case 4:
		Ux  = dss04(Xl,Xu,SIZE,U);
		//BC at x = 0
		Ux[0]=0.0;
		//BC at x = 1
		Ux[SIZE-1]=0.0;
		Uxx = dss04(Xl,Xu,SIZE,Ux);
		break;
	case 6:
		Ux  = dss06(Xl,Xu,SIZE,U);
		//BC at x = 0
		Ux[0]=0.0;
		//BC at x = 1
		Ux[SIZE-1]=0.0;
		Uxx = dss06(Xl,Xu,SIZE,Ux);
		break;
	case 8:
		Ux  = dss08(Xl,Xu,SIZE,U);
		//BC at x = 0
		Ux[0]=0.0;
		//BC at x = 1
		Ux[SIZE-1]=0.0;
		Uxx = dss08(Xl,Xu,SIZE,Ux);
		break;
	case 10:
		Ux  = dss10(Xl,Xu,SIZE,U);
		//BC at x = 0
		Ux[0]=0.0;
		//BC at x = 1
		Ux[SIZE-1]=0.0;Uxx = dss10(Xl,Xu,SIZE,Ux);
		break;
	default:
		exit(0);
	}

	for (int i=0; i<SIZE ; i++)
	    ut[i]=VIS*Uxx[i]-U[i]*Ux[i];
}

int main(int /* argc */ , char** /* argv */) {

    using namespace std;
    using namespace boost::numeric::odeint;

    array_type u(SIZE);
    array_type X(SIZE);

    double xl = 0.0;
    double xu = 1.0;
    double dx = (xu-xl)/(SIZE-1);

    for (int i=0; i<SIZE; i++) {
    	X[i] = xl+(i)*dx;
    	u[i]=phi(X[i],0.0);
    }

    void (*FuncPtr) (const array_type&, array_type&, const double);
    if (BC == 1) // Dirichlet
    	FuncPtr = burgers_solution;
    if (BC == 2) // Newman
    	FuncPtr = burgers_solution_BC_Newman;

    string fn= "/home/julio/temp/Burgers"+ to_string(CASE) +".csv";

    std::ofstream out;
    out.open(fn); // append instead of overwrite
    out <<"t,";
    for (auto i : X)
    	out << i <<",";

    out << "\n";
    out.close();

    //size_t steps = integrate( FuncPtr, u, 0.0 , 1.0 , 0.1, streaming_observer( cout ));
    //size_t steps = integrate( FuncPtr, u, 0.0 , 1.0 , 0.1, output_observer_append( fn ));
    size_t steps = integrate( FuncPtr, u, 0.0 , 1.0, 0.1, output_observer_append_error( fn ));

    runge_kutta4< state_type > stepper;

//	//[ harm_iterator_const_step]
//	std::for_each( make_const_step_time_iterator_begin( stepper , FuncPtr, u , 0.0 , 0.30, 0.005 ) ,
//				   make_const_step_time_iterator_end( stepper , FuncPtr, u ) ,
//				   []( std::pair< const state_type & , const double & > x ) {
//	    				string fn_= "/home/julio/temp/Burgers.csv";
//
//						//if (x.second == MYTIME)
//							std::ofstream fout(fn_, std::ios_base::app);
//							fout.open(fn_);
//							fout<<x.second<<",";
//							for (int k=0; k<x.first.size(); k++) {
//								double x_v = k*(1.0/200.0);
//								double temp = phi(x_v,x.second);
//								double error = (temp-x.first[k]);
//								cout << x.second << " "<<x_v<<" "<< x.first[k]<<" "<<error<< "\n";
//								fout <<x.first[k] << ",";
//								cout<<error<<"\n";
//							}
//							fout<<"\n";
//							//fout.close();
//	} );
//	//]



    //cout << "!!!Hello World!!! with " <<steps<< endl; // prints !!!Hello World!!!


	return 0;
}

double phi(double x, double t) {

	/* The one-dimensional Burgers' equation is

	   ut = mu*uxx - u*ux

	 where mu is a positive constant.

	 An analytical solution to the PDE is

	             0.1*exp(-a) + 0.5*exp(-b) + exp(-c)
	   ua(x,t) = -----------------------------------
	                 exp(-a) + exp(-b) + exp(-c)

	 where

	   a = (0.05/mu)*(x - 0.5 + 4.95*t)

	   b = (0.25/mu)*(x - 0.5 + 0.75*t)

	   c = (0.5/mu)*(x - 0.375)

	*/

	// Function phi computes the exact solution of Burgers’ equation for comparison with
	// the numerical solution. It is also used to define the initial and boundary conditions
	// for the numerical solution.

	//
	// Analytical solution
	//Euler, Navier Stokes, and Burgers Equations 107
	double a = (0.05/VIS)*(x-0.5+4.95*t);
	double b = (0.25/VIS)*(x-0.5+0.75*t);
	double c=( 0.5/VIS)*(x-0.375);

	return (0.1*exp(-a) + 0.5*exp(-b) + exp(-c))/(exp(-a) + exp(-b) + exp(-c));
}

double analytical_solution (double x, double time) {

	return (exp(-(M_PI*M_PI/4.0)*time) *  sin(M_PI_2*x));
}

double analytical_solution_green (double time, double x) {

	return (exp(-x*x/(4.0*time))/(2.0*sqrt(M_PI*time)));

}

ostream& operator<<(ostream& os, const array_type& a)
{
	int count=0;
	os<<"Array type:\n";
	for (auto i : a) {
		cout<<i<<"\n";
	count++;
	}

		cout<<"\n";

    return os;
}

