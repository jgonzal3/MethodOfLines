//============================================================================
// Name        : GreenPDE.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include "dssCPP.h"

using namespace std;

#define MYTIME 0.5
#define ERROR 1

#include <boost/numeric/odeint.hpp>

int invoke(int x, int y, int (*func)(int, int));
double analytical_solution (double, double);
double analytical_solution_green (double , double);
ostream& operator<<(ostream& os, const array_type& a);

struct output_observer_append
{
    string filename_;
    size_t count_;
    output_observer_append( const string &filename ) : filename_( filename ) , count_( 0 ) { }

    void operator()( const array_type &x , double t )
    {
    	int n=x.size();
    	std::ofstream fout(filename_, std::ios_base::app);
        //fout<<t<<",";
        for( int i=0 ; i<n-1 ; ++i ){
        	fout <<x[i] << ",";
        }
        fout<<x[n-1];
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

        if (t == MYTIME) {
        	m_out << t << ' '<<"\n";
        	for (size_t i=0; i<x.size(); i++) {
        		double x_v = -10.0+0.2*i;
        		double temp = analytical_solution_green(MYTIME,x_v);
        		if (ERROR == 0)
        			m_out<<"x="<<x_v<<", x[i]="<<x[i]<<" : "<< (temp-x[i]) <<"\n";
        		else
        			m_out<< (temp-x[i]) <<"\n";
        	}
        }
    }
};

void system_of_equation_v44 ( const array_type &u , array_type &ut , const double /* t */ )
{
	int n = 101;
	array_type temp = u;
	array_type ux(n,0);
	int nl = 1;
	int nu = 1;
	array_type uxx =  dss44(-10.0, 10.0, n, temp, ux, nl, nu);

	ut = uxx;
}


void system_of_equation_v50 ( const array_type &u , array_type &ut , const double /* t */ )
{
	int n = 101;
	array_type temp = u;
	array_type ux(n,0);
	int nl = 1;
	int nu = 1;
	array_type uxx =  dss50(-10.0, 10.0, n, temp, ux, nl, nu);

	ut = uxx;
}

void system_of_equation_green ( const array_type &u , array_type &ut , const double /* t */ )
{
	int n=101;
	array_type temp = u;
	temp[0] = 0.0;
	temp[n-1] = 0.0;
	array_type ux =  dss08(-10.0, 10.0, n, temp);
	array_type uxx = dss08(-10.0, 10.0, n, ux);

	ut = uxx;

	ut[0] = 0.0;
	ut[n-1]=0.0;

}

int main(int /* argc */ , char** /* argv */) {

    using namespace std;
    using namespace boost::numeric::odeint;

    long n = 101;
    double xl=-10.0;
    double xu= 10.0;
    double dx=(xu-xl)/(n-1);
    array_type x(n);
    array_type u(n);

    for (int k=0; k<n; k++) {
    	x[k]=xl+k*dx;
    }

	// Initial condition

    int n2=(n-1)/2+1;
    for (int i = 0; i<n; i++) {
      if (i == n2-1)
    	  u[i]=25.0*dx;
      else
    	  u[i]=0.0;
    }

    int mf = 3;

    void (*FuncPtr) (const array_type&, array_type&, const double);
    if (mf == 4)
    	FuncPtr = system_of_equation_v50;
    else if (mf == 3)
    	FuncPtr = system_of_equation_v44;
    else if (mf == 2)
    	FuncPtr = system_of_equation_green;
    else
    	FuncPtr = mol_green;

    string fn= "/home/julio/temp/green.csv";
   	//size_t steps = integrate( FuncPtr, u, 0.0 , 0.5 , 0.1, streaming_observer( cout ));
    size_t steps = integrate( FuncPtr, u, 0.0 , 2.0 , 0.1, output_observer_append( fn ));

    runge_kutta4< state_type > stepper;

	//[ harm_iterator_const_step]
//	std::for_each( make_const_step_time_iterator_begin( stepper , FuncPtr, u , 0.01 , 0.5, 0.01 ) ,
//				   make_const_step_time_iterator_end( stepper , FuncPtr, u ) ,
//				   []( std::pair< const state_type & , const double & > x ) {
//						//if (x.second == MYTIME)
//							for (size_t k=0; k<x.first.size(); k++) {
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

