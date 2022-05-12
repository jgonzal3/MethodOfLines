//============================================================================
// Name        : MOL.cpp
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



#include <boost/numeric/odeint.hpp>

int invoke(int x, int y, int (*func)(int, int));
double analytical_solution (double x, double time);

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
        for( int i=1 ; i<=n-2 ; ++i ){
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
//        m_out << t;
//        for (size_t i=1; i<=x.size(); i++)
//        	m_out << "\t" << x[i];
//        m_out << "\n";
    	double temp = analytical_solution(0.5,t);
        m_out << t << ' ';
        //m_out<<abs(x[10] - temp) <<"\n";
        m_out<<x[10]<<","<<temp<<"\n";

    }
};

void system_of_equation_v42 ( const array_type &u , array_type &ut , const double /* t */ )
{
	int n = 21;
	array_type temp = u;
	array_type ux(n,0);
	temp[0] = 0.0;
	int nl = 1;
	int nu = 2;
	array_type uxx =  dss42(0.0, 1.0, n, temp, ux, nl, nu);

	ut = uxx;
}

void system_of_equation_v44 ( const array_type &u , array_type &ut , const double /* t */ )
{
	int n = 21;
	array_type temp = u;
	array_type ux(n,0);
	temp[0] = 0.0;
	int nl = 1;
	int nu = 2;
	array_type uxx =  dss44(0.0, 1.0, n, temp, ux, nl, nu);

	ut = uxx;
}

void system_of_equation_v46 ( const array_type &u , array_type &ut , const double /* t */ )
{
	int n = 21;
	array_type temp = u;
	array_type ux(n,0);
	temp[0] = 0.0;
	int nl = 1;
	int nu = 2;
	array_type uxx =  dss46(0.0, 1.0, n, temp, ux, nl, nu);

	ut = uxx;
}

void system_of_equation_v48 ( const array_type &u , array_type &ut , const double /* t */ )
{
	int n = 21;
	array_type temp = u;
	array_type ux(n,0);
	temp[0] = 0.0;
	int nl = 1;
	int nu = 2;
	array_type uxx =  dss48(0.0, 1.0, n, temp, ux, nl, nu);

	ut = uxx;
}

void system_of_equation_v50 ( const array_type &u , array_type &ut , const double /* t */ )
{
	int n = 21;
	array_type temp = u;
	array_type ux(n,0);
	temp[0] = 0.0;
	int nl = 1;
	int nu = 2;
	array_type uxx =  dss50(0.0, 1.0, n, temp, ux, nl, nu);

	ut = uxx;
}
void system_of_equation ( const array_type &u , array_type &ut , const double /* t */ )
{
	int n = u.size();
	n=21;
	array_type temp = u;
	temp[0] = 0.0;
	array_type ux =  dss02(0.0, 1.0, n, temp);
	ux[n-1] = 0.0;
	array_type uxx = dss02(0.0, 1.0, n, ux);

	ut = uxx;
}

void system_of_equation_v04 ( const array_type &u , array_type &ut , const double /* t */ )
{
	int n = u.size();
	n=21;
	array_type temp = u;
	temp[0] = 0.0;
	array_type ux =  dss04(0.0, 1.0, n, temp);
	ux[n-1] = 0.0;
	array_type uxx = dss04(0.0, 1.0, n, ux);

	ut = uxx;
}

void system_of_equation_v06 ( const array_type &u , array_type &ut , const double /* t */ )
{
	int n = u.size();
	n=21;
	array_type temp = u;
	temp[0] = 0.0;
	array_type ux =  dss06(0.0, 1.0, n, temp);
	ux[n-1] = 0.0;
	array_type uxx = dss06(0.0, 1.0, n, ux);

	ut = uxx;
}

void system_of_equation_v08 ( const array_type &u , array_type &ut , const double /* t */ )
{
	int n = u.size();
	n=21;
	array_type temp = u;
	temp[0] = 0.0;
	array_type ux =  dss08(0.0, 1.0, n, temp);
	ux[n-1] = 0.0;
	array_type uxx = dss08(0.0, 1.0, n, ux);

	ut = uxx;
}

void system_of_equation_v10 ( const array_type &u , array_type &ut , const double /* t */ )
{
	int n = u.size();
	n=21;
	array_type temp = u;
	temp[0] = 0.0;
	array_type ux =  dss10(0.0, 1.0, n, temp);
	ux[n-1] = 0.0;
	array_type uxx = dss10(0.0, 1.0, n, ux);

	ut = uxx;
}
int main(int /* argc */ , char** /* argv */) {

    using namespace std;
    using namespace boost::numeric::odeint;

	int n=21;

	// Initial condition
    array_type u(n);
    for (int i=0; i<n; i++) {
    	u[i] = sin(M_PI_2*i/(n-1));;
    }

    void (*FuncPtr) (const array_type&, array_type&, const double);
    //FuncPtr = system_of_equation_mol;
	//FuncPtr = system_of_equation;
	//FuncPtr = system_of_equation_v42;
	FuncPtr = mol;
	//FuncPtr = system_of_equation_v06;
	//FuncPtr = system_of_equation_v48;
	//FuncPtr = system_of_equation_v50;

   	size_t steps = integrate( FuncPtr, u, 0.0 , 2.5 , 0.1, streaming_observer( cout ));

    cout << "!!!Hello World!!! with " <<steps<< endl; // prints !!!Hello World!!!


	return 0;
}

double analytical_solution (double x, double time) {

	return (exp(-(M_PI*M_PI/4.0)*time) *  sin(M_PI_2*x));
}



