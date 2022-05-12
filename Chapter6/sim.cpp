//============================================================================
// Name        : sim.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
using namespace std;
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <fstream>

#include "dssCPP.h"
#define SIZE 201
#define CASE 1

using namespace std;
#include <boost/numeric/odeint.hpp>

double phi(double t, double x);
ostream& operator<<(ostream& os, const array_type& a);
array_type simp(double xl, double xu, int n, array_type u);
array_type operator*(const array_type& op1, const array_type& op2);

struct output_observer_append
{
    string filename_;
    size_t count_;
    output_observer_append( const string &filename ) : filename_( filename ) , count_( 0 ) { }

    void operator()( const array_type &x , double t )
    {
    	auto start = x.begin();
    	auto end   = x.begin() + SIZE;
    	array_type w2v2(SIZE);
    	array_type w(SIZE);
    	array_type v(SIZE);
    	copy(start+SIZE, end+SIZE-1, w.begin());
    	copy(start, end, v.begin());

    	for (int k=0; k<SIZE; k++)
    		w2v2[k] = v[k]*v[k] + w[k]*w[k];

    	std::ofstream fout(filename_, std::ios_base::app);
    	fout<<t<<",";
        for( int i=0 ; i<SIZE ; ++i ){
       		fout <<w2v2[i] << ",";
        }
        fout<<w2v2[SIZE-1];
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
//    	auto start = x.begin();
//    	auto end   = x.begin() + SIZE;
    	array_type w2v2(SIZE);
//    	array_type w(SIZE);
//    	array_type v(SIZE);
//    	copy(start+SIZE, end+SIZE-1, w.begin());
//    	copy(start, end, v.begin());

    	for (int k=0; k<SIZE; k++)
    		w2v2[k] = x[k]*x[k] + x[k+SIZE]*x[k+SIZE];

    	std::ofstream fout(filename_, std::ios_base::app);
    	fout<<t<<",";
        for( int i=0 ; i<SIZE ; ++i ){
        	double s_x = -10 + i*(50.0/200.0);
        	double real_u = phi(s_x,t);
        	cout<<w2v2[i]-real_u<<"\n";
       		fout <<w2v2[i]-real_u << ",";
        }
        fout<<w2v2[SIZE-1];
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
void schordingeer_solution ( const array_type &u , array_type &ut , const double t/* t */ )
{
	array_type v(SIZE,0),   vx(SIZE,0);
	array_type w(SIZE),   wx(SIZE,0);
	array_type vxx(SIZE), wxx(SIZE);
	array_type vt(SIZE),  wt(SIZE);
	array_type v2w2(SIZE);

	double Xl=-10.0;
	double Xu= 40.0;
	int nl=1;
	int nu=1;
	double q=1.0;

	for (int i=0;i<SIZE;i++) {
		v[i] = u[i];
		w[i] = u[i+SIZE];
	}

	//v[0 ]     = 0.0;   removed because no impact
	//v[SIZE-1] = 0.0;   removed because no impact
	vxx=dss44(Xl,Xu,SIZE,v,vx,nl,nu);

	//w[0  ] = 0.0;     removed because no impact
	//w[SIZE-1] = 0.0;  removed because no impact
	wxx=dss44(Xl,Xu,SIZE,w,wx,nl,nu);

	// ODEs at the boundaries
	vt[0] = 0.0;
	wt[0] = 0.0;

	vt[SIZE-1] = 0.0;
	wt[SIZE-1] = 0.0;

	// ODEs at the interior points
	//The Cubic Schrodinger Equation (133)
	for (int i=1; i<SIZE-2; i++) {
		// v**2 + w**2
		v2w2[i] = v[i]*v[i] + w[i]*w[i];
		// vt
		vt[i] = -wxx[i] - q*v2w2[i]*w[i];
		// wt
		wt[i] =  vxx[i] + q*v2w2[i]*v[i];
	}

	// Two vectors to one vector
	for (int k=0; k<SIZE; k++) {
		ut[k] =vt[k];
		ut[k+SIZE] = wt[k];
	}

}
double phi(double x, double t) {

// Function ua computes the analytical solution to the CSE

// The following exact solution is |u(x,t)|^2 { (√2 sech(x − t))^2 } which is compared with the numerical solution
	double sech = 2.0/(exp(x-t) + exp(-(x-t)));
	double sech2 = sech*sech;
	return (2.0*sech2);
	return (2.0*(2.0/(exp(x-t)+exp(-(x-t))))*(2.0/(exp(x-t)+exp(-(x-t)))));
}

int main() {
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!

	using namespace std;
    using namespace boost::numeric::odeint;

	array_type x(SIZE), v(SIZE), w(SIZE);
	array_type u(2*SIZE);
	double sch;

	double rt2 = sqrt(2.0);

    double xl = -10.0;
    double xu = 40.0;
    double dx = (xu-xl)/(SIZE-1);

    // IC over the spatial grid
	for (int i=0; i<SIZE; i++) {
		x[i] = xl + i*dx;
		sch=2.0/(exp(x[i])+exp(-x[i]));
		v[i] = rt2*cos(0.5*x[i])*sch;
		w[i] = rt2*sin(0.5*x[i])*sch;
	}
	 //
	 // Two vectors to one vector
	for (int k=0; k<SIZE; k++) {
		u[k  ]  = v[k];
		u[k+SIZE]  = w[k];
	}

	//cout<<u;
    void (*FuncPtr) (const array_type&, array_type&, const double);
    FuncPtr = schordingeer_solution;

    string fn= "/home/julio/temp/schodinger"+ to_string(CASE) +".csv";

    std::ofstream out;
    out.open(fn); // append instead of overwrite
    out <<"t,";
    for (int i=0; i<SIZE; i++)
    	out<<x[i]<<",";
    out<<"\n";
    //out.close();


    //size_t steps = integrate( FuncPtr, u, 0.0 , 1.0 , 0.1, streaming_observer( cout ));
    //size_t steps = integrate( FuncPtr, u, 0.0 , 30.0 , 0.1, output_observer_append( fn ));
    //size_t steps = integrate( FuncPtr, u, 0.0 , 30.0, 0.1, output_observer_append_error( fn ));

    runge_kutta4< state_type > stepper1;
    adams_bashforth< 8 , state_type > stepper;

    const double dt = 0.0001;
    array_type W(SIZE);
	array_type w2v2(SIZE);
	array_type V(SIZE);

	int count = 0;
	for( double t=0.0 ; t<30.0 ; t+= dt ) {
        stepper.do_step( FuncPtr, u , t , dt);
        auto start = u.begin();
        auto end   = u.begin() + SIZE;
        copy(start+SIZE, end+SIZE-1, W.begin());
        copy(start, end, V.begin());
        w2v2 = V*W;
        if (count%10000 == 0){
			for (auto i : w2v2){
				out<<i<<",";
			}
			out<<"\n";
        }
		count++;
	}


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

}

// File: simp.m
// From Chapter 6 - 'The Cubic Schrodinger Equation' of the book:
// William E Schiesser and Graham W Griffiths (2009).
// A Compendium of Partial Differential Equation Models
// - Method of Lines Analysis with Matlab,
// Cambridge University Press (ISBN-13: 9780521519861).
array_type simp(double xl, double xu, int n, array_type u) {

	// Function simp computes three integral invariants by Simpson's rule

	//	   / ∞                       i=n-2
	//    / u(x,t)dx  ≈ (h/3)*[|u1|2 + SUM(4|u |^2 + 2|u   |^2) + |u |^2 ]
	// −∞/                           i=2      i         i+1         n

	double q = 1.0;

	array_type v(n), w(n);
	array_type vx(n), wx(n);
	array_type uabs(n);
	array_type uint(3,0);
	array_type uxabs(n);
	array_type INT(n);

	//  One vector to two vectors
	for (int i=0;i<n;i++) {
		v[i] = u[i];
		w[i] = u[i+n];
	}

	for (int k=1;k<4; k++) {
		double h = (xu - xl)/(n-1);

	//   I1
		if( k == 1) {
//		 		 / ∞
//		u1(t) = / (|u(x, t)|^2) dx
//			 −∞/
			uabs[0  ] = u[0]*u[0] + v[0]*v[0];
			uabs[n-1] = u[n-1]*u[n-1] + v[n-1]*v[n-1];
			uint[0  ] = uabs[0] - uabs[n-1];
			for (int i = 2; i<n; i=i+2) {
				uabs[i-1]= u[i-1]*u[i-1]+v[i-1]*v[i-1];
				uabs[i  ]= u[i  ]*u[i  ]+v[i  ]*v[i  ];
				uint[0]  = uint[0] + 4.0*uabs[i-1]+2.0*uabs[i];
			}
			uint[0] = h/3.0*uint[0];
		}
		if ( k == 2) {
//		 		      /∞
//		u2(t) = 2 *  / (vw  − wv ) dx
//			      −∞/	   x	 x
			vx=dss04(xl,xu,n,v);
		    wx=dss04(xl,xu,n,w);
		    INT[0  ] = v[0  ]*wx[0  ]-w[0  ]*vx[0  ];
		    INT[n-1] = v[n-1]*wx[n-1]-w[n-1]*vx[n-1];
		    uint[1 ] = INT[0] - INT[n-1];
		    for (int i = 2; i < n; i=i+2) {
		       INT[i-1]=v[i-1]*wx[i-1]-w[i-1]*vx[i-1];
		       INT[i  ]=v[i  ]*wx[i  ]-w[i  ]*vx[i  ];
		       uint[1]=uint[1]+4.0*INT[i-1]+2.0*INT[i];
		    }
		    uint[1]=h/3.0*uint[1];
		}
		if ( k == 3) {
//		|u|^2 = |(u + i w )|^2 = u ^2 + w ^2
//				   x     x        x      x
//		|u|^4 = (|u|^2)^2
//		 		      ∞
//		u3(t) = 2 * int (|u |^2 − 0.5*q|u|^4) dx
//		  	       −∞	   x

			uxabs[0  ]=vx[0  ]*vx[0  ]+wx[0  ]*wx[0  ];
			uxabs[n-1]=vx[n-1]*vx[n-1]+wx[n-1]*wx[n-1];
			uint[2]=uxabs[0]-0.5*q*uabs[0]*uabs[0]-(uxabs[n-1]-0.5*q*uabs[n-1]*uabs[n-1]);
			for (int i=2; i<n; i=i+2) {
				uxabs[i-1]=vx[i-1]*vx[i-1]+wx[i-1]*wx[i-1];
				uxabs[i]  =vx[i  ]*vx[i  ]+wx[i  ]*wx[i  ];
				uint[2]=uint[2] + 4.0*(uxabs[i-1] - 0.5*q*uabs[i-1]*uabs[i-1]) + 2.0*(uxabs[i]-0.5*q*uabs[i]*uabs[i]);
			}
			uint[2]=2.0*h/3.0*uint[2];
		}
	}

	return uint;
}
//

array_type operator*(const array_type& op1, const array_type& op2)
{
	array_type op3(SIZE);


	for (int k=0; k<SIZE; k++) {
		op3[k] = op1[k]*op1[k] + op2[k]*op2[k];
	}

    return op3;
}

ostream& operator<<(ostream& os, const array_type& a)
{
	int count=0;
	//os<<"Array type:\n";
	for (auto i : a) {
		cout<<i<<",";
	count++;
	}

		cout<<"\n";

    return os;
}

//	std::vector<double> u = {0.000036425172583,
//	  -0.076916409156505,
//	   0.809573181179732,
//	   0.049157678262127,
//	  -0.039697367676954,
//	  -0.005881510298376,
//	   0.015282106043724,
//	  -0.001296668619999,
//	  -0.006483537132089,
//	   0.000000000000000,
//	   0.000123135842523,
//	  -0.013491832721215,
//	   0.220744304808428,
//	   0.037166522931547,
//	   0.007582288246685,
//	  -0.029626014314798,
//	  -0.000587146946925,
//	   0.003465279303277,
//	   0.001626050411134,
//	   0.000000000000000};

//I1 =     4.9926
//I2 =     0.0160
//I3 =    -4.4904


