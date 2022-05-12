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
#define CASE 2
#define SIZE 301
#define SAMPLE_SIZE 100

using namespace std;
#include <boost/numeric/odeint.hpp>

double ua(double x, double t, double c);
ostream& operator<<(ostream& os, const array_type& a);
array_type operator*(const array_type& op1, const array_type& op2);
array_type simp(double xl, double xu, int n, array_type u);
void (*FuncPtr) (const array_type&, array_type&, const double);

struct output_observer_append
{
    string filename_;
    size_t count_;
    output_observer_append( const string &filename ) : filename_( filename ) , count_( 0 ) { }

    void operator()( const array_type &x , double t )
    {
        ++count_;
    	std::ofstream fout(filename_, std::ios_base::app);
    	 if (count_%SAMPLE_SIZE == 0) {
     		 cout<<simp(-30, 70, SIZE, x);
     		 fout<<t<<",";
    		 for( int i=0 ; i<SIZE ; ++i)
    			 fout <<x[i] << ",";
    		 fout<<x[SIZE-1]<<"\n";
    	 }
    }
};

struct output_observer_append_error
{
    string filename_;
    size_t count_;
    output_observer_append_error( const string &filename ) : filename_( filename ) , count_( 0 ) { }

    void operator()( const array_type &x , double t )
    {

        ++count_;
    	double real_u,s_x;
    	std::ofstream fout(filename_, std::ios_base::app);
    	if (count_%SAMPLE_SIZE == 0) {
    		cout<<simp(-30, 70, SIZE, x)<<"\n";
    		fout<<t<<",";
    		for( int i=0 ; i<SIZE ; ++i ){
    			s_x = -10 + i*(50.0/200.0);
    			real_u = ua(s_x,t, 1.0);
    			//cout<<x[i]-real_u<<"\n";
    			fout <<x[i]-real_u << ",";
    		}
    		fout<<x[SIZE-1]-ua(70,30, 1.0)<<"\n";
    	}

    }
};

struct streaming_observer
{
    std::ostream& m_out;
    mutable size_t count_;

    streaming_observer( std::ostream &out ) : m_out( out ) , count_( 0 ){ }

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
        ++(this->count_);
		m_out << t << ' '<<"\n";
		for (auto i : x)
			m_out<<t<<","<<i<<endl;
    }
};

void kortwert_solution(const array_type &u , array_type &ut , const double t/* t */){

	array_type ux = dss04(-30.0, 70.0, SIZE, u);
	array_type uxxx = dss007(-30.0, 70.0, SIZE, u);


	for (int i=0; i<SIZE; i++)
	    ut[i]=-uxxx[i] - 6.0*u[i]*ux[i];

}

double ua(double x, double t, double c) {

	double sech = 2.0/(exp(0.5*sqrt(c)*(x-c*t)) + exp(-0.5*sqrt(c)*(x-c*t)));
	double sech2  = sech*sech;
	return (0.5*c*sech2);
}

int main() {
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!

	using namespace std;
    using namespace boost::numeric::odeint;

	array_type x(SIZE);
	array_type u(SIZE);

    string fn= "/home/julio/temp/korteweg"+ to_string(CASE) +".csv";

    std::ofstream out;
    out.open(fn); // append instead of overwrite
    out <<"t,";
    for (int i=0; i<SIZE; i++)
    	out<<x[i]<<",";
    out<<"\n";
    out.close();

	double c, c1, c2;
    double xl = -30.0;
    double xu = 70.0;
    double dx = (xu-xl)/(SIZE-1);

    int ncase = 2;

    if (ncase == 1) {
    	c=1.0;
    } else {
    	c1=2.0;
    	c2=0.5;
    }

    double expm, expp,	pulse1,	pulse2;

	// Case 1 - Single pulse
	if(ncase==1) {
	    for (int i=0; i<SIZE; i++ ) {
	      x[i]=-30.0+i*dx;
	      u[i]=ua(x[i],0.0, 1.0);
	    }
	}
	// Case 2 - two pulses
	if(ncase==2){
		for (int i=0; i<SIZE; i++ ) {
			x[i]=-30.0+i*dx;
			expm=exp(-1.0/2.0*sqrt(c1)*(x[i]+15.0));
			expp=exp( 1.0/2.0*sqrt(c1)*(x[i]+15.0));
			pulse1=(1.0/2.0)*c1*(2.0/(expp+expm))*(2.0/(expp+expm));

			expm=exp(-1.0/2.0*sqrt(c2)*(x[i]-15.0));
			expp=exp( 1.0/2.0*sqrt(c2)*(x[i]-15.0));
			pulse2=(1.0/2.0)*c2*(2.0/(expp+expm))*(2.0/(expp+expm));

			u[i]=pulse1+pulse2;
		}
	}

    FuncPtr = kortwert_solution;




    size_t steps = integrate( FuncPtr, u, 0.0 , 1.0 , 0.1, streaming_observer( cout ));
    //size_t steps = integrate( FuncPtr, u, 0.0 , 30.0 , 0.1, output_observer_append( fn ));
    //size_t steps = integrate( FuncPtr, u, 0.0 , 30.0, 0.1, output_observer_append_error( fn ));

//    runge_kutta4< state_type > stepper1;
//    adams_bashforth< 8 , state_type > stepper;
//
//    const double dt = 0.0001;
//    array_type W(SIZE);
//	array_type w2v2(SIZE);
//	array_type V(SIZE);
//
//	int count = 0;
//	for( double t=0.0 ; t<30.0 ; t+= dt ) {
//        stepper.do_step( FuncPtr, u , t , dt);
//        auto start = u.begin();
//        auto end   = u.begin() + SIZE;
//        copy(start+SIZE, end+SIZE-1, W.begin());
//        copy(start, end, V.begin());
//        w2v2 = V*W;
//        if (count%10000 == 0){
//			for (auto i : w2v2){
//				out<<i<<",";
//			}
//			out<<"\n";
//        }
//		count++;
//	}




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



    cout << "!!!Hello World!!! with " <<steps<< endl; // prints !!!Hello World!!!

}
array_type simp(double xl, double xu, int n, array_type u) {

	// Function simp computes three integral invariants by Simpson's rule

	//	   / ∞                       i=n-2
	//    / u(x,t)dx  ≈ (h/3)*[ u  + SUM( 4u + 2u   ) + u ]
	// −∞/                       1    i=2   i    i+1     n

	array_type v(n), ux(n), ux2(n);
	array_type uabs(n);
	array_type uint(3,0);
	array_type uxabs(n);
	array_type INT(n);
	double h = (xu - xl)/(n-1);

	for (int k=1;k<4; k++) {
	// From https://en.wikipedia.org/wiki/Simpson%27s_rule
	// I1 conservation of mass


		if( k == 1) {
			//		 		 / ∞
			//		u1(t) = /  u(x, t) dx
			//			 −∞/

			uint[0  ] = u[0]+ u[n-1];

			for (int i = 1; i<n/2; i++)
					uint[0]  = uint[0] + 4.0*u[2*i-1];

			for (int j = 1; j<=n/2-1; j++)
				uint[0]  = uint[0] + 2.0*u[2*j];

			uint[0] = h/3.0*uint[0];
		}
		if( k == 2) {
			//		 		      / ∞
			//		u2(t) = 0.5* /  u(x, t)^2 dx
			//			      −∞/

			uint[1] = u[0]*u[0] + u[n-1]*u[n-1];

			for (int i = 1; i<n/2; i++)
					uint[1]  = uint[1] + 4.0*u[2*i-1]*u[2*i-1];

			for (int j = 1; j<=n/2-1; j++)
				uint[1]  = uint[1] + 2.0*u[2*j]*u[2*j];

			uint[1] = 0.5*h/3.0*uint[1];
		}
		if( k == 3) {
			//		 		 / ∞
			//		u3(t) = / 2*u(x, t)^3 - ux(x,t)^2 dx
			//			 −∞/
			ux=dss04(xl,xu,n,u);
			ux2 = ux*ux;
			uint[2] = 2.0*u[0]*u[0]*u[0]-ux2[0]+ 2.0*u[n-1]*u[n-1]*u[n-1]-ux2[0];

			for (int i = 1; i<n; i++)
				uint[2]  = uint[2] + 4.0*(2.0*u[2*i-1]*u[2*i-1]*u[2*i-1] - ux2[2*i-1]);

			for (int j = 1; j<=n/2-1; j++)
				uint[2]  = uint[2] + 2.0*(2.0*u[2*j]*u[2*j]*u[2*j] - ux2[2*j]);

			uint[2] = h/3.0*uint[2];
		}
	}

	return uint;
}

array_type operator*(const array_type& op1, const array_type& op2)
{
	array_type op3(SIZE);

	size_t p1 = op1.size();
	size_t p2 = op2.size();

	try{
		assert (p1 == p2);
		for (int k=0; k<SIZE; k++) {
			op3[k] = op1[k]*op2[k];
		}
	} catch (size_t p1) {
		cout<<p1;
		cout<<"\n";
	} catch (size_t p2) {
		cout<<p2;
		exit(0);
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

//		 		 / ∞
//		u1(t) = /  u(x, t) dx
//			 −∞/

	// From https://en.wikipedia.org/wiki/Simpson%27s_rule
// used only if n is multiple of three.
//if( k == 1) {
////		 		 / ∞
////		u1(t) = /  u(x, t) dx
////			 −∞/
//
//	uint[0  ] = u[0]+ u[n-1];
//	for (int i = 1; i<n; i++) {
//		if (i%3 != 0)
//			uint[0]  = uint[0] + 3.0*u[i];
//	}
//	for (int j = 1; j<=n/3-1; j++) {
//		uint[0]  = uint[0] + 2.0*u[3*j];
//	}
//	uint[0] = 3*h/8.0*uint[0];
//}
//if( k == 3) {
//		 		      / ∞
//		u1(t) = 0.5* /  u(x, t)^2 dx
//			      −∞/
//	ux=dss04(xl,xu,n,u);
//	ux2 = ux*ux;
//	uint[2] =           2.0*u[0  ]*u[0  ]*u[0  ] - ux2[0];
//	uint[2] = uint[2] + 2.0*u[n-1]*u[n-1]*u[n-1] - ux2[0];
//	for (int i = 1; i<n; i++) {
//		if (i%3 != 0)
//			uint[2]  = uint[2] + 3.0*(2.0*u[i]*u[i]*u[i] - ux2[i]);
//	}
//	for (int j = 1; j<=n/3-1; j++) {
//		uint[2]  = uint[2] + 2.0*(2.0*u[3*j]*u[3*j]*u[3*j] - ux2[3*j]);
//	}
//	uint[2] = 3.0*h/8.0*uint[2];
//}

/*
//array_type un = {0.0000, 0.0034, 0.0034, 0.0000, 0.0275, 0.0275, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000};
//
// un = {0.000,    0.000, 0.000,   -0.0045,0.0011,    0.0045,-0.0057,    0.0027,-0.0028,
// 0.0041,-0.0006,   -0.0058,0.0082,   -0.0043,0.0019,   -0.0034,0.0012,    0.0052,-0.0098,
// 0.0073,-0.0012,    0.0005, 0.0007,   -0.0047, 0.009,   -0.0081, 0.0022,    0.0026,-0.0046,
// 0.0062,-0.007,    0.0036, 0.0007,   -0.0022, 0.0026,   -0.0052, 0.0076,   -0.0032,-0.0039,
// 0.0065,-0.0026,    0.0012,-0.0052,    0.0044, 0.0039,   -0.0105, 0.0061,    0.0028,-0.0032,
//-0.0022, 0.0017,    0.006,-0.0092,    0.002, 0.0054,   -0.0039, -0.0008,
//-0.0011,    0.0073,   -0.0055,   -0.0035,    0.0098,   -0.0049,   -0.0028,    0.0012,    0.0037,
//-0.0031,   -0.005,    0.0091,   -0.002,   -0.0056,   0.0034,    0.0035,   -0.0023,   -0.0062,    0.0087,
// 0.0012,   -0.0101,    0.0056,    0.0052,   -0.0062,   -0.0032,   0.006,    0.0025,   -0.009,
// 0.0022,    0.0082,   -0.0057,   -0.0048,    0.0059,    0.0039,
//-0.0079,   -0.002,    0.0104,   -0.004,   -0.0068,    0.0054,   0.005,   -0.0051,   -0.0055,
// 0.0096,    0.0014,   -0.0111,    0.0043,    0.0065,   -0.0054,   -0.0066,    0.0082,    0.005,
//-0.0117,    0.0026,    0.0089,   -0.0055,   -0.006,    0.0053,    0.0068,   -0.0101,   -0.0013,
// 0.0107,   -0.0048,   -0.0048,    0.0023,    0.0062,   -0.0051,   -0.0054,    0.0092,   -0.0015,
//-0.0044,    0.0011,    0.0035,   -0.0007,   -0.0052,    0.0057,   -0.001,   -0.0018,    0.0022,
// 0.0002,    0.0009,   -0.0012,    0.0043,   -0.0028,    0.000,    0.0056,   -0.0015,    0.0001,
// 0.0034,    0.0074,   -0.0031,    0.0041,    0.0141,    0.003,    0.0071,    0.0135,    0.0206,    0.0123,    0.0207,    0.0391,    0.0318,    0.0437,    0.0586,    0.0725,    0.080,
// 0.0981,    0.1303,    0.1351,    0.1634,    0.1926,    0.207,    0.2242,    0.2394,    0.2559,
// 0.2418,    0.2427,    0.2416,    0.2084,    0.1949,    0.1743,    0.1538,    0.1242,    0.104,    0.0993,    0.0685,    0.0585,    0.0514,
// 0.0416,     0.0325,    0.0185,    0.0282,    0.0158,    0.0087,    0.0114,    0.0103,    0.0103,   -0.0044,    0.008,    0.0083,   -0.0036,
// 0.0011,    0.0039,    0.0079,   -0.0084,    0.0009,    0.0102,   -0.006,   -0.0032,    0.0031,
// 0.0066,   -0.0064,   -0.0041,    0.0114,   -0.0027,   -0.0054,    0.0028,    0.008,   -0.0006,
//-0.0027,    0.0149,    0.0137,    0.011,    0.0273,    0.0508,    0.0714,    0.1031,    0.1722,
// 0.267,    0.3778,    0.5347,    0.7295,    0.8866,    0.9782,    0.9953,    0.8879,    0.7078,    0.5414,    0.3871,    0.2561,
// 0.1617,    0.1179,    0.0754,    0.0328,    0.0319,    0.0208,    0.0135,   -0.0031,    0.0054,
// 0.0155,   -0.0115,    0.0021,    0.0040,    0.0049,   -0.0054,   -0.0071,    0.0173,   -0.0102,
//-0.004,    0.006,    0.0025,   -0.0014,   -0.0103,    0.0156,   -0.0044,   -0.0097,    0.0085,    0.0014,   -0.002,
//-0.0081,    0.0117,    0.0006,   -0.0123,    0.0082,    0.0048,   -0.0052,   -0.0064,    0.0092,    0.0019,
//-0.0115,    0.0041,    0.008,   -0.0059,   -0.0071,    0.010,    0.0011,   -0.0086,    0.0016,
// 0.008,   -0.0034,   -0.0092,    0.0111,    0.0001,   -0.0086,    0.0029,    0.0043,    0.0006,   -0.0111,    0.0104,
// 0.0025,   -0.0107,    0.0053,    0.000,    0.000,    0.000};
//
// simp(xl,xu,301,un);*/
