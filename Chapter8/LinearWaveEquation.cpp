#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <boost/numeric/odeint.hpp>

#include "dssCPP.h"
#define CASE 4
#define SIZE 401
#define SAMPLE_SIZE 1
#define X_L 0.0
#define X_U 100.0
#define NDSS -1
#define N_L 2
#define N_U 2

using namespace std;

double ua(double x, double t, double c);
ostream& operator<<(ostream& os, const array_type& a);
array_type operator*(const array_type& op1, const array_type& op2);
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
    	 //if (count_%SAMPLE_SIZE == 0) {
     		 fout<<t<<",";
    		 for( int i=0 ; i<SIZE-1 ; ++i)
    			 fout <<x[i] << ",";
    		 fout<<x[SIZE-1]<<"\n";
    	// }
    }
};

struct output_observer_append_error
{
    string filename_;
    size_t count_;
    output_observer_append_error( const string &filename ) : filename_( filename ) , count_( 0 ) { }

    void operator()( const array_type &x , double t ) {
    	++count_;
    	double real_u,s_x;
    	std::ofstream fout(filename_, std::ios_base::app);
    	if (count_%SAMPLE_SIZE == 0) {
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


void linear_solution(const array_type &u , array_type &ut , const double t/* t */){

	array_type u1(SIZE);  // array that represents u
	array_type u2(SIZE);  // array that represents ut
	array_type u1xx(SIZE);  // array that represents u
	array_type u1t(SIZE);  // array that represents u
	array_type u2t(SIZE);  // array that represents ut

	double dx = (X_U-X_L)/(SIZE-1);
	double dx2 = dx*dx;
	double cs = 1.0;

	for (int k=0; k<SIZE; k++) {
		u1[k] = u[k];
		u2[k] = u[k+SIZE];
	}

	// Calculate u1xx with ux(1) = ux(n) = 0 BCs
	  for (int i=0; i<SIZE; i++) {
	     if (i == 0)
	    	 u1xx[i]=2.0*(u1[i+1]-u1[i])/dx2;
	     else if (i == SIZE-1)
	    	 u1xx[i]=2.0*(u1[i-1]-u1[i])/dx2;
	     else
	    	 u1xx[i]=(u1[i+1]-2.0*u1[i]+u1[i-1])/dx2;
	  }

	// PDE
	  for (int i=0; i<SIZE; i++) {
	    u1t[i]=u2[i];
	    u2t[i]=cs*u1xx[i];
	  }

//	  Two vectors to one vector
	  for (int i=0; i<SIZE; i++) {
	    ut[i]  =u1t[i];
	    ut[i+SIZE]=u2t[i];
	  }
}

array_type ua(array_type x, double t, double c) {

	double a=2.0;
	array_type uanal(SIZE);
	array_type pulse1(SIZE), pulse2(SIZE);

	switch (CASE){

	case 1: //Rectangular pulse

		for (int i=0; i<SIZE; i++) {
			if (x[i]<(45.0-c*t))
				pulse1[i]=0.0;
			else if (x[i]>(55.0-c*t))
				pulse1[i]=0.0;
			else if (x[i]>=(45.0-c*t) && x[i]<=(50.0-c*t))
				pulse1[i]=0.5;
			else if (x[i]>=(50.0-c*t) && x[i]<=(55.0-c*t))
				pulse1[i]=0.5;

			if (x[i]<(45.0+c*t))
				pulse2[i]=0.0;
			else if (x[i]>(55.0+c*t))
				pulse2[i]=0.0;
			else if (x[i]>=(45.0+c*t) && x[i]<=(50.0+c*t))
				pulse2[i]=0.5;
			else if (x[i]>=(50.0+c*t) && x[i]<=(55.0+c*t))
				pulse2[i]=0.5;

			uanal[i]=pulse1[i]+pulse2[i];
		}
		break;
	case 2:  //Triangular pulse

		for (int i=0; i<SIZE; i++) {
			if (x[i]<(45.0-c*t))
				pulse1[i]=0.0;
			else if (x[i]>(55.0-c*t))
				pulse1[i]=0.0;
			else if (x[i]>=(45.0-c*t) && x[i]<=(50.0-c*t))
				pulse1[i] = (x[i]-(45.0-c*t))/10.0;
			else if (x[i]>=(50.0-c*t) && x[i]<=(55.0-c*t))
				pulse1[i]=((55.0-c*t)-x[i])/10.0;

			if (x[i]<(45.0+c*t))
				pulse2[i]=0.0;
			else if (x[i]>(55.0+c*t))
				pulse2[i]=0.0;
			else if (x[i]>=(45.0+c*t) && x[i]<=(50.0+c*t))
				pulse2[i]=(x[i]-(45.0+c*t))/10.0;
			else if (x[i]>=(50.0+c*t) && x[i]<=(55.0+c*t))
				pulse2[i]=((55.0+c*t)-x[i])/10.0;

			uanal[i]=pulse1[i]+pulse2[i];
		}
		break;
	case 3: // Sine pulse

		for (int i=0; i<SIZE; i++) {
			if (x[i]<(45.0-c*t))
				pulse1[i]=0.0;
			else if (x[i]>(55.0-c*t))
				pulse1[i]=0.0;
			el0se if (x[i]>=(45.0-c*t) && x[i]<=(50.0-c*t))
				pulse1[i]=(1.0/2.0)*sin(M_PI*(x[i]-(45.0-c*t))/10.0);
			else if (x[i]>=(50.0-c*t) && x[i]<=(55.0-c*t))
				pulse1[i]=(1.0/2.0)*sin(M_PI*((55.0-c*t)-x[i])/10.0);

			if (x[i]<(45.0+c*t))
				pulse2[i]=0.0;
			else if (x[i]>(55.0+c*t))
				pulse2[i]=0.0;
			else if (x[i]>=(45.0+c*t) && x[i]<=(50.0+c*t))
				pulse2[i]=(1.0/2.0)*sin(M_PI*(x[i]-(45.0+c*t))/10.0);
			else if (x[i]>=(50.0+c*t) && x[i]<=(55.0+c*t))
				pulse2[i]=(1.0/2.0)*sin(M_PI*((55.0+c*t)-x[i])/10.0);

			uanal[i]=pulse1[i]+pulse2[i];
		}
	break;
	case 4: // Gausian pulse
		a=0.2;
		for (int i=0; i<SIZE; i++) {
			if (x[i]<(45.0-c*t))
				pulse1[i]=0.0;
			else if (x[i] > (55.0-c*t))
				pulse1[i]=0.0;
			else if (x[i]>=(45.0-c*t) && x[i]<=(55.0-c*t))
				pulse1[i]=(1.0/2.0)*exp(-a*(x[i]-(50.0-c*t))*(x[i]-(50.0-c*t))/10.0);

			if (x[i]<(45.0+c*t))
				pulse2[i]=0.0;
			else if (x[i]>(55.0+c*t))
				pulse2[i]=0.0;
			else if (x[i]>=(45.0+c*t) && x[i]<=(55.0+c*t))
				pulse2[i]=(1.0/2.0)*exp(-a*(x[i]-(50.0+c*t))*(x[i]-(50.0+c*t))/10.0);

			uanal[i]=pulse1[i]+pulse2[i];
		}
		break;
	default:
		break;
	}
	return (uanal);
}

void linear_solution_highOrder(const array_type &u , array_type &ut , const double t/* t */){
	array_type u1(SIZE);  // array that represents u
	array_type u2(SIZE);  // array that represents ut
	array_type u1xx(SIZE);  // array that represents u
	array_type u1x(SIZE);  // array that represents u
	array_type u1t(SIZE);  // array that represents u
	array_type u2t(SIZE);  // array that represents ut

	double dx = (X_U-X_L)/(SIZE-1);
	double dx2 = dx*dx;
	double cs = 1.0;

	for (int k=0; k<SIZE; k++) {
		u1[k] = u[k];
		u2[k] = u[k+SIZE];
	}
	if (NDSS== 2)
		u1x=dss02(X_L,X_U,SIZE,u1); // second order
	else if (NDSS== 4)
		u1x=dss04(X_L,X_U,SIZE,u1); // fourth order
	else if(NDSS== 6)
		u1x=dss06(X_L,X_U,SIZE,u1); // sixth order
	else if(NDSS== 8)
		u1x=dss08(X_L,X_U,SIZE,u1); // eighth order
	else if(NDSS==10)
		u1x=dss10(X_L,X_U,SIZE,u1); // tenth order


 	// BC at x = 0
	// The condition is not given but will be enforced because u(0,t)=0
	u1x[0] = 0.0;

	// BC at x = 1
	u1x[SIZE-1] = 0.0;

	// Calculate u1xx
	if (NDSS== 2) u1xx=dss02(X_L,X_U,SIZE,u1x); // second order
	else if(NDSS== 4)
		u1xx=dss04(X_L,X_U,SIZE,u1x); // fourth order
	else if(NDSS== 6)
		u1xx=dss06(X_L,X_U,SIZE,u1x); // sixth order
	else if(NDSS== 8)
		u1xx=dss08(X_L,X_U,SIZE,u1x); // eighth order
	else if(NDSS==10)
		u1xx=dss10(X_L,X_U,SIZE,u1x); // tenth order

	// PDE

	for (int i=0; i<SIZE; i++) {
		u1t[i] = u2[i];
		u2t[i] = cs*u1xx[i];
	}

//	  Two vectors to one vector
	for (int i=0; i<SIZE; i++) {
		ut[i]  = u1t[i];
		ut[i+SIZE] = u2t[i];
	}

}

void linear_solution_highOrder_second(const array_type &u , array_type &ut , const double t/* t */){
	array_type u1(SIZE);  // array that represents u
	array_type u2(SIZE);  // array that represents ut
	array_type u1xx(SIZE);  // array that represents u
	array_type u1x(SIZE,0);  // array that represents u
	array_type u1t(SIZE);  // array that represents u
	array_type u2t(SIZE);  // array that represents ut

	double dx = (X_U-X_L)/(SIZE-1);
	double dx2 = dx*dx;
	double cs = 1.0;

	for (int k=0; k<SIZE; k++) {
		u1[k] = u[k];
		u2[k] = u[k+SIZE];
	}

 	// BC at x = 0
	// The condition is not given but will be enforced because u(0,t)=0
	u1x[0] = 0.0;

	// BC at x = 1
	u1x[SIZE-1] = 0.0;

	if (NDSS==42)
		u1xx=dss42(X_L,X_U,SIZE,u1,u1x,N_L,N_U);
	// second order
	else if(NDSS==44)
		u1xx=dss44(X_L,X_U,SIZE,u1,u1x,N_L,N_U);
	// fourth order
	else if(NDSS==46)
		u1xx=dss46(X_L,X_U,SIZE,u1,u1x,N_L,N_U);
	// sixth order
	else if(NDSS==48)
		u1xx=dss48(X_L,X_U,SIZE,u1,u1x,N_L,N_U);
	// eighth order
	else if(NDSS==50)
		u1xx=dss50(X_L,X_U,SIZE,u1,u1x,N_L,N_U);
	// tenth order

	// PDE
	for (int i=0; i<SIZE; i++) {
		u1t[i] = u2[i];
		u2t[i] = cs*u1xx[i];
	}

//	  Two vectors to one vector
	for (int i=0; i<SIZE; i++) {
		ut[i]  = u1t[i];
		ut[i+SIZE] = u2t[i];
	}

}

int main(int argc, char **argv) {
	std::cout << "Hello world\n";

    using namespace boost::numeric::odeint;

	array_type x(SIZE);
	array_type u(2*SIZE,0);

    string fn= "/home/julio/temp/linearEquation"+ to_string(CASE) + "_"+to_string(NDSS) +".csv";

    std::ofstream out;
    out.open(fn); // append instead of overwrite
    out <<"t,";

	double dx=(X_U-X_L)/(SIZE - 1);

	// Initialise the spatial grid
	for (int i=0; i<SIZE; i++) {
		x[i]=i*dx;
		out<<x[i]<<",";
	}

    out<<"\n,";
	for (int i=0; i<SIZE-1; i++) {
		out<<i+1<<",";
	}
	out<<"\n";
    out.close();

	// Initialise the solution

	array_type u1=ua(x, 0, 1.0);

	for (int i=0; i<SIZE; i++) {
	    u[i]      = u1[i];
	    u[i+SIZE] = 0.0;
	}

    //
//	if (NDSS < 11)
//		FuncPtr = linear_solution_highOrder;
//	else if (NDSS > 40)
//		FuncPtr =  linear_solution_highOrder_second;
//	else if (NDSS == -1)
		FuncPtr = linear_solution;

    //size_t steps = integrate( FuncPtr, u, 0.0 , 30.0 , 0.1, streaming_observer( cout ));
    size_t steps = integrate( FuncPtr, u, 0.0 , 30.0 , 0.1, output_observer_append( fn ));
    //size_t steps = integrate( FuncPtr, u, 0.0 , 30.0, 0.1, output_observer_append_error( fn ));

    cout << "!!!Hello World!!! with " <<steps<< endl; // prints !!!Hello World!!!

	return 0;
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
	for (auto i : a) {
		cout<<i<<"\n";
	count++;
	}
    return os;
}

