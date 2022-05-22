#include <iostream>
#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <boost/numeric/odeint.hpp>
#include "CONSTANTS.hpp"

#include <boost/numeric/odeint/iterator/adaptive_iterator.hpp>
#include <boost/numeric/odeint/iterator/adaptive_time_iterator.hpp>

#include "dssCPP.h"
#define NDSS 2
#define INTEGRATOR 7
/*
 *  0: traditional integrate
 *  1: bulirsch_stoer;
 *  2: adams_bashforth_moulton [own initialised/initialised
 *  3: runge_kutta_fehlberg78 with observer;
 *  4: runge_kutta_fehlberg78 with local output;
 *  5: runge_kutta_fehlberg78 with single step with local output
 *  6: bulirsch_stoer : stepper
	   runge_kutta_fehlberg78: error_stepper_type
	   observer output
	7: controlled_adams_bashforth_moulton (not a good option)
*/

using namespace std;

double ua(double, double);
ostream& operator<<(ostream& os, const array_type& a);
void (*FuncPtr) (const array_type&, array_type&, const double);
ostream& operator<<(ostream& os, const matrix_type &);
ostream& operator<<(ostream& os, const matrix3D_type &);
double simp( double, double ,array_type );

struct output_observer_append
{
    string filename_;
    size_t count_;
    output_observer_append( const string &filename ) : filename_( filename ) , count_( 0 ) { }

    void operator()( const array_type &x , double t )
    {

    	double dr = r0/(NR-1);
    	array_type r2(NR, 0.0);
     	array_type x_(NR,0.0);
     	for (int j=0; j<NR; j++){
    		r2[j] = (j*dr)*(j*dr);
    	}

    	for (int i=0; i<NR; i++)
    	      x_[i] = 4.0*M_PI*r2[i]*x[i];

    	cout<<simp(0,1,x_)/((4.0/3.0)*M_PI*r2[NR-1])<<endl;
    	//cout<<simp(0,1,x_)/((4.0/3.0)*M_PI*r2[NR-1])-simp(0,1,x)<<endl;//cout<<x;
        ++count_;
    	std::ofstream fout(filename_, std::ios_base::app);
     	fout<<t<<",";
  		for (int i=0; i<NR; i++) {
    			fout <<x[i] << ",";
    		}
    	fout<<"\n";
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
    	//if (count_//SAMPLE_SIZE == 0) {
    		fout<<t<<",";
    		for( int i=0 ; i<2*NZ*NR; ++i ){
    			s_x = 0.0 + i*(1.0/50.0);
    			real_u = ua(s_x,t);
    			cout<<x[i]-real_u<<"\n";
    			fout<<x[i]-real_u << ",";
    		}
    		fout<<"\n";
    	//}
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

double ua(double x, double t) {

	return 0.0;

}

struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};

void difussion_pde_spherical_solution (const array_type &u , array_type &ut , const double t/* t */){

	// 1D to 2D matrices
    // Grid in axial direction
	array_type r(NR);
	array_type f(NR);
	array_type z(NZ);
	array_type ur(NR);
	array_type urr(NR);

	double fs;

	// Grid in radial direction
	double dr = r0/(NR-1);

	for (int j=0; j<NR; j++){
		r[j] = j*dr;
		f[j]=exp(-(r[j]*r[j]/STD*STD));
	}

	double drs = dr*dr;

	// Step through the grid points in r and z
	for (auto i=0; i<NR; i++) {

	//
		// ur
		if(i == 0) {
			//ur[i] = 0.0;
			ur[i] = 4.0*(u[i+1]-u[i])/drs;  // 2*urr
		}
		else if (i == NR-1) {
			//ur[i] = (u[i-1]-u[i])/dr;
			ur[i] = 0.0;
		}
		else {
		    ur[i] = (1.0/r[i])*(u[i+1]- u[i-1])/dr;
		}

		// urr
		if(i == 0) {
			urr[i]=2.0*(u[i+1]-u[i])/drs;
		}
		else if( i == NR-1) {
			urr[i]=2.0*(u[i-1]-u[i])/drs;
		}
		else {
			urr[i]=(u[i+1]-2.0*u[i]+u[i-1])/drs;
		}
		//

		// PDEs

	    if(t<=tau)
	      fs=f[i];
	    else
	      fs=0.0;

//	    This is another way of saying the in the case of r == 0, the application of L'hopital resuls as
//		two times the second derivative plus the term associated with the second derivative itself.
//		This resuls as 3 times the second derivative.
//	    if (i == 0)
//	    	ut[i] = 3*D*urr[i]  + fs ;
//	   else
	    	ut[i]=D*(urr[i] + ur[i]) + fs;

		}
}

double simp( double xl,double xu,array_type u) {
//%
//% Function simp computes an integral numerically by Simpson’s rule

	double uint;

	double h=(xu-xl)/(NR - 1);
	uint = u[0]-u[NR];
	for (int i=2; i<NR; i=i+2)
		uint  = uint + 4.0*u[i-1] + 2.0*u[i];

	uint = h/3.0*uint;

	return uint;

}

void difussion_pde_solution_highOrder (const array_type &u , array_type &ut , const double t/* t */){

	array_type ur(NR, 0.0);
	array_type urr(NR, 0.0);
	array_type f(NR, 0.0);
	array_type r(NR, 0.0);
	double fs;

	double dr = r0/(NR-1);

	for (int j=0; j<NR; j++) {
		r[j] = j*dr;
		f[j]=exp(-(r[j]*r[j]/STD*STD));
	}

	// Spatial grid in r
	for (int i=0; i<NR; i++) {
		// ur
		ur = dss04(0.0,r0,NR,u);
		//After this call, BCs (14.10a) and (14.10b) are imposed to ensure the correct (zero) values of ur at r = 0 and r0.
		ur[NR-1] = 0.0;
		ur[0] = 0.0;
		// urr
		urr=dss04(0.0,r0,NR,ur);
	}

	//
	// PDE
	for (int i=0; i<NR; i++) {
		if(t<=tau)
			fs=f[i];
		else
			fs=0.0;

		if(i == 0)
			ut[i]=D*3.0*urr[i]+fs;
		else
			ut[i]=D*(urr[i]+2.0/r[i]*ur[i])+fs;
	}

}

int main(int argc, char **argv) {
	std::cout << "Hello world\n";

    using namespace boost::numeric::odeint;

	array_type u(NR, 0.0);

    string fn= "/home/julio/temp/SIMPLE_DIFFUSION_SPHE"+to_string(NDSS)+"_integrator_"+to_string(INTEGRATOR) +".csv";

    std::ofstream out;
    out.open(fn); // append instead of overwrite
    out <<"t,";

	if (NDSS > 1)
		FuncPtr = difussion_pde_solution_highOrder;
	else
		FuncPtr = difussion_pde_spherical_solution;

	if (INTEGRATOR == 0) {

    //size_t steps = integrate( FuncPtr, u, 0.0 , 20.0, 0.1, streaming_observer( cout ));
    size_t steps = integrate( FuncPtr, u, 0.0 , 2.5,  0.1, output_observer_append( fn ));
    //size_t steps = integrate( FuncPtr, u, 0.0 , 20.0, 0.1, output_observer_append_error( fn ));
	cout << "!!!Hello World!!! with " <<steps<< endl; // prints !!!Hello World!!!

   //cout << "!!!Hello World!!! with " <<steps<< endl; // prints !!!Hello World!!!
	} else if (INTEGRATOR == 1) {
    typedef bulirsch_stoer< state_type > stepper_type;
    stepper_type stepper( 1E-9 , 1E-9 , 1.0 , 0.0 );
    double dt = 0.1;
//    stepper.adjust_size( u );
//    //stepper.try_step( FuncPtr, u , t , dt );
//    std::cout << "starting integration..." << std::endl;
//
	size_t steps = integrate_adaptive( stepper , FuncPtr , u , 0.0 , 2.5, dt, output_observer_append(fn));
	cout << "!!!Hello World!!! with " <<steps<< endl; // prints !!!Hello World!!!
//
	} else if (INTEGRATOR == 2) {

		int self = 1;

//    // adams_bashforth_moulton stepper example
//
        double t = 0.0 , dt = 0.01;

//        //[ multistep_detail_example
		for( double t=0.0 ; t<2.5 ; t+= dt ) {

			if (self == 1) {
				adams_bashforth_moulton< 5, state_type > abm;
				abm.initialize( FuncPtr , u , t , dt );
				abm.do_step( FuncPtr , u , t , dt);
				cout<<u<<endl;
			} else {
//        //[ multistep_detail_own_stepper_initialization
				adams_bashforth_moulton<2, state_type > abm;
				abm.initialize( runge_kutta_fehlberg78< state_type >() , FuncPtr , u , t , dt );
				abm.do_step( FuncPtr , u , t , dt);
				cout<<t<<","<<u<<endl;
			}
        }


	} else if (INTEGRATOR == 3) {
		runge_kutta_fehlberg78< state_type > stepper;
		typedef runge_kutta_fehlberg78< state_type > error_stepper_type;
		typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;

		//[integrate_adapt_full
		double abs_err = 1.0e-12 , rel_err = 1.0e-12 , a_x = 1.0 , a_dxdt = 1.0;
		controlled_stepper_type controlled_stepper(
			default_error_checker< double , range_algebra , default_operations >( abs_err , rel_err , a_x , a_dxdt ) );
		size_t steps = integrate_adaptive( controlled_stepper , FuncPtr, u, 0.0 , 2.50 , 0.001, output_observer_append (fn));
		cout << "!!!Hello World!!! with " <<steps<< endl; // prints !!!Hello World!!!

	} else if (INTEGRATOR == 4) {
		runge_kutta_fehlberg78< state_type > stepper;
		typedef runge_kutta_fehlberg78< state_type > error_stepper_type;
		typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
		vector<state_type> x_vec;

		vector<double> times;
		//[integrate_adapt_full
		double abs_err = 1.0e-10 , rel_err = 1.0e-6 , a_x = 1.0 , a_dxdt = 1.0;
		controlled_stepper_type controlled_stepper(
				default_error_checker< double , range_algebra , default_operations >( abs_err , rel_err , a_x , a_dxdt ) );
		size_t steps = integrate_adaptive( controlled_stepper , FuncPtr, u, 0.0 , 2.50 , 0.01, push_back_state_and_time( x_vec , times ));
		for( size_t i=0; i<=steps; i++ )
		     cout << times[i] << '\t' << x_vec[i][0] << '\t' << x_vec[i][1] << '\t' << x_vec[i][2] <<'\n';
		cout << "!!!Hello World!!! with " <<steps<< endl; // prints !!!Hello World!!!
	} else if (INTEGRATOR == 5) {

		double dr = r0/(NR-1);
		array_type r2(NR, 0.0);
		array_type x_(NR,0.0);
		runge_kutta_fehlberg78< state_type > stepper;

		for (int j=0; j<NR; j++){
			r2[j] = (j*dr)*(j*dr);
		}

		const double dt = 0.002;
		for( double t=0.0 ; t<2.5 ; t+= dt ) {
			stepper.do_step( FuncPtr , u , t , dt);
			for (int i=0; i<NR; i++)
				  x_[i] = 4.0*M_PI*r2[i]*u[i];

			cout<<simp(0,1,x_)/((4.0/3.0)*M_PI*r2[NR-1])<<endl;
		}
		size_t steps = int(1/dt);
		cout << "!!!Hello World!!! with " <<steps<< endl; // prints !!!Hello World!!!
	} else if (INTEGRATOR == 6) {
		typedef bulirsch_stoer< state_type > stepper_type;
		stepper_type stepper( 1E-9 , 1E-9 , 1.0 , 0.0 );
		typedef runge_kutta_fehlberg78< state_type > error_stepper_type;
		typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
		vector<state_type> x_vec;

		vector<double> times;
		//[integrate_adapt_full
		double abs_err = 1.0e-10 , rel_err = 1.0e-6 , a_x = 1.0 , a_dxdt = 1.0;
		controlled_stepper_type controlled_stepper(
				default_error_checker< double , range_algebra , default_operations >( abs_err , rel_err , a_x , a_dxdt ) );
		size_t steps = integrate_adaptive( controlled_stepper , FuncPtr, u, 0.0 , 2.50 , 0.01, output_observer_append (fn));
		cout << "!!!Hello World!!! with " <<steps<< endl; // prints !!!Hello World!!!
	}else if (INTEGRATOR == 7) {

        double t = 0.0 , dt = 0.01;

//        //[ multistep_detail_example
		for( double t=0.0 ; t<2.5 ; t+= dt ) {
			controlled_adams_bashforth_moulton<adaptive_adams_bashforth_moulton<5, state_type> > s9;
			adams_bashforth_moulton< 5, state_type > abm;
			s9.try_step( FuncPtr , u , t , dt);
			cout<<t<<","<<u;
        }
	}

	return 0;
}

ostream& operator<<(ostream& os, const array_type& a)
{
	int count=0;
	for (auto i : a) {
		cout<<i<<",";
	count++;
	}
	cout<<"\n";

    return os;
}

ostream& operator<<(ostream& os, const matrix_type & mat)
//void printMatrix(matrix_type vect )
{
	cout<<"Calling print function \n";

    for (vector<double> vect1D : mat)
    {
        for (double x : vect1D)
        {
            cout << x << " ";
        }
        cout << endl;
    }

    return os;

}


// From https://iq.opengenus.org/3d-vectors-in-cpp/
ostream& operator<<(ostream& os, const matrix3D_type & mat) {


	  for(int k=0;k<mat.size();k++)
	  {
	    for(int i=0;i<mat[k].size();i++)
	    {
	      for(int j=0;j<mat[k][i].size();j++)
	      {
	        cout<<mat[i][j][k]<<" ";
	      }
	      cout<<endl;
	    }
	    cout<<endl;
	  }
//
//
//	for(vector<vector<double>> i : mat) {
//
//		for(vector<double> j: i) {
//
//			for(double k: j) {
//				cout<<k<<" ";
//			}
//			cout<<endl;
//		}
//		cout<<endl;
//  	}

    return os;

}


