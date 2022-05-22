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
#define INTEGRATOR 6
#define TF 1.0
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
matrix_type operator*(double, const matrix_type&);

struct output_observer_append
{
    string filename_;
    size_t count_;
    output_observer_append( const string &filename ) : filename_( filename ) , count_( 0 ) { }

    void operator()( const array_type &x , double t )
    {

    	++count_;
    	std::ofstream fout(filename_, std::ios_base::app);
  		fout<<t<<",\n";
  		for (int i=0; i<NR; i++) {
  			for (auto j=0; j<NTH; j++) {
    			fout <<x[i*NTH + j] << ",";
    			cout <<x[i*NTH + j] << ",";
  			}
  	    	fout<<"\n";
  	    	cout<<"\n";
  		}
    	fout<<"\n";
    	cout<<"\n";
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
    		for( int i=0 ; i<2*NTH*NR; ++i ){
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
	array_type z(NTH);
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
////
//// Function simp computes an integral numerically by Simpson’s rule

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

array_type extract_column(matrix_type a, int pos) {
	array_type ret(NTH,0.0);

	for (int k=0; k<NTH; k++)
		ret[k] = a[k][pos];
	return (ret);

}
void difussion_3D_pde_solution_highOrder (const array_type &u , array_type &ut , const double t/* t */){

	array_type r(NR, 0.0);
	array_type th(NTH, 0.0);
	array_type y0(NTH,0.0);
	matrix_type f(NR,y0);
	matrix_type fs(NR,y0);
	matrix_type ut_(NR, y0);
	matrix_type u_(NR, y0);

	array_type u1d(NR,0.0);
	array_type ur1d(NR,0.0);
	array_type urr1d(NR,0.0);
	array_type uth1d(NTH,0.0);
	array_type uthth1d(NTH,0.0);
	matrix_type urr(NR, y0);
	matrix_type ur(NR, y0);
	matrix_type uth(NTH, y0);
	matrix_type uthth(NTH, y0);

	// Radial grid
	double dr = r0/(NR-1);
	double dth=th0/(NTH-1);

	for (int i=0; i<NR; i++)
		r[i]=i*dr;

	// Angular grid
	for (int i=0; i<NTH; i++)
		th[i]=i*dth;

	//
	// Inhomogeneous term
	for (int i=0; i<NR; i++) {
		for (int j=0; j<NTH; j++) {
			double a = r[i]*sin(th[j])/STD_pi2;
			double b = r[i]*cos(th[j])/STD_0;
			f[i][j]=exp(-(a*a) -(b*b));

		}
	}

//	//
//	// 1D to 2D
	for (int i=0; i<NR; i++) {
		for (int j=0; j<NTH; j++) {
			u_[i][j]=u[i*NTH+j];
		}
	}

//	//
//	// Spatial grids in r and theta
	for (int i=0; i<NR; i++) {
		for (int j=0; j<NTH; j++) {

	//	// ur
			u1d = extract_column(u_,j);
			// u1d=u(:,j);
			ur1d=dss04(0.0,r0,NR,u1d);
			ur1d[NR-1]=0.0;
			ur1d[0 ]=0.0;

			for (int m=0; m<NR; m++)
				ur[m][j] = ur1d[m];


			// urr
			urr1d=dss04(0.0,r0,NR,ur1d);


			for (int m=0; m<NR; m++)
				urr[m][j] = urr1d[m];
			//
			// uth


			u1d = u_[i];


			uth1d=dss04(0.0,th0,NTH,u1d);
			uth1d[NTH-1] = 0.0;
			uth1d[0]     = 0.0;
			uth[i] = uth1d;
			//
			// uthth
			uthth1d=dss04(0.0,th0,NTH,uth1d);
			uthth[i] = uthth1d;

		}
	}

	// Inhomogeneous term
	if (t <= tau)
		fs=f;
	else
	    fs=0.0*f;

	// PDE
	for (int i=0; i<NR; i++) {
		for (int j=0; j<NTH; j++) {
			//   r=0
			if (i == 0)
				// PDE for r = 0
				ut_[i][j]=D*(urr[i][j]+2.0*urr[i][j])+fs[i][j];
			// r ~= 0
			else {  //if(i~=1)
				// PDEs for r ~= 0
				if (j == 0 || j == NTH-1)
					// th = 0, pi/2
					ut_[i][j]=D*(urr[i][j]+(2.0/r[i])*ur[i][j] + (2.0/(r[i]*r[i]))*uthth[i][j])+fs[i][j];
				else
					// th ~= 0, pi/2
					ut_[i][j]=D*(urr[i][j]+(2.0/r[i])*ur[i][j] + (1.0/(r[i]*r[i]))*(uthth[i][j]+(cos(th[j])/sin(th[j]))*uth[i][j]))+fs[i][j];
			}
		}
	}
	  //
	  // 2D to 1D
	for (int i=0; i<NR; i++) {
		for (int j=0; j<NTH; j++) {
			ut[i*NTH + j]=ut_[i][j];
		}
	}
}

//
int main(int argc, char **argv) {
	std::cout << "Hello world\n";

    using namespace boost::numeric::odeint;

    // Initial condition
	array_type u(NR*NTH, 0.0);
	array_type y0(NTH,0.0);
	matrix_type u_c(NR,y0);

//	for (int k=0; k<NR*NTH; k++)
//		u[k] = k;

    string fn= "/home/julio/temp/SIMPLE_DIFFUSION_SPHE_3D_"+to_string(NDSS)+"_integrator_"+to_string(INTEGRATOR) +".csv";

    std::ofstream out;
    out.open(fn); // append instead of overwrite
    out <<"t,";

    FuncPtr = difussion_3D_pde_solution_highOrder;

	if (INTEGRATOR == 0) {

    //size_t steps = integrate( FuncPtr, u, 0.0 , 20.0, 0.1, streaming_observer( cout ));
    size_t steps = integrate( FuncPtr, u, 0.0 , 1.0,  0.1, output_observer_append( fn ));
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

		runge_kutta_fehlberg78< state_type > stepper;

		double eps = 1e-8;
		const double dt = 0.01;
		for( double t=0.0 ; t<1.0 ; t+= dt ) {
			stepper.do_step( FuncPtr , u , t , dt);

			if ( (t> 0 && t < dt+eps) || (t> 0.25 && t < 0.25+eps)|| (t> 0.5 && t < 0.5+eps)
					|| (t> 0.75 && t < 0.75+eps)|| (t> 0.99 && t < 1.99+eps) || t >= 1.0 ) {
				cout<<t<<endl;
				for (int i=0; i<NR; i++) {
					for (auto j=0; j<NTH; j++) {
						cout <<u[i*NTH + j] << ",";
					}
					cout<<"\n";
				}
				cout<<"\n";
			}
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
		size_t steps = integrate_adaptive( controlled_stepper , FuncPtr, u, 0.0 , TF , 0.01, output_observer_append (fn));
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

matrix_type operator*(double value, const matrix_type& rhs) {

	matrix_type result = rhs;

	for(int i=0;i<rhs.size();i++) {
		transform(result[i].begin(),
				result[i].end(),
				result[i].begin(),
				[value](auto &c){ return c*value; });
	}

	//cout<<result;

	return (result);
}
