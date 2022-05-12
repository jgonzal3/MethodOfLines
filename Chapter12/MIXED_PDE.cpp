#include <iostream>
#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <boost/numeric/odeint.hpp>

#include "dssCPP.h"
#define CASE 4
#define SIZE 50
#define SAMPLE_SIZE 1

#define X_L 0.0
#define Y_L 0.0
#define Z_L 0.0

#define X_U 1.0
#define Y_U 1.0
#define Z_U 1.0

#define NDSS 44

using namespace std;

double ua(double, double);
ostream& operator<<(ostream& os, const array_type& a);
array_type operator*(const array_type& op1, const array_type& op2);
void (*FuncPtr) (const array_type&, array_type&, const double);
ostream& operator<<(ostream& os, const matrix_type &);
ostream& operator<<(ostream& os, const matrix3D_type &);

struct output_observer_append
{
    string filename_;
    size_t count_;
    output_observer_append( const string &filename ) : filename_( filename ) , count_( 0 ) { }

    void operator()( const array_type &x , double t )
    {
        ++count_;
    	std::ofstream fout(filename_, std::ios_base::app);
     	fout<<t<<",";
    	for( int i=0 ; i<SIZE ; ++i) {
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
    	//if (count_%SAMPLE_SIZE == 0) {
    		fout<<t<<",";
    		for( int i=0 ; i<SIZE; ++i ){
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

	double L = (X_U - X_L);
	double L2 = L*L;
	double PI2 = M_PI*M_PI;

    double a=1.0/(1.0 - PI2/L2);
    return (sin(M_PI*x/L)*exp(a*t));

}

void mixed_pde_solution_CM (const array_type &u , array_type &ut , const double t/* t */){

	double dx=(X_U-X_L)/(SIZE);
	double dx2 = dx*dx;

	array_type temp = u;
	// Convert u() from the ODE solver in a matrix
    Eigen::MatrixXd b = Eigen::Map<Eigen::Matrix<double, SIZE, 1> >(temp.data());

	Eigen::MatrixXd cm = Eigen::MatrixXd::Zero(SIZE, SIZE);
    Eigen::MatrixXd x;


	// Create the coefficient matrix in Eigen

	for (int i=0; i<SIZE; i++) {
		for (int j=0; j<SIZE; j++) {
			if(i == j) {
				cm(i,j) = (1.0-2.0/dx2);
			}
			else if (abs(i-j) == 1)
			{
				cm(i,j) = 1.0/dx2;
			}
		}
	}

//  Solve the system of equations
	if (fabs(cm.determinant()) < 1e-9) {
		cout<<cm.determinant();
		cout<<"Exiting as determinant of A is less than eps: ";
		exit(0);
	}
	else {
    	x = cm.colPivHouseholderQr().solve(b);
		//x = cm.fullPivLu().solve(b);
	}

	// If you need to cast one matrix type into another, use the cast operator
	// cast a matrix of floats with doubles vector
	// Eigen::MatrixXf f = d.cast <float> ();
	array_type vec(x.data(), x.data() + x.rows() * x.cols());
	ut = vec;

}

int main(int argc, char **argv) {
	std::cout << "Hello world\n";

    using namespace boost::numeric::odeint;

	array_type x(SIZE);
	array_type u(SIZE);
	//array_type ut(SIZE);
	array_type y0(SIZE,0.0);
	matrix_type cm(SIZE, y0);


    string fn= "/home/julio/temp/Mixed_PDE_"+ to_string(CASE) + "_"+to_string(NDSS) +".csv";

    std::ofstream out;
    out.open(fn); // append instead of overwrite
    out <<"t,";

	double dx=(X_U-X_L)/(SIZE);
	double L = (X_U - X_L);

	// IC over the spatial grid
	for (int i=0; i<SIZE; i++) {
		// Initialise the spatial grid as uniform grid
		x[i]=X_L + i*dx;
		// Initial condition
		u[i] = sin(M_PI*x[i]/L);
		out<<x[i]<<",";
	}

    out<<"\n";
    out.close();



/*	Equation (12.1) is first order in t and second order in x (through
 *                    3    2
 * the mixed partial ∂ u/∂x ∂t). It therefore requires one initial condition
 * (IC) and two boundary conditions (BCs). The IC is taken as
 *
	u(x, t = 0) = sin(πx/L) (12.2)

	and the two BCs as

	u(x = 0, t) = u(x = L, t) = 0

	*/



    //
	FuncPtr = mixed_pde_solution_CM;
	//mixed_pde_solution_CM (u, ut, 1.0);

	//cout<<ut;

    //size_t steps = integrate( FuncPtr, y, 0.0 , 20.0, 0.1, streaming_observer( cout ));
    size_t steps = integrate( FuncPtr, u, 0.0 , 20.0 , 0.1, output_observer_append( fn ));
    //size_t steps = integrate( FuncPtr, u, 0.0 , 20.0, 0.1, output_observer_append_error( fn ));

   //cout << "!!!Hello World!!! with " <<steps<< endl; // prints !!!Hello World!!!

	return 0;
}

array_type operator*(const array_type& op1, const array_type& op2)
{
	array_type op3(SIZE);

	size_t p1 = op1.size();
	size_t p2 = op2.size();
	array_type y0(SIZE,0.0);

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

