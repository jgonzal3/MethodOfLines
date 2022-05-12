#include <iostream>
#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <boost/numeric/odeint.hpp>

#include "dssCPP.h"
#define CASE 4
#define SIZE 11
#define SAMPLE_SIZE 1
#define X_L 0.0
#define Y_L 0.0
#define Y_U 1.0
#define X_U 1.0
#define NDSS 44

using namespace std;

double ua(double, double, int);
ostream& operator<<(ostream& os, const array_type& a);
array_type operator*(const array_type& op1, const array_type& op2);
void (*FuncPtr) (const array_type&, array_type&, const double);
ostream& operator<<(ostream& os, const matrix_type &);

struct output_observer_append
{
    string filename_;
    size_t count_;
    output_observer_append( const string &filename ) : filename_( filename ) , count_( 0 ) { }

    void operator()( const array_type &x , double t )
    {
        ++count_;
    	std::ofstream fout(filename_, std::ios_base::app);
    	for( int i=0 ; i<SIZE ; ++i) {
     	fout<<t<<",";
    		for( int i=0 ; i<SIZE-1 ; ++i) {
    			fout <<x[i] << ",";
    		    }
			fout<<x[SIZE-1]<<"\n";
    	}

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

array_type extract_column(matrix_type a, int pos) {
	array_type ret(SIZE,0.0);

	for (int k=0; k<SIZE; k++)
		ret[k] = a[k][pos];
	return (ret);

}

double ua(double x, double t, int a) {

	return 0.0;

}

void laplace_solution(const array_type &u , array_type &ut , const double t/* t */){

	/* The routines that produces good results are DSS2 and DSS4.
	 * DSS6 and DSS8 are fine at the initial but as they approach 1, not fine
	 * DSS10 not find at all
	 */

	// Initially zero derivatives in x, y, t
	array_type y0(SIZE,0.0);
	array_type u1d(SIZE,0.0);
	array_type ux1d(SIZE,0.0);
	array_type uxx1d(SIZE,0.0);
	array_type uy1d(SIZE,0.0);
	array_type uyy1d(SIZE,0.0);
	matrix_type uxx(SIZE, y0);
	matrix_type uyy(SIZE, y0);
	matrix_type u_(SIZE, y0);
	matrix_type ut_(SIZE, y0);

	//
	// 1D to 2D matrix conversion
	for (int i = 0; i<SIZE; i++) {
		for (int j = 0; j<SIZE; j++) {
			u_[i][j] = u[i*SIZE+j];
		}
	}


	for (int j=0; j<SIZE; j++) {
		// ux

		u1d = extract_column(u_,j);
		if (NDSS == 2) 	ux1d=dss02(X_L,X_U,SIZE,u1d); // second order
		if (NDSS == 4) 	ux1d=dss04(X_L,X_U,SIZE,u1d); // fourth order
		if (NDSS == 6)  ux1d=dss06(X_L,X_U,SIZE,u1d); // sixth order
		if (NDSS == 8)  ux1d=dss08(X_L,X_U,SIZE,u1d); // eighth order
		if (NDSS == 10)	ux1d=dss10(X_L,X_U,SIZE,u1d); // tenth order
		u1d = extract_column(u_,j);

		/*
		For the BCs in x, we use the homogeneous Neumann conditions
		∂u(x = 0, y, t)
		--------------- = 0.0     (10.4)
		      ∂x
		*/

		u1d[0] = 0.0;
		/*
		∂u(x = 1, y, t)
		--------------- = 0     (10.5)
		       ∂x
		*/

		ux1d[SIZE-1] = 0.0;
	//
	//	// Calculate uxx1d
		if (NDSS == 2) 	uxx1d=dss02(X_L,X_U,SIZE,ux1d); // second order
		if (NDSS == 4)	uxx1d=dss04(X_L,X_U,SIZE,ux1d); // fourth order
		if (NDSS == 6)	uxx1d=dss06(X_L,X_U,SIZE,ux1d); // sixth order
		if (NDSS == 8)	uxx1d=dss08(X_L,X_U,SIZE,ux1d); // eighth order
		if (NDSS == 10)	uxx1d=dss10(X_L,X_U,SIZE,ux1d); // tenth order

		for (int m=0; m<SIZE; m++)
			uxx[j][m] = uxx1d[m];


	}

	for (int i=0; i<SIZE; i++) {
			// uy

		/*
		 * For the BCs in y, we use the Dirichlet conditions
		 *
		 *  u(x, y = 0, t) = 0 		(10.6)
		 *  u(x, y = 1, t) = 1 		(10.7)
		 *
		 */


		u1d = u_[i];
		if (NDSS == 2) 	uy1d=dss02(X_L,X_U,SIZE,u1d); // second order
		if (NDSS == 4) 	uy1d=dss04(X_L,X_U,SIZE,u1d); // fourth order
		if (NDSS == 6)  uy1d=dss06(X_L,X_U,SIZE,u1d); // sixth order
		if (NDSS == 8)  uy1d=dss08(X_L,X_U,SIZE,u1d); // eighth order
		if (NDSS == 10)	uy1d=dss10(X_L,X_U,SIZE,u1d); // tenth order

		// Calculate uyy1d

		if (NDSS == 2) 	uyy1d=dss02(X_L,X_U,SIZE,uy1d); // second order
		if (NDSS == 4)	uyy1d=dss04(X_L,X_U,SIZE,uy1d); // fourth order
		if (NDSS == 6)	uyy1d=dss06(X_L,X_U,SIZE,uy1d); // sixth order
		if (NDSS == 8)	uyy1d=dss08(X_L,X_U,SIZE,uy1d); // eighth order
		if (NDSS == 10)	uyy1d=dss10(X_L,X_U,SIZE,uy1d); // tenth order

		uyy[i] = uyy1d;

	}

	for (int i = 0; i<SIZE; i++) {
		for (int j = 0; j<SIZE; j++) {
			ut_[i][j] = uxx[i][j] + uyy[i][j];
		}
	}


	for (int m=0; m<SIZE; m++) {
		ut_[m][0] = 0.0;
		ut_[m][SIZE-1] = 0.0;
	}

	//
	// 2D to 1D matrix conversion
	for (int i=0; i<SIZE; i++) {
		for (int j=0; j<SIZE; j++){
			ut[i*SIZE+j]=ut_[i][j];
		}
	}

}

void laplace_solution_highOrder(const array_type &u , array_type &ut , const double t/* t */){

	int n_l, n_u;

	array_type y0(SIZE,0.0);
	array_type u1d(SIZE,0.0);
	array_type ux1d(SIZE,0.0);
	array_type uxx1d(SIZE,0.0);
	array_type uy1d(SIZE,0.0);
	array_type uyy1d(SIZE,0.0);
	matrix_type uxx(SIZE, y0);
	matrix_type uyy(SIZE, y0);
	matrix_type u_(SIZE, y0);
	matrix_type ut_(SIZE, y0);
	//
	// 1D to 2D matrix conversion
	for (int i = 0; i<SIZE; i++) {
		for (int j = 0; j<SIZE; j++) {
			u_[i][j] = u[i*SIZE+j];
		}
	}



	for (int j=0; j<SIZE; j++) {
		u1d = extract_column(u_,j);


		// Since BCs (10.4) and (10.5) are Neumann,
		// nl = nu = 2 are inputs to dssXX corresponding to x = 0, 1, respectively.
		n_l = 2;
		n_u = 2;

		// BC at x = 0
		// The condition is not given but will be enforced because u(0,t)=0
		ux1d[0] = 0.0;

		// BC at x = 1
		ux1d[SIZE-1] = 0.0;

		if (NDSS == 42)	uxx1d=dss42(X_L,X_U,SIZE,u1d,ux1d,n_l,n_u);
		// second order
		if(NDSS == 44)	uxx1d=dss44(X_L,X_U,SIZE,u1d,ux1d,n_l,n_u);
		// fourth order
		if(NDSS == 46)	uxx1d=dss46(X_L,X_U,SIZE,u1d,ux1d,n_l,n_u);
		// sixth order
		if(NDSS == 48)	uxx1d=dss48(X_L,X_U,SIZE,u1d,ux1d,n_l,n_u);
		// eighth order
		if(NDSS == 50)	uxx1d=dss50(X_L,X_U,SIZE,u1d,ux1d,n_l,n_u);
		// tenth order

		for (int m=0; m<SIZE; m++)
			uxx[j][m] = uxx1d[m];

	}

	for (int i=0; i<SIZE; i++) {
		u1d = u_[i];
		// Since BCs (10.4) and (10.5) are Dirichlet, then at y = 0, 1 nl = nu = 1.

		n_l = 1;
		n_u = 1;
		if (NDSS == 42) uyy1d=dss42(Y_L,Y_U,SIZE,u1d,uy1d,n_l,n_u);
		// second order
		if(NDSS == 44)  uyy1d=dss44(Y_L,Y_U,SIZE,u1d,uy1d,n_l,n_u);
		// fourth order
		if(NDSS == 46)	uyy1d=dss46(Y_L,Y_U,SIZE,u1d,uy1d,n_l,n_u);
		// sixth order
		if(NDSS == 48)	uyy1d=dss48(X_L,X_U,SIZE,u1d,uy1d,n_l,n_u);
		// eighth order
		if(NDSS == 50)	uyy1d=dss50(X_L,X_U,SIZE,u1d,uy1d,n_l,n_u);
		// tenth order

		uyy[i] = uyy1d;

	}

	for (int i = 0; i<SIZE; i++) {
		for (int j = 0; j<SIZE; j++) {
			ut_[i][j] = uxx[i][j] + uyy[i][j];
		}
	}

	for (int m=0; m<SIZE; m++) {
		ut_[m][0] = 0.0;
		ut_[m][SIZE-1] = 0.0;
	}

	//
	// 2D to 1D matrix conversion
	for (int i=0; i<SIZE; i++) {
		for (int j=0; j<SIZE; j++){
			ut[i*SIZE+j]=ut_[i][j];
		}
	}
}

void laplace_solution_mol (const array_type &u , array_type &ut , const double t/* t */){

	//
	// Function pde_1 computes the temporal derivative in the
	// pseudo transient solution of Laplace’s equation by explicit
	// finite differences
	//

	double dx= (X_U - X_L )/(SIZE-1);
	double dx2 = dx*dx;

	double dy= (Y_U - Y_L)/(SIZE-1);
	double dy2 = dy*dy;

	//
	// Initially zero derivatives in x, y, t
	array_type y0(SIZE,0.0);
	matrix_type uxx(SIZE, y0);
	matrix_type uyy(SIZE, y0);
	matrix_type ut0(SIZE, y0);
	matrix_type u_(SIZE, y0);
	matrix_type ut_(SIZE, y0);

	//
	// 1D to 2D matrix conversion
	for (int i = 0; i<SIZE; i++) {
		for (int j = 0; j<SIZE; j++) {
			u_[i][j] = u[i*SIZE+j];
		}
	}

	//
	// PDE

	for (int i=0; i<SIZE; i++) {
		for (int j=0; j<SIZE; j++) {
			//
			// uxx
			if( i == 0 )
				uxx[i][j] = 2.0*( u_[i+1][j] - u_[i][j])/dx2;
			else if( i == SIZE-1 )
				uxx[i][j] = 2.0*( u_[i-1][j] - u_[i][j])/dx2;
			else
				uxx[i][j] = ( u_[i+1][j] - 2.0*u_[i][j] + u_[i-1][j])/dx2;

			//
			// uyy
			if ( j != 0 && j != SIZE-1)
				uyy[i][j]=(u_[i][j+1]-2.0*u_[i][j]+u_[i][j-1])/dy2;

			//
			// ut = uxx + uyy
			if(j == 0 || j == SIZE-1)
				ut_[i][j]=0.0;
			else
				ut_[i][j]=uxx[i][j]+uyy[i][j];

		}
	}

	//
	// 2D to 1D matrix conversion
	for (int i=0; i<SIZE; i++) {
		for (int j=0; j<SIZE; j++){
			ut[i*SIZE+j]=ut_[i][j];
		}
	}
}


int main(int argc, char **argv) {
	std::cout << "Hello world\n";

    using namespace boost::numeric::odeint;

	array_type x(SIZE);
	array_type y0(SIZE,0.0);
	array_type y(SIZE*SIZE,0.0);
	matrix_type u(SIZE, y0);
	matrix_type u_(SIZE, y0);
	array_type u1(SIZE,0); // Et
	array_type u2(SIZE,0); // Ett



    string fn= "/home/julio/temp/laplaceEquation"+ to_string(CASE) + "_"+to_string(NDSS) +".csv";

    std::ofstream out;
    out.open(fn); // append instead of overwrite
    out <<"t,";

	double dx=(X_U-X_L)/(SIZE - 1);

	// Initialise the spatial grid
	for (int i=0; i<SIZE; i++) {
		x[i]=i*dx;
		out<<x[i]<<",";
	}
    out<<"\n";
    out.close();

	// Initialise the solution
	int nx=SIZE;
	int ny=SIZE;

 	for (int i=0; i<nx; i++) {
		for (int j=0; j<ny; j++ ) {
			if(j == 0)
				u[i][j] =  0.0;
			else if(j == ny-1)
				u[i][j]=1.0;
			else
				u[i][j]=0.5;
		}
	}

	// 2D to 1D matrix conversion

 	for (int i=0; i<nx; i++) {
 		for (int j=0; j<ny; j++ ) {
 			y[i*ny+j]=u[i][j];
 		}
 	}

//	// 1D to 2D matrix conversion
//	for (int i = 0; i<SIZE; i++) {
//		for (int j = 0; j<SIZE; j++) {
//			u_[i][j] = y[i*SIZE+j];
//		}
//	}


    //
	if (NDSS < 11)
		FuncPtr = laplace_solution;
	else if (NDSS > 40)
		FuncPtr =  laplace_solution_highOrder;
	else
		FuncPtr = laplace_solution_mol;

    //size_t steps = integrate( FuncPtr, u, 0.0 , 0.15 , 0.1, streaming_observer( cout ));
    size_t steps = integrate( FuncPtr, y, 0.0 , 0.25 , 0.1, output_observer_append( fn ));
    //size_t steps = integrate( FuncPtr, u, 0.0 , 0.15, 0.1, output_observer_append_error( fn ));

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
