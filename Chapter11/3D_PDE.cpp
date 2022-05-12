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
#define Z_L 0.0

#define X_U 1.0
#define Y_U 1.0
#define Z_U 1.0

#define NDSS 44

using namespace std;

double ua(double, double, int);
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
    		for( int j=0 ; j<SIZE ; ++j) {
    	     	for (int k=0; k<SIZE; k++) {
    				fout <<x[k*SIZE*SIZE + j*SIZE +i] << ",";
    			}
    			//fout<<x[SIZE-1]<<"\n";
    	     }
			    //fout<<x[SIZE-1]<<"\n";
    			fout<<",";
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

array_type extract_column_i(const matrix3D_type& a, int j, int k) {
	array_type ret(SIZE,0.0);

	for (int i=0; i<SIZE; i++) {
		ret[i] = a[i][j][k];
	}
	return (ret);

}

array_type extract_column_j(const matrix3D_type &a, int i, int k) {
	array_type ret(SIZE,0.0);

	for (int j=0; j<SIZE; j++) {
		ret[j] = a[i][j][k];
	}
	return (ret);

}

array_type extract_column_k(const matrix3D_type &a, int i, int j) {
	return (a[i][j]);

}

void insert_column_i(matrix3D_type &a, array_type b, int j, int k) {

	for (int i=0; i<SIZE; i++) {
		a[i][j][k] = b[i];
	}

}

void insert_column_j(matrix3D_type& a, array_type b, int i, int k) {

	for (int j=0; j<SIZE; j++) {
		a[i][j][k]= b[j];
	}

}

void insert_column_k(matrix3D_type& a, array_type b, int i, int j) {
	a[i][j] = b;
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
	array_type u1d(SIZE,0.0);
	array_type ux1d(SIZE,0.0);
	array_type uxx1d(SIZE,0.0);
	array_type uy1d(SIZE,0.0);
	array_type uyy1d(SIZE,0.0);
	array_type uz1d(SIZE,0.0);
	array_type uzz1d(SIZE,0.0);
	// Initially zero derivatives in x, y, textract

	array_type y0(SIZE,0.0);
	matrix_type u0(SIZE, y0);
	matrix3D_type uxx(SIZE, u0);
	matrix3D_type uyy(SIZE, u0);
	matrix3D_type uzz(SIZE, u0);
	matrix3D_type ut0(SIZE, u0);
	matrix3D_type ut_(SIZE, u0);
	matrix3D_type u_(SIZE, u0);

	//
	// 1D to 3D matrix conversion
	for (int i = 0; i<SIZE; i++)
		for (int j = 0; j<SIZE; j++)
			for (int k = 0; k<SIZE; k++)
				u_[i][j][k] = u[i*SIZE*SIZE +j*SIZE + k];

	for (int j=0; j<SIZE; j++) {
		for (int k=0; k<SIZE; k++) {
			// ux
			u1d=extract_column_i(u_,j,k);
			if (NDSS == 2) 	ux1d=dss02(X_L,X_U,SIZE,u1d); // second order
			if (NDSS == 4) 	ux1d=dss04(X_L,X_U,SIZE,u1d); // fourth order
			if (NDSS == 6)  ux1d=dss06(X_L,X_U,SIZE,u1d); // sixth order
			if (NDSS == 8)  ux1d=dss08(X_L,X_U,SIZE,u1d); // eighth order
			if (NDSS == 10)	ux1d=dss10(X_L,X_U,SIZE,u1d); // tenth order

			//uxx

			// Calculate uxx1d
			if (NDSS == 2) 	uxx1d=dss02(X_L,X_U,SIZE,ux1d); // second order
			if (NDSS == 4)	uxx1d=dss04(X_L,X_U,SIZE,ux1d); // fourth order
			if (NDSS == 6)	uxx1d=dss06(X_L,X_U,SIZE,ux1d); // sixth order
			if (NDSS == 8)	uxx1d=dss08(X_L,X_U,SIZE,ux1d); // eighth order
			if (NDSS == 10)	uxx1d=dss10(X_L,X_U,SIZE,ux1d); // tenth order

			insert_column_i(uxx, uxx1d,j,k);
		}
	}

	for (int i=0; i<SIZE; i++) {
		for (int k=0; k<SIZE; k++) {

			// uy
			u1d=extract_column_j(u_,i,k);
			if (NDSS == 2) 	uy1d=dss02(Y_L,Y_U,SIZE,u1d); // second order
			if (NDSS == 4) 	uy1d=dss04(Y_L,Y_U,SIZE,u1d); // fourth order
			if (NDSS == 6)  uy1d=dss06(Y_L,Y_U,SIZE,u1d); // sixth order
			if (NDSS == 8)  uy1d=dss08(Y_L,Y_U,SIZE,u1d); // eighth order
			if (NDSS == 10)	uy1d=dss10(Y_L,Y_U,SIZE,u1d); // tenth order

			// Calculate uyy1d
			uy1d[0] = 0.0;
			uy1d[SIZE-1] = 0.0;

			if (NDSS == 2) 	uyy1d=dss02(Y_L,Y_U,SIZE,uy1d); // second order
			if (NDSS == 4)	uyy1d=dss04(Y_L,Y_U,SIZE,uy1d); // fourth order
			if (NDSS == 6)	uyy1d=dss06(Y_L,Y_U,SIZE,uy1d); // sixth order
			if (NDSS == 8)	uyy1d=dss08(Y_L,Y_U,SIZE,uy1d); // eighth order
			if (NDSS == 10)	uyy1d=dss10(Y_L,Y_U,SIZE,uy1d); // tenth order

			insert_column_j(uyy, uyy1d,i,k);

		}
	}

	for (int i=0; i<SIZE; i++) {
		for (int j=0; j<SIZE; j++) {

			u1d=extract_column_k(u_,i,j);

			if (NDSS == 2) 	uz1d=dss02(Z_L,Z_U,SIZE,u1d); // second order
			if (NDSS == 4) 	uz1d=dss04(Z_L,Z_U,SIZE,u1d); // fourth order
			if (NDSS == 6)  uz1d=dss06(Z_L,Z_U,SIZE,u1d); // sixth order
			if (NDSS == 8)  uz1d=dss08(Z_L,Z_U,SIZE,u1d); // eighth order
			if (NDSS == 10)	uz1d=dss10(Z_L,Z_U,SIZE,u1d); // tenth order

			// Calculate uzz1d

			uz1d[0] =1.0-u1d[0];
			uz1d[SIZE-1]=1.0-u1d[SIZE-1];

			if (NDSS == 2) 	uzz1d=dss02(Z_L,Z_U,SIZE,uz1d); // second order
			if (NDSS == 4)	uzz1d=dss04(Z_L,Z_U,SIZE,uz1d); // fourth order
			if (NDSS == 6)	uzz1d=dss06(Z_L,Z_U,SIZE,uz1d); // sixth order
			if (NDSS == 8)	uzz1d=dss08(Z_L,Z_U,SIZE,uz1d); // eighth order
			if (NDSS == 10)	uzz1d=dss10(Z_L,Z_U,SIZE,uz1d); // tenth order

			insert_column_k(uzz, uzz1d,i,j);
		}
	}

	for (int i = 0; i<SIZE; i++)
		for (int j = 0; j<SIZE; j++)
			for (int k = 0; k<SIZE; k++)
				ut_[i][j][k] = uxx[i][j][k] + uyy[i][j][k] + uzz[i][j][k];

	for (int j=0; j<SIZE-1; j++){
		for (int k = 0; k<SIZE-1; k++){
			ut_[0][j][k] = 0.0;
			ut_[SIZE-1][j][k] = 0.0;
		}
	}

	// 3D to 1D matrix conversion
	for (int i=0; i<SIZE; i++) {
		for (int j=0; j<SIZE; j++) {
			for (int k=0; k<SIZE; k++ ) {
				ut[i*SIZE*SIZE +j*SIZE + k]=ut_[i][j][k];
			}
		}
	}
}


void laplace_solution_highOrder(const array_type &u , array_type &ut , const double t/* t */){

	int n_l, n_u;

	/* The routines that produces good results are DSS2 and DSS4.
	 * DSS6 and DSS8 are fine at the initial but as they approach 1, not fine
	 * DSS10 not find at all
	 */

	// Initially zero derivatives in x, y, t
	array_type u1d(SIZE,0.0);
	array_type ux1d(SIZE,0.0);
	array_type uxx1d(SIZE,0.0);
	array_type uy1d(SIZE,0.0);
	array_type uyy1d(SIZE,0.0);
	array_type uz1d(SIZE,0.0);
	array_type uzz1d(SIZE,0.0);
	// Initially zero derivatives in x, y, textract

	array_type y0(SIZE,0.0);
	matrix_type u0(SIZE, y0);
	matrix3D_type uxx(SIZE, u0);
	matrix3D_type uyy(SIZE, u0);
	matrix3D_type uzz(SIZE, u0);
	matrix3D_type ut_(SIZE, u0);
	matrix3D_type u_(SIZE, u0);

	//
	// 1D to 3D matrix conversion
	for (int i = 0; i<SIZE; i++)
		for (int j = 0; j<SIZE; j++)
			for (int k = 0; k<SIZE; k++)
				u_[i][j][k] = u[i*SIZE*SIZE +j*SIZE + k];

	for (int j=0; j<SIZE; j++) {
		for (int k=0; k<SIZE; k++) {
			// Since BCs (10.4) and (10.5) are Neumann,
			// nl = nu = 2 are inputs to dssXX corresponding to x = 0, 1, respectively.
			n_l = 1;
			n_u = 1;
			u1d=extract_column_i(u_,j,k);

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

			insert_column_i(uxx, uxx1d,j,k);

		}
	}

	for (int i=0; i<SIZE; i++) {
		for (int k=0; k<SIZE; k++) {
			// uy
			u1d=extract_column_j(u_,i,k);

			n_l = 2;
			n_u = 2;
			if (NDSS == 42) uyy1d=dss42(Y_L,Y_U,SIZE,u1d,uy1d,n_l,n_u);
			// second order
			if(NDSS == 44)  uyy1d=dss44(Y_L,Y_U,SIZE,u1d,uy1d,n_l,n_u);
			// fourth order
			if(NDSS == 46)	uyy1d=dss46(Y_L,Y_U,SIZE,u1d,uy1d,n_l,n_u);
			// sixth order
			if(NDSS == 48)	uyy1d=dss48(Y_L,Y_U,SIZE,u1d,uy1d,n_l,n_u);
			// eighth order
			if(NDSS == 50)	uyy1d=dss50(Y_L,Y_U,SIZE,u1d,uy1d,n_l,n_u);
			// tenth order
			insert_column_j(uyy, uyy1d,i,k);
		}
	}

	for (int i=0; i<SIZE; i++) {
		for (int j=0; j<SIZE; j++) {

			u1d=extract_column_k(u_,i,j);
			n_l = 2;
			n_u = 2;

			uz1d[0] =1.0-u1d[0];
			uz1d[SIZE-1]=1.0-u1d[SIZE-1];

			if (NDSS == 42) uzz1d=dss42(Z_L,Z_U,SIZE,u1d,uz1d,n_l,n_u);
			// second order
			if(NDSS == 44)  uzz1d=dss44(Z_L,Z_U,SIZE,u1d,uz1d,n_l,n_u);
			// fourth order
			if(NDSS == 46)	uzz1d=dss46(Z_L,Z_U,SIZE,u1d,uz1d,n_l,n_u);
			// sixth order
			if(NDSS == 48)	uzz1d=dss48(Z_L,Z_U,SIZE,u1d,uz1d,n_l,n_u);
			// eighth order
			if(NDSS == 50)	uzz1d=dss50(Z_L,Z_U,SIZE,u1d,uz1d,n_l,n_u);
			// tenth order
			insert_column_k(uzz, uzz1d,i,j);

		}
	}

	for (int i = 0; i<SIZE; i++)
		for (int j = 0; j<SIZE; j++)
			for (int k = 0; k<SIZE; k++)
				ut_[i][j][k] = uxx[i][j][k] + uyy[i][j][k] + uzz[i][j][k];

	for (int j=0; j<SIZE-1; j++){
		for (int k = 0; k<SIZE-1; k++){
			ut_[0][j][k] = 0.0;
			ut_[SIZE-1][j][k] = 0.0;
		}
	}

	// 3D to 1D matrix conversion
	for (int i=0; i<SIZE; i++) {
		for (int j=0; j<SIZE; j++) {
			for (int k=0; k<SIZE; k++ ) {
				ut[i*SIZE*SIZE +j*SIZE + k]=ut_[i][j][k];
			}
		}
	}

}

void laplace3D_solution_mol (const array_type &u , array_type &ut , const double t/* t */){

	//
	// Function pde_1 computes the temporal derivative in the
	// pseudo transient solution of Laplace’s equation in 3D by explicit
	// finite differences
	//

	double dx= (X_U - X_L )/(SIZE-1);
	double dx2 = dx*dx;

	double dy= (Y_U - Y_L)/(SIZE-1);
	double dy2 = dy*dy;

	double dz= (Z_U - Z_L)/(SIZE-1);
	double dz2 = dz*dz;
	//
	// Initially zero derivatives in x, y, t
	array_type y0(SIZE,0.0);
	matrix_type u0(SIZE, y0);
	matrix3D_type uxx(SIZE, u0);
	matrix3D_type uyy(SIZE, u0);
	matrix3D_type uzz(SIZE, u0);
	matrix3D_type ut0(SIZE, u0);
	matrix3D_type ut_(SIZE, u0);
	matrix3D_type u_(SIZE, u0);

	//
	// 1D to 2D matrix conversion
	for (int i = 0; i<SIZE; i++)
		for (int j = 0; j<SIZE; j++)
			for (int k = 0; k<SIZE; k++)
				u_[i][j][k] = u[i*SIZE*SIZE +j*SIZE + k];
//
	// PDE

	for (int i=0; i<SIZE; i++) {
		for (int j=0; j<SIZE; j++) {
			for (int k=0;k<SIZE; k++) {
			//

			// uxx
				if (i != 0 && i != SIZE-1){
					uxx[i][j][k] = (u_[i+1][j][k] - 2.0*u_[i][j][k] + u_[i-1][j][k])/dx2;
				}

			// uyy
				if( j == 0 )
					uyy[i][j][k] = 2.0*( u_[i][j+1][k] - u_[i][j][k])/dy2;
				else if( j == SIZE-1 )
					uyy[i][j][k] = 2.0*( u_[i][j-1][k] - u_[i][j][k])/dy2;
				else
					uyy[i][j][k] = ( u_[i][j+1][k] - 2.0*u_[i][j][k] + u_[i][j-1][k])/dy2;

			// uzz
				double u_00, u_11;

				if (k == 0) {
					u_00 = u_[i][j][k+1] - 2.0*dz*(1.0-u_[i][j][k]);
					uzz[i][j][k] = (u_[i][j][k+1] - 2.0*u_[i][j][k] + u_00)/dz2;
				} else if (k == SIZE-1) {
					u_11 = u_[i][j][k-1] + 2.0*dz*(1.0-u_[i][j][k]);
					uzz[i][j][k] = (u_11 - 2.0*u_[i][j][k] + u_[i][j][k-1])/dz2;
				}
				else {
					uzz[i][j][k] = (u_[i][j][k+1] - 2.0*u_[i][j][k] + u_[i][j][k-1])/dz2;
				}


				//
				// ut = uxx + uyy + uzz
				if (i == 0 || i == SIZE-1)
	//				For the BCs in x, we use the non-homogeneous or in-homogeneous (nonzero) Dirichlet conditions
	//				u(x = 0, y, z, t) = u(x = 1, y, z, t) = 1. In other words, since the boundary values are
	//				constant, their derivatives in t are zero.
					ut_[i][j][k] = 0.0;
				else
					ut_[i][j][k] = uxx[i][j][k] + uyy[i][j][k] + uzz[i][j][k];

			}
		}
	}


	//
	// 3D to 1D matrix conversion
	for (int i=0; i<SIZE; i++) {
		for (int j=0; j<SIZE; j++) {
			for (int k=0; k<SIZE; k++ ) {
				ut[i*SIZE*SIZE +j*SIZE + k]=ut_[i][j][k];
			}
		}
	}
}


int main(int argc, char **argv) {
	std::cout << "Hello world\n";

    using namespace boost::numeric::odeint;

	array_type x(SIZE);
	array_type y0(SIZE,0.0);
	matrix_type u0(SIZE, y0);
	matrix3D_type u(SIZE,u0);

	array_type y(SIZE*SIZE*SIZE,0.0);


    string fn= "/home/julio/temp/3Dimensional_PDE_"+ to_string(CASE) + "_"+to_string(NDSS) +".csv";

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
	int nz=SIZE;

	int c=0;
 	for (int i=0; i<nx; i++) {
		for (int j=0; j<ny; j++ ) {
			for (int k=0; k<nz; k++) {
//				u[i][j][k] = c++;
				if (i == 0)
					u[i][j][k] = 1.0;
				else if (i == nx-1)
					u[i][j][k] = 1.0;
				else
					u[i][j][k] = 0.0;
			}
		}
	}

	// 3D to 1D matrix conversion
 	for (int i=0; i<nx; i++) {
 		for (int j=0; j<ny; j++ ) {
 	 		for (int k=0; k<nz; k++ ) {
 	 			y[i*ny*nz +j*nz + k]=u[i][j][k];
 	 		}
 		}
 	}

    //
	if (NDSS < 11)
		FuncPtr = laplace_solution;
	else if (NDSS > 40)
		FuncPtr =  laplace_solution_highOrder;
	else
		FuncPtr = laplace3D_solution_mol;

    //size_t steps = integrate( FuncPtr, y, 0.0 , 1.0, 0.1, streaming_observer( cout ));
    size_t steps = integrate( FuncPtr, y, 0.0 , 1.0 , 0.1, output_observer_append( fn ));
    //size_t steps = integrate( FuncPtr, y, 0.0 , 1.0, 0.1, output_observer_append_error( fn ));

   cout << "!!!Hello World!!! with " <<steps<< endl; // prints !!!Hello World!!!

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
