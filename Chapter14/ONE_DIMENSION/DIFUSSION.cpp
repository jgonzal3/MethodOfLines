#include <iostream>
#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <boost/numeric/odeint.hpp>
#include "CONSTANTS.hpp"

#include "dssCPP.h"
#define NDSS 6

using namespace std;

double ua(double, double);
ostream& operator<<(ostream& os, const array_type& a);
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
    	cout<<t<<"\n";
    	cout<<x;
        ++count_;
    	std::ofstream fout(filename_, std::ios_base::app);
        fout<<"t_Vr-->,z_V,0,1,2,3,4,5,6\n";
        for( int j=0 ; j<NZ ; ++j) {
     	fout<<t<<","<<(j+1)*5<<",";
   		for (int i=0; i<NR; i++) {
    			fout <<x[j*NR + i] << ",";
    		}
    		fout<<"\n";
    		}
		fout<<"\n";

        for( int j=0 ; j<NZ ; ++j) {
     	fout<<t<<","<<(j+1)*5<<"," ;
   		for (int i=0; i<NR; i++) {
    			fout<<x[j*NR + i + NR*NZ] << ",";
    		}
    		fout<<"\n";
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

void difussion_pde_solution (const array_type &u , array_type &ut , const double t/* t */){

	// 1D to 2D matrices
    // Grid in axial direction
	array_type y0(NR);
	array_type r(NR);
	array_type z(NZ);

	matrix_type ca(NZ,y0);
	matrix_type Tk(NZ,y0);
	matrix_type Tkr(NZ,y0);
	matrix_type Tkz(NZ,y0);
	matrix_type car(NZ,y0);
	matrix_type caz(NZ,y0);
	matrix_type Tkrr(NZ,y0);
	matrix_type Tkzz(NZ,y0);
	matrix_type carr(NZ,y0);
	matrix_type cazz(NZ,y0);
	matrix_type cat(NZ,y0);
	matrix_type Tkt(NZ,y0);

	array_type y(2*NR*NZ);
//
	double dz = (1.0*ZL)/NZ;
	for (int i=1; i<=NZ; i++)
		z[i] = i*dz;

	// Grid in radial direction
	double dr = r0/(NR-1);

	for (int j=0; j<NR; j++)
		r[j] = j*dr;

	double drs = dr*dr;

	for (auto i = 0; i<NZ; i++) {
		for (auto j = 0; j<NR; j++) {
			int ij = i*NR + j;
			ca[i][j] = u[ij];
			Tk[i][j] = u[ij+NR*NZ];
		}
	}


	// Step through the grid points in r and z
	for (auto i=0; i<NZ; i++) {
		for (auto j=0; j<NR; j++) {
	//
		// (1/r)*car, (1/r)*Tkr
		if(j==0) {
			car[i][j]=2.0*(ca[i][j+1]-ca[i][j])/drs;
			Tkr[i][j]=2.0*(Tk[i][j+1]-Tk[i][j])/drs;
		}
		else if (j == NR-1) {
			car[i][j]=0.0;
			Tkr[i][j]=(1.0/r[j])*(h/k)*(Tw-Tk[i][j]);
		}
		else {
			car[i][j]=(1.0/r[j])*(ca[i][j+1]-ca[i][j-1])/(2.0*dr);
			Tkr[i][j]=(1.0/r[j])*(Tk[i][j+1]-Tk[i][j-1])/(2.0*dr);
		}

		// carr, Tkrr
		if(j == 0) {
			carr[i][j]=2.0*(ca[i][j+1]-ca[i][j])/drs;
			Tkrr[i][j]=2.0*(Tk[i][j+1]-Tk[i][j])/drs;
		}
		else if(j == NR-1) {
			carr[i][j]=2.0*(ca[i][j-1]-ca[i][j])/drs;
			double Tkf=Tk[i][j-1]+2.0*dr*(h/k)*(Tw-Tk[i][j]);
			Tkrr[i][j]=(Tkf-2.0*Tk[i][j]+Tk[i][j-1])/drs;
		}
		else {
			carr[i][j]=(ca[i][j+1]-2.0*ca[i][j]+ca[i][j-1])/drs;
			Tkrr[i][j]=(Tk[i][j+1]-2.0*Tk[i][j]+Tk[i][j-1])/drs;
		}
		//
		// caz, Tkz
		if (i == 0) {
			caz[i][j]=(ca[i][j]-cae)/dz;
			Tkz[i][j]=(Tk[i][j]-Tke)/dz;
		}
		else {
			caz[i][j]=(ca[i][j]-ca[i-1][j])/dz;
			Tkz[i][j]=(Tk[i][j]-Tk[i-1][j])/dz;
		}

		// PDEs
		double rk=rk0*exp(-E/(R*Tk[i][j]))*(ca[i][j]*ca[i][j]);
		cat[i][j]=Dc*(carr[i][j]+car[i][j])-v*caz[i][j]-rk;
		Tkt[i][j]=Dt*(Tkrr[i][j]+Tkr[i][j])-v*Tkz[i][j]-(dH*rk)/(rho*Cp);

		}
	}

	//2D to 1D matrices
	for (int i=0; i<NZ; i++) {
		for (int j=0; j<NR; j++) {
			int ij= i*NR+j;
			ut[ij]=cat[i][j];
			ut[ij+NR*NZ]=Tkt[i][j];
		}
	}

	//cout<<ut;

}

array_type extract_column(matrix_type a, int pos) {
	array_type ret(NR,0.0);

	for (int k=0; k<NR; k++)
		ret[k] = a[k][pos];
	return (ret);

}

void difussion_pde_solution_highOrder (const array_type &u , array_type &ut , const double t/* t */){

	array_type y0(NR);
	array_type r(NR);
	array_type z(NZ);

	array_type ca1d(NZ);
	array_type Tk1d(NZ);
	array_type car1d(NZ);
	array_type Tkr1d(NZ);
	array_type carr1d(NZ);
	array_type Tkrr1d(NZ);

	matrix_type ca(NZ,y0);
	matrix_type Tk(NZ,y0);
	matrix_type Tkr(NZ,y0);
	matrix_type Tkz(NZ,y0);
	matrix_type car(NZ,y0);
	matrix_type caz(NZ,y0);
	matrix_type Tkrr(NZ,y0);
	matrix_type Tkzz(NZ,y0);
	matrix_type carr(NZ,y0);
	matrix_type cazz(NZ,y0);
	matrix_type cat(NZ,y0);
	matrix_type Tkt(NZ,y0);

	array_type y(2*NR*NZ);
//
	double dz = (1.0*ZL)/NZ;
	for (int i=1; i<=NZ; i++)
		z[i] = i*dz;

	// Grid in radial direction
	double dr = r0/(NR-1);

	for (int j=0; j<NR; j++)
		r[j] = j*dr;

	double drs = dr*dr;

	for (auto i = 0; i<NZ; i++) {
		for (auto j = 0; j<NR; j++) {
			int ij = i*NR + j;
			ca[i][j] = u[ij];
			Tk[i][j] = u[ij+NR*NZ];
		}
	}

	//
	// Step through the grid points in r and z
	for (int i=0; i<NZ; i++) {
		ca1d = ca[i];
		Tk1d = Tk[i];

//		//
//		// car, Tkr
		car1d=dss08(0.0,r0,NR,ca1d);
		car[i]= car1d;
		car[i][0]= 0.0;
		car[i][NR-1]=0.0;
//
		Tkr1d=dss08(0.0,r0,NR,Tk1d);
		Tkr[i] = Tkr1d;
		Tkr[i][0]= 0.0;
		Tkr[i][NR-1] = (h/k)*(Tw-Tk[i][NR-1]);

		// carr, Tkrr
//
		car1d[0   ] = 0.0;
		car1d[NR-1] = 0.0;
//
		double nl = 2.0;
		double nu = 2.0;
//
		carr1d=dss48(0.0,r0,NR,ca1d,car1d,nl,nu);
		carr[i] = carr1d;

		Tkr1d[0]=0.0;
		Tkr1d[NR-1] = (h/k)*(Tw-Tk1d[NR-1]);

		nl=2;
		nu=2;
		Tkrr1d=dss48(0.0,r0,NR,Tk1d,Tkr1d,nl,nu);

		Tkrr[i] = Tkrr1d;

		for (int j=0; j<NR; j++) {
//			//
//			// (1/r)*car, (1/r)*Tkr
			if(j != 0) {
				car[i][j]=(1.0/r[j])*car[i][j];
				Tkr[i][j]=(1.0/r[j])*Tkr[i][j];
			}
			if(j == 0)
				carr[i][j] = 2.0*carr[i][j];
			if(j == 0) {
				Tkrr[i][j] = 2.0*Tkrr[i][j];
			}
//
//			//
//			// caz, Tkz
			if (i == 0) {
				caz[i][j]=(ca[i][j]-cae)/dz;
				Tkz[i][j]=(Tk[i][j]-Tke)/dz;
			}
			else {
				caz[i][j]=(ca[i][j]-ca[i-1][j])/dz;
				Tkz[i][j]=(Tk[i][j]-Tk[i-1][j])/dz;
			}

//
//			//
//
//			// PDEs
			double rk=rk0*exp(-E/(R*Tk[i][j]))*(ca[i][j]*ca[i][j]);
			cat[i][j]=Dc*(carr[i][j]+car[i][j])-v*caz[i][j]-rk;
			Tkt[i][j]=Dt*(Tkrr[i][j]+Tkr[i][j])-v*Tkz[i][j]-(dH*rk)/(rho*Cp);
		}
	}


	//2D to 1D matrices
	for (int i=0; i<NZ; i++) {
		for (int j=0; j<NR; j++) {
			int ij= i*NR+j;
			ut[ij]=cat[i][j];
			ut[ij+NR*NZ]=Tkt[i][j];
		}
	}
	//

	cout<<ut;
}

int main(int argc, char **argv) {
	std::cout << "Hello world\n";

    using namespace boost::numeric::odeint;

	array_type u(2*NR*NZ);
	array_type y(2*NR*NZ);
	array_type y0(NR,0.0);
	matrix_type ca(NZ, y0);
	matrix_type Tk(NZ, y0);


    string fn= "/home/julio/temp/Difussion_PDE_"+to_string(NDSS) +".csv";

    std::ofstream out;
    out.open(fn); // append instead of overwrite
    out <<"t,";
	// Initial condition
	for (int i=0; i<NZ; i++) {
		for (int j=0; j<NR; j++) {
			ca[i][j] = ca0;
			Tk[i][j] = Tk0;
			u[i*NR + j        ] = ca[i][j];
			u[i*NR + j + NZ*NR] = Tk[i][j];
		}
	}


/*	Equation (12.1) is first order in t and second order in x (through
 *                    3    2
 * the mixed partial ∂ u/∂x ∂t). It therefore requires one initial condition
 * (IC) and two boundary conditions (BCs). The IC is taken as
 *
	u(x, t = 0) = sin(πx/L) (12.2)

	and the two BCs as

	u(x = 0, t) = u(x = L, t) = 0

	*/
	/*
	 *                                                       2
	    C ∂Tk/∂t = − 1/ρr ∂(rqh)/∂r − vC ∂Tk/∂z - (dHkr/ρ) c
	     p								p 					a

                   	  2    2                                        2
	    ∂Tk/∂t 	= D (∂ Tk/∂ r + 1/r ∂Tk/∂r) − v ∂Tk∂z − dHkr/(ρC )c
                   t                                            p  a

	    kr  = k0 exp(−E/(RTk))
	*/
	if (NDSS > 1)
		FuncPtr = difussion_pde_solution_highOrder;
	else
		FuncPtr = difussion_pde_solution;

    //size_t steps = integrate( FuncPtr, u, 0.0 , 20.0, 0.1, streaming_observer( cout ));
    size_t steps = integrate( FuncPtr, u, 0.0 , 200.0 , 0.1, output_observer_append( fn ));
    //size_t steps = integrate( FuncPtr, u, 0.0 , 20.0, 0.1, output_observer_append_error( fn ));

   cout << "!!!Hello World!!! with " <<steps<< endl; // prints !!!Hello World!!!
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


