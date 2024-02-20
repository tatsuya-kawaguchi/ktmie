//
//
//  ktmie   Mie code by using recurr equations.
//
//  Created by Tatsuya Kawaguchi
//
//  compile:
//    c++ -O2 ktmie.cpp -std=c++11

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <vector>
using std::vector;

#include <fstream>
using std::ofstream;

#include <complex>
using std::complex;

#include <string>
using std::string;

#include <iomanip>
#include <cmath>
#include <algorithm>




int ktmie(
    const double x,                 //  size parameter
    const complex<double> m,        //  refractive index
    const int N_ANG,                //  number of angle
    vector<complex<double> >& S1,   //  return value
    vector<complex<double> >& S2    //  return value
);


int main(int argc, char* argv[])
{
    //////////////////////////////////////////////////////
    //  size parameter
    //////////////////////////////////////////////////////
    const double x = 100.0;

    //////////////////////////////////////////////////////
    //  refractive index of both sphere and environment.
    //////////////////////////////////////////////////////
#if 1

    //  water droplet in air
    const complex<double> m1 = complex<double>(1.33, 0.0);
    const double m2 = 1.00;

#elif 0

    //  gas bubble in water
    const complex<double> m1 = complex<double>(1.00, 0.0);
    const double m2 = 1.33;

#else 

    //  absorbing sphere
    const complex<double> m1 = complex<double>(1.33, -0.0);
    const double m2 = 1.00;

#endif

    //////////////////////////////////////////////////////
    //  number of anguler division
    //
    //  e.g.
    //
    //  181 = 1.0 deg/div
    //  1801 = 0.1 deg/div
    //  18001 = 0.01 deg/div
    //////////////////////////////////////////////////////
    const int N_ANG = 1801;



    //////////////////////////////////
    //
    //  parameters of physical values
    //
    const complex<double> m = m1/m2;



    vector<complex<double> > S1;
    vector<complex<double> > S2;


    ktmie(x, m, N_ANG, S1, S2);


    if(1) {
        string fname = "ktmie.csv";
        ofstream ofs(fname);
        if(!ofs) {
            cerr << "ERROR: fail to open " << fname << endl;
            return 1;
        }
        ofs << "x=," << x << "\n";
        ofs << "m1=," << m1.real() << "," << m1.imag() << "\n";
        ofs << "m2=," << m2 << "\n";
        ofs << "Nun of Angles=," << N_ANG << "\n";
        ofs << "\n";
        
        //  write data
        ofs << "theta,i1,i2,log10(i1),log10(i2)\n\n";
        for(int i=0; i<N_ANG; i++) {
            ofs << (double)i / (N_ANG-1) * 180 << "," // in degree
                << norm(S1[i]) << "," 
                << norm(S2[i]) << ","
                << log10(norm(S1[i])) << ","
                << log10(norm(S2[i])) << "\n";
        }
        ofs << "\n";

        ofs.close();
    }

    return 0;
}






int ktmie(
    const double x,                 //  size parameter
    const complex<double> m,        //  refractive index
    const int N_ANG,                //  number of angle
    vector<complex<double> >& S1,   //  return value
    vector<complex<double> >& S2    //  return value
)
{

	/////////////////////////////////
	//
	//  parameters for computation
	//
	//	definition of the imaginary number, $i$
	const complex<double> Z_ONE = complex<double>(0.0, 1.0);



	//  store theta(rad)
	vector<double> theta(N_ANG);
	for(int i=0; i<theta.size(); i++) {
		theta[i] = i * M_PI / (N_ANG-1);  //  radian 
	}

	//  relative complex refractive index
	const complex<double> y = m * x;

	//  iteration limit number
	const int n_end = x+4*pow(x,1.0/3.0)+2 + 50;

	//  if the delta S becomes very small, end the iteration
	const double EPS = 1.0e-10;

	//  Array of S1 and S2
	S1.resize(theta.size(), complex<double>(0.0, 0.0));
	S2.resize(theta.size(), complex<double>(0.0, 0.0));



	/////////////////////////////////////////////////
	//
	// Riccati-Bessel functions and relating valiables
	// for recurrence calcuration of a_n and b_n
	//
	//  psi_n_x,  psi_n_y; --- recurrence
	//  dpsi_n_x, dpsi_n_y --- could be calc by psi
	//  zeta_n_x        --- could be obtained from psi and chi
	//  dzeta_n_x;      --- could be obtained from zeta = psi, chi
	//
	//  Consequently, psi, chi must be calcurated by recurrence formula
	//
	/////////////////////////////////////////////////

	complex<double> psi_n_x;	//	== S_n(x)
	complex<double> psi_n1_x;	//	psi_{n-1}
	complex<double> psi_n2_x;	//	psi_{n-2}

	complex<double> psi_n_y;	//	== S_n(y)
	complex<double> psi_n1_y;	//	psi_{n-1}
	complex<double> psi_n2_y;	//	psi_{n-2}

	complex<double> chi_n_x;	//	== C_n(x)
	complex<double> chi_n1_x;	//	chi_{n-1}
	complex<double> chi_n2_x;	//	chi_{n-2}


	//	set initial values
	psi_n2_x = cos(x);	//	psi_{n-2}
	psi_n1_x = sin(x);	//	psi_{n-1}

	psi_n2_y = cos(y);	//	psi_{n-2}
	psi_n1_y = sin(y);	//	psi_{n-1}

	chi_n2_x = -sin(x);	//	chi_{n-2}
	chi_n1_x =  cos(x);	//	chi_{n-1}


    //////////////////////////////////////////////////////
    //
    //  angular functions
    //
    vector<double> pi_n, pi_n1, pi_n2;    //  pi_n(cos theta)
    //  tau_n is obtained from pi_n and pi_n1

	//  set initial values
	pi_n2.resize(theta.size(), 0.0);
	pi_n1.resize(theta.size(), 1.0);


	/////////////////////////////////
	//  master loop
	/////////////////////////////////
	for(int n=1; n<=n_end; n++) {

		cerr << "n=" << n << " of " << n_end << endl;


		////////////////////////////////////////////////////////////////////
		//  (1) set array of angular function
		////////////////////////////////////////////////////////////////////

		if( n == 1 ) {
			pi_n = pi_n1;
		} else {
			for(int i=0; i<theta.size(); i++) {
				const double x = cos(theta[i]); //  x is not a size parameter here, but a cos(theta).
				pi_n[i] = (2.0*n-1.0)/(n-1.0) * x * pi_n1[i] - n/(n-1.0) * pi_n2[i];
			}
		}

		////////////////////////////////////////////////////////////////////
		//  (2) set coefficients, a_n and b_n
		////////////////////////////////////////////////////////////////////

		//  For the calculation of coefficients, a_n and b_n, the following terms are required.
		//
		//  \psi_n(x)       psi_n_x
		//  \psi_n(y)       psi_n_y
		//  \psi'_n(x)      dpsi_n_x
		//  \psi'_n(y)      dpsi_n_y
		//  \zeta_n(x)      zeta_n_x
		//  \zeta'_n(x)     dzeta_n_x

		//	(2-1) calc n'th values of recurrence formura
		psi_n_x = (2.0*n-1.0)/(x) * psi_n1_x - psi_n2_x;	//	psi_1, psi_2, ...
		psi_n_y = (2.0*n-1.0)/(y) * psi_n1_y - psi_n2_y;	//	psi_1, psi_2, ...
		chi_n_x = (2.0*n-1.0)/(x) * chi_n1_x - chi_n2_x;	//	chi_1, chi_2, ...
		
		// (2-2) calc other functions
		const complex<double> zeta_n1_x = psi_n1_x + Z_ONE * chi_n1_x;
		const complex<double> zeta_n_x  = psi_n_x  + Z_ONE * chi_n_x ;
		
		const complex<double> dpsi_n_x  = psi_n1_x  - (double)n / x * psi_n_x;
		const complex<double> dpsi_n_y  = psi_n1_y  - (double)n / y * psi_n_y;
		const complex<double> dzeta_n_x = zeta_n1_x - (double)n / x * zeta_n_x;
		
		
		//  (2-3) calc a_n and b_n
		//  Definition of the coefficients is appeared in the book by van de Hulst, p.123
		const complex<double> a_n
		 = (    dpsi_n_y * psi_n_x - m * psi_n_y * dpsi_n_x) / (    dpsi_n_y * zeta_n_x - m * psi_n_y * dzeta_n_x);
		
		const complex<double> b_n
		 = (m * dpsi_n_y * psi_n_x -     psi_n_y * dpsi_n_x) / (m * dpsi_n_y * zeta_n_x -     psi_n_y * dzeta_n_x);
		
		////////////////////////////////////////////////////////////////////
		//  (3) calc n'th term of the S, in p.125
		////////////////////////////////////////////////////////////////////

		//  the n'th S
		vector<complex<double> > S1_n(S1.size(), 0.0);
		vector<complex<double> > S2_n(S2.size(), 0.0);
		
		for(int i=0; i<theta.size(); i++) {
			const double x = cos(theta[i]);
			const double tau_n = n * x * pi_n[i] - (n+1.0)*pi_n1[i];
		
			S1_n[i] = (2.0*n+1.0)/(n*(n+1.0)) * (a_n * pi_n[i] + b_n * tau_n );
			S2_n[i] = (2.0*n+1.0)/(n*(n+1.0)) * (b_n * pi_n[i] + a_n * tau_n );
		}
		
		//////////////////////////////////////////////////////////////////
		//  check the n'th term of S1 and S2
		//  if the absolute values become small, terminate the loop of n
		//////////////////////////////////////////////////////////////////
		{
			//  define unnamed function for the comparison of the complex values
			auto complex_compare = [](complex<double> a, complex<double> b)
			{
				return abs(a) < abs(b);
			};
	
			const complex<double> max_S1_n = *std::max_element(S1_n.begin(), S1_n.end(), complex_compare);
			const complex<double> max_S2_n = *std::max_element(S2_n.begin(), S2_n.end(), complex_compare);
	
			if(fmax(abs(max_S1_n), abs(max_S2_n)) < EPS) {
				cerr << "break loop of n" << endl;
				break;
			}
		}
	
		// update S1 and S2
		for(int i=0; i<theta.size(); i++) {
			S1[i] += S1_n[i];
			S2[i] += S2_n[i];
		}
		
		////////////////////////////////////////////////////////////////////
		// (4-1) update Riccati-Bessel functions for the next loop
		////////////////////////////////////////////////////////////////////
		{
			psi_n2_x = psi_n1_x;
			psi_n1_x = psi_n_x;
		
			psi_n2_y = psi_n1_y;
			psi_n1_y = psi_n_y;
		
			chi_n2_x = chi_n1_x;
			chi_n1_x = chi_n_x;
		}
		
		////////////////////////////////////////////////////////////////////
		//  (4-2) update angular functions for the next loop
		////////////////////////////////////////////////////////////////////
		if( n >= 2 ) {
			pi_n2 = pi_n1;
			pi_n1 = pi_n;
		}
		
		if(0) { 
			//  write data
			cout << "n=," << n << ",,";
			for(int i=0; i<theta.size(); i++) {
//				cout << theta[i]*180/M_PI << "," << log10(norm(S1[i])) << "," << log10(norm(S2[i])) << "\n";
				cout << log10(norm(S1[i])) << ",";
//				cout << log10(norm(S2[i])) << ",";
			}
			cout << "\n";
		}
	
	}
	
    return 0;
}
