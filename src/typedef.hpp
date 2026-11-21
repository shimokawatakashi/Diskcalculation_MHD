#ifndef TYPEDEF_HPP_
#define TYPEDEF_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <fstream>
#include <cmath>
#include <complex>
#include <time.h>

typedef unsigned int UnsignedInt;
typedef std::vector<unsigned int> UnsignedIntVector;
typedef std::vector<int>          IntVector;
typedef std::vector<std::vector<int> > IntMatrix;
typedef std::vector<double>       DoubleVector;
typedef std::vector<std::vector<double > > DoubleMatrix;
typedef std::string String;
typedef std::vector<std::string>  StringVector;
typedef std::vector<std::complex<int> > IcomplexVector;
typedef std::complex<double> DoubleComplex;
typedef std::vector<std::complex<double> > DcomplexVector;
typedef std::vector<DcomplexVector > DcomplecMatrix;

#define FAILED_TO_OPEN(file_stream, file_name) if(!file_stream) {std::cerr << "##ERROR : Failed to open " << file_name << std::endl; exit(1); }
#define CHECK(str) std::cout << str << std::endl;

#define SQR(x) ( (x)*(x) )
#define CUB(x) ( (x)*(x)*(x) )
#define SIGN(x) ( (x >= 0.0 ? 1.0 : -1.0))

struct Mag_field
{
    DoubleMatrix Alpha_RZ_;
    DoubleVector Alpha_;
    DoubleVector B_;
    DoubleVector beta_mid;
    DoubleMatrix beta_RZ;
    double Delta_time = 0.0;
};

struct Dust_Single
{
    double _dust_inner_density_ ;
    double _dust_size_ ;
    double _dust_cross_section_ ;
    double _dust_mass_;
    double _dust_number_density_;
};

struct Dust_MRN
{
    int _number_bin_ ;
    double _dust_inner_density_ ;
    DoubleVector _dust_size_ ;
    DoubleVector _dust_cross_section_ ;
    DoubleVector _dust_mass_;
    DoubleVector _dust_number_density_;
};

struct Resistivity
{
    double _eta_O_;
    double _eta_A_;
    DoubleVector _etaO_mean_;
    DoubleVector _etaA_mean_;
};

#endif /* TYPEDEF_HPP_ */