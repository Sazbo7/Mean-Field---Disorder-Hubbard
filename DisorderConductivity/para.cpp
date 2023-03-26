
#include "susc.h"
#include "src/localCondCalc.h"

extern "C" {

int compute_paramagnetic_susceptibility(
  // input
  double temperature,
  double frequency,
  double broadening,
  int num_bond,
  int const * rowsites,
  int const * colsites,
  void const * v_amplitudes,
  double const * displacements_x,
  double const * displacements_y,
  int num_eigen,
  int b1, // number of sites in direction 1
  int b2, // number of sites in direction 2
  double const * eigenvalues,
  void const* v_eigenvectors_up,
  void const* v_eigenvectors_dn,
  // output
  void * v_susceptibility_x_up,
  void * v_susceptibility_x_dn,
  void * v_susceptibility_y_up,
  void * v_susceptibility_y_dn
)
{
    using namespace std;
    using namespace localcondcalc;

    complex<double> const * amplitudes
    = reinterpret_cast<complex<double> const*>(v_amplitudes);
    complex<double> const * eigenvectors_up
    = reinterpret_cast<complex<double> const*> (v_eigenvectors_up);
    complex<double> const * eigenvectors_dn
    = reinterpret_cast<complex<double> const*> (v_eigenvectors_dn);
    complex<double> * susceptibility_x_up
    = reinterpret_cast<complex<double> *> (v_susceptibility_x_up);
    complex<double> * susceptibility_x_dn
    = reinterpret_cast<complex<double> *> (v_susceptibility_x_dn);
    complex<double> * susceptibility_y_up
    = reinterpret_cast<complex<double> *> (v_susceptibility_y_up);
    complex<double> * susceptibility_y_dn
    = reinterpret_cast<complex<double> *> (v_susceptibility_y_dn);


    vector<double> nrgs( eigenvalues , eigenvalues + num_eigen );
    vector<complex<double> > eigVecsUp( eigenvectors_up , eigenvectors_up + num_eigen*b1*b2 );
    vector<complex<double> > eigVecsDn( eigenvectors_dn , eigenvectors_dn + num_eigen*b1*b2 );

    //Convert bond list
    vector<tuple<int,int,int,int> > bondList;
    for (int ii = 0 ; ii < num_bond ; ++ii )
    {
        auto r = rowsites[ii];
        auto c = colsites[ii];
        auto dx = static_cast<int>(round(displacements_x[ii]));
        auto dy = static_cast<int>(round(displacements_y[ii]));
        bondList.push_back(make_tuple(r, c, dx, dy));
    }

    auto numStates = num_eigen;

    // Make lind coefs with mu = 0.0
    auto staticLindhards = getStaticLindhards(nrgs, temperature, 0.0);

    std::cout << "Full # of Lindhard Coefficients = " << numStates * (numStates-1) / 2 << std::endl;
    std::cout << "Reduced # of Lindhard Coefficients = " << staticLindhards.size() << std::endl;

    auto dynamicLindhards = getDynamicLindhards(nrgs,
                                                staticLindhards,
                                                frequency,
                                                broadening,
                                                temperature,
                                                0.0);

    std::cout << "Computing Up" << std::endl;
    {
        std::cout << "Computing Current Elements" << std::endl;
        auto currentElementsUp = getCurrentElements(bondList, staticLindhards, eigVecsUp, numStates);
        std::cout << "Computing Susceptibility" << std::endl;
#pragma omp parallel for
        for (int ii = 0 ; ii < num_bond ; ++ii )
        {
            auto susUp = getBondFieldSusceptibility(bondList, dynamicLindhards, currentElementsUp, numStates, ii);
            susceptibility_x_up[ii] = get<0>(susUp);
            susceptibility_y_up[ii] = get<1>(susUp);
        }
    }

    std::cout << "Computing Dn" << std::endl;
    {
        std::cout << "Computing Current Elements" << std::endl;
        auto currentElementsDn = getCurrentElements(bondList, staticLindhards, eigVecsDn, numStates);
        std::cout << "Computing Susceptibility" << std::endl;
#pragma omp parallel for
        for (int ii = 0 ; ii < num_bond ; ++ii )
        {
            auto susDn = getBondFieldSusceptibility(bondList, dynamicLindhards, currentElementsDn, numStates, ii);
            susceptibility_x_dn[ii] = get<0>(susDn);
            susceptibility_y_dn[ii] = get<1>(susDn);
        }
    }

    return 0;
}



int compute_paramagnetic_susceptibility_naive(
  // input
  double temperature,
  double frequency,
  double broadening,
  int num_bond,
  int const * rowsites,
  int const * colsites,
  void const * v_amplitudes,
  double const * displacements_x,
  double const * displacements_y,
  int num_eigen,
  int b1, // number of sites in direction 1
  int b2, // number of sites in direction 2
  double const * eigenvalues,
  void const* v_eigenvectors_up,
  void const* v_eigenvectors_dn,
  // output
  void * v_susceptibility_x_up,
  void * v_susceptibility_x_dn,
  void * v_susceptibility_y_up,
  void * v_susceptibility_y_dn
)
{
    using namespace std;
    using namespace localcondcalc;

    complex<double> const * amplitudes
    = reinterpret_cast<complex<double> const*>(v_amplitudes);
    complex<double> const * eigenvectors_up
    = reinterpret_cast<complex<double> const*> (v_eigenvectors_up);
    complex<double> const * eigenvectors_dn
    = reinterpret_cast<complex<double> const*> (v_eigenvectors_dn);
    complex<double> * susceptibility_x_up
    = reinterpret_cast<complex<double> *> (v_susceptibility_x_up);
    complex<double> * susceptibility_x_dn
    = reinterpret_cast<complex<double> *> (v_susceptibility_x_dn);
    complex<double> * susceptibility_y_up
    = reinterpret_cast<complex<double> *> (v_susceptibility_y_up);
    complex<double> * susceptibility_y_dn
    = reinterpret_cast<complex<double> *> (v_susceptibility_y_dn);


    vector<double> nrgs( eigenvalues , eigenvalues + num_eigen );
    vector<complex<double> > eigValsUp( eigenvectors_up , eigenvectors_up + num_eigen*b1*b2 );
    vector<complex<double> > eigValsDn( eigenvectors_dn , eigenvectors_dn + num_eigen*b1*b2 );

    //Convert bond list
    vector<tuple<int,int,int,int> > bondList;
    for ( int ii = 0 ; ii < num_bond ; ii++ )
    {
        tuple<int,int,int,int> tempTuple( rowsites[ii] , colsites[ii] ,
                                          (int) round(displacements_x[ii]) , (int) round(displacements_y[ii]) );
        bondList.push_back(tempTuple);
    }

    // Make lind coefs with mu = 0.0
    auto lindCoefs
    = naive::getLindCoefs( nrgs , frequency , broadening , temperature , 0.0 );

    // Calculate
    for ( int ii = 0 ; ii < num_bond ; ii++ )
    {
        tuple< complex<double>, complex<double> > susUp
        = naive::getBondCorrFunc( bondList , eigValsUp , lindCoefs , get<0>(bondList[ii]) , get<1>(bondList[ii]) );
        tuple< complex<double>, complex<double> > susDn
        = naive::getBondCorrFunc( bondList , eigValsDn , lindCoefs , get<0>(bondList[ii]) , get<1>(bondList[ii]) );
        susceptibility_x_up[ii] = get<0>(susUp);
        susceptibility_x_dn[ii] = get<0>(susDn);
        susceptibility_y_up[ii] = get<1>(susUp);
        susceptibility_y_dn[ii] = get<1>(susDn);
    }

    return 0;
}

}
