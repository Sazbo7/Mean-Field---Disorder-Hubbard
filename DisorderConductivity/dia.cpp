
#include <algorithm>
#include <complex>
#include <vector>

#include "susc.h"

typedef std::complex<double> complex128_t;

//! Default tolerance for double type
static const double kDoubleTolerance = std::sqrt(std::numeric_limits<double>::epsilon());

inline
std::function<double(double)>
makeFermiDirac(double temperature,
               double mu=0.0,
               double tolerance=kDoubleTolerance)
{
    if (temperature < tolerance) { // zero temperature
        auto fermiDirac = [=](double nrg) -> double {
            if (nrg < -tolerance + mu) {
                return 1.0;
            } else if (nrg <= tolerance + mu) {
                return 0.5;
            } else {
                return 0.0;
            }
        };
        return fermiDirac;
    } else { // nonzero temperature
        double inverse_temperature = 1.0 / temperature;
        auto fermiDirac = [=](double nrg) -> double {
            return 1.0 / (1.0 + exp(inverse_temperature * (nrg - mu)));
        };
        return fermiDirac;
    }
}

extern "C" {

int compute_diamagnetic_susceptibility(
  // input
  double temperature,
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
  complex128_t const * amplitudes
    = reinterpret_cast<complex128_t const*>(v_amplitudes); 
  complex128_t const * eigenvectors_up
    = reinterpret_cast<complex128_t const*> (v_eigenvectors_up); 
  complex128_t const * eigenvectors_dn
    = reinterpret_cast<complex128_t const*> (v_eigenvectors_dn);
  complex128_t * susceptibility_x_up
    = reinterpret_cast<complex128_t *> (v_susceptibility_x_up);
  complex128_t * susceptibility_x_dn
    = reinterpret_cast<complex128_t *> (v_susceptibility_x_dn);
  complex128_t * susceptibility_y_up
    = reinterpret_cast<complex128_t *> (v_susceptibility_y_up);
  complex128_t * susceptibility_y_dn
    = reinterpret_cast<complex128_t *> (v_susceptibility_y_dn);


  std::vector<double> fermi(num_eigen);
  std::transform(eigenvalues, eigenvalues + num_eigen,
                  fermi.begin(), makeFermiDirac(temperature));

  int num_site = b1 * b2;

  auto psi_up = [&](int idx_site, int idx_eigen) -> complex128_t {
    return eigenvectors_up[idx_eigen * num_site + idx_site];
  };
  auto psi_dn = [&](int idx_site, int idx_eigen) -> complex128_t {
    return eigenvectors_dn[idx_eigen * num_site + idx_site];
  };

  auto psi_up_conj = [&](int idx_site, int idx_eigen) -> complex128_t {
    return std::conj(eigenvectors_up[idx_eigen * num_site + idx_site]);
  };
  auto psi_dn_conj = [&](int idx_site, int idx_eigen) -> complex128_t {
    return std::conj(eigenvectors_dn[idx_eigen * num_site + idx_site]);
  };


#pragma omp parallel for
  for (int idx_bond = 0 ; idx_bond < num_bond ; ++idx_bond)
  {
    int i = rowsites[idx_bond];
    int j = colsites[idx_bond];
    double dx = displacements_x[idx_bond];
    double dy = displacements_y[idx_bond];

    complex128_t bondsusc_up = 0;
    complex128_t bondsusc_dn = 0;
    auto ampl = amplitudes[idx_bond];

    for (int idx_eigen = 0 ; idx_eigen < num_eigen ; ++idx_eigen) {
      bondsusc_up += 2 * std::real(ampl * psi_up_conj(i, idx_eigen) * psi_up(j, idx_eigen))
                     * fermi[idx_eigen];
      bondsusc_dn += 2 * std::real(ampl * psi_dn_conj(i, idx_eigen) * psi_dn(j, idx_eigen))
                     * fermi[idx_eigen];
    }
    susceptibility_x_up[idx_bond] = bondsusc_up * dx;
    susceptibility_y_up[idx_bond] = bondsusc_up * dy;
    susceptibility_x_dn[idx_bond] = bondsusc_dn * dx;
    susceptibility_y_dn[idx_bond] = bondsusc_dn * dy;
  }
  return 0;
}

}


