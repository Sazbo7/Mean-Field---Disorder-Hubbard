#ifndef _SUSC_H_
#define _SUSC_H_

extern "C" {

int compute_diamagnetic_susceptibility(
  // input
  double temperature,
  int num_bond,
  int const * rowsites,
  int const * colsites,
  void const * amplitudes,
  double const * displacements_x,
  double const * displacements_y,
  int num_eigen,
  int b1, // number of sites in direction 1
  int b2, // number of sites in direction 2
  double const * eigenvalues,
  void const* eigenvectors_up,
  void const* eigenvectors_dn,
  // output
  void * susceptibility_x_up,
  void * susceptibility_x_dn,
  void * susceptibility_y_up,
  void * susceptibility_y_dn
);

int compute_paramagnetic_susceptibility(
  // input
  double temperature,
  double frequency,
  double broadening,
  int num_bond,
  int const * rowsites,
  int const * colsites,
  void const * amplitudes,
  double const * displacements_x,
  double const * displacements_y,
  int num_eigen,
  int b1, // number of sites in direction 1
  int b2, // number of sites in direction 2
  double const * eigenvalues,
  void const* eigenvectors_up,
  void const* eigenvectors_dn,
  // output
  void * susceptibility_x_up,
  void * susceptibility_x_dn,
  void * susceptibility_y_up,
  void * susceptibility_y_dn
);


int compute_paramagnetic_susceptibility_naive(
  // input
  double temperature,
  double frequency,
  double broadening,
  int num_bond,
  int const * rowsites,
  int const * colsites,
  void const * amplitudes,
  double const * displacements_x,
  double const * displacements_y,
  int num_eigen,
  int b1, // number of sites in direction 1
  int b2, // number of sites in direction 2
  double const * eigenvalues,
  void const* eigenvectors_up,
  void const* eigenvectors_dn,
  // output
  void * susceptibility_x_up,
  void * susceptibility_x_dn,
  void * susceptibility_y_up,
  void * susceptibility_y_dn
);

}

#endif _SUSC_H_