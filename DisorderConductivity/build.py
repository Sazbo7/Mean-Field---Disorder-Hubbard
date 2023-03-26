import cffi
import os, sys
import platform

cdef = r'''
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
'''

def main():
    with open("dia.cpp", "r") as inputfile:
        source_dia = inputfile.read()
    with open("para.cpp", "r") as inputfile:
        source_para = inputfile.read()

    operating_system = platform.system()
    if operating_system == 'Windows':
        extra_compile_args = []
    elif operating_system == 'Linux':
        extra_compile_args = ['-std=c++11', '-fopenmp']
        extra_compile_args += ['-DBENCHMARK=1']
    elif operating_system == 'Darwin':
        extra_compile_args = ['-std=c++11']
        extra_compile_args += ['-DBENCHMARK=1']
    else:
        raise NotImplementedError("Support for {} not implemented.".format(operating_system))

    ffibuilder = cffi.FFI()
    ffibuilder.cdef(cdef)
    ffibuilder.set_source("c_susc",
                          source_dia + source_para,
                          extra_compile_args=extra_compile_args,
                          source_extension=".cpp")

    ffibuilder.compile(verbose=True)

if __name__ == '__main__':
    main()
