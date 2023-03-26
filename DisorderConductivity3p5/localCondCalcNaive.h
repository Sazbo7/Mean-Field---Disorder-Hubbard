#pragma once

namespace localcondcalc {

using namespace std;

namespace naive {

/*! Get Lindhard Coefficients

  @param nrgIn energies of eigenstates (array of size numStates)
  @param omega real part of frequency
  @param broadening imaginary part of frequency
  @param temp temperature
  @param mu chemical potential

  @return matrix whose (n,m) coefficient is
          (f_m - f_n) / ( omega + i eta + E_m - E_n)
*/
inline
vector <complex<double>> getLindCoefs(const vector<double> &nrgIn,
                                      double omega, double broadening, double temp, double mu) {
    complex<double> imNrg(omega, broadening);
    int numStates = nrgIn.size();

    vector<double> fermiIn(numStates, 0.0);
    double tolerance = ::sqrt(std::numeric_limits<double>::epsilon());

    if (temp < tolerance) { // zero temperature
        auto fermi = [=](double nrg) -> double {
          if (nrg < -tolerance + mu) {
              return 1.0;
          } else if (nrg <= tolerance + mu) {
              return 0.5;
          } else {
              return 0.0;
          }
        };
        for (int nn = 0; nn < numStates; ++nn) {
            fermiIn[nn] = fermi(nrgIn[nn]);
        }
    } else { // nonzero temperature
        double inv_temp = 1.0 / temp;
        auto fermi = [=](double nrg) -> double {
          return 1.0 / (1.0 + exp(inv_temp * (nrg - mu)));
        };
        for (int nn = 0; nn < numStates; ++nn) {
            fermiIn[nn] = fermi(nrgIn[nn]);
        }
    }

    vector <complex<double>> output;
    output.reserve(numStates * numStates);
    for (int nn = 0; nn < numStates; ++nn) {
        for (int mm = 0; mm < numStates; ++mm) {
            output.push_back((fermiIn[mm] - fermiIn[nn])
                                 / (imNrg + nrgIn[mm] - nrgIn[nn]));
        }
    }

    return output;
}

/*! Get Bond-Bond susceptibility Lambda(i,j,k,l)

  @param statesIn state vectors (matrix of size (numStates x num_sites))
  @param lindCoefIn Lindhard coefficient
                     (square matrix of size (numStates x numStates))
  @param siteI
  @param siteJ
  @param siteK
  @param siteL
  @return Lambda(i,j,k,l)
*/
complex<double> getLocalSuscept(const vector <complex<double>> &statesIn,
                                const vector <complex<double>> &lindCoefIn,
                                int siteI, int siteJ, int siteK, int siteL) {
    const complex<double> im(0.0, 1.0);
    complex<double> output(0.0, 0.0);
    int numStates = sqrt(lindCoefIn.size());
    int num_sites = statesIn.size() / numStates;

#if DEBUG
    std::cout << "siteI = " << siteI << std::endl;
    std::cout << "siteJ = " << siteJ << std::endl;
    std::cout << "siteK = " << siteK << std::endl;
    std::cout << "siteL = " << siteL << std::endl;

    std::cout << "stateIn.size() = " << statesIn.size() << std::endl;
    std::cout << "lindCoefIn.size() = " << lindCoefIn.size() << std::endl;
    std::cout << "numStates = " << numStates << std::endl;
#endif

    auto psi = [&](int n, int i) -> complex<double> {
      return statesIn[n * num_sites + i];
    };
    auto psiC = [&](int n, int i) -> complex<double> {
      return conj(statesIn[n * num_sites + i]);
    };
    auto lindFunc = [&](int n, int m) -> complex<double> {
      return lindCoefIn[n * numStates + m];
    };

    for (int nn = 0; nn < numStates; ++nn) {
        for (int mm = 0; mm < numStates; ++mm) {
            auto Jij = im * (psiC(nn, siteI) * psi(mm, siteJ) - psiC(nn, siteJ) * psi(mm, siteI));
            auto Jkl = im * (psiC(mm, siteK) * psi(nn, siteL) - psiC(mm, siteL) * psi(nn, siteK));
            output += Jij * Jkl * lindFunc(nn, mm);
        }
    }
    output /= static_cast<double>(-numStates);
    return output;
}

complex<double> getLocalSuscept2(const vector <complex<double>> &statesIn,
                                 const vector <complex<double>> &lindCoefIn,
                                 int siteI, int siteJ, int siteK, int siteL) {
    const complex<double> im(0.0, 1.0);
    complex<double> output(0.0, 0.0);
    int numStates = sqrt(lindCoefIn.size());
    int num_sites = statesIn.size() / numStates;

#if 1
    std::cout << "siteI = " << siteI << std::endl;
    std::cout << "siteJ = " << siteJ << std::endl;
    std::cout << "siteK = " << siteK << std::endl;
    std::cout << "siteL = " << siteL << std::endl;

    std::cout << "stateIn.size() = " << statesIn.size() << std::endl;
    std::cout << "lindCoefIn.size() = " << lindCoefIn.size() << std::endl;
    std::cout << "numStates = " << numStates << std::endl;
#endif

    auto psi = [&](int n, int i) -> complex<double> {
      return statesIn[n * num_sites + i];
    };
    auto psiC = [&](int n, int i) -> complex<double> {
      return conj(statesIn[n * num_sites + i]);
    };
    auto lindFunc = [&](int n, int m) -> complex<double> {
      return lindCoefIn[n * numStates + m];
    };

    for (int nn = 0; nn < numStates; ++nn) {
        for (int mm = nn + 1; mm < numStates; ++mm) {
            auto Jij = psiC(nn, siteI) * psi(mm, siteJ) - psiC(nn, siteJ) * psi(mm, siteI);
            auto Jkl = psiC(nn, siteK) * psi(mm, siteL) - psiC(nn, siteL) * psi(mm, siteK);
            auto JJ = Jij * conj(Jkl);
            output += JJ * lindFunc(nn, mm);
            output += conj(JJ) * lindFunc(mm, nn);
        }
    }
    output /= static_cast<double>(numStates);
    return output;
}

/*! Get bond correlation function

  Lambda(i, j, mu)

  @param bondList list of bonds.
         elements of bondlist are tuples with
         (siteI,siteJ,xBondDirection,yBondDirection)
  @param statesIn a square matrix of state vectors
  @param lindCoefIn Lindhard coefficients
  @param siteI
  @param siteJ
  @return (Lambda(i,j,x), Lambda(i,j,y))
*/
inline
tuple <complex<double>, complex<double>>
getBondCorrFunc(const vector <tuple<int, int, int, int>> &bondList,
                const vector <complex<double>> &statesIn,
                const vector <complex<double>> &lindCoefIn,
                int siteI, int siteJ) {
    complex<double> outputX(0.0, 0.0), outputY(0.0, 0.0);

    for (auto const &bond : bondList) {
        auto siteK = get<0>(bond);
        auto siteL = get<1>(bond);
        auto vecX = get<2>(bond);
        auto vecY = get<3>(bond);
        double vX = static_cast<double>(vecX);
        double vY = static_cast<double>(vecY);

        auto localSuscept = getLocalSuscept(statesIn, lindCoefIn,
                                            siteI, siteJ, siteK, siteL);
        outputX += localSuscept * vX;
        outputY += localSuscept * vY;
    }
    return make_tuple(outputX, outputY);
}

inline
tuple <complex<double>, complex<double>>
getBondCorrFunc2(const vector <tuple<int, int, int, int>> &bondList,
                 const vector <complex<double>> &statesIn,
                 const vector <complex<double>> &lindCoefIn,
                 int siteI, int siteJ) {
    complex<double> outputX(0.0, 0.0), outputY(0.0, 0.0);

    for (auto const &bond : bondList) {
        auto siteK = get<0>(bond);
        auto siteL = get<1>(bond);
        auto vecX = get<2>(bond);
        auto vecY = get<3>(bond);
        double vX = static_cast<double>(vecX);
        double vY = static_cast<double>(vecY);

        auto localSuscept = getLocalSuscept2(statesIn, lindCoefIn,
                                             siteI, siteJ, siteK, siteL);
        outputX += localSuscept * vX;
        outputY += localSuscept * vY;
    }
    return make_tuple(outputX, outputY);
}

} // namespace naive

} // namespace localcondcalc