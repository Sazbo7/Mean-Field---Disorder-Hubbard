#pragma once
//===========================================================================//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//===========================================================================//
//
//  localCondCalc.h
//
//  Programmers:    Tim McCormick
//
//  Date:       June 11, 2017
//
//  Purpose:
//
//  Major Modifications:
//
//===========================================================================//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//===========================================================================//

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//                            Include Statements
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

#include <cassert>
#include <cmath>
#include <complex>
#include <limits>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <functional>
#include <algorithm>
#include <chrono>

#define SHOWDEBUG(x) { std::cout << (#x) << " = " << (x) << std::endl; }


namespace aux{
template<std::size_t...> struct seq{};

template<std::size_t N, std::size_t... Is>
struct gen_seq : gen_seq<N-1, N-1, Is...>{};

template<std::size_t... Is>
struct gen_seq<0, Is...> : seq<Is...>{};

template<class Ch, class Tr, class Tuple, std::size_t... Is>
void print_tuple(std::basic_ostream<Ch,Tr>& os, Tuple const& t, seq<Is...>){
  using swallow = int[];
  (void)swallow{0, (void(os << (Is == 0? "" : ", ") << std::get<Is>(t)), 0)...};
}
} // aux::

template<class Ch, class Tr, class... Args>
auto operator<<(std::basic_ostream<Ch, Tr>& os, std::tuple<Args...> const& t)
    -> std::basic_ostream<Ch, Tr>&
{
  os << "(";
  aux::print_tuple(os, t, aux::gen_seq<sizeof...(Args)>());
  return os << ")";
}


namespace localcondcalc {

using namespace std;

template <typename T>
void vecPrint( const vector<T>& input )
{
    for(auto const & v : input)
    {
        cout << v << " ";
    }
    cout << endl;
}


/*! Zip two vectors of type T to a vector of complex<T>
 *
 * @tparam T
 * @param reIn
 * @param imIn
 * @return
 */
template <typename T>
std::vector<std::complex<T>>
convertToCompVect(const std::vector<T>& real, const std::vector<T>& imag )
{
    std::size_t len = std::min(real.size(), imag.size());
    std::vector< std::complex<T> > output;
    output.reserve(len);
    for (std::size_t ii = 0 ; ii < len ; ++ii ) {
        output.push_back(std::complex<T>(real[ii], imag[ii]));
    }
    return output;
}


/*! Zip two vectors of type T to a vector of complex<T>
 *
 * @tparam T
 * @param reIn
 * @param imIn
 * @return
 */
template <typename T>
std::vector<std::complex<T>>
zipComplex(std::vector<T> const & real, std::vector<T> const & imag )
{
    std::size_t len = std::min(real.size(), imag.size());
    std::vector< std::complex<T> > output;
    output.reserve(len);
    for (std::size_t i = 0 ; i < len ; ++i ) {
        output.push_back(std::complex<T>(real[i], imag[i]));
    }
    return output;
}


//! Default tolerance for double type
static const double kDoubleTolerance = std::sqrt(std::numeric_limits<double>::epsilon());


/*! Make a Fermi-Dirac distribution function
 *
 * @param temperature
 * @param mu
 * @param tolerance
 * @return Fermi-Dirac distribution function f
 */
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



/*! Compute static Lindhard coefficients
 *
 * Lindhard(n,m) = ( f(m) - f(n) ) / ( E(m) - E(n) )
 *
 * |  n  |  m  |  Lindhard(n,m)  |
 * |:---:|:---:|:---------------:|
 * |  0  |  3  |  Lindhard(0,3)  |
 * |  0  |  4  |  Lindhard(0,4)  |
 * | ... | ... |      ...        |
 *
 * Keep the ones with `| Lindhard(n,m) | > tolerance`
 * Also, throw away pairs with too similar energies even though Lindhard is large.
 * This is because we are only interested in the dynamic Lindhard function,
 * with non-zero imaginary frequency.
 *
 * @param eigenEnergies Energies of the eigenstates
 * @param temperature Temperature of the system (for Fermi-Dirac function)
 * @param mu Chemical potential
 * @param tolerance Tolerance. Defaults to kDoubleTolerance
 * @return  Table of (n, m, Lindhard(n,m) )
 */
inline
vector<tuple<int, int, double>>
getStaticLindhards(
    const vector<double>& eigenEnergies,
    double temperature,
    double mu=0.0,
    double tolerance=kDoubleTolerance)
{
    ssize_t numStates = eigenEnergies.size();
    vector<double> fermiValues(numStates, 0.0);
    transform(eigenEnergies.begin(),
              eigenEnergies.end(), 
              fermiValues.begin(),
              makeFermiDirac(temperature, mu, tolerance));
    vector<tuple<int, int, double>> output;
    output.reserve(numStates*(numStates-1)/2);
    for (int nn = 0 ; nn < numStates ; ++nn ) {
        for (int mm = nn+1 ; mm < numStates ; ++mm ) {
            if (abs(eigenEnergies[mm] - eigenEnergies[nn]) < tolerance) {
                continue;
            }
            auto val = (fermiValues[mm] - fermiValues[nn])
                       / (eigenEnergies[mm] - eigenEnergies[nn]);
            if (abs(val) > tolerance) {
                output.push_back(make_tuple(nn, mm, val));
            }
        }
    }
    output.shrink_to_fit();
    return output;
}


/*! Compute dynamic Lindhard function
 *
 * Compute dynmic Lindhard functions for eigenstate pairs whose
 * static Lindhard exists (i.e. non-zero)
 *
 *       f(m) - f(n)            f(n) - f(m)
 * ( ------------------ , ------------------ )
 *   ω+iη + E(m) - E(n)    ω+iη + E(n) - E(m)
 * = (Lindhard^+(n,m) , Lindhard^-(n,m))
 * Here, n < m.
 *
 * @param eigenEnergies
 * @param staticLindhards
 * @param frequency
 * @param broadening
 * @param temperature
 * @param mu
 * @param tolerance
 * @return List of
 */
inline
vector<tuple<complex<double>, complex<double>>>
getDynamicLindhards(
    const vector<double>& eigenEnergies,
    const vector<tuple<int, int, double>> staticLindhards,
    double frequency,
    double broadening,
    double temperature,
    double mu=0.0,
    double tolerance=kDoubleTolerance)
{
    complex<double> omega(frequency, broadening);
    vector<double> fermiValues(eigenEnergies.size());
    transform(eigenEnergies.begin(),
              eigenEnergies.end(),
              fermiValues.begin(),
              makeFermiDirac(temperature, mu, tolerance));

    vector<tuple<complex<double>, complex<double>>> output;
    output.reserve(staticLindhards.size());

    auto computeDynamicLindhard
      = [&](const tuple<int, int, double> &sl)
          -> tuple<complex<double>, complex<double>> {
        auto nn = get<0>(sl);
        auto mm = get<1>(sl);
        auto val1 = (fermiValues[mm] - fermiValues[nn])
            / (omega + eigenEnergies[mm] - eigenEnergies[nn]);
        auto val2 = (fermiValues[nn] - fermiValues[mm])
            / (omega + eigenEnergies[nn] - eigenEnergies[mm]);
        return make_tuple(val1, val2);
    };

    transform(staticLindhards.begin(),
              staticLindhards.end(),
              back_inserter(output),
              computeDynamicLindhard);

    return output;
}

/*! Compute current elements <n| J(ij) |m>
 *
 * Compute only for eigenstate pairs with nonzero static Lindhard.
 *
 * @param bonds
 * @param staticLindhards
 * @param states
 * @param numStates
 * @return
 */
inline
vector<complex<double>>
getCurrentElements(
    const vector<tuple<int, int, int, int>>& bonds,
    const vector<tuple<int, int, double>>& staticLindhards,
    const vector<complex<double>>& states,
    int numStates)
{
    #ifdef BENCHMARK
    auto startTime = chrono::steady_clock::now();
    #endif
    static const complex<double> I(0.0, 1.0);

    ssize_t numBonds = bonds.size();
    ssize_t numLindhards = staticLindhards.size();
    ssize_t numSites = states.size() / numStates;

    assert(numStates > 0);
    assert(numStates <= states.size());
    assert(states.size() % numStates == 0);
    assert(numLindhards <= numStates * (numStates-1) / 2);

    auto psi = [&](ssize_t n, ssize_t i) -> complex<double> {
        assert(n >= 0 && i >= 0);
        assert(i < numSites);
        assert(n < numStates);
        return states[n*numSites + i];
    };
    auto psiC = [&](ssize_t n, ssize_t i) -> complex<double> {
        assert(n >= 0 && i >= 0);
        assert(i < numSites);
        assert(n < numStates);
        return conj( states[n*numSites + i] );
    };

    vector<complex<double>> out(numBonds * numLindhards);

#pragma omp parallel for
    for (ssize_t idxBond = 0 ; idxBond < numBonds ; ++idxBond) {
        auto idxSite1 = get<0>(bonds[idxBond]);
        auto idxSite2 = get<1>(bonds[idxBond]);

        for (ssize_t idxLindhard = 0 ; idxLindhard < numLindhards ; ++idxLindhard) {
            auto const & lc = staticLindhards[idxLindhard];
            ssize_t idxState1 = get<0>(lc);
            ssize_t idxState2 = get<1>(lc);
            auto Jij = (psiC(idxState1, idxSite1)*psi(idxState2, idxSite2)
                        - psiC(idxState1, idxSite2)*psi(idxState2, idxSite1)) * I;
            out[idxBond * numLindhards + idxLindhard] = Jij;
        }
    }
    #ifdef BENCHMARK
    auto endTime = chrono::steady_clock::now();
    cout << "getCurrentElements : ";
    cout << chrono::duration<double, milli>(endTime - startTime).count()
         << " ms"
         << std::endl;
    #endif
    return out;
}

/*! Compute Λ(ij; kl)
 *
 * \f[
 *   \Lambda(ij; kl)
 *     = \sum_{n<m}
 *       \langle n \vert J(ij) \vert m \rangle
 *       \langle n \vert J(kl) \vert m \rangle^*
 *       Lindhard^+(n,m)
 *       +
 *       \langle n \vert J(ij) \vert m \rangle^*
 *       \langle n \vert J(kl) \vert m \rangle
 *       Lindhard^-(n,m)
 * \f]
 *
 * @param dynamicLindhards List of Lindhard(n,m)^+ and Lindhard(n,m)^-
 * @param currentElements  List of <n|J(ij)|m>
 * @param numStates Number of states.
 * @param idxBond1 Bond "ij"
 * @param idxBond2 Bond "kl"
 * @return
 */
inline
complex<double>
getBondBondSusceptibility(
    const vector<tuple<complex<double>, complex<double>>> & dynamicLindhards,
    const vector<complex<double>>& currentElements,
    int numStates,
    int idxBond1,
    int idxBond2)
{
    ssize_t numLindhards = dynamicLindhards.size();
    auto bondBlockBegin1 = idxBond1 * numLindhards;
    auto bondBlockBegin2 = idxBond2 * numLindhards;

    assert(currentElements.size() >= numLindhards);
    assert(currentElements.size() % numLindhards == 0);

    complex<double> output(0.0, 0.0);

    for (ssize_t idxLindhard = 0 ; idxLindhard < numLindhards ; ++idxLindhard) {
        auto const & dl = dynamicLindhards[idxLindhard];
        auto lcVal1 = get<0>(dl);
        auto lcVal2 = get<1>(dl);

        auto Jij = currentElements[bondBlockBegin1 + idxLindhard];
        auto Jkl = currentElements[bondBlockBegin2 + idxLindhard];
        auto JJ = Jij * conj(Jkl);
        output += JJ * lcVal1 + conj(JJ) * lcVal2;
    }

    output /= static_cast<double>(-numStates);
    return output;
}


/*! Compute Λ(ij; x) and Λ(ij; y)
 *
 * @param bonds List of bonds (i,j, dx, dy)
 * @param dynamicLindhards List of Lindhard^+ and Lindhard^-
 * @param currentElements List of <n|J(ij)|m>
 * @param numStates Number of eigenstates
 * @param idxBond1 Bond "ij"
 * @return Tuple ( Λ(ij; x), Λ(ij; y) )
 */
inline
tuple<complex<double>, complex<double>>
getBondFieldSusceptibility(
    const vector<tuple<int, int, int, int>>& bonds,
    const vector<tuple<complex<double>, complex<double>>> & dynamicLindhards,
    const vector<complex<double>>& currentElements,
    int numStates,
    int idxBond1)
{
    #ifdef BENCHMARK
    auto startTime = std::chrono::steady_clock::now();
    #endif

    complex<double> outputX(0.0, 0.0), outputY(0.0, 0.0);
    int numBonds = bonds.size();
    for (int idxBond2 = 0 ; idxBond2 < numBonds ; ++idxBond2) {
        auto bond2 = bonds[idxBond2];
        auto dX = static_cast<double>(get<2>(bond2));
        auto dY = static_cast<double>(get<3>(bond2));

        auto bbs = getBondBondSusceptibility(dynamicLindhards,
                                             currentElements,
                                             numStates,
                                             idxBond1,
                                             idxBond2);
        outputX += bbs * dX;
        outputY += bbs * dY;
    }

    #ifdef BENCHMARK
    auto endTime = chrono::steady_clock::now();
    cout << "getBondFieldSusceptibility : ";
    cout << chrono::duration<double, milli>(endTime - startTime).count()
         << " ms"
         << std::endl;
    #endif

    return make_tuple(outputX, outputY);
}


} //namespace localcondcalc


#include "localCondCalcNaive.h"
