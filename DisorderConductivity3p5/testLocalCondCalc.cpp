
#include <vector>
#include <string>
#include <tuple>
#include <fstream>
#include <random>
#include <iterator>
#include <algorithm>
#include <iomanip>

#include "localCondCalc.h"

using namespace std;

// Imports a single array of doubles
vector<double> importVector( const string& nameIn )
{
    ifstream inFile( nameIn.c_str() );
    vector<double> output;
    double temp;
    while ( inFile >> temp )
    {
        output.push_back( temp );
    }
    inFile.close();
    return output;
}


// generate a bond list with tuples with
//   (siteI,siteJ,xBondDirection,yBondDirection)
//
// The convention is
// 
// Coordinates = (Row,Column) = (X,Y)
//
// +---------------------------------> Y
// |   Column  0      1      2
// | Row
// |   0       0 ->   3 ->   6 ->
// |           |      |      |
// |           V      V      V
// |
// |   1       1 ->   4 ->   7 ->  
// |           |      |      |
// |           V      V      V
// |
// |   2       2 ->   5 ->   8 ->
// |           |      |      |
// |           V      V      V
// V
// X
//
vector< tuple<int,int,int,int> > genBondList( int num_sites )
{
    int n1 = sqrt( num_sites );
    int n2 = n1;
    assert(fabs(num_sites - n1*n1) < sqrt(std::numeric_limits<double>::epsilon()));

    vector< tuple<int,int,int,int> > output;
    output.reserve(num_sites*2);
    
    // loop through number of sites
    for ( int ii = 0 ; ii < num_sites ; ++ii )
    {
        // right
        if ( ii >= num_sites-n1 )
        {
            tuple<int,int,int,int> tempRi( ii , ii-(n1-1)*n1 , 0 , 1 );
            output.push_back(tempRi);
        }
        else
        {
            tuple<int,int,int,int> tempRi( ii , ii+n1 ,  0 , 1 );
            output.push_back(tempRi);
        }

        // down
        if ( ii%n1 == n1-1 )
        {
            tuple<int,int,int,int> tempUp( ii , ii-n1+1 , 1 , 0 );
            output.push_back(tempUp);
        }
        else
        {
            tuple<int,int,int,int> tempUp( ii , ii+1 , 1 , 0 );
            output.push_back(tempUp);
        }
    }
    return output;
}


int testNaive ()
{
    using namespace localcondcalc;
    using namespace localcondcalc::naive;

    complex<double> imI(0.0,1.0);

    string filename1  = "statesTestRe1.txt";
    string filename2 = "statesTestIm1.txt";

    // Imaginary and real flipped from previous
    vector<double> nrgTest = importVector( "nrgTest1.txt" );

    if (nrgTest.empty()) {
        cout << "Empty file" << endl;
        exit(1);
    }

    vector<double> tester1 = importVector( filename1 );
    vector<double> tester2 = importVector( filename2 );

    if (tester1.empty() || tester2.empty()) {
        cout << "Empty files" << endl;
        exit(1);
    } else if (tester1.size() != tester2.size()) {
        cout << "Length mismatch" << endl;
        exit(1);
    }

    vector< complex<double> > compTest = convertToCompVect( tester1 , tester2 );

    vector< tuple<int,int,int,int> > bondTest = genBondList( 9 );

    vector<complex<double> > lindTest = getLindCoefs( nrgTest , 0.0001, 0.000001 , 0.001, 0.0);

    cout << "lindTest = " << lindTest[0] << endl;

    // test bond list
    /*
    vector<tuple<int,int> > testTuple;
    for ( int tt = 0 ; tt < bondTest.size() ; tt++ )
    {
        cout << get<0>(bondTest[tt]) << " " << get<1>(bondTest[tt]) << " "
             << get<2>(bondTest[tt]) << " " << get<3>(bondTest[tt]) << endl;
    }*/

    //cout << compTest[0][1] << endl;
    cout << getBondCorrFunc( bondTest , compTest , lindTest , 0 , 1) << endl;
    //cout << getLocalSuscept( compTest , nrgTest , 0 , 1 , 3 ,4, 0.001 , 0.00001, 0.1 , 0.0 ) << endl;

    return 0;
}


int testOptimized()
{
    using namespace std;
    using namespace localcondcalc;

    wcout.precision(std::numeric_limits<double>::max_digits10);

    complex<double> imI(0.0,1.0);
    int b1 = 3;
    int b2 = 3;

    int numSites = b1 * b2;
    int numStates = 2 * b1 * b2; // == Hilbert space dimension
    int numBonds = 2 * numSites;

    double frequency = 0.01;
    double broadening = 0.001;
    double temperature = 0.01;
    double mu = 0.0;

    mt19937 randomGenerator(1);
    uniform_real_distribution<double> uniformDistribution(-1.0, 1.0);

    vector<double> eigenenergies(numStates);
    generate(eigenenergies.begin(), eigenenergies.end(),
             [&]() -> double { return uniformDistribution(randomGenerator); });

    vector<complex<double>> eigenvectors(numStates * numSites); // Yes. Sites!.
    generate(eigenvectors.begin(), eigenvectors.end(),
             [&]() -> std::complex<double> {
                 auto r = uniformDistribution(randomGenerator);
                 auto i = uniformDistribution(randomGenerator);
                 return complex<double>(r, i);
             });

    auto bonds = genBondList( numSites );

    auto lindCoefs = naive::getLindCoefs(eigenenergies, frequency, broadening, temperature, mu );

    auto staticLindhards = getStaticLindhards(eigenenergies, temperature, mu);
    auto dynamicLindhards = getDynamicLindhards(eigenenergies,
                                                staticLindhards,
                                                frequency,
                                                broadening,
                                                temperature,
                                                mu);

    auto currentElements = getCurrentElements(bonds, staticLindhards, eigenvectors, numStates);

    // Test if new algorithm gives same results as the naive one.

    for (int idxBond1 = 0 ; idxBond1 < numBonds ; ++idxBond1) {
        for (int idxBond2 = 0 ; idxBond2 < numBonds ; ++idxBond2) {
            auto localSuscept = getBondBondSusceptibility(dynamicLindhards,
                                                          currentElements,
                                                          numStates,
                                                          idxBond1, idxBond2);

            auto localSuscept2 = naive::getLocalSuscept(eigenvectors,
                                                        lindCoefs,
                                                        get<0>(bonds[idxBond1]),
                                                        get<1>(bonds[idxBond1]),
                                                        get<0>(bonds[idxBond2]),
                                                        get<1>(bonds[idxBond2]));

            wcout << "localSuscept (" <<  idxBond1 << ", " << idxBond2 << ") = ";
            wcout << localSuscept << endl;
            wcout << "localSuscept2(" <<  idxBond1 << ", " << idxBond2 << ") = ";
            wcout << localSuscept2 << endl;
            wcout << endl;

        }
    }

    for (int idxBond1 = 0 ; idxBond1 < numBonds ; ++idxBond1) {
        auto bondCorrFunc = getBondFieldSusceptibility(bonds,
                                                       dynamicLindhards,
                                                       currentElements,
                                                       numStates,
                                                       idxBond1);

        auto bondCorrFunc2 = naive::getBondCorrFunc(bonds,
                                                    eigenvectors,
                                                    lindCoefs,
                                                    get<0>(bonds[idxBond1]),
                                                    get<1>(bonds[idxBond1]));
        wcout << "bondCorrFunc (" << idxBond1 << ") = " << bondCorrFunc << endl;
        wcout << "bondCorrFunc2(" << idxBond1 << ") = " << bondCorrFunc2 << endl;
    }
    return 0;
    
}

int main(int argc, char** argv)
{
    return testOptimized();
}

