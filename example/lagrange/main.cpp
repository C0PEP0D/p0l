// std includes
#include <vector>
#include <iostream>
#include <cmath> // cos
// thirdparty includes
#include <Eigen/Dense>
// lib includes
#include "p0l/lagrange.h"
#include "p0l/polynome.h"
#include "p0l/to_string.h"

using TypeScalar = double;
// Space
const unsigned int DIM = 2;
using tSpaceVector = Eigen::Matrix<TypeScalar, DIM, 1>;
template<typename ...Args>
using tView = Eigen::Map<Args...>;
// Data points
const unsigned int N = 4;
const double DX = 0.1 * 2.0 * M_PI;

double f(const double x) {
    return std::cos(x);
}

int main() { 
	// Build axis
    std::vector<double> x(N);
    std::generate(x.begin(), x.end(), [x = -DX] () mutable { return x += DX; });
    // Data init
    std::vector<double> y(N);
    for(unsigned int index = 0; index < N; index++) {
        y[index] = f(x[index]);
    }
    // Interpolation
    std::vector<double> a = p0l::lagrange::interpolation(x.data(), x.size(), y.data(), 0.0);
    std::cout << "Interpolation polynome A: " << p0l::polynome::toString(a.data(), a.size()) << std::endl;
    std::cout << "Evaluation: A(0 DX) = " << p0l::polynome::evaluation(a.data(), a.size(), 0 * DX) << "   f(0 DX) = " << f(0 * DX) << std::endl;
    std::cout << "Evaluation: A(0.5 DX) = " << p0l::polynome::evaluation(a.data(), a.size(), 0.5 * DX) << "   f(0.5 DX) = " << f(0.5 * DX) << std::endl;
    std::cout << "Evaluation: A(1 DX) = " << p0l::polynome::evaluation(a.data(), a.size(), 1 * DX) << "   f(1 DX) = " << f(1 * DX) << std::endl;
    std::cout << "Evaluation: A(1.5 DX) = " << p0l::polynome::evaluation(a.data(), a.size(), 1.5 * DX) << "   f(1.5 DX) = " << f(1.5 * DX) << std::endl;
    std::cout << "Evaluation: A(2 DX) = " << p0l::polynome::evaluation(a.data(), a.size(), 2 * DX) << "   f(2 DX) = " << f(2 * DX) << std::endl;
    std::cout << "Evaluation: A(2.5 DX) = " << p0l::polynome::evaluation(a.data(), a.size(), 2.5 * DX) << "   f(2.5 DX) = " << f(2.5 * DX) << std::endl;
    std::cout << "Evaluation: A(3 DX) = " << p0l::polynome::evaluation(a.data(), a.size(), 3 * DX) << "   f(3 DX) = " << f(3 * DX) << std::endl;
}
