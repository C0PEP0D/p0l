// std includes
#include <vector>
#include <iostream>
#include <random> // random_device, default_random_engine, uniform_real_distribution
#include <memory> // shared_ptr
#include <cmath> // sin, cos
// thirdparty includes
#include <Eigen/Dense>
// lib includes
#include "m0sh/non_uniform.h"
#include "p0l/polynome.h"
#include "p0l/lagrange.h"
#include "p0l/to_string.h"

using TypeScalar = double;
// space
const unsigned int DIM = 2;
using tSpaceVector = Eigen::Matrix<TypeScalar, DIM, 1>;
template<typename ...Args>
using tView = Eigen::Map<Args...>;
// mesh
using tMesh = m0sh::NonUniform<tSpaceVector, tView, DIM>;
// mesh parameters
const std::size_t n = 101;
const double dx = 2.0 * M_PI / (n - 1);
// random setup
std::random_device r;
std::default_random_engine e(r());
std::uniform_real_distribution<double> uniform(0.0, n * dx); // 4 * n * dx
// interpolation parameter
const unsigned int order = 2;

double f(const double x, const double y) {
    return std::cos(x) + std::sin(y);
}

void test(const std::vector<double>* pAxes, const double* pZValues) {
    // select x
    const tSpaceVector x = {uniform(e), uniform(e)};
    // interpolate
    std::vector<double> a = p0l::lagrange::interpolationMeshPoint<tMesh>(pAxes, pZValues, 0.0, x.data(), order, true);
    std::vector<unsigned int> aDegrees = {order, order}; // TODO: better
    std::cout << "\nInterpolation Polynome:\n" << p0l::polynome::multivariate::toString(a.data(), aDegrees) << "\n" << std::endl;
    // compute error
    const double analy = f(x[0], x[1]);
    double interp = p0l::polynome::multivariate::evaluation(a.data(), aDegrees, x.data());
    double error = std::abs(analy - interp);
    double relative = error / analy;
    // print
    std::cout << "Interpolation order " << order << " x : (" << x.transpose() << ") Analy : " << analy << " Result : " << interp << " | error = " << error << " | relative = " << relative << std::endl;
}

int main() { 
	// build axis
    std::vector<double> axis(n);
    std::generate(axis.begin(), axis.end(), [x = -dx] () mutable { return x += dx; });
    // build grid
    std::vector<std::vector<double>> axes(DIM, axis);
    std::vector<unsigned int> nbPointsPerAxis = tMesh::nbPointsPerAxis(axes.data());
    // data init
    const std::vector<unsigned int> nbCellsPerAxis = tMesh::nbCellsPerAxis(nbPointsPerAxis.data());
    std::vector<double> zValues(tMesh::nb(nbCellsPerAxis.data()));
    for(unsigned int index = 0; index < zValues.size(); index++) {
        tSpaceVector x = tMesh::positionPointPeriodic(axes.data(), index);
        zValues[index] = f(x[0], x[1]);
    }
    // test
    test(axes.data(), zValues.data());
    test(axes.data(), zValues.data());
    test(axes.data(), zValues.data());
    test(axes.data(), zValues.data());
    test(axes.data(), zValues.data());
    test(axes.data(), zValues.data());
}
