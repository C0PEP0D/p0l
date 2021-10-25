// std includes
#include <vector>
#include <iostream>
#include <random> // random_device, default_random_engine, uniform_real_distribution
#include <memory> // shared_ptr
#include <cmath> // sin, cos
// thirdparty includes
#include <Eigen/Dense>
// lib includes
#include "m0sh/regular.h"
#include "m0sh/structured_sub.h"
#include "p0l/interpolation.h"

using TypeScalar = double;
// Space
const unsigned int DIM = 2;
using TypeVector = Eigen::Matrix<TypeScalar, DIM, 1>;
template<typename ...Args>
using TypeRef = Eigen::Ref<Args...>;
// Mesh
template<typename ...Args>
using TypeContainer = std::vector<Args...>;
using TypeMeshStructured = m0sh::Structured<TypeVector, TypeRef, TypeContainer>;
using TypeMeshStructuredSub = m0sh::StructuredSub<TypeVector, TypeRef, TypeContainer>;
using TypeMeshStructuredRegular = m0sh::Regular<TypeVector, TypeRef, TypeContainer>;
// Data
const std::size_t np = 5;
const std::size_t n = 100;
const double l = 2 * M_PI;

double f(const double x, const double y) {
    return std::cos(x) + std::sin(y);
}

void print(const std::shared_ptr<TypeMeshStructured>& sMesh, const TypeContainer<TypeScalar> q, std::uniform_real_distribution<TypeScalar>& uniform, std::default_random_engine& e) {
    const TypeVector x = {uniform(e), uniform(e)};
    const double analy = f(x[0], x[1]);
    double interp = p0l::lagrangeMesh<TypeMeshStructured, TypeContainer, double, TypeVector, TypeRef, TypeMeshStructuredSub>(sMesh, q, x, np);
    double error = std::abs(analy - interp);
    double relative = error / analy;
    std::cout << "Interpolation using " << np << " grid points (" << std::pow(np, DIM) << " points)." << " x : \n" << x << "\n Analy : " << analy << " Result : " << interp << " | error = " << error << " | relative = " << relative << std::endl;
}

int main() { 
    // Init
    std::shared_ptr<TypeMeshStructured> sMesh = std::make_shared<TypeMeshStructuredRegular>(TypeContainer<std::size_t>(DIM, n), TypeContainer<TypeScalar>(DIM, l), TypeVector::Constant(0.0));
    TypeContainer<TypeScalar> q(sMesh->size());
    for(std::size_t cellIndex = 0; cellIndex < q.size(); cellIndex++) {
        TypeVector x = sMesh->x(cellIndex);
        q[cellIndex] = f(x[0], x[1]);
    }
    // Random setup
    std::random_device r;
    std::default_random_engine e(r());
    std::uniform_real_distribution<TypeScalar> uniform(0, l);
    // Print
    print(sMesh, q, uniform, e);
    print(sMesh, q, uniform, e);
    print(sMesh, q, uniform, e);
    print(sMesh, q, uniform, e);
    print(sMesh, q, uniform, e);
    print(sMesh, q, uniform, e);
}
