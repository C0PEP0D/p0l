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
#include "v0l/bin/file_data.h"

using TypeScalar = double;
// Space
const unsigned int DIM = 3;
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
const std::size_t order = 4; // interpolation order

double f(const double x, const double y) {
    return std::cos(x) + std::sin(y);
}

int main() { 
    // Init
    v0l::FileData<float> data("../data/v.vtk", 0);
    // build length
    std::vector<double> lengths(data.meta.spacing.size());
    TypeVector origin;
    for(unsigned int i = 0; i < lengths.size(); i++) {
        lengths[i] = data.meta.spacing[i] * data.meta.dimensions[i];
        origin[i] = data.meta.origin[i];
    }
    // build create mesh
    std::shared_ptr<TypeMeshStructured> sMesh = std::make_shared<TypeMeshStructuredRegular>(data.meta.dimensions, lengths, origin);
    // interpolate
    std::cout << "interpolated value: " << p0l::lagrangeMesh<TypeMeshStructured, v0l::FileData, float, TypeVector, TypeRef, TypeMeshStructuredSub>(sMesh, data, TypeVector::Random() * 0.5, order + 1) << std::endl;
}
