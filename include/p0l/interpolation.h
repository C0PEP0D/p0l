#ifndef P0L_INTERPOLATION_H
#define P0L_INTERPOLATION_H
#pragma once

#include <cmath>
#include <cassert>
#include <algorithm>
#include <functional>
#include <numeric>

namespace p0l {

// Lagrange 1D

// Source : https://en.wikipedia.org/wiki/Lagrange_polynomial
template<template<typename...> class TypeContainer>
double lagrangeBasis(const TypeContainer<double>& xs, const std::size_t& i, const double& x) {
    // debug
    if (i >= xs.size()) {
        std::cout << "lagrangeBasis i is to big: i:" << i << " size:" << xs.size() << std::endl;
    }
    //if (x < xs.front() || x > xs.back()) {
    //    std::cout << "x is not in range: x:" << x << " front:" << xs.front() << " back:" << xs.back() << std::endl;
    //}
    // Init result
    double y = 1.0;
    // Compute basis
    if (x != xs[i]) {
        for (std::size_t j = 0; j < xs.size(); j++) 
        { 
            if (j != i) {
                y *= (x - xs[j]) / (xs[i] - xs[j]); 
            }
        } 
    }
    // Return result
    return y;
}

// Source : https://en.wikipedia.org/wiki/Lagrange_polynomial
// Source : https://www.geeksforgeeks.org/lagranges-interpolation/
template<template<typename...> class TypeAxisContainer, template<typename...> class TypeDataContainer, typename TypeInterpolated>
TypeInterpolated lagrange(const TypeAxisContainer<double>& xs, const TypeDataContainer<TypeInterpolated>& ys, const double& x) {
    // Init result
    TypeInterpolated result = ys[0] * lagrangeBasis(xs, 0, x);
    // Compute each term
    for (std::size_t i = 1; i < xs.size(); i++) 
    { 
        result += ys[i] * lagrangeBasis(xs, i, x); 
    }
    // Retrun result
    return result;
}

// Lagrange Vector

template<template<typename...> class TypeVectorContainer, typename TypeVector, template<typename...> class TypeDataContainer, typename TypeInterpolated, template<typename...> class TypeRef>
TypeInterpolated lagrangeVector(const TypeVectorContainer<TypeVector>& xs, const TypeDataContainer<TypeInterpolated>& ys, const TypeRef<const TypeVector>& x) {
    // Init result
    TypeInterpolated result = ys[0];
    for (std::size_t j = 1; j < ys.size(); j++) {
        for (std::size_t k = 0; k < x.size(); k++) { 
            if (xs[0][k] != xs[j][k]) {
                result *= (x[k] - xs[j][k]) / (xs[0][k] - xs[j][k]); 
            }
        }
    }
    // Compute each term
    for (std::size_t i = 1; i < ys.size(); i++) {
        TypeInterpolated term = ys[i];
        for (std::size_t j = 0; j < ys.size(); j++) {
            if (j != i) {
                for (std::size_t k = 0; k < x.size(); k++) { 
                    if (xs[i][k] != xs[j][k]) {
                        term *= (x[k] - xs[j][k]) / (xs[i][k] - xs[j][k]); 
                    }
                }
            }
        }
        // Add term to result
        result += term;
    }
    // Retrun result
    return result;
}

// Lagrange Mesh

template<typename TypeMesh, template<typename ...> class TypeDataContainer, typename TypeInterpolated, typename TypeVector, template<typename...> class TypeRef>
TypeInterpolated lagrangeMesh(const TypeMesh& mesh, const TypeDataContainer<TypeInterpolated>& ys, const TypeRef<const TypeVector>& x) {
    // Compute mesh grid
    std::vector<std::vector<double>> grid(x.size());
    for(std::size_t i = 0; i < mesh.n.size(); i++) {
        grid[i].resize(mesh.n[i]);
        std::vector<int> ijk = {0, 0, 0};
        for(std::size_t j = 0; j < mesh.n[i]; j++) {
            ijk[i] = j;
            grid[i][j] = mesh.x(ijk)[i];
        }
    }
    // Get indexs
    std::vector<std::size_t> indexs = mesh.indexs();
    // Init result
    TypeInterpolated result = ys[indexs.back()];
    std::vector<int> ijk = mesh.ijk(indexs.back());
    for(std::size_t k = 0; k < x.size(); k++) {
        result *= lagrangeBasis(grid[k], ijk[k], x[k]); 
    }
    indexs.pop_back();
    // Compute each term
    for (const auto& cellIndex : indexs) { 
        TypeInterpolated term = ys[cellIndex];
        // Get ijk
        std::vector<int> ijk = mesh.ijk(cellIndex);
        // Compute lagrange interpolation
        for(std::size_t k = 0; k < x.size(); k++) {
            term *= lagrangeBasis(grid[k], ijk[k], x[k]); 
        }
        result += term;
    }
    // Retrun result
    return result;
}

template<typename TypeMesh, template<typename ...> class TypeContainer, typename TypeInterpolated, typename TypeVector, template<typename...> class TypeRef, typename TypeSubMesh>
TypeInterpolated lagrangeMesh(const std::shared_ptr<TypeMesh>& sMesh, const TypeContainer<TypeInterpolated>& ys, const TypeRef<const TypeVector>& x, const std::size_t& n, const bool periodic = true) {
    // Compute offset
    std::vector<int> offset = sMesh->ijk(x);
    for(auto& o : offset) {
        o -= n/2;
    }
    // If not periodic then be sure you are inside
    if(not periodic) {
        for(std::size_t i = 0; i < offset.size(); i++) {
            if(offset[i] < 0) {
                offset[i] = 0;
            }
            else if(offset[i] > sMesh->n[i] - n) {
                offset[i] = sMesh->n[i] - n;
            }
        }
    }
    // Build subMesh and compute
    return lagrangeMesh<TypeSubMesh, TypeContainer, TypeInterpolated, TypeVector, TypeRef>(TypeSubMesh(std::vector<std::size_t>(x.size(), n), offset, sMesh), ys, x);
}

}

#endif