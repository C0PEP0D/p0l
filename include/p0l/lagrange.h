#ifndef P0L_LAGRANGE_H
#define P0L_LAGRANGE_H
#pragma once

// #include <cmath>
// #include <cassert>
// #include <algorithm>
// #include <functional>
// #include <numeric>
#include <vector>

#include "p0l/polynome.h"

namespace p0l {

namespace lagrange {
    // Source : https://en.wikipedia.org/wiki/Lagrange_polynomial
    std::vector<double> basis(const double* pXValues, const unsigned int nbXValues, const unsigned int index) {
        std::vector<double> output = polynome::monome(0, 1.0, 0.0);
        for (unsigned int i = 0; i < nbXValues; ++i) {
            if (i != index) {
                // build new term
                std::vector<double> term = polynome::monome(1, 1.0, 0.0);
                polynome::substract(term, pXValues[i]);
                polynome::multiply(term,  1.0/(pXValues[index] - pXValues[i]));
                // multiply to output
                polynome::multiply(output, term.data(), term.size(), 0.0); 
            }
        }
        return output;
    }

    // Source : https://en.wikipedia.org/wiki/Lagrange_polynomial
    // Source : https://www.geeksforgeeks.org/lagranges-interpolation/
    template<typename tInterpolated>
    std::vector<tInterpolated> interpolation(const double* pXValues, const unsigned int nbXValues, const tInterpolated* pYValues, const tInterpolated& zero) {
        std::vector<tInterpolated> output = polynome::monome(0, zero, zero);
        for(unsigned int i = 0; i < nbXValues; ++i) {
            // build new term
            std::vector<double> term = basis(pXValues, nbXValues, i);
            polynome::multiply(term, pYValues[i]);
            // add to output
            polynome::add(output, term.data(), term.size(), zero);
        }
        return output;
    }

    template<typename tMesh>
    std::vector<unsigned int> interpolationMeshPointDegrees(const std::vector<double>* pMeshAxesPoints) {
        // get number of points per axis
        const std::vector<unsigned int> nbPointsPerAxis = tMesh::nbPointsPerAxis(pMeshAxesPoints);
        // compute degrees
        std::vector<unsigned int> output = nbPointsPerAxis;
        for (unsigned int index = 0; index < output.size(); ++index) {
            output[index] -= 1;
        }
        return output;
    }

    template<typename tMesh, typename tInterpolated>
    std::vector<tInterpolated> interpolationMeshPoint(const std::vector<double>* pMeshAxesPoints, const tInterpolated* pYValues, const tInterpolated& zero) {
        // indexs
        const std::vector<unsigned int> nbPointsPerAxis = tMesh::nbPointsPerAxis(pMeshAxesPoints);
        const std::vector<unsigned int> indexs = tMesh::indexs(nbPointsPerAxis.data());
        // init
        std::vector<unsigned int> outputDegrees(tMesh::Dim, 0);
        std::vector<tInterpolated> output = polynome::multivariate::monome(outputDegrees, zero, zero);
        for (const auto& pointIndex : indexs) {
            // init term
            std::vector<unsigned int> termDegrees(tMesh::Dim, 0);
            std::vector<tInterpolated> term = polynome::multivariate::monome(termDegrees, pYValues[pointIndex], zero);
            // get ijk corresponding to pointIndex
            std::vector<int> ijk = tMesh::ijkPoint(pMeshAxesPoints, pointIndex);
            for(std::size_t k = 0; k < tMesh::Dim; k++) {
                std::vector<double> _basis = basis(pMeshAxesPoints[k].data(), pMeshAxesPoints[k].size(), ijk[k]);
                // compute basis degrees
                std::vector<unsigned int> basisDegrees(tMesh::Dim, 0);
                basisDegrees[k] = _basis.size() - 1;
                // operate on polynome
                polynome::multivariate::multiply(term, termDegrees, _basis.data(), basisDegrees, zero);
            }
            polynome::multivariate::add(output, outputDegrees, term.data(), termDegrees, zero);
        }
        return output;
    }

    template<typename tMesh, typename tInterpolated>
    std::vector<tInterpolated> interpolationMeshPoint(const std::vector<double>* pMeshAxesPoints, const tInterpolated* pYValues, const tInterpolated& zero, const double* pX, const unsigned int order, const bool isPeriodic = true) {
        const unsigned int numberOfPoints = order + 1;
        // compute offset
        std::vector<int> offset = tMesh::ijkCell(pMeshAxesPoints, pX);
        for(int& o : offset) {
            o -= numberOfPoints/2;
        }
        if(not isPeriodic) {
            // If not periodic then be sure you are inside
            for(std::size_t i = 0; i < offset.size(); i++) {
                if(offset[i] < 0) {
                    offset[i] = 0;
                }
                else if(offset[i] > pMeshAxesPoints[i].size() - numberOfPoints) {
                    offset[i] = pMeshAxesPoints[i].size() - numberOfPoints;
                }
            }
        }
        // compute sub mesh
        std::vector<std::vector<double>> subMeshAxesPoints(tMesh::Dim, std::vector<double>(numberOfPoints));
        for(unsigned int i = 0; i < tMesh::Dim; i++) {
            const double lengthAxis = pMeshAxesPoints[i][pMeshAxesPoints[i].size() - 1] - pMeshAxesPoints[i][0];
            for(unsigned int j = 0; j < numberOfPoints; j++) {
                std::div_t div = std::div((int)(offset[i] + j), (int)(pMeshAxesPoints[i].size() - 1));
                if (div.rem < 0) {
                    div.rem += pMeshAxesPoints[i].size() - 1;
                    div.quot -= 1;
                }
                subMeshAxesPoints[i][j] = pMeshAxesPoints[i][div.rem] + lengthAxis * div.quot;
            }
        }
        // compute sub values
        std::vector<tInterpolated> subYValues(std::pow(numberOfPoints, tMesh::Dim));
        std::vector<unsigned int> subNbPointsPerAxis = tMesh::nbPointsPerAxis(subMeshAxesPoints.data());
        std::vector<unsigned int> nbPointsPerAxis = tMesh::nbPointsPerAxis(pMeshAxesPoints);
        std::vector<unsigned int> subIndexs;
        if (isPeriodic) {
            subIndexs = tMesh::subIndexsPointsPeriodic(nbPointsPerAxis.data(), subNbPointsPerAxis.data(), offset.data());
        } else {
            subIndexs = tMesh::subIndexsPoints(nbPointsPerAxis.data(), subNbPointsPerAxis.data(), offset.data());
        }
        for (unsigned int index = 0; index < subYValues.size(); index++) {
            subYValues[index] = pYValues[subIndexs[index]];
        }
        // build subMesh and compute
        return interpolationMeshPoint<tMesh, tInterpolated>(subMeshAxesPoints.data(), subYValues.data(), zero);
    }
}

}

#endif
