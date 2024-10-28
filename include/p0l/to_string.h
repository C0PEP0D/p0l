#ifndef P0L_TO_STRING_H
#define P0L_TO_STRING_H
#pragma once

#include <sstream>
#include <string>

namespace p0l {

namespace polynome {

    std::string toString(const double* pPolynome, const unsigned int sizePolynome) {
        std::stringstream output;
        if (sizePolynome > 0) {
            output << pPolynome[0];
            if (sizePolynome > 1) {
                if (pPolynome[1] == 1.0) {
                    output << "  +  x ";
                } else if (pPolynome[1] != 0.0) {
                    output << "  +  " << pPolynome[1] << " x ";
                }
                for (unsigned int index = 2; index < sizePolynome; ++index) {
                    if (pPolynome[index] == 1.0) {
                        output << "  +  x" << index << " ";
                    } else if (pPolynome[index] != 0.0) {
                        output << "  +  " << pPolynome[index] << " x" << index << " ";
                    }
                }
            }
        } else {
            output << 0.0;
        }
        return output.str();
    }

    namespace multivariate {
        std::string toString(const double* pPolynome, const std::vector<unsigned int>& degrees) {
            std::stringstream output;
            const unsigned int sizePolynome = p0l::polynome::multivariate::size(degrees);
            if (sizePolynome > 0) {
                output << pPolynome[0];
                if (sizePolynome > 1) {
                    for (unsigned int index = 1; index < sizePolynome; ++index) {
                        std::vector<unsigned int> ijk = p0l::polynome::multivariate::ijk(degrees, index);
                        if (pPolynome[index] != 0.0) {
                            if (pPolynome[index] == 1.0) {
                                output << "  +  ";
                                
                            } else {
                                output << "  +  " << pPolynome[index];
                            }
                            for (unsigned int variableIndex = 0; variableIndex < ijk.size(); ++variableIndex) {
                                if (ijk[variableIndex] == 1) {
                                    output << " " << char(int('a') + ((int('x') - int('a')) + variableIndex) % 26);
                                }
                                else if (ijk[variableIndex] != 0) {
                                    output << " " << char(int('x') + variableIndex) << ijk[variableIndex];
                                }
                            }
                        }
                    }
                }
            } else {
                output << 0.0;
            }
            return output.str();
        }  
    }
}

}

#endif
