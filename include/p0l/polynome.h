#ifndef P0L_POLYNOME_H
#define P0L_POLYNOME_H
#pragma once

#include <cassert>
#include <vector>

// TODO: SUM, DIFFERENCE, PRODUCT, DERIVATIVE AND PRIMITIVE (INSTANCING FUNCTIONS)

namespace p0l {

namespace polynome {

	template<typename tPolynome>
	std::vector<tPolynome> monome(const unsigned int degree, const tPolynome value, const tPolynome& zero) {
		std::vector<tPolynome> output(degree+1, zero);
		output[degree] = value;
		return output;
	}

	template<typename tA, typename tB>
	void add(std::vector<tA>& a, const tB* pB, const unsigned int sizeB, const tA& zero) {
		if (a.size() < sizeB) {
			a.resize(sizeB, zero);
		}
		for(unsigned int index = 0; index < sizeB; ++index) {
			a[index] += pB[index];
		}
	}

	template<typename tPolynome>
	void add(std::vector<tPolynome>& a, const tPolynome value) {
		a[0] += value;
	}

	template<typename tA, typename tB>
	void substract(std::vector<tA>& a, const tB* pB, const unsigned int sizeB, const tA& zero) {
		if (a.size() < sizeB) {
			a.resize(sizeB, zero);
		}
		for(unsigned int index = 0; index < sizeB; ++index) {
			a[index] -= pB[index];
		}
	}

	template<typename tPolynome>
	void substract(std::vector<tPolynome>& a, const tPolynome value) {
		a[0] -= value;
	}

	template<typename tA, typename tB>
	void multiply(std::vector<tA>& a, const tB* pB, const unsigned int sizeB, const tA& zero) { // TODO: improve, at the moment the implementation is very naïve
		std::vector<tA> c(a.size() + sizeB - 1, zero);
		for(unsigned int i = 0; i < a.size(); ++i) {
			for(unsigned int j = 0; j < sizeB; ++j) {
				c[i+j] += a[i] * pB[j];
			}
		}
		a = c;
	}

	template<typename tPolynome>
	void multiply(std::vector<tPolynome>& a, const double scalar) {
		for(unsigned int index = 0; index < a.size(); ++index) {
			a[index] *= scalar;
		}
	}
	template<typename tA>
	void raise(std::vector<tA>& a, const unsigned int power) {
		std::vector<tA> c = a;
		for (unsigned int index = 1; index < power; ++index) {
			multiply(c, a);
		}
		a = c;
	}

	template<typename tA, typename tB>
	void compose(std::vector<tA>& a, const tB* pB, const unsigned int sizeB, const tA& zero) {
		// deal with 0 and 1
		std::vector<tA> c(1, zero); // c = a[0]
		add(c, a[0]);
		std::vector<tB> d(pB, pB + sizeB); // d = b
		multiply(d, a[1]); // d = a[1] * b
		add(c, d.data(), d.size()); // c = a[0] + a[1] * b
		for(unsigned int i = 2; i < a.size(); ++i) {
			d = std::vector<tB>(pB, pB + sizeB); // d = b
			raise(d, i); // d = b^i
			multiply(d, a[i]); // d = a[i] * b^i
			add(c, d);// c = a[0] + a[1] * b + ... + a[i] * b^i
		}
		a = c;
	}

	template<typename tPolynome>
	void differentiate(std::vector<tPolynome>& a) {
		if (a.size() > 1) {
			for(unsigned int index = 0; index < a.size() - 1; ++index) {
				a[index] = (index+1) * a[index+1];
			}
			a.pop_back();
		} else {
			a[0] *= 0.0;
		}
	}

	template<typename tPolynome>
	void differentiate(std::vector<tPolynome>& a, const unsigned int order) {
		for(unsigned int index = 0; index < order; ++index) {
			differentiate(a);
		}
	}

	template<typename tPolynome>
	void integrate(std::vector<tPolynome>& a, const tPolynome& constant) {
		a.resize(a.size() + 1);
		for(unsigned int index = a.size() - 1; index > 0; --index) {
			a[index] = a[index-1] / index;
		}
		a[0] = constant;
	}

	template<typename tPolynome>
	void integrate(std::vector<tPolynome>& a, const tPolynome& constant, const unsigned int order) {
		for(unsigned int index = 0; index < order; ++index) {
			integrate(a, constant);
		}
	}

	template<typename tPolynome>
	tPolynome evaluation(const tPolynome* pA, const unsigned int sizeA, const double x) {
		tPolynome output = pA[0];
		double xPower = x;
		for(unsigned int index = 1; index < sizeA; ++index) {
			output += pA[index] * xPower;
			xPower *= x;
		}
		return output;
	}

	namespace multivariate {
		bool is_valid_ijk(const std::vector<unsigned int>& degrees, const std::vector<unsigned int>& _ijk) {
			for(unsigned int k = 0; k < degrees.size(); k++) {
				if(_ijk[k] > degrees[k]) {
					return false;
				}
			}
			return true;
		};
	
		std::vector<unsigned int> ijk(const std::vector<unsigned int>& degrees, const unsigned int _index) {
			// Compute nDiv
			std::vector<unsigned int> nDiv(degrees.size(), 1);
			for(unsigned int k = 1; k < degrees.size(); k++) {
				nDiv[k] = nDiv[k-1] * (degrees[k-1] + 1);
			}
			// Compute ijk
			std::ldiv_t dv = {0, (long int)_index};
			std::vector<unsigned int> ijk_(degrees.size());
			for(unsigned int k = 0; k < degrees.size(); k++) {
				dv = std::div(dv.rem, nDiv[degrees.size() - k - 1]);
				ijk_[degrees.size() - k - 1] = dv.quot;
			}
			return ijk_;
		};

		unsigned int index(const std::vector<unsigned int>& degrees, const std::vector<unsigned int>& _ijk) {
			if (not _ijk.empty()) {
				// Compute index
				unsigned int nTot = 1;
				unsigned int _index = _ijk[0];
				assert((_ijk[0] < degrees[0] + 1) && "ijk is not in range");
				for(unsigned int i = 1; i < _ijk.size(); i++) {
					assert((_ijk[i] < degrees[i] + 1) && "ijk is not in range");
					nTot *= (degrees[i-1] + 1);
					_index += nTot * _ijk[i];
				}
				return _index;
			} else {
				return 0;
			}
		};

		unsigned int size(const std::vector<unsigned int>& degrees) {
			unsigned int output = 1;
			for(unsigned int _index = 0; _index < degrees.size(); ++_index) {
				output *= (degrees[_index] + 1);
			}
			return output;
		};

		template<typename tPolynome>
		std::vector<tPolynome> monome(const std::vector<unsigned int>& degrees, const tPolynome& value, const tPolynome& zero) {
			const unsigned int _size = size(degrees);
			std::vector<tPolynome> output(_size, zero);
			output[_size - 1] = value;
			return output;
		}

		template<typename tA, typename tB>
		void add(std::vector<tA>& a, std::vector<unsigned int>& degreesA, const tB* pB, const std::vector<unsigned int>& degreesB, const tA& zero) {
			// create c
			std::vector<unsigned int> degreesC(std::max(degreesA.size(), degreesB.size()));
			unsigned int sizeC = 1;
			for(unsigned int _index = 0; _index < degreesC.size(); ++_index) {
				if (_index < degreesA.size() && _index < degreesB.size()) {
					degreesC[_index] = std::max(degreesA[_index], degreesB[_index]);
				} else if (_index < degreesA.size()) {
					degreesC[_index] = degreesA[_index];
				} else {
					degreesC[_index] = degreesB[_index];
				}
				sizeC *= (degreesC[_index] + 1);
			}
			std::vector<tA> c(sizeC, zero);
			// perform addition
			for(unsigned int indexC = 0; indexC < c.size(); ++indexC) {
				std::vector<unsigned int> _ijk = ijk(degreesC, indexC);
				if (is_valid_ijk(degreesA, _ijk)) {
					const unsigned int indexA = index(degreesA, _ijk);
					if(is_valid_ijk(degreesB, _ijk)) {
						const unsigned int indexB = index(degreesB, _ijk);
						c[indexC] = a[indexA] + pB[indexB];
					} else {
						c[indexC] = a[indexA];
					}
				} else if (is_valid_ijk(degreesB, _ijk)) {
					const unsigned int indexB = index(degreesB, _ijk);
					c[indexC] = pB[indexB];
				}
			}
			// store c in a
			a = c;
			degreesA = degreesC;
		}

		template<typename tPolynome>
		void add(std::vector<tPolynome>& a, std::vector<unsigned int>& degreesA, const tPolynome value) {
			a[0] += value;
		}

		template<typename tA, typename tB>
		void substract(std::vector<tA>& a, std::vector<unsigned int>& degreesA, const tB* pB, const std::vector<unsigned int>& degreesB, const tA& zero) {
			// create c
			std::vector<unsigned int> degreesC(std::max(degreesA.size(), degreesB.size()));
			unsigned int sizeC = 1;
			for(unsigned int _index = 0; _index < degreesC.size(); ++_index) {
				if (_index < degreesA.size() && _index < degreesB.size()) {
					degreesC[_index] = std::max(degreesA[_index], degreesB[_index]);
				} else if (_index < degreesA.size()) {
					degreesC[_index] = degreesA[_index];
				} else {
					degreesC[_index] = degreesB[_index];
				}
				sizeC *= (degreesC[_index] + 1);
			}
			std::vector<tA> c(sizeC, zero);
			   // perform substraction
			for(unsigned int indexC = 0; indexC < c.size(); ++indexC) {
				std::vector<unsigned int> _ijk = ijk(degreesC, indexC);
				if (is_valid_ijk(degreesA, _ijk)) {
					const unsigned int indexA = index(degreesA, _ijk);
					if(is_valid_ijk(degreesB, _ijk)) {
						const unsigned int indexB = index(degreesB, _ijk);
						c[indexC] = a[indexA] - pB[indexB];
					} else {
						c[indexC] = a[indexA];
					}
				} else if (is_valid_ijk(degreesB, _ijk)) {
					const unsigned int indexB = index(degreesB, _ijk);
					c[indexC] = -pB[indexB];
				}
			}
			// store c in a
			a = c;
			degreesA = degreesC;
		}
	
		template<typename tPolynome>
		void substract(std::vector<tPolynome>& a, std::vector<unsigned int>& degreesA, const tPolynome value) {
			a[0] -= value;
		}

		template<typename tA, typename tB>
		void multiply(std::vector<tA>& a, std::vector<unsigned int>& degreesA, const tB* pB, const std::vector<unsigned int>& degreesB, const tA& zero) {
			// create c
			std::vector<unsigned int> degreesC(std::max(degreesA.size(), degreesB.size()));
			unsigned int sizeC = 1;
			for(unsigned int _index = 0; _index < degreesC.size(); ++_index) {
				if (_index < degreesA.size() && _index < degreesB.size()) {
					degreesC[_index] = degreesA[_index] + degreesB[_index];
				} else if (_index < degreesA.size()) {
					degreesC[_index] = degreesA[_index];
				} else {
					degreesC[_index] = degreesB[_index];
				}
				sizeC *= (degreesC[_index] + 1);
			}
			std::vector<tA> c(sizeC, zero);
			// perform multiplication
			for(unsigned int indexA = 0; indexA < a.size(); ++indexA) {
				const std::vector<unsigned int> ijkA = ijk(degreesA, indexA);
				for(unsigned int indexB = 0; indexB < size(degreesB); ++indexB) {
					const std::vector<unsigned int> ijkB = ijk(degreesB, indexB);
					// compute ijkC
					std::vector<unsigned int> ijkC(degreesC.size());
					for (unsigned int _index = 0; _index < ijkC.size(); ++_index) {
						if (_index < ijkA.size()) {
							if (_index < ijkB.size()) {
								ijkC[_index] = ijkA[_index] + ijkB[_index];   
							} else {
								ijkC[_index] = ijkA[_index];
							}
						} else {
							ijkC[_index] = ijkB[_index];
						}
					}
					// update c
					c[index(degreesC, ijkC)] += a[indexA] * pB[indexB];
				}
			}
			// store c in a
			a = c;
			degreesA = degreesC;
		}
	
		template<typename tPolynome>
		void multiply(std::vector<tPolynome>& a, std::vector<unsigned int>& degreesA, const double scalar) {
			for(unsigned int _index = 0; _index < a.size(); ++_index) {
				a[_index] *= scalar;
			}
		}

		// TODO: deal with case degree 0
		template<typename tPolynome>
		void differentiate(std::vector<tPolynome>& a, std::vector<unsigned int>& degreesA, const unsigned int variableIndex) {
			// create b
			std::vector<unsigned int> degreesB = degreesA;
			degreesB[variableIndex] -= 1;
			std::vector<tPolynome> b(size(degreesB));
			// differentiate
			for(unsigned int indexB = 0; indexB < b.size(); ++indexB) {
				const std::vector<unsigned int> ijkB = ijk(degreesB, indexB);
				std::vector<unsigned int> ijkA = ijkB;
				ijkA[variableIndex] += 1;
				b[indexB] = ijkA[variableIndex] * a[index(degreesA, ijkA)];
			}
			// store b in a
			a = b;
			degreesA = degreesB;
		}
	
		template<typename tPolynome>
		void differentiate(std::vector<tPolynome>& a, std::vector<unsigned int>& degreesA, const unsigned int variableIndex, const unsigned int order) {
			for(unsigned int _index = 0; _index < order; ++_index) {
				differentiate(a, degreesA, variableIndex);
			}
		}
	
		template<typename tPolynome>
		void integrate(std::vector<tPolynome>& a, std::vector<unsigned int>& degreesA, const unsigned int variableIndex, const tPolynome& constant) {
			// create b
			std::vector<unsigned int> degreesB = degreesA;
			degreesB[variableIndex] += 1;
			std::vector<tPolynome> b(size(degreesB));
			// integrate
			for(unsigned int indexA = 0; indexA < a.size(); ++indexA) {
				const std::vector<unsigned int> ijkA = ijk(degreesA, indexA);
				std::vector<unsigned int> ijkB = ijkA;
				ijkB[variableIndex] += 1;
				b[index(degreesB, ijkB)] = a[indexA] / ijkB[variableIndex];
			}
			// zero
			std::vector<unsigned int> ijkB(degreesB.size(), 0);
			for(unsigned int _index = 0; _index < ijkB[variableIndex]; ++_index) {
				ijkB[variableIndex] = _index;
				b[index(degreesB, ijkB)] = constant;
			}
			// store b in a
			a = b;
			degreesA = degreesB;
		}
	
		template<typename tPolynome>
		void integrate(std::vector<tPolynome>& a, std::vector<unsigned int>& degreesA, const unsigned int variableIndex, const tPolynome& constant, const unsigned int order) {
			for(unsigned int _index = 0; _index < order; ++_index) {
				integrate(a, degreesA, variableIndex, constant);
			}
		}

		template<typename tPolynome>
		tPolynome evaluation(const tPolynome* pA, const std::vector<unsigned int>& degreesA, const double* pX) {
			tPolynome output = pA[0];
			for (unsigned int indexA = 1; indexA < size(degreesA); ++indexA) {
				tPolynome term = pA[indexA];
				const std::vector<unsigned int> ijkA = ijk(degreesA, indexA);
				for (unsigned int _index = 0; _index < ijkA.size(); ++_index) {
					term *= std::pow(pX[_index], ijkA[_index]);
				}
				output += term;
			}
			return output;
		}
	}
}

}

#endif
