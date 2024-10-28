// std includes
#include <vector>
#include <iostream>
#include <random> // random_device, default_random_engine, uniform_real_distribution
#include <memory> // shared_ptr
#include <cmath> // sin, cos
// lib includes
#include "p0l/polynome.h"
#include "p0l/to_string.h"

int main() {
    std::cout << "monovariate polynomes:\n";
	// add
    std::vector<double> a({3.0, 4.0, 10.0});
    std::vector<double> b({0.0, 1.0, 0.0, 4.0});
    // print
    std::cout << "A: ";
    std::cout << p0l::polynome::toString(a.data(), a.size()) << "\n";
    std::cout << "B: ";
    std::cout << p0l::polynome::toString(b.data(), b.size()) << "\n";
    // add
    p0l::polynome::add(a, b.data(), b.size(), 0.0);
    std::cout << "A + B: ";
    std::cout << p0l::polynome::toString(a.data(), a.size()) << "\n";
    std::cout << "(A + B)(2) = " << p0l::polynome::evaluation(a.data(), a.size(), 2.0) << "\n";
    // separator
    std::cout << "\n";
    // substract
    std::vector<double> c({3.0, 4.0, 10.0});
    std::vector<double> d({0.0, 1.0, 0.0, 4.0});
    std::cout << "C: ";
    std::cout << p0l::polynome::toString(c.data(), c.size()) << "\n";
    std::cout << "D: ";
    std::cout << p0l::polynome::toString(d.data(), d.size()) << "\n";
    // substract
    p0l::polynome::substract(c, d.data(), d.size(), 0.0);
    std::cout << "C - D: ";
    std::cout << p0l::polynome::toString(c.data(), c.size()) << "\n";
    std::cout << "(C - D)(2) = " << p0l::polynome::evaluation(c.data(), c.size(), 2.0) << "\n";
    // separator
    std::cout << "\n";
    // multiply
    std::vector<double> e({1.0, 1.0, 1.0});
    std::vector<double> f({1.0, 1.0, 1.0});
    std::cout << "E: ";
    std::cout << p0l::polynome::toString(e.data(), e.size()) << "\n";
    std::cout << "F: ";
    std::cout << p0l::polynome::toString(f.data(), f.size()) << "\n";
    // multiply
    p0l::polynome::multiply(e, f.data(), f.size(), 0.0);
    std::cout << "E * F: ";
    std::cout << p0l::polynome::toString(e.data(), e.size()) << "\n";
    std::cout << "(E * F)(2) = " << p0l::polynome::evaluation(e.data(), e.size(), 2.0) << "\n";
    // separator
    std::cout << "\n";
    // calculus
    std::vector<double> g({1.0, 1.0, 1.0});
    std::vector<double> h({1.0, 1.0, 1.0});
    std::cout << "G: ";
    std::cout << p0l::polynome::toString(g.data(), g.size()) << "\n";
    std::cout << "h: ";
    std::cout << p0l::polynome::toString(h.data(), h.size()) << "\n";
    // calculus
    p0l::polynome::differentiate(g);
    std::cout << "g = G': ";
    std::cout << p0l::polynome::toString(g.data(), g.size()) << "\n";
    p0l::polynome::integrate(h, 0.0);
    std::cout << "H = int h: ";
    std::cout << p0l::polynome::toString(h.data(), h.size()) << "\n";

    std::cout << "\n\nmultivariate polynomes:\n";
    // x
    std::vector<double> x({2.0, 3.0});
    // add
    std::vector<unsigned int> degreesAA = {1, 2};
    std::vector<double> aa = p0l::polynome::multivariate::monome(degreesAA, 2.0, 0.0);
    aa[p0l::polynome::multivariate::index(degreesAA, {0, 1})] = 1.0;
    aa[p0l::polynome::multivariate::index(degreesAA, {1, 0})] = 1.0;
    std::vector<unsigned int> degreesBB = {2, 0};
    std::vector<double> bb = p0l::polynome::multivariate::monome(degreesBB, 2.0, 0.0);
    bb[p0l::polynome::multivariate::index(degreesBB, {1, 0})] = 1.0;
    // print
    std::cout << "A: ";
    std::cout << p0l::polynome::multivariate::toString(aa.data(), degreesAA) << "\n";
    std::cout << "B: ";
    std::cout << p0l::polynome::multivariate::toString(bb.data(), degreesBB) << "\n";
    // add
    p0l::polynome::multivariate::add(aa, degreesAA, bb.data(), degreesBB, 0.0);
    std::cout << "A + B: ";
    std::cout << p0l::polynome::multivariate::toString(aa.data(), degreesAA) << "\n";
    std::cout << "(A + B)(2, 3) = " << p0l::polynome::multivariate::evaluation(aa.data(), degreesAA, x.data()) << "\n";
    // separator
    std::cout << "\n";
    // substract
    std::vector<unsigned int> degreesCC = {1, 2};
    std::vector<double> cc = p0l::polynome::multivariate::monome(degreesCC, 2.0, 0.0);
    cc[p0l::polynome::multivariate::index(degreesCC, {0, 1})] = 1.0;
    cc[p0l::polynome::multivariate::index(degreesCC, {1, 0})] = 1.0;
    std::vector<unsigned int> degreesDD = {2, 0};
    std::vector<double> dd = p0l::polynome::multivariate::monome(degreesDD, 2.0, 0.0);
    dd[p0l::polynome::multivariate::index(degreesDD, {1, 0})] = 1.0;
    // print
    std::cout << "C: ";
    std::cout << p0l::polynome::multivariate::toString(cc.data(), degreesCC) << "\n";
    std::cout << "D: ";
    std::cout << p0l::polynome::multivariate::toString(dd.data(), degreesDD) << "\n";
    // substract
    p0l::polynome::multivariate::substract(cc, degreesCC, dd.data(), degreesDD, 0.0);
    std::cout << "C - D: ";
    std::cout << p0l::polynome::multivariate::toString(cc.data(), degreesCC) << "\n";
    std::cout << "(C - D)(2, 3) = " << p0l::polynome::multivariate::evaluation(cc.data(), degreesCC, x.data()) << "\n";
    // separator
    std::cout << "\n";
    // multiply
    std::vector<unsigned int> degreesEE = {1, 2};
    std::vector<double> ee = p0l::polynome::multivariate::monome(degreesEE, 2.0, 0.0);
    ee[p0l::polynome::multivariate::index(degreesEE, {0, 1})] = 1.0;
    ee[p0l::polynome::multivariate::index(degreesEE, {1, 0})] = 1.0;
    std::vector<unsigned int> degreesFF = {2, 0};
    std::vector<double> ff = p0l::polynome::multivariate::monome(degreesFF, 2.0, 0.0);
    ff[p0l::polynome::multivariate::index(degreesFF, {1, 0})] = 1.0;
    // print
    std::cout << "E: ";
    std::cout << p0l::polynome::multivariate::toString(ee.data(), degreesEE) << "\n";
    std::cout << "F: ";
    std::cout << p0l::polynome::multivariate::toString(ff.data(), degreesFF) << "\n";
    // multiply
    p0l::polynome::multivariate::multiply(ee, degreesEE, ff.data(), degreesFF, 0.0);
    std::cout << "E * F: ";
    std::cout << p0l::polynome::multivariate::toString(ee.data(), degreesEE) << "\n";
    std::cout << "(E * F)(2, 3) = " << p0l::polynome::multivariate::evaluation(ee.data(), degreesEE, x.data()) << "\n";
    // separator
    std::cout << "\n";
    // calculus
    std::vector<unsigned int> degreesGG = {1, 2};
    std::vector<double> gg = p0l::polynome::multivariate::monome(degreesGG, 2.0, 0.0);
    gg[p0l::polynome::multivariate::index(degreesGG, {0, 1})] = 1.0;
    gg[p0l::polynome::multivariate::index(degreesGG, {1, 0})] = 1.0;
    std::vector<unsigned int> degreesHH = {2, 2};
    std::vector<double> hh = p0l::polynome::multivariate::monome(degreesHH, 2.0, 0.0);
    hh[p0l::polynome::multivariate::index(degreesHH, {1, 0})] = 1.0;
    // print
    std::cout << "G: ";
    std::cout << p0l::polynome::multivariate::toString(gg.data(), degreesGG) << "\n";
    std::cout << "H: ";
    std::cout << p0l::polynome::multivariate::toString(hh.data(), degreesHH) << "\n";
    // calculus
    p0l::polynome::multivariate::differentiate(gg, degreesGG, 0);
    std::cout << "dG/dx: ";
    std::cout << p0l::polynome::multivariate::toString(gg.data(), degreesGG) << "\n";
    p0l::polynome::multivariate::integrate(hh, degreesHH, 1, 0.0);
    std::cout << "int H dy: ";
    std::cout << p0l::polynome::multivariate::toString(hh.data(), degreesHH) << "\n";
    // end
    std::cout << "\n";
}
