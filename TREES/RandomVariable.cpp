//
//  RandomVariable.cpp
//  TREES
//
//  Created by alex on 6/10/22.
//

//
//  RandomVariable.cpp
//  PseudoRandomNumberGenerator_Extension
//
//  Created by alex on 6/4/22.
//

#include <ctime>
#include <iostream>
#include <cmath>
#include "RandomVariable.hpp"
// CPP file needs to define the functions.

RandomVariable::RandomVariable(void) {
    seed = (int) time (NULL);
    //std::cout << "Seed value: " << seed << std::endl;
    std::cout << "Default constructor. The seed equals " << seed << std::endl;
}

RandomVariable::RandomVariable(int x) {
    seed = x;
    std::cout << "Alternate constructor. The seed equals " << seed << std::endl;
}

double RandomVariable::uniformRv(void) {
    int hi = seed / 127773;
    int lo = seed % 127773;
    int test = 16807 * lo - 2836 * hi;
    
    if (test > 0)
        seed = test;
    else
        seed = test + 2147483647;
    // end if
    
    return (double)(seed) / (double) 2147483647;
}

double RandomVariable::uniformRv(double lower, double upper) {
    return ( lower + uniformRv() * (upper - lower) );
}

double RandomVariable::exponentialRv(double lambda) {
    return -log(uniformRv()) / lambda;
}
