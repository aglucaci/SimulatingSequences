//
//  RandomVariable.hpp
//  TREES
//
//  Created by alex on 6/10/22.
//

//
//  RandomVariable.hpp
//  PseudoRandomNumberGenerator_Extension
//
//  Created by alex on 6/4/22.
//


#pragma once

#ifndef RandomVariable_hpp
#define RandomVariable_hpp

#include <stdio.h>

// Helper file gives a skeleton of the class file.

class RandomVariable {
    public:
        RandomVariable(void);
        RandomVariable(int x);
        double uniformRv(void);
        double uniformRv(double lower, double upper);
        double exponentialRv(double lambda);
    
    protected:
        int seed;
};




#endif /* RandomVariable_hpp */




//#endif /* RandomVariable_hpp */
