//
//  Alignment.hpp
//  TREES
//
//  Created by alex on 6/13/22.
//

#ifndef Alignment_hpp
#define Alignment_hpp

#include <stdio.h>
#include <vector>
#include <string>

class RandomVariable;
class Tree;

class Alignment {
    public:
        Alignment(Tree* t, RandomVariable* rv, int ns, double ep[6], double bfp[4]);
        ~Alignment(void);
        void print(void);
    
    
    private:
        void initializeRateMatrix(double ep[6], double bfp[4]);
        void scaleRateMatrix(double bfp[4]);
        int numTaxa;
        int numSites;
        double q[4][4];
        int** matrix;
        std::vector<std::string> names;
        char convertIndexToNucleotide(int x);
};

#endif /* Alignment_hpp */
