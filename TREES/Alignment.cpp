//
//  Alignment.cpp
//  TREES
//
//  Created by alex on 6/13/22.
//

#include <iomanip>
#include <iostream>
#include "Alignment.hpp"
#include "Tree.hpp"
#include "Node.hpp"
#include "RandomVariable.hpp"


//#define DEBUG_RATE_MATRIX
#undef DEBUG_RATE_MATRIX

//class RandomVariable;
//class Node;

Alignment::Alignment(Tree* t, RandomVariable* rv, int ns, double ep[6], double bfp[4]) {
    // set the instance variable
    numSites = ns;
    // initialize and scale the rate matrix
    initializeRateMatrix(ep, bfp);
    scaleRateMatrix(bfp);
    // get the traversal order of the nodes for the tree
    std::vector<Node*> traversalSequence = t->getTraversalOrder();
    // allocate enough memory to hold the DNA sequences for each node in the tree
    int** nodeSequences = new int*[traversalSequence.size()];
    nodeSequences[0] = new int[traversalSequence.size() * numSites];
    for (int i=1; i<traversalSequence.size(); i++)
        nodeSequences[i] = nodeSequences[i-1] + numSites;
    for (int i=0; i<traversalSequence.size(); i++)
        for (int j=0; j<numSites; j++)
            nodeSequences[i][j] = 0;
    // loop over the sites
    for (int c=0; c<numSites; c++)
    {
        // traverse the nodes in preorder, simulating along each branch
        for (int n=(int)traversalSequence.size()-1; n>=0; n--) {
            // pass
            Node* p = traversalSequence[n];
            if (p->getAnc() == NULL) {
                // we are at the root of the tree
                double u = rv->uniformRv(), sum = 0.0;
                for (int i=0; i<4; i++) {
                    sum += bfp[i];
                    if (u < sum) {
                        nodeSequences[p->getIndex()][c] = i;
                        break;
                    }
                }
            }
            else {
            // we are at an internal node or a tip
                int curNuc = nodeSequences[p->getAnc()->getIndex()][c];
                double duration = p->getBranchLength();
                double t = 0.0;
                while (t < duration) {
                    double rate = -q[curNuc][curNuc];
                    t += rv->exponentialRv(rate);
                    if (t < duration) {
                        // we have a change...what kind is it?
                        double u = rv->uniformRv(), sum = 0.0;
                        for (int j=0; j<4; j++) {
                            if (j != curNuc) {
                                sum += q[curNuc][j] / rate;
                                if (u < sum) {
                                    curNuc = j;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                nodeSequences[p->getIndex()][c] = curNuc;
            }
        }
    }
            

    // allocate a matrix for just the tip nodes and fill it in
    numTaxa = 0;
    for (int i=0; i<traversalSequence.size(); i++) {
        if (traversalSequence[i]->getLft() == NULL)
            numTaxa++;
    }

    matrix = new int*[numTaxa];
    matrix[0] = new int[numTaxa * numSites];
    for (int i=1; i<numTaxa; i++)
        matrix[i] = matrix[i-1] + numSites;

    for (int i=0; i<numTaxa; i++)
        for (int j=0; j<numSites; j++)
            matrix[i][j] = 0;

    for (int n=0; n<traversalSequence.size(); n++) {
        Node* p = traversalSequence[n];
        if (p->getLft() == NULL) {
            names.push_back(p->getName());
            int idx = p->getIndex();
            for (int c=0; c<numSites; c++)
                matrix[idx][c] = nodeSequences[idx][c];
        }
    }
    // free the memory that was allocated
    delete [] nodeSequences[0];
    delete [] nodeSequences;
}

void Alignment::initializeRateMatrix(double ep[6], double bfp[4]) {
    // the exchangeability parameters in ep are in the
    // order A <-> C, A <-> G, A <-> T, C <-> G, C <-> T,
    // and G <-> T. The base frequency, or stationary
    // probability, parameters are in the order A, C, G, and T.
    // set the off-diagonal components
    for (int i=0, k=0; i<4; i++) {
        for (int j=i+1; j<4; j++) {
            q[i][j] = ep[k] * bfp[j];
            q[j][i] = ep[k] * bfp[i];
            k++;
        }
    }
    
    // set the diagonal components
    for (int i=0; i<4; i++) {
        double sum = 0.0;
        for (int j=0; j<4; j++) {
            if (i != j)
                sum += q[i][j];
        }
        q[i][i] = -sum;
    }
    
    # if defined(DEBUG_RATE_MATRIX)
    std::cout << "Rate matrix before scaling:" << std::endl;
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            if (q[i][j] > 0.0)
                std::cout << std::fixed << std::setprecision(3) << " " << q[i][j] << " ";
            else
                std::cout << std::fixed << std::setprecision(3) << q[i][j] << " ";
        }
        std::cout << std::endl;
    }
    # endif
}


void Alignment::scaleRateMatrix(double bfp[4]) {
    // rescale the rate matrix such that the average
    // rate of substitution is one
    double averageRate = 0.0;
    for (int i=0; i<4; i++)
        averageRate += -(bfp[i] * q[i][i]);
    // end for
    
    double scaler = 1.0 / averageRate;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            q[i][j] *= scaler;
        // end inner for
    // end for
    
    # if defined(DEBUG_RATE_MATRIX)
    std::cout << "Rate matrix after scaling:" << std::endl;
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            if (q[i][j] > 0.0)
                std::cout << std::fixed << std::setprecision(3) << " " << q[i][j] << " ";
            else
                std::cout << std::fixed << std::setprecision(3) << q[i][j] << " ";
        }
        std::cout << std::endl;
    }
    
    # endif
}


Alignment::~Alignment(void) {
    delete [] matrix[0];
    delete [] matrix;
    
}

void Alignment::print(void) {
    // get length of longest name
    int length = 0;
    for (int i=0; i<names.size(); i++) {
        if (names[i].size() > length)
            length = (int)names[i].size();
    }
    
    // print the alignment
    for (int i=0; i<numTaxa; i++) {
        std::cout << " " << names[i];
    
        for (int j=0; j<length-names[i].size(); j++)
            std::cout << " ";
        std::cout << " ";
    
        for (int j=0; j<numSites; j++)
            //std::cout << matrix[i][j];
            std::cout << convertIndexToNucleotide( matrix[i][j] );
            
        std::cout << std::endl;
       
    }
}

char Alignment::convertIndexToNucleotide(int x) {
    char nucs[4] = { 'A', 'C', 'G', 'T' };
    return nucs[x];
}
