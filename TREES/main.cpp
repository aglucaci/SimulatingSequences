//
//  main.cpp
//  TREES
//
//  Created by alex on 6/8/22.
//

#include <iostream>
//#include "Node.hpp"
#include <string>
#include "RandomVariable.hpp"
#include "Tree.hpp"

//#include <iostream>
//#include <string>
#include "Alignment.hpp"
//#include "RandomVariable.hpp"
//#include "Tree.hpp"

//#include "RandomVariable.hpp"
//#include "Tree.hpp"


int main(int argc, char* argv[]) {
    //RandomVariable rv;
    //Tree t(3.0, 1.0, 1.0, &rv);
    //Tree t(3.0, 1.0, 1.0, &rv);
    
    //Tree t(5.0, 1.0, 1.0, &rv);
    //t.listNodes();
    //std::cout << t.getNewickRepresentation() << std::endl;
    
    //std::string myTree = "(Taxon_I:0.3,((Taxon_II:0.1,Taxon_III:0.1):0.1, Taxon_IV:0.2):0.1);";
    //Tree t(myTree);
    //t.listNodes();
    //std::cout << t.getNewickRepresentation() << std::endl;
    
    
    // instantiate the random number generator
    RandomVariable rv;
    
    // set the parameters of the simulation
    int numSites = 100;
    
    std::string myTree = "(Taxon_I:0.3,((Taxon_II:0.1,Taxon_III:0.1):0.1,Taxon_IV:0.2):0.1);";
    
    double exchangabilityParms[6] = { 1.0, 5.0, 1.0, 1.0, 5.0, 1.0 };
    
    double baseFrequencyParms[4] = { 0.4, 0.3, 0.2, 0.1 };
    
    // instantiate the tree and alignment objects
    Tree t(myTree);
    
    Alignment myAlignment(&t, &rv, numSites, exchangabilityParms, baseFrequencyParms);
    
    myAlignment.print();
    
    return 0;
}
