//
//  Tree.hpp
//  TREES
//
//  Created by alex on 6/8/22.
//

#ifndef Tree_hpp
#define Tree_hpp

#include <stdio.h>
#include <vector>
#include <set>

//#include <set>
#include <sstream>
#include <string>
//#include <vector>

class Node;
class RandomVariable;

class Tree {
    public:
        Tree(double lambda, double mu, double duration,
        RandomVariable* rv);
        ~Tree(void);
        std::vector<Node*>& getTraversalOrder(void) { return postOrderSequence; }
        void listNodes(void);
        std::string getNewickRepresentation(void);
        void setNumExtant(int x) { numExtant = x; }
        long getNumExtant(void) { return numExtant; }
        Tree(std::string newickStr);
    
    protected:
        Tree(void) { }
        void initializeTraversalOrder(void);
        void passDown(Node* p);
        Node* root;
        void writeTree(Node* p, std::stringstream& ss);
        std::vector<Node*> nodes;
        std::vector<Node*> postOrderSequence;
        Node* chooseNodeFromSet(std::set<Node*>& s, RandomVariable* rv);
        long numExtant;
        std::vector<std::string> parseNewickString(std::string ns);
};
#endif /* Tree_hpp */
