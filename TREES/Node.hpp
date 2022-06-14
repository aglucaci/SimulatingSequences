//
//  Node.hpp
//  TREES
//
//  Created by alex on 6/8/22.
//

#ifndef Node_hpp
#define Node_hpp

#include <stdio.h>
#include <string>

class Node {
    public:
        Node(void);
        Node* getLft(void) { return left; }
        Node* getRht(void) { return right; }
        Node* getAnc(void) { return ancestor; }
        int getIndex(void) { return index; }
        std::string getName(void) { return name; }
        double getBranchLength(void) { return branchLength; }
        void setLft(Node* p) { left = p; }
        void setRht(Node* p) { right = p; }
        void setAnc(Node* p) { ancestor = p; }
        void setIndex(int x) { index = x; }
        void setName(std::string s) { name = s; }
        void setBranchLength(double x) { branchLength = x; }
        void print(void);
        void setTime(double x) { t = x; }
        double getTime(void) { return t; }
        //void set_numExtant(int x) { numExtant = x; }
        //int get_numExtant(void) { return numExtant; }

    protected:
        Node* left;
        Node* right;
        Node* ancestor;
        int index;
        std::string name;
        double branchLength;
        double t;
        //int numExtant;
    
};
    
    
#endif /* Node_hpp */
