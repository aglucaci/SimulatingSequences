//
//  Tree.cpp
//  TREES
//
//  Created by alex on 6/8/22.
//

#include <iostream>
#include "Node.hpp"
#include "Tree.hpp"
#include "RandomVariable.hpp"
//#include <iostream>
#include <string>
#include <set>
#include <cmath>
#include <iomanip>

//#define DEBUG_NEWICK_PARSER
#undef DEBUG_NEWICK_PARSER

Tree::Tree(std::string newickStr) {
    // break the string into its component parts
    std::vector<std::string> tokens = parseNewickString(newickStr);
    
    
    // build up the tree from the parsed Newick string
    bool readingBranchLength = false;
    Node* p = NULL;
    for (int i=0; i<tokens.size(); i++) {
        std::string token = tokens[i];
        //std::cout << token << std::endl;
        if (token == "(") {
            // new node
            Node* newNode = new Node;
            nodes.push_back(newNode);
            if (p == NULL)
                root = newNode;
            else {
                newNode->setAnc(p);
                if (p->getLft() == NULL) {
                    p->setLft(newNode);
                }
                else {
                    p->setRht(newNode);
                    }
                }
                p = newNode;
            }
            else if (token == ")" || token == ",") {
                // move down one node
                if (p->getAnc() == NULL) {
                    std::cout << "Error: We cannot find an expected ancestor" << std::endl;
                    exit(1);
                }
                p = p->getAnc();
            }
            else if (token == ":")
            {
                // begin reading a branch length
                readingBranchLength = true;
            }
            else if (token == ";") {
                // finished!
                if (p != root)
                {
                    std::cout << "Error: We expect to finish at the root node" << std::endl;
                    exit(1);
                }
            }
            else
            {
                // we have a taxon name or a branch length
                if (readingBranchLength == true) {
                    double x = stod(token);
                    p->setBranchLength(x);
                    readingBranchLength = false;
                }
                else
                {
                    Node* newNode = new Node;
                    nodes.push_back(newNode);
                    newNode->setAnc(p);
                    if (p->getLft() == NULL)
                        p->setLft(newNode);
                    else
                        p->setRht(newNode);
                    newNode->setName(token);
                    p = newNode;
                }
            }
        }
    
    // index the nodes
    int ndeIdx = 0;
    for (int i=0; i<nodes.size(); i++) {
        if (nodes[i]->getLft() == NULL)
        nodes[i]->setIndex(ndeIdx++);
    }
    for (int i=0; i<nodes.size(); i++) {
        if (nodes[i]->getLft() != NULL)
        nodes[i]->setIndex(ndeIdx++);
    }
    
    // initialzie the postorder traversal sequence
    initializeTraversalOrder();
}
    

std::vector<std::string> Tree::parseNewickString(std::string ns) {
    std::vector<std::string> tks;
    for (int i=0; i<ns.size(); i++) {
        char c = ns[i];
        if (c == '(' || c == ')' || c == ',' || c == ':' || c == ';') {
            std::string tempStr;
            tempStr = c;
            tks.push_back(tempStr);
        }
        else {
            int j = i;
            std::string tempStr = "";
            while ( !(c == '(' || c == ')' || c == ',' || c == ':' || c == ';') ) {
                tempStr += c;
                j++;
                c = ns[j];
                }
                i = j-1;
                tks.push_back(tempStr);
            }
        }
    
    # if defined(DEBUG_NEWICK_PARSER)
        std::cout << "The Newick string, broken into its parts:" << std::endl;
        for (int i=0; i<tks.size(); i++)
        std::cout << " tks[" << i << "] = \"" << tks[i] << "\"" << std::endl;
    # endif
    return tks;
}


Tree::Tree(double lambda, double mu, double duration, RandomVariable* rv) {
    // generate a birth-death with parameter lambda, mu, and duration
    std::cout << "Calling the birth-death tree constructor" << std::endl;
    
    // initialize the single lineage, adding
    // the descendant to a list of active nodes
    
    nodes.push_back( new Node );
    nodes.push_back( new Node );
    
    nodes[0]->setLft(nodes[1]);
    nodes[1]->setAnc(nodes[0]);
    
    std::set<Node*> activeNodes;
    activeNodes.insert( nodes[1] );
    root = nodes[0];
    
    // generate the full tree
    double t = 0.0;
    
    while (t < duration) {
        // increment t using the exponential distribution
        double rate = activeNodes.size() * (lambda + mu);
        t += rv->exponentialRv(rate);
        // if t is still less than duration, go ahead do the speciation
        // or extinction thing on a randomly selected active lineage
        if (t < duration) {
            // choose a node
            //Node* p = NULL;
            Node* p = chooseNodeFromSet(activeNodes, rv);
            p -> setTime(t);
            // choose a type of event
            double u = rv->uniformRv();
            
            if ( u < lambda / (lambda+mu) ) {
                // speciation event
                Node* newLft = new Node; // 1. allocate new left and
                Node* newRht = new Node; // right nodes
                nodes.push_back(newLft); // 2. add the new nodes to the tree
                nodes.push_back(newRht);
                newLft->setAnc(p); // 3. set the ancestor of both new
                newRht->setAnc(p); // nodes to be p
                p->setLft(newLft); // 4. set the left and right values of
                p->setRht(newRht); // p to be the new nodes
                activeNodes.erase(p); // 5. modify the list of active nodes
                activeNodes.insert(newLft);
                activeNodes.insert(newRht);
                
            }
            else {
                // extinction event
                // extinction event
                activeNodes.erase(p); // poor p ... he’s dead
            }
        }
    }
    
    // clean up
    //long numExtant;
    numExtant = activeNodes.size();
    
    for (Node* nde : activeNodes) {
        nde -> setTime(duration);
    }
    
    initializeTraversalOrder();
    
    // set the index variable and assign branch lengths from the node times
    int nodeIdx = 0;
    
    for (int i=0; i < postOrderSequence.size(); i++) {
        //Node* p = postOrderTraversalOrder()[i];
        Node* p = postOrderSequence[i];
        if (p->getLft() == NULL && p->getRht() == NULL) {
            p->setIndex(nodeIdx++);
            p->setName( std::to_string(nodeIdx) );
        }
        if (p->getAnc() != NULL)
            p->setBranchLength( p -> getTime() - p -> getAnc() -> getTime() );
        // end if
    }
    
    //for (int i=0; i < postOrderSequence.size(); i++) {
    for (int i=0; i < postOrderSequence.size(); i++) {
        //Node* p = postOrderTraversalOrder[i];
        Node* p = postOrderSequence[i];
        if ( !(p->getLft() == NULL && p->getRht() == NULL) )
            p->setIndex(nodeIdx++);
        // end if
    }
}

Tree::~Tree(void) {
    std::cout << "Calling the Tree destructor with " << nodes.size() << " nodes" << std::endl;
    for (int i=0; i<nodes.size(); i++)
        delete nodes[i];
}

void Tree::listNodes(void) {
    for (int i=0; i<nodes.size(); i++) {
        nodes[i]->print();
    }
}

void Tree::initializeTraversalOrder(void) {
    postOrderSequence.clear();
    passDown(root);
}

void Tree::passDown(Node* p) {
    if (p != NULL)
    {
        passDown(p->getLft());
        passDown(p->getRht());
        postOrderSequence.push_back(p);
    }
    //root = nodes[4];
    //initializeTraversalOrder();
}

Node* Tree::chooseNodeFromSet(std::set<Node*>& s, RandomVariable* rv) {
    int whichNode = (int)(s.size() * rv->uniformRv());
    int i = 0;
    for (Node* nde : s) { // for node in s, the node set
        if (whichNode == i)
            return nde;
        // end if
        i++;
    }
    return NULL;
}


std::string Tree::getNewickRepresentation(void) {
    std::stringstream ss;
    if (root->getLft() != NULL && root->getRht() != NULL)
        writeTree(root, ss);
    else
        writeTree(root->getLft(), ss);
    // end if
    std::string newick = ss.str();
    return newick;
}

void Tree::writeTree(Node* p, std::stringstream& ss) {
    if (p != NULL) {
        if (p->getLft() == NULL) {
            ss << p->getName() << ":" << std::fixed << std::setprecision(5) <<
            p->getBranchLength();
        }
        else {
            ss << "(";
            writeTree (p->getLft(), ss);
            ss << ",";
            writeTree (p->getRht(), ss);
            if (p->getAnc() == NULL)
                ss << ")";
            else
                ss << "):" << std::fixed << std::setprecision(5) << p->getBranchLength();
        }
    }
}
