//
//  Node.cpp
//  TREES
//
//  Created by alex on 6/8/22.
//

#include "Node.hpp"
#include <iostream>

Node::Node(void) {
    left = NULL;
    right = NULL;
    ancestor = NULL;
    index = 0;
    name = "";
    branchLength = 0.0;
    t = 0.0;
}

void Node::print(void) {
    std::cout << "Node " << index << " (" << this << ")" << std::endl;
    std::cout << " Lft: " << left << std::endl;
    std::cout << " Rht: " << right << std::endl;
    std::cout << " Anc: " << ancestor << std::endl;
    std::cout << " Name: \"" << name << "\"" << std::endl;
    std::cout << " Brlen: " << branchLength << std::endl;
    std::cout << " Time: " << t << std::endl;
}


