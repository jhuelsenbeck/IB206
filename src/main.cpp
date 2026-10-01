#include <iostream>
#include <string>
#include "Alignment.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"



int main(int argc, const char * argv[]) {

    RandomVariable rng;

    std::string newickString = "(A:0.3,(D:0.2,(B:0.1,C:0.1):0.1):0.1);";
    Tree t(newickString);
    t.print();
    
    Alignment alignment(&rng, &t, 10000);
    alignment.print();
    
//    Tree bdTree(&rng, 2.0, 1, 7.0);
//    bdTree.print();
//    std::string newickStr = bdTree.getNewickString();
//    std::cout << newickStr << std::endl;
    
    
    return EXIT_SUCCESS;
}
