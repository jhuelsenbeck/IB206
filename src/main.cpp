#include <iostream>
#include <string>
#include "Tree.hpp"


int main(int argc, const char * argv[]) {

    std::string newickString = "((A,B),(C,D));";
    Tree t(newickString);
    t.print();
    
    return EXIT_SUCCESS;
}
