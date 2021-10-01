#include <iostream>
#include "node.hpp"
#include "tree_manip.hpp"
#include "libhmsbeagle/beagle.h"

using namespace strom;

const double Node::_smallest_edge_length = 1.0e-12;

int main(int argc, const char *argv[]) {
    std::cout << "Starting..." << std::endl;
    TreeManip tm;
    tm.createTestTree();
    std::cout << "\nFinished!" << std::endl;
    std::cout << tm.makeNewick(3, false) << std::endl;
    std::cout << tm.makeNewick(3, true) << std::endl;
    return 0;
}
