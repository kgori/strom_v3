#define CLI_PARSER_CLI11
#include "strom.hpp"
#include "tree_summary.hpp"
#include <iostream>

using namespace strom;

std::string Strom::_program_name = "strom";
unsigned Strom::_major_version = 1;
unsigned Strom::_minor_version = 0;
const double Node::_smallest_edge_length = 1.0e-12;

int main(int argc, const char *argv[]) {
    Strom strom;

    try {
        strom.processCommandLineOptions(argc, argv);
        strom.run();
    } catch (std::exception &x) {
        std::cerr << "Exception: " << x.what() << std::endl;
        std::cerr << "Aborted." << std::endl;
    } catch (...) {
        std::cerr << "Exception of unknown type!\n";
    }

    return 0;
}
