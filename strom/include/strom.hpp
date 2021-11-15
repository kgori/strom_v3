//
// Created by Kevin Gori on 10/10/2021.
//

#pragma once
#include "tree_summary.hpp"

#include <CLI11.hpp>
#include <fmt/core.h>

#include <fstream>
#include <iostream>

inline bool exists(const std::string &name) {
    std::ifstream f(name.c_str());
    return f.good();
}

namespace strom {
class Strom {
public:
    Strom();

    void clear();

    void processCommandLineOptions(int argc, const char *argv[]);

    void run();

private:
    std::string _data_file_name;
    std::string _tree_file_name;

    TreeSummary::SharedPtr _tree_summary;

    static std::string _program_name;
    static unsigned _major_version;
    static unsigned _minor_version;
};

inline Strom::Strom() {
    clear();
}

inline void Strom::clear() {
    _data_file_name = "";
    _tree_file_name = "";
    _tree_summary = nullptr;
}

inline void Strom::processCommandLineOptions(int argc, const char *argv[]) {
    CLI::App app{"strom"};
    app.add_option("datafile", _data_file_name);
    app.add_option("treefile", _tree_file_name);

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &err) {
        std::cout << err.what() << std::endl;
        app.exit(err);
    }
}

inline void Strom::run() {
    std::cout << "Starting..." << std::endl;

    try {
        // Create new TreeSummary
        _tree_summary = std::make_shared<TreeSummary>();

        // Read the user-specified tree file
        _tree_summary->readTreefile(_tree_file_name, 0);
        Tree::SharedPtr tree = _tree_summary->getTree(0);

        // Summarise the trees read
        _tree_summary->showSummary();
    } catch (XStrom &x) {
        std::cerr << "Strom encountered a problem:\n " << x.what() << std::endl;
    }

    std::cout << "Finished!" << std::endl;
}

} // namespace strom