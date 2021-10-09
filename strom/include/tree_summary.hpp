//
// Created by Kevin Gori on 08/10/2021.
//

#pragma once
#include <algorithm>
#include <cassert>
#include <fmt/core.h>
#include <fstream>
#include <map>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/view/reverse.hpp>
#include <set>
#include <vector>

#include "ncl/nxsmultiformat.h"

#include "split.hpp"
#include "tree_manip.hpp"
#include "xstrom.hpp"

namespace strom {

class TreeSummary {
public:
    void readTreefile(const std::string &filename, unsigned skip);

    void showSummary() const;

    typename Tree::SharedPtr getTree(unsigned index);

    std::string getNewick(unsigned index);

    void clear();

private:
    Split::treemap_t _treeIDs;
    std::vector<std::string> _newicks;

public:
    typedef std::shared_ptr<TreeSummary> SharedPtr;
};

inline typename Tree::SharedPtr TreeSummary::getTree(unsigned int index) {
    if (index > _newicks.size()) {
        throw XStrom("getTree called with index greater than number of trees");
    }

    TreeManip tm;

    tm.buildFromNewick(_newicks[index], false, false);

    return tm.getTree();
}

inline std::string TreeSummary::getNewick(unsigned int index) {
    if (index > _newicks.size()) {
        throw XStrom("getNewick called with index greater than number of trees");
    }
    return _newicks[index];
}

inline void TreeSummary::clear() {
    _newicks.clear();
    _treeIDs.clear();
}

inline void TreeSummary::readTreefile(const std::string &filename, unsigned int skip) {
    TreeManip tm;
    Split::treeid_t splitset;

    // See http://phylo.bio.ku.edu/ncldocs/v2.1/funcdocs/index.html for NCL documentation
    MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
    try {
        nexusReader.ReadFilepath(filename.c_str(), MultiFormatReader::NEXUS_FORMAT);
    } catch (...) {
        nexusReader.DeleteBlocksFromFactories();
        throw;
    }

    unsigned numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
    for (int i = 0; i < numTaxaBlocks; ++i) {
        clear();
        NxsTaxaBlock *taxaBlock = nexusReader.GetTaxaBlock(i);
        std::string taxaBlockTitle = taxaBlock->GetTitle();

        const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(taxaBlock);
        for (unsigned j = 0; j < nTreesBlocks; ++j) {
            const NxsTreesBlock *treesBlock = nexusReader.GetTreesBlock(taxaBlock, j);
            unsigned nTrees = treesBlock->GetNumTrees();

            if (skip < nTrees) {
                for (unsigned t = skip; t < nTrees; ++t) {
                    const NxsFullTreeDescription &d = treesBlock->GetFullTreeDescription(t);

                    // store the newick tree description
                    const std::string &newick = d.GetNewick();
                    _newicks.push_back(newick);
                    auto tree_index = static_cast<unsigned>(_newicks.size()) - 1;

                    // build the tree
                    tm.buildFromNewick(newick, false, false);

                    // store set of splits
                    splitset.clear();
                    tm.storeSplits(splitset);

                    auto iter = _treeIDs.lower_bound(splitset);

                    if (iter == _treeIDs.end() || iter->first != splitset) {
                        // splitset key not found in map, so need to create an entry
                        std::vector<unsigned> v(1, tree_index);
                        _treeIDs.insert(iter, Split::treemap_t::value_type(splitset, v));
                    } else {
                        // splitset key was found in map, so need to add this tree's index to vector
                        iter->second.push_back(tree_index);
                    }
                }// trees loop
            }    // skip loop
        }        //TREES block loop
    }            // TAXA block loop
    nexusReader.DeleteBlocksFromFactories();
}

inline void TreeSummary::showSummary() const {
    // Produce some output to show that it works
    fmt::print(FMT_STRING("\nRead {:d} trees from file\n"), _newicks.size());

    // Show all unique topologies with a list of trees that have that topology
    // Also create a map that can be used to sort topologies by their sample frequency
    typedef std::pair<unsigned, unsigned> sorted_pair_t;
    std::vector<sorted_pair_t> sorted;
    int t = 0;
    for (auto &key_value_pair : _treeIDs) {
        unsigned topology = ++t;
        auto ntrees = static_cast<unsigned>(key_value_pair.second.size());
        sorted.emplace_back(ntrees, topology);
        fmt::print(FMT_STRING("Topology {:d} seen in these {:d} trees:\n {}\n"),
                   topology,
                   ntrees,
                   fmt::join(key_value_pair.second.begin(), key_value_pair.second.end(), " "));
    }

    // Show sorted histogram data
    ranges::sort(sorted);
    fmt::print("\nTopologies sorted by sample frequency:\n");
    fmt::print(FMT_STRING("{:^20s} {:^20s}\n"), "topology", "frequency");
    for (auto &ntrees_topol_pair : ranges::views::reverse(sorted)) {
        unsigned n = ntrees_topol_pair.first;
        unsigned t = ntrees_topol_pair.second;
        fmt::print(FMT_STRING("{:^20d} {:^20d}\n"), t, n);
    }
}

}// namespace strom