//
// Created by Kevin Gori on 13/10/2021.
//

#pragma once

#include <cmath>
#include <fmt/core.h>
#include <limits>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/split.hpp>
#include <range/v3/view/transform.hpp>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include "datatype.hpp"
#include "genetic_code.hpp"
#include "xstrom.hpp"

namespace strom {
class Partition {
public:
    typedef std::match_results<std::string::const_iterator>::const_reference regex_match_t;
    typedef std::tuple<unsigned, unsigned, unsigned, unsigned> subset_range_t;
    typedef std::vector<subset_range_t> partition_t;
    typedef std::vector<DataType> datatype_vect_t;
    typedef std::vector<unsigned> subset_sizes_vect_t;
    typedef std::vector<std::string> subset_names_vect_t;
    typedef std::shared_ptr<Partition> SharedPtr;

    Partition();

    [[nodiscard]] unsigned getNumSites() const;

    [[nodiscard]] unsigned getNumSubsets() const;

    [[nodiscard]] std::string getSubsetName(unsigned subset) const;

    [[nodiscard]] const partition_t &getSubsetRangeVect() const;

    [[nodiscard]] unsigned findSubsetByName(const std::string &subset_name) const;

    [[nodiscard]] unsigned findSubsetForSite(unsigned site_index) const;

    [[nodiscard]] bool siteInSubset(unsigned site_index, unsigned subset_index) const;

    [[nodiscard]] DataType getDataTypeForSubset(unsigned subset_index) const;

    [[nodiscard]] const datatype_vect_t &getSubsetDataTypes() const;

    [[nodiscard]] unsigned numSitesInSubset(unsigned subset_index) const;

    [[nodiscard]] subset_sizes_vect_t calcSubsetSizes() const;

    void defaultPartition(unsigned nsites = std::numeric_limits<unsigned>::max());

    void parseSubsetDefinition(std::string &s);

    void finalize(unsigned nsites);

    void clear();

private:
    int extractIntFromRegexMatch(regex_match_t s, unsigned min_value);

    void addSubsetRange(unsigned subset_index, std::string range_definition);

    void addSubset(unsigned subset_index, std::string subset_definition);

    unsigned _num_sites;
    unsigned _num_subsets;
    subset_names_vect_t _subset_names;
    partition_t _subset_ranges;
    datatype_vect_t _subset_data_types;

    const unsigned _infinity;
};

inline Partition::Partition() : _infinity(std::numeric_limits<unsigned>::max()) {
    clear();
}

unsigned Partition::getNumSites() const {
    return _num_sites;
}

unsigned Partition::getNumSubsets() const {
    return _num_subsets;
}

std::string Partition::getSubsetName(unsigned subset) const {
    assert(subset < _num_subsets);
    return _subset_names[subset];
}

const Partition::partition_t &Partition::getSubsetRangeVect() const {
    return _subset_ranges;
}

DataType Partition::getDataTypeForSubset(unsigned int subset_index) const {
    assert(subset_index < _subset_data_types.size());
    return _subset_data_types[subset_index];
}

unsigned Partition::findSubsetByName(const std::string &subset_name) const {
    auto iter = std::find(_subset_names.begin(), _subset_names.end(), subset_name);
    if (iter == _subset_names.end()) {
        throw XStrom(fmt::format(FMT_STRING("Specified subset name \"{:s}\" not found in partition"), subset_name));
    }
    return static_cast<unsigned>(std::distance(_subset_names.begin(), iter));
}

unsigned Partition::findSubsetForSite(unsigned int site_index) const {
    for (auto &t : _subset_ranges) {
        auto [begin_site, end_site, stride, site_subset] = t;
        bool inside_range = site_index >= begin_site && site_index <= end_site;
        if (inside_range && (site_index - begin_site) % stride == 0) {
            return site_subset;
        }
    }
    throw XStrom(fmt::format(FMT_STRING("Site {:d} not found in any subset of partition"), (site_index + 1)));
}

bool Partition::siteInSubset(unsigned int site_index, unsigned int subset_index) const {
    unsigned which_subset = findSubsetForSite(site_index);
    return which_subset == subset_index;
}

const Partition::datatype_vect_t &Partition::getSubsetDataTypes() const {
    return _subset_data_types;
}

unsigned Partition::numSitesInSubset(unsigned int subset_index) const {
    unsigned nsites = 0;
    for (auto &t : _subset_ranges) {
        auto [begin_site, end_site, stride, site_subset] = t;
        if (site_subset == subset_index) {
            unsigned n = end_site - begin_site + 1;
            nsites += static_cast<unsigned>(floor(n / stride) + (n % stride == 0 ? 0 : 1));
        }
    }
    return nsites;
}

Partition::subset_sizes_vect_t Partition::calcSubsetSizes() const {
    assert(_num_sites > 0);
    std::vector<unsigned> nsites_vect(_num_subsets, 0);
    for (auto &t : _subset_ranges) {
        auto [begin_site, end_site, stride, site_subset] = t;
        unsigned hull = end_site - begin_site + 1;
        auto n = static_cast<unsigned>(floor(hull / stride) + (hull % stride) == 0 ? 0 : 1);
        nsites_vect[site_subset] += n;
    }
    return nsites_vect;
}

void Partition::defaultPartition(unsigned int nsites) {
    clear();
    _num_sites = nsites;
    _num_subsets = 1;
    _subset_ranges[0] = std::make_tuple(1, nsites, 1, 0);
}

// Made split (using ranges) into a free function to simplify uses elsewhere
std::vector<std::string> split(std::string_view s, const char delim) {
    auto v = s |
             ranges::views::split(delim) |
             ranges::views::transform([](auto rng) {
                 return rng | ranges::to<std::string>;
             }) |
             ranges::to<std::vector>;
    return v;
}

void Partition::parseSubsetDefinition(std::string &s) {
    // first separate part before colon (stored in v[0]) from part after colon (stored in v[1])
    auto v = split(s, ':');

    if (v.size() != 2) {
        throw XStrom("Expecting exactly one colon in partition subset definition");
    }

    std::string before_colon = v[0];
    std::string subset_definition = v[1];

    // now see if before_colon contains a data type specification in square brackets
    const char *pattern_string = R"((.+?)\s*(\[(\S+?)\])*)";
    std::regex re(pattern_string);
    std::smatch match_obj;
    bool matched = std::regex_match(before_colon, match_obj, re);
    if (!matched) {
        throw XStrom(fmt::format(
                FMT_STRING("Could not interpret \"{:s}\" as a subset label with optional data type in square brackets"),
                before_colon));
    }

    // match_obj always yields two strings that can be indexed using operator[]
    // match_obj[0] equals the entire subset label/type string (e.g. "rbcL[codon,standard]")
    // match_obj[1] equals the subset label (e.g. "rbcL")

    // Two more elements will exist if the user has specified a data type for this partition subset
    // match_obj[2] equals data type inside square brackets (e.g. "[codon,standard]")
    // match_obj[3] equals data type only (e.g. "codon,standard")

    std::string subset_name = match_obj[1].str();
    DataType dt;// nucleotide by default
    std::string datatype = "nucleotide";
    if (match_obj.size() == 4 && match_obj[3].length() > 0) {
        datatype = match_obj[3].str();
        std::transform(datatype.begin(), datatype.end(), datatype.begin(),
                       [](unsigned char c) { return std::tolower(c); });

        // Check for comma plus genetic code in case of codon
        std::regex re(R"(codon\s*,\s*(\S+))");
        std::smatch m;
        if (std::regex_match(datatype, m, re)) {
            dt.setCodon();
            std::string genetic_code_name = m[1].str();
            dt.setGeneticCodeFromName(genetic_code_name);
        } else if (datatype == "codon") {
            dt.setCodon();// assumes standard
        } else if (datatype == "protein") {
            dt.setProtein();
        } else if (datatype == "nucleotide") {
            dt.setNucleotide();
        } else if (datatype == "standard") {
            dt.setStandard();
        } else {
            throw XStrom(fmt::format(FMT_STRING(
                                             "Datatype \"{:s}\" specified for subset(s) \"{:s}\" is invalid: must be either nucleotide, codon, protein, or standard"),
                                     datatype, subset_name));
        }
    }

    // Remove default subset if there is one
    unsigned end_site = std::get<1>(_subset_ranges[0]);
    if (_num_subsets == 1 && end_site == _infinity) {
        _subset_names.clear();
        _subset_data_types.clear();
        _subset_ranges.clear();
    } else if (subset_name == "default") {
        throw XStrom("Cannot specifiy \"default\" partition subset after already defining other subsets");
    }
    _subset_names.push_back(subset_name);
    _subset_data_types.push_back(dt);
    _num_subsets = static_cast<unsigned>(_subset_names.size());
    addSubset(_num_subsets - 1, subset_definition);

    std::cout << fmt::format(FMT_STRING("Partition subset {:s} comprises sites {:s} and has type {:s}"),
                             subset_name, subset_definition, datatype)
              << std::endl;
}

void Partition::finalize(unsigned int nsites) {
    if (_num_sites == 0) {
        defaultPartition(nsites);
        return;
    }

    // First sanity check:
    //   nsites is the number of sites read in from the alignment file
    //   _num_sites is the maximum site index specified in any partition subset.
    //   These two numbers should be the same.
    if (_num_sites != nsites) {
        throw XStrom(fmt::format(FMT_STRING("Number of sites specified by the partition ({:d}) does not match the actual number of sites ({:d})"),
                                 _num_sites, nsites));
    }

    // Second sanity check: ensure no sites were left out of all partition subsets
    // Third sanity check: ensure that no sites were included in more than one partition subset
    std::vector<int> tmp(nsites, -1);
    for (auto &t : _subset_ranges) {
        auto [begin_site, end_site, stride, site_subset] = t;
        for (unsigned s = begin_site; s <= end_site; s += stride) {
            if (tmp[s - 1] != -1) {
                throw XStrom("Some sites were included in more than on partition subset");
            } else {
                tmp[s - 1] = static_cast<int>(site_subset);
            }
        }
    }
    if (std::find(tmp.begin(), tmp.end(), -1) != tmp.end()) {
        throw XStrom("Some sites were not included in any partition subset");
    }
    tmp.clear();
}

void Partition::clear() {
    _num_sites = 0;
    _num_subsets = 1;
    _subset_data_types.clear();
    _subset_data_types.push_back(DataType());
    _subset_names.clear();
    _subset_names.push_back("default");
    _subset_ranges.clear();
    _subset_ranges.push_back(std::make_tuple(1, _infinity, 1, 0));
}

int Partition::extractIntFromRegexMatch(std::match_results<std::string::const_iterator>::const_reference s,
                                        unsigned int min_value) {
    int int_value = min_value;
    if (s.length() > 0) {
        std::string str_value = s.str();
        try {
            int_value = std::stoi(str_value);
        } catch (std::invalid_argument) {
            throw XStrom(fmt::format(FMT_STRING("Could not interpret \"{:s}\" as a number in partition subset definition"),
                                     str_value));
        }

        // sanity check
        if (int_value < static_cast<int>(min_value)) {
            throw XStrom(fmt::format(FMT_STRING("Value specified in partition subset definition ({:d}) is lower than minimum value ({:d})"),
                                     int_value, min_value));
        }
    }

    return int_value;
}

void Partition::addSubsetRange(unsigned int subset_index, std::string range_definition) {
    // match patterns like these: "1-.\3" "1-1000" "1001-."
    const char *pattern_string = R"((\d+)\s*(-\s*([0-9.]+)(\\\s*(\d+))*)*)";
    std::regex re(pattern_string);
    std::smatch match_obj;
    bool matched = std::regex_match(range_definition, match_obj, re);
    if (!matched) {
        throw XStrom(fmt::format(FMT_STRING("Could not interpret \"{:s}\" as a range of site indices"),
                                 range_definition));
    }

    // match_obj always yields 6 strings that can be indexed using operator[]
    // match_obj[0] equals the entire site range (e.g. "1-.\3")
    // match_obj[1] equals the beginning site index (e.g. "1")
    // match_obj[2] equals everything after the beginning site index (e.g. "-.\3")
    // match_obj[3] equals "" or the ending site index or placeholder (e.g. ".")
    // match_obj[4] equals "" or everything after ending site index (e.g. "\3")
    // match_obj[5] equals "" or the step value (e.g. "3")
    int ibegin = extractIntFromRegexMatch(match_obj[1], 1);
    int iend = extractIntFromRegexMatch(match_obj[3], ibegin);
    int istep = extractIntFromRegexMatch(match_obj[5], 1);

    // record the triplet
    _subset_ranges.push_back(std::make_tuple(ibegin, iend, istep, subset_index));

    // determine last site in subset
    unsigned last_site_in_subset = iend - ((iend - ibegin) % istep);
    if (last_site_in_subset > _num_sites) {
        _num_sites = last_site_in_subset;
    }
}

void Partition::addSubset(unsigned int subset_index, std::string subset_definition) {
    std::vector<std::string> parts = split(subset_definition, ',');
    for (const std::string& subset_component : parts) {
        addSubsetRange(subset_index, subset_component);
    }
}

}// namespace strom