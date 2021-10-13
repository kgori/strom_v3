//
// Created by Kevin Gori on 12/10/2021.
//

#pragma once

#include "xstrom.hpp"
#include <algorithm>
#include <cctype>

void tolowercase(std::string &s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](auto c) { return std::tolower(c); });
}

namespace strom {

class Data;

class Model;

class QMatrix;

class GeneticCode {
    friend class Data;

    friend class Model;

    friend class QMatrix;

public:
    typedef std::map<int, int> genetic_code_map_t;
    typedef std::map<char, unsigned> amino_acid_map_t;
    typedef std::vector<unsigned> amino_acid_vect_t;
    typedef std::vector<std::string> codon_vect_t;
    typedef std::vector<char> amino_acid_symbol_vect_t;
    typedef std::map<std::string, std::string> genetic_code_definitions_t;
    typedef std::vector<std::string> genetic_code_names_t;

    GeneticCode();

    explicit GeneticCode(std::string name);

    ~GeneticCode() = default;

    [[nodiscard]] std::string getGeneticCodeName() const;

    void useGeneticCode(const std::string &name);

    [[nodiscard]] unsigned getNumNonStopCodons() const;

    [[nodiscard]] int getStateCode(int triplet_index) const;

    [[nodiscard]] char getAminoAcidAbbrev(unsigned aa_index) const;

    void copyCodons(codon_vect_t &codon_vect) const;

    void copyAminoAcids(amino_acid_vect_t &aa_vect) const;

    [[nodiscard]] static genetic_code_names_t getRecognizedGeneticCodeNames();

    static bool isRecognizedGeneticCodeName(const std::string &name);

    static void ensureGeneticCodeNameIsValid(const std::string &name);

private:
    void buildGeneticCodeTranslators();

    std::string _genetic_code_name;
    genetic_code_map_t _genetic_code_map;
    amino_acid_map_t _amino_acid_map;

    amino_acid_vect_t _amino_acids;
    codon_vect_t _codons;

    const amino_acid_symbol_vect_t _all_amino_acids = {'A', 'C', 'D', 'E', 'F',
                                                       'G', 'H', 'I', 'K', 'L',
                                                       'M', 'N', 'P', 'Q', 'R',
                                                       'S', 'T', 'V', 'W', 'Y'};

    const codon_vect_t _all_codons = {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT",
                                      "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
                                      "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT",
                                      "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
                                      "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
                                      "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
                                      "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT",
                                      "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};

    static genetic_code_definitions_t _definitions;

public:
    typedef std::shared_ptr<GeneticCode> SharedPtr;
};

inline GeneticCode::GeneticCode() {
    useGeneticCode("standard");
}

inline GeneticCode::GeneticCode(std::string name) {
    useGeneticCode(name);
}

inline std::string GeneticCode::getGeneticCodeName() const {
    return _genetic_code_name;
}

inline void GeneticCode::useGeneticCode(const std::string &name) {
    _genetic_code_name = name;
    buildGeneticCodeTranslators();
}

inline unsigned GeneticCode::getNumNonStopCodons() const {
    return static_cast<unsigned>(_codons.size());
}

inline int GeneticCode::getStateCode(int triplet_index) const {
    return _genetic_code_map.at(triplet_index);
}

inline char GeneticCode::getAminoAcidAbbrev(unsigned int aa_index) const {
    return _all_amino_acids[aa_index];
}

inline void GeneticCode::copyCodons(codon_vect_t &codon_vect) const {
    codon_vect.resize(_codons.size());
    std::copy(_codons.begin(), _codons.end(), codon_vect.begin());
}

inline void GeneticCode::copyAminoAcids(amino_acid_vect_t &aa_vect) const {
    aa_vect.resize(_amino_acids.size());
    std::copy(_amino_acids.begin(), _amino_acids.end(), aa_vect.begin());
}

inline void GeneticCode::buildGeneticCodeTranslators() {
    _amino_acid_map.clear();
    for (unsigned i = 0; i < _all_amino_acids.size(); ++i) {
        char aa = _all_amino_acids[i];
        _amino_acid_map[aa] = i;
    }

    ensureGeneticCodeNameIsValid(_genetic_code_name);
    std::string gcode_aa = _definitions[_genetic_code_name];// e.g. "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"

    int k = 0;
    int state_code = 0;
    _codons.clear();
    _amino_acids.clear();
    _genetic_code_map.clear();
    for (char ch : gcode_aa) {
        if (ch != '*') {
            _genetic_code_map[k] = state_code++;
            _codons.push_back(_all_codons[k]);
            _amino_acids.push_back(_amino_acid_map[ch]);
        }
        ++k;
    }
}

inline genetic_code_names_t GeneticCode::getRecognizedGeneticCodeNames() {
    genetic_code_names_t names;
    for (const auto &[k, _] : _definitions) {
        names.push_back(k);
    }
    return names;
}

inline bool GeneticCode::isRecognizedGeneticCodeName(const std::string &name) {
    std::string lcname = name;
    tolowercase(lcname);
    return _definitions.find(lcname) != _definitions.end();
}

inline void GeneticCode::ensureGeneticCodeNameIsValid(const std::string &name) {
    if (!isRecognizedGeneticCodeName(name)) {
        auto valid_genetic_code_names = getRecognizedGeneticCodeNames();
        std::cout << "Recognized genetic codes:\n";
        for (std::string &name : valid_genetic_code_names) {
            std::cout << "  " << name << '\n';
        }
        std::cout << std::endl;
        throw XStrom(fmt::format(FMT_STRING("{:s} is not a recognized genetic code"), name));
    }
}
}// namespace strom
