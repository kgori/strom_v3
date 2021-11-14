//
// Created by Kevin Gori on 12/10/2021.
//

#pragma once

#include "genetic_code.hpp"
#include <fmt/core.h>

namespace strom {

class DataType {
public:
    DataType();
    ~DataType() = default;

    void setNucleotide();

    void setCodon();

    void setProtein();

    void setStandard();

    [[nodiscard]] bool isNucleotide() const;

    [[nodiscard]] bool isCodon() const;

    [[nodiscard]] bool isProtein() const;

    [[nodiscard]] bool isStandard() const;

    void setStandardNumStates(unsigned nstates);

    void setGeneticCodeFromName(const std::string &genetic_code_name);

    void setGeneticCode(GeneticCode::SharedPtr gcode);

    [[nodiscard]] unsigned getDataType() const;

    [[nodiscard]] unsigned getNumStates() const;

    [[nodiscard]] std::string getDataTypeAsString() const;

    [[nodiscard]] const GeneticCode::SharedPtr getGeneticCode() const;

    static std::string translateDataTypeToString(unsigned datatype);

private:
    unsigned _datatype{};
    unsigned _num_states{};
    GeneticCode::SharedPtr _genetic_code;

    enum DataTypes
    {
        NUCLEOTIDE,
        CODON,
        PROTEIN,
        STANDARD
    };
};

inline DataType::DataType() {
    setNucleotide();
}

inline void DataType::setNucleotide() {
    _datatype = DataTypes::NUCLEOTIDE;
    _num_states = 4;
    _genetic_code = nullptr;
}

inline void DataType::setCodon() {
    _datatype = DataTypes::CODON;
    _genetic_code = std::make_shared<GeneticCode>("standard");
    _num_states = _genetic_code->getNumNonStopCodons();
}

inline void DataType::setProtein() {
    _datatype = DataTypes::PROTEIN;
    _num_states = 20;
    _genetic_code = nullptr;
}

inline void DataType::setStandard() {
    _datatype = DataTypes::STANDARD;
    _num_states = 2;
    _genetic_code = nullptr;
}

inline bool DataType::isNucleotide() const {
    return _datatype == DataTypes::NUCLEOTIDE;
}

inline bool DataType::isCodon() const {
    return _datatype == DataTypes::CODON;
}

inline bool DataType::isProtein() const {
    return _datatype == DataTypes::PROTEIN;
}

inline bool DataType::isStandard() const {
    return _datatype == DataTypes::STANDARD;
}

inline void DataType::setGeneticCodeFromName(const std::string &genetic_code_name) {
    assert(isCodon());
    _genetic_code = std::make_shared<GeneticCode>(genetic_code_name);
}

inline void DataType::setGeneticCode(GeneticCode::SharedPtr gcode) {
    assert(isCodon());
    assert(gcode);
    _genetic_code = gcode;
}

inline void DataType::setStandardNumStates(unsigned int nstates) {
    _datatype = DataTypes::STANDARD;
    _num_states = nstates;
    _genetic_code = nullptr;
}

inline unsigned DataType::getDataType() const {
    return _datatype;
}

inline unsigned DataType::getNumStates() const {
    return _num_states;
}

inline const GeneticCode::SharedPtr DataType::getGeneticCode() const {
    assert(isCodon());
    return _genetic_code;
}

inline std::string DataType::getDataTypeAsString() const {
    std::string s = translateDataTypeToString(_datatype);
    if (isCodon()) {
        s += fmt::format(FMT_STRING(",{:s}"), _genetic_code->getGeneticCodeName());
    }
    return s;
}

inline std::string DataType::translateDataTypeToString(unsigned int datatype) {
    assert(datatype == DataTypes::NUCLEOTIDE || datatype == DataTypes::CODON ||
           datatype == DataTypes::PROTEIN || datatype == DataTypes::STANDARD);
    if (datatype == DataTypes::NUCLEOTIDE) {
        return std::string("nucleotide");
    } else if (datatype == DataTypes::CODON) {
        return std::string("codon");
    } else if (datatype == DataTypes::PROTEIN) {
        return std::string("protein");
    } else {
        return std::string("standard");
    }
}

}// namespace strom