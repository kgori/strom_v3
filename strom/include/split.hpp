//
// Created by Kevin Gori on 06/10/2021.
//

#pragma once

#include <cassert>
#include <climits>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <vector>

namespace strom {

class Split {
public:
    Split();

    bool operator==(const Split &other) const;

    bool operator<(const Split &other) const;

    void clear();

    void resize(unsigned nleaves);

    typedef unsigned long split_unit_t;
    typedef std::vector<split_unit_t> split_t;
    typedef std::set<Split> treeid_t;
    typedef std::map<treeid_t, std::vector<unsigned>> treemap_t;
    typedef std::tuple<unsigned, unsigned, unsigned> split_metrics_t;

    [[nodiscard]] split_unit_t getBits(unsigned unit_index) const;

    [[nodiscard]] bool getBitAt(unsigned leaf_index) const;

    void setBitAt(unsigned leaf_index);

    void addSplit(const Split &other);

    [[nodiscard]] bool isEquivalent(const Split &other) const;

    [[nodiscard]] bool isCompatibleWith(const Split &other) const;

    [[nodiscard]] bool conflictsWith(const Split &other) const;

    [[nodiscard]] std::string createPatternRepresentation() const;

    [[nodiscard]] split_metrics_t getSplitMetrics() const;

private:
    split_unit_t _mask;
    split_t _bits;
    unsigned _bits_per_unit;
    unsigned _nleaves;

public:
    typedef std::shared_ptr<Split> SharedPtr;
};

inline Split::Split() {
    _mask = 0L;
    _nleaves = 0;
    _bits_per_unit = (CHAR_BIT) * sizeof(Split::split_unit_t);

    clear();
}

inline void Split::clear() {
    for (auto &u : _bits) {
        u = 0L;
    }
}

inline bool Split::operator==(const Split &other) const {
    return (_bits == other._bits);
}

inline bool Split::operator<(const Split &other) const {
    assert(_bits.size() == other._bits.size());
    return (_bits < other._bits);
}

inline void Split::resize(unsigned int nleaves) {
    _nleaves = nleaves;
    unsigned nunits = 1 + ((nleaves - 1) / _bits_per_unit);
    _bits.resize(nunits);

    // Create mask used to select only those bits used in final unit
    unsigned num_unused_bits = nunits * _bits_per_unit - _nleaves;
    unsigned num_used_bits = _bits_per_unit - num_unused_bits;
    _mask = 0L;
    split_unit_t unity = 1;
    for (unsigned i = 0; i < num_used_bits; ++i) {
        _mask |= (unity << i);
    }

    clear();
}

inline void Split::setBitAt(unsigned int leaf_index) {
    unsigned unit_index = leaf_index / _bits_per_unit;
    unsigned bit_index = leaf_index - unit_index * _bits_per_unit;
    split_unit_t unity = 1;
    split_unit_t bit_to_set = unity << bit_index;
    _bits[unit_index] |= bit_to_set;
}

inline Split::split_unit_t Split::getBits(unsigned int unit_index) const {
    assert(unit_index < _bits.size());
    return _bits[unit_index];
}

inline bool Split::getBitAt(unsigned int leaf_index) const {
    unsigned unit_index = leaf_index / _bits_per_unit;
    unsigned bit_index = leaf_index - unit_index * _bits_per_unit;
    split_unit_t unity = 1;
    split_unit_t bit_to_get = unity << bit_index;
    return static_cast<bool>(_bits[unit_index] & bit_to_get);
}

inline void Split::addSplit(const Split &other) {
    auto nunits = static_cast<unsigned>(_bits.size());
    assert(nunits == other._bits.size());
    for (unsigned i = 0; i < nunits; ++i) {
        _bits[i] |= other._bits[i];
    }
}

inline std::string Split::createPatternRepresentation() const {
    std::string s;
    unsigned int ntax_added = 0;
    split_unit_t unity = 1;
    split_unit_t zero = 0;
    for (unsigned long bitunit : _bits) {
        for (unsigned j = 0; j < _bits_per_unit; ++j) {
            split_unit_t bitmask = (unity << j);
            bool bit_is_set = (bitunit & bitmask) > zero;
            s += bit_is_set ? '*' : '-';
            if (++ntax_added == _nleaves) {
                break;
            }
        }
    }
    return s;
}

inline bool Split::isEquivalent(const Split &other) const {
    auto nunits = static_cast<unsigned>(_bits.size());
    assert(nunits > 0);

    // polarity 1 means that the root is on the same side in both splits
    // polarity 2 means the roots are inverted
    unsigned volatile polarity = 0;
    for (unsigned i = 0; i < nunits; ++i) {
        split_unit_t a = _bits[i];
        split_unit_t b = other._bits[i];
        bool a_equals_b = (a == b);
        bool a_equals_inverse_b = (a == ~b);
        if (i == nunits - 1) {
            a_equals_inverse_b = (a == (~b & _mask));
        }
        bool ok = (a_equals_b || a_equals_inverse_b);
        if (ok) {
            if (polarity == 0) {
                // First unit examined gives the polarity
                if (a_equals_b) {
                    polarity = 1;
                } else {
                    polarity = 2;
                }
            } else {
                // Polarity was already determined by the first unit in the loop
                if ((polarity == 1 && !a_equals_b) || (polarity == 2 && !a_equals_inverse_b)) {
                    return false;
                }
            }
        } else {
            return false;
        }
    }
    return true;
}

inline bool Split::isCompatibleWith(const Split &other) const {
    for (unsigned i = 0; i < _bits.size(); ++i) {
        split_unit_t a = _bits[i];
        split_unit_t b = other._bits[i];
        split_unit_t a_and_b = (a & b);
        bool equals_a = (a_and_b == a);
        bool equals_b = (a_and_b == b);
        if (a_and_b && !(equals_a || equals_b)) {
            return false;
        }
    }
    return true;
}

inline bool Split::conflictsWith(const Split &other) const {
    return !isCompatibleWith(other);
}

}// namespace strom
