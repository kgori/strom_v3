//
// Created by Kevin Gori on 04/11/2019.
//

#pragma once

#include "node.hpp"
#include <iostream>
#include <memory>

namespace strom {

class TreeManip;
//class Likelihood;
//class Updater;

class Tree {

    friend class TreeManip;
    //friend class Likelihood;
    //friend class Updater;

public:
    Tree();

    [[nodiscard]] bool isRooted() const;

    [[nodiscard]] unsigned numLeaves() const;

    [[nodiscard]] unsigned numInternals() const;

    [[nodiscard]] unsigned numNodes() const;

private:
    void clear();

    bool _is_rooted;
    Node *_root;
    unsigned _nleaves;
    unsigned _ninternals;
    Node::PtrVector _preorder;
    Node::PtrVector _levelorder;
    Node::Vector _nodes;

public:
    typedef std::shared_ptr<Tree> SharedPtr;
};

inline Tree::Tree() {
    clear();
}

inline void Tree::clear() {
    _is_rooted = false;
    _root = nullptr;
    _nleaves = 0;
    _ninternals = 0;
    _nodes.clear();
    _preorder.clear();
    _levelorder.clear();
}

inline bool Tree::isRooted() const {
    return _is_rooted;
}

inline unsigned Tree::numLeaves() const {
    return _nleaves;
}

inline unsigned Tree::numInternals() const {
    return _ninternals;
}

inline unsigned Tree::numNodes() const {
    return (unsigned) _nodes.size();
}

}// namespace strom
