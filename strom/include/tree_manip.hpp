//
// Created by Kevin Gori on 27/07/2021.
//

#pragma once

#include "tree.hpp"
#include "xstrom.hpp"
#include <cassert>
#include <fmt/core.h>
#include <memory>
#include <queue>
#include <range/v3/view/reverse.hpp>
#include <regex>
#include <set>
#include <stack>

namespace strom {

class TreeManip {
public:
    TreeManip();

    explicit TreeManip(Tree::SharedPtr t);

    ~TreeManip();

    void setTree(Tree::SharedPtr t);

    Tree::SharedPtr getTree();

    [[nodiscard]] double calcTreeLength() const;

    [[nodiscard]] unsigned countEdges() const;

    void scaleAllEdgeLengths(double scaler);

    void createTestTree();

    [[nodiscard]] std::string makeNewick(unsigned precision, bool use_names = false) const;

    void buildFromNewick(const std::string &newick, bool rooted, bool allow_polytomies);

    void rerootAtNodeNumber(int node_number);

    void clear();

private:
    Node *findNextPreorder(Node *nd);

    void refreshPreorder();

    void refreshLevelOrder();

    void renumberInternals();

    void rerootAtNode(Node *prospective_root);

    void extractNodeNumberFromName(Node *nd, std::set<unsigned> &used);

    void extractEdgeLen(Node *nd, const std::string &edge_length_string);

    [[nodiscard]] unsigned countNewickLeaves(const std::string &newick) const;

    void stripOutNexusComments(std::string &newick) const;

    bool canHaveSibling(Node *nd, bool rooted, bool allow_polytomies);

    Tree::SharedPtr _tree;

public:
    typedef std::shared_ptr<TreeManip> SharedPtr;
};

inline TreeManip::TreeManip() {
    std::cout << "Constructing a TreeManip" << std::endl;
    clear();
}

inline TreeManip::TreeManip(Tree::SharedPtr t) {
    std::cout << "Constructing a TreeManip with a supplied tree" << std::endl;
    clear();
    setTree(t);
}

inline TreeManip::~TreeManip() {
    std::cout << "Deleting a TreeManip" << std::endl;
}

inline void TreeManip::clear() {
    _tree.reset();
}

inline void TreeManip::setTree(Tree::SharedPtr t) {
    assert(t);
    _tree = t;
}

inline Tree::SharedPtr TreeManip::getTree() {
    return _tree;
}

inline double TreeManip::calcTreeLength() const {
    double tree_length = 0;
    for (auto nd : _tree->_preorder) {
        tree_length += nd->_edge_length;
    }
    return tree_length;
}

inline unsigned TreeManip::countEdges() const {
    return static_cast<unsigned>(_tree->_preorder.size());
}

inline void TreeManip::scaleAllEdgeLengths(double scaler) {
    for (auto nd : _tree->_preorder) {
        nd->setEdgeLength(nd->getEdgeLength() * scaler);
    }
}

inline void TreeManip::createTestTree() {
    clear();
    _tree = std::make_shared<Tree>();
    _tree->_nodes.resize(6);

    Node *root_node = &_tree->_nodes[0];
    Node *first_internal = &_tree->_nodes[1];
    Node *second_internal = &_tree->_nodes[2];
    Node *first_leaf = &_tree->_nodes[3];
    Node *second_leaf = &_tree->_nodes[4];
    Node *third_leaf = &_tree->_nodes[5];

    // Here is the structure of the tree (numbers in
    // parentheses are node numbers, other numbers
    // are edge lengths):
    //
    // first_leaf (0)   second_leaf (1)   third_leaf (2)
    //      \              /                  /
    //       \ 0.1        / 0.1              /
    //        \          /                  /
    //     second_internal (3)             / 0.2
    //             \                      /
    //              \ 0.1                /
    //               \                  /
    //                first_internal (4)
    //                        |
    //                        | 0.1
    //                        |
    //                    root_node (5)
    //

    root_node->_parent = nullptr;
    root_node->_left_child = first_internal;
    root_node->_right_sib = nullptr;
    root_node->_number = 5;
    root_node->_name = "root_node";
    root_node->_edge_length = 0.0;

    first_internal->_parent = root_node;
    first_internal->_left_child = second_internal;
    first_internal->_right_sib = nullptr;
    first_internal->_number = 4;
    first_internal->_name = "first_internal_node";
    first_internal->_edge_length = 0.1;

    second_internal->_parent = first_internal;
    second_internal->_left_child = first_leaf;
    second_internal->_right_sib = third_leaf;
    second_internal->_number = 3;
    second_internal->_name = "second_internal_node";
    second_internal->_edge_length = 0.1;

    first_leaf->_parent = second_internal;
    first_leaf->_left_child = nullptr;
    first_leaf->_right_sib = second_leaf;
    first_leaf->_number = 0;
    first_leaf->_name = "first_leaf";
    first_leaf->_edge_length = 0.1;

    second_leaf->_parent = second_internal;
    second_leaf->_left_child = nullptr;
    second_leaf->_right_sib = nullptr;
    second_leaf->_number = 1;
    second_leaf->_name = "second_leaf";
    second_leaf->_edge_length = 0.1;

    third_leaf->_parent = first_internal;
    third_leaf->_left_child = nullptr;
    third_leaf->_right_sib = nullptr;
    third_leaf->_number = 2;
    third_leaf->_name = "third_leaf";
    third_leaf->_edge_length = 0.2;

    _tree->_is_rooted = true;
    _tree->_root = root_node;
    _tree->_nleaves = 3;

    _tree->_preorder.push_back(first_internal);
    _tree->_preorder.push_back(second_internal);
    _tree->_preorder.push_back(first_leaf);
    _tree->_preorder.push_back(second_leaf);
    _tree->_preorder.push_back(third_leaf);

    _tree->_levelorder.push_back(first_internal);
    _tree->_levelorder.push_back(second_internal);
    _tree->_levelorder.push_back(third_leaf);
    _tree->_levelorder.push_back(first_leaf);
    _tree->_levelorder.push_back(second_leaf);
}

inline std::string TreeManip::makeNewick(unsigned precision, bool use_names) const {
    std::string newick;
    const auto tip_node_name_format = fmt::format("{{:s}}:{{:.{:d}f}}", precision);
    const auto tip_node_number_format = fmt::format("{{:d}}:{{:.{:d}f}}", precision);
    const auto internal_node_format = fmt::format("):{{:.{:d}f}}", precision);
    std::stack<Node *> node_stack;

    Node *root_tip = (_tree->_is_rooted ? nullptr : _tree->_root);

    for (auto nd : _tree->_preorder) {
        if (nd->_left_child) {
            newick += "(";
            node_stack.push(nd);
            if (root_tip) {
                if (use_names) {
                    newick += fmt::format(tip_node_name_format, root_tip->_name, nd->_edge_length);
                } else {
                    newick += fmt::format(tip_node_number_format, root_tip->_number + 1, nd->_edge_length);
                }

                newick += ",";
                root_tip = nullptr;
            }
        } else {
            if (use_names) {
                newick += fmt::format(tip_node_name_format, nd->_name, nd->_edge_length);
            } else {
                newick += fmt::format(tip_node_number_format, nd->_number + 1, nd->_edge_length);
            }
            if (nd->_right_sib) {
                newick += ",";
            } else {
                Node *popped = (node_stack.empty() ? nullptr : node_stack.top());
                while (popped && !popped->_right_sib) {
                    node_stack.pop();
                    if (node_stack.empty()) {
                        newick += ")";
                        popped = nullptr;
                    } else {
                        newick += fmt::format(internal_node_format, popped->_edge_length);
                        popped = node_stack.top();
                    }
                }
                if (popped && popped->_right_sib) {
                    node_stack.pop();
                    newick += fmt::format(internal_node_format, popped->_edge_length);
                    newick += ",";
                }
            }
        }
    }
    return newick + ";";
}

inline void TreeManip::extractNodeNumberFromName(Node *nd, std::set<unsigned> &used) {
    assert(nd);
    unsigned x = 0;
    try {
        x = std::stoi(nd->_name);
    } catch (std::invalid_argument &) {
        // node name could not be converted to an integer value
        throw XStrom(fmt::format(FMT_STRING("node name {:s} not interpretable as a positive integer"), nd->_name));
    }

    // attempt to insert x into the set of already used node numbers
    std::pair<std::set<unsigned>::iterator, bool> insert_result = used.insert(x);
    if (insert_result.second) {
        // insertion was made, so x has not already been used
        nd->_number = x - 1;
    } else {
        throw XStrom(fmt::format(FMT_STRING("leaf number {:d} used more than once"), x));
    }
}

inline void TreeManip::extractEdgeLen(Node *nd, const std::string &edge_length_string) {
    assert(nd);
    double d = 0.0;
    try {
        d = std::stod(edge_length_string);
    } catch (std::invalid_argument &) {
        throw XStrom(fmt::format(FMT_STRING("{:s} is not interpretatble as a floating point number"), edge_length_string));
    }
    nd->setEdgeLength(d);
}

inline unsigned TreeManip::countNewickLeaves(const std::string &newick) const {
    std::regex taxonexpr("[(,]\\s*(\\d+|\\S+?|['].+?['])\\s*(?=[,):])");
    std::sregex_iterator m1(newick.begin(), newick.end(), taxonexpr);
    std::sregex_iterator m2;
    return static_cast<unsigned>(std::distance(m1, m2));
}

inline void TreeManip::stripOutNexusComments(std::string &newick) const {
    std::regex commentexpr("\\[.*?\\]");
    newick = std::regex_replace(newick, commentexpr, std::string(""));
}

inline Node *TreeManip::findNextPreorder(Node *nd) {
    assert(nd);
    Node *next = nullptr;
    if (!nd->_left_child && !nd->_right_sib) {
        // nd has no children and no siblings, so the next preorder node must be the right sibling
        // of the first ancestral node that has a right sibling
        Node *anc = nd->_parent;
        while (anc && !anc->_right_sib) {
            anc = anc->_parent;
        }
        if (anc) {
            // One was found
            next = anc->_right_sib;
        } else {
            // nd is the last preorder node in the tree
            next = nullptr;
        }
    } else if (nd->_right_sib && !nd->_left_child) {
        // nd has no children, but it does have a right sibling
        next = nd->_right_sib;
    } else if (nd->_left_child) {
        next = nd->_left_child;
    }
    return next;
}

/*
 * Rebuild the preorder vector of Node pointers.
 * This must be called any time the structure of the tree is changed
 */
inline void TreeManip::refreshPreorder() {
    // Create a vector of Nodes in preorder sequence
    _tree->_preorder.clear();
    _tree->_preorder.reserve(_tree->_nodes.size() - 1);// _preorder does not include the root node

    if (!_tree->_root) {
        return;
    }

    Node *first_preorder = _tree->_root->_left_child;

    // sanity check: the first preorder node should be the only child of the root
    assert(first_preorder->_right_sib == nullptr);

    Node *nd = first_preorder;
    _tree->_preorder.push_back(nd);

    while (true) {
        nd = findNextPreorder(nd);
        if (nd) {
            _tree->_preorder.push_back(nd);
        } else {
            break;
        }
    }
}

inline void TreeManip::refreshLevelOrder() {
    if (!_tree->_root) {
        return;
    }

    std::queue<Node *> q;

    _tree->_levelorder.clear();
    _tree->_levelorder.reserve(_tree->_nodes.size() - 1);

    Node *nd = _tree->_root->_left_child;

    assert(nd->_right_sib == nullptr);

    q.push(nd);

    while (!q.empty()) {
        nd = q.front();
        q.pop();

        _tree->_levelorder.push_back(nd);

        Node *child = nd->_left_child;
        if (child) {
            q.push(child);
            child = child->_right_sib;
            while (child) {
                q.push(child);
                child = child->_right_sib;
            }
        }
    }
}

inline void TreeManip::renumberInternals() {
    assert(!_tree->_preorder.empty());

    // Renumber internal nodes in preorder sequence
    unsigned curr_internal = _tree->_nleaves;
    for (auto nd : ranges::reverse_view(_tree->_preorder)) {
        if (nd->_left_child) {
            // nd is an internal node
            nd->_number = curr_internal++;
        }
    }

    // Root node is not included in the preorder vector, so if the root is
    // internal, renumber it here
    if (_tree->_is_rooted) {
        _tree->_root->_number = curr_internal++;
    }

    _tree->_ninternals = curr_internal - _tree->_nleaves;

    // If the tree has polytomies, then there are Nodes in the tree that
    // have not yet been numbered. Their _number is currently set to -1,
    // so we find nodes with _number equal to -1, and renumber them here
    for (auto &nd : _tree->_nodes) {
        if (nd._number == -1) {
            nd._number = curr_internal++;
        }
    }
}

inline bool TreeManip::canHaveSibling(Node *nd, bool rooted, bool allow_polytomies) {
    assert(nd);

    if (!nd->_parent) {
        // trying to give root node a sibling is not valid
        return false;
    }

    if (allow_polytomies) {
        // A tree that allows polytomies can always add siblings (provided the node
        // is not the root node)
        return true;
    }

    bool nd_can_have_sibling = true;
    if (nd != nd->_parent->_left_child) {
        if (nd->_parent->_parent) {
            // trying to give a sibling to a sibling of nd, and nd's parent is not the root
            nd_can_have_sibling = false;
        } else {
            if (rooted || nd != nd->_parent->_left_child->_right_sib) {
                // root node has exactly two children in a rooted tree, or trying to give a root node more than three children
                nd_can_have_sibling = false;
            }
        }
    }
    return nd_can_have_sibling;
}

inline void TreeManip::rerootAtNodeNumber(int node_number) {
    // Locate the node with the given node_number
    Node *nd = nullptr;
    for (auto &curr : _tree->_nodes) {
        if (curr._number == node_number) {
            nd = &curr;
            break;
        }
    }

    if (!nd) {
        throw XStrom(fmt::format(FMT_STRING("No node found with node number {:d}"), node_number));
    }

    if (nd != _tree->_root) {
        if (nd->_left_child) {
            throw XStrom(fmt::format(FMT_STRING("Cannot currently root trees at internal nodes (e.g. node {:d}"), nd->_number));
        }
        rerootAtNode(nd);
    }
}

inline void TreeManip::rerootAtNode(Node *prospective_root) {
    Node *a = prospective_root;
    Node *b = prospective_root->_parent;
    Node *c = nullptr;
    Node *d = nullptr;
    Node *p = nullptr;
    a->_parent = nullptr;
    double tmp_edgelen = 0.0;
    double prev_edgelen = a->getEdgeLength();

    while (b) {
        // Prune node a from b
        if (a == b->_left_child) {
            if (a->_right_sib) {
                b->_left_child = a->_right_sib;
                a->_right_sib = nullptr;
            } else {
                b->_left_child = nullptr;
            }
        } else {
            c = b->_left_child;
            while (c->_right_sib != a) {
                c = c->_right_sib;
            }
            d = a->_right_sib;
            c->_right_sib = d;
            a->_right_sib = nullptr;
        }

        // Graft node b onto node a (but don't unhook b from its parent yet)
        if (a->_left_child) {
            c = a->_left_child;
            while (c->_right_sib) {
                c = c->_right_sib;
            }
            c->_right_sib = b;
        } else {
            a->_left_child = b;
        }

        // Rotate
        p = a;
        a = b;
        b = b->_parent;
        a->_parent = p;

        // Swap nd's edge length with its new parent's edge length
        tmp_edgelen = a->getEdgeLength();
        a->setEdgeLength(prev_edgelen);
        prev_edgelen = tmp_edgelen;
    }
    prospective_root->setEdgeLength(0.0);
    _tree->_root = prospective_root;
    refreshPreorder();
    refreshLevelOrder();
}

inline void TreeManip::buildFromNewick(const std::string &newick, bool rooted, bool allow_polytomies) {
    _tree = std::make_shared<Tree>();
    _tree->_is_rooted = rooted;

    std::set<unsigned> used;// Ensure no two leaf nodes have the same number
    unsigned curr_leaf = 0;
    unsigned num_edge_lengths = 0;
    unsigned curr_node_index = 0;

    // Remove comments from newick string
    std::string commentless_newick = newick;
    stripOutNexusComments(commentless_newick);

    // Resize the nodes vector
    _tree->_nleaves = countNewickLeaves(commentless_newick);
    if (_tree->_nleaves < 4) {
        throw XStrom("Expecting newick tree description to have at least four leaves");
    }
    unsigned max_nodes = 2 * _tree->_nleaves - (rooted ? 0 : 2);
    _tree->_nodes.resize(max_nodes);

    for (auto &nd : _tree->_nodes) {
        nd._number = -1;
    }

    try {
        // Root node
        Node *nd = &_tree->_nodes[curr_node_index];
        _tree->_root = nd;

        if (_tree->_is_rooted) {
            nd = &_tree->_nodes[++curr_node_index];
            nd->_parent = &_tree->_nodes[curr_node_index - 1];
            nd->_parent->_left_child = nd;
        }

        // Define some flags for keeping track of operations
        enum
        {
            Prev_Tok_LParen = 0x01,// previous token was '('
            Prev_Tok_RParen = 0x02,// previous token was ')'
            Prev_Tok_Colon = 0x04, // previous token was ':'
            Prev_Tok_Comma = 0x08, // previous token was ','
            Prev_Tok_Name = 0x10,  // previous token was a node name
            Prev_Tok_Edgelen = 0x20// previous token was an edge length
        };
        auto previous = Prev_Tok_LParen;

        // Some useful flag combinations
        auto LParen_Valid = (Prev_Tok_LParen | Prev_Tok_Comma);
        auto RParen_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_Edgelen);
        auto Comma_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_Edgelen);
        auto Colon_Valid = (Prev_Tok_RParen | Prev_Tok_Name);
        auto Name_Valid = (Prev_Tok_LParen | Prev_Tok_RParen | Prev_Tok_Comma);

        // Set to true while reading an edge length
        bool inside_edge_length = false;
        std::string edge_length_str;
        unsigned edge_length_position = 0;

        // Set to true while reading a node name surrounded by single quotes
        bool inside_quoted_name = false;

        // Set to true while reading a node name, no quotes
        bool inside_unquoted_name = false;

        // Set to start of each node name
        unsigned node_name_position = 0;

        // Loop through the characters in the newick
        unsigned position_in_string = 0;
        for (auto ch : commentless_newick) {
            position_in_string++;

            if (inside_quoted_name) {
                if (ch == '\'') {
                    inside_quoted_name = false;
                    node_name_position = 0;
                    if (!nd->_left_child) {
                        extractNodeNumberFromName(nd, used);
                        curr_leaf++;
                    }
                    previous = Prev_Tok_Name;
                } else if (iswspace(ch)) {
                    nd->_name += ' ';
                } else {
                    nd->_name += ch;
                }
                continue;
            } else if (inside_unquoted_name) {
                if (ch == '(') {
                    throw XStrom(fmt::format(FMT_STRING("Unexpected left parenthesis inside node name at position {:d} in tree description"), node_name_position));
                }

                if (iswspace(ch) || ch == ':' || ch == ',' || ch == ')') {
                    inside_unquoted_name = false;

                    // Expect a node name only after a left paren, a comma, or a right paren
                    if (!(previous & Name_Valid)) {
                        throw XStrom(fmt::format(FMT_STRING("Unexpected node name ({:s}) at position {:d} in tree description"), nd->_name, node_name_position));
                    }

                    if (!nd->_left_child) {
                        extractNodeNumberFromName(nd, used);
                        curr_leaf++;
                    }
                    previous = Prev_Tok_Name;
                } else {
                    nd->_name += ch;
                    continue;
                }
            } else if (inside_edge_length) {
                if (ch == ',' || ch == ')' || iswspace(ch)) {
                    inside_edge_length = false;
                    edge_length_position = 0;
                    extractEdgeLen(nd, edge_length_str);
                    ++num_edge_lengths;
                    previous = Prev_Tok_Edgelen;
                } else {
                    // Floating point and scientific notation parsing
                    bool valid = (ch == 'e' || ch == 'E' || ch == '.' || ch == '-' || ch == '+' || isdigit(ch));
                    if (!valid) {
                        throw XStrom(fmt::format(FMT_STRING("Invalid branch length character {:c} at position {:d} in tree_description"),
                                                 ch, position_in_string));
                    }
                    edge_length_str += ch;
                    continue;
                }
            }

            if (iswspace(ch)) {
                continue;
            }

            switch (ch) {
                case ';':
                    break;

                case ')':
                    // If nd is bottommost node, expecting left paren or semicolon, but not right paren
                    if (!nd->_parent) {
                        throw XStrom(fmt::format(FMT_STRING("Too many right parentheses at position {:d} in tree description"), position_in_string));
                    }

                    // Expect right paren only after an edge length, a node name, or another right paren
                    if (!(previous & RParen_Valid)) {
                        throw XStrom(fmt::format(FMT_STRING("Unexpected right parenthesisat position {:d} in tree description"), position_in_string));
                    }
                    // Go down a level
                    nd = nd->_parent;
                    if (!nd->_left_child->_right_sib) {
                        throw XStrom(fmt::format(FMT_STRING("Internal node has only one child at position {:d} in tree description"), position_in_string));
                    }
                    previous = Prev_Tok_RParen;
                    break;

                case ':':
                    // Expect colon only after a node name or another right paren
                    if (!(previous & Colon_Valid)) {
                        throw XStrom(fmt::format(FMT_STRING("Unexpected colon at position {:d} in tree description"), position_in_string));
                    }
                    previous = Prev_Tok_Colon;
                    break;

                case ',':
                    // Expect comma only after an edge length, a node name, or a right paren
                    if (!nd->_parent || !(previous & Comma_Valid)) {
                        throw XStrom(fmt::format(FMT_STRING("Unexpected comma at position {:d} in tree description"), position_in_string));
                    }

                    // Check for polytomies
                    if (!canHaveSibling(nd, rooted, allow_polytomies)) {
                        throw XStrom(fmt::format(FMT_STRING("Polytomy found in the following tree description but polytomies prohibited:\n{:s}"), newick));
                    }

                    // Create the sibling
                    curr_node_index++;
                    if (curr_node_index == _tree->_nodes.size()) {
                        throw XStrom(fmt::format(FMT_STRING("Too many nodes specified by tree description ({:d} nodes allocated for {:d} leaves)"), _tree->_nodes.size(), _tree->_nleaves));
                    }
                    nd->_right_sib = &_tree->_nodes[curr_node_index];
                    nd->_right_sib->_parent = nd->_parent;
                    nd = nd->_right_sib;
                    previous = Prev_Tok_Comma;
                    break;

                case '(':
                    // Expect left paren only after a comma or another left paren
                    if (!(previous & LParen_Valid)) {
                        throw XStrom(fmt::format(FMT_STRING("Not expecting left parenthesis at position {:d} in tree description"), position_in_string));
                    }
                    // Create new node above and to the left of the current node
                    assert(!nd->_left_child);
                    curr_node_index++;
                    if (curr_node_index == _tree->_nodes.size()) {
                        throw XStrom(fmt::format(FMT_STRING("malformed tree description (more than {:d} nodes specified)"), _tree->_nodes.size()));
                    }
                    nd->_left_child = &_tree->_nodes[curr_node_index];
                    nd->_left_child->_parent = nd;
                    nd = nd->_left_child;
                    previous = Prev_Tok_LParen;
                    break;

                case '\'':
                    // Encountered an apostrophe, which always indicates the start of a
                    // node name (but note that node names do not have to be quoted)

                    // Expect node name only after a left paren (child's name), a comma (sib's name)
                    // or a right paren (parent's name)
                    if (!(previous & Name_Valid)) {
                        throw XStrom(fmt::format(FMT_STRING("Not expecting node name at position {:d} in tree description"), position_in_string));
                    }

                    // Get the rest of the name
                    nd->_name.clear();

                    inside_quoted_name = true;
                    node_name_position = position_in_string;

                    break;

                default:
                    // Get here if ch is not one of ();:,'

                    // Expecting either an edge length or an unquoted node name
                    if (previous == Prev_Tok_Colon) {
                        // Edge length expected (e.g. "235", "0.12345", "1.7e-3")
                        inside_edge_length = true;
                        edge_length_position = position_in_string;
                        edge_length_str = ch;
                    } else {
                        // Get the node name
                        nd->_name = ch;

                        inside_unquoted_name = true;
                        node_name_position = position_in_string;
                    }
            }// end of switch statement
        }

        if (inside_unquoted_name) {
            throw XStrom(fmt::format(FMT_STRING("Tree description ended before end of node name starting at position {:d} was found"), node_name_position));
        }
        if (inside_quoted_name) {
            throw XStrom(fmt::format(FMT_STRING("Expecting single quote to mark the end of node name at position {:d} in tree description"), node_name_position));
        }
        if (inside_edge_length) {
            throw XStrom(fmt::format(FMT_STRING("Tree description ended before end of edge length starting at position {:d} was found"), edge_length_position));
        }

        if (_tree->_is_rooted) {
            refreshPreorder();
            refreshLevelOrder();
        } else {
            // Root at leaf number 0
            rerootAtNodeNumber(0);
        }
        renumberInternals();
    } catch (XStrom &x) {
        clear();
        throw x;
    }
}

}// namespace strom