//
// Created by Jonathan Ebrahimian on 7/16/20.
//
/*
AUTHORS
 Original Python:
    Leland McInnes <leland.mcinnes@gmail.com>
 C++ translation & optimization:
    Jonathan Ebrahimian <jebrahimian@mail.smu.edu>
    Connor Meehan <connor.gw.meehan@gmail.com>
    Stephen Meehan <swmeehan@stanford.edu>

ALGORITHMS
   1 Dong, Wei, Charikar, Moses, Li, Kai; 
     Efficient K-Nearest Neighbor Graph Construction for Generic Similarity Measures; 
     https://www.cs.princeton.edu/cass/papers/www11.pdf 
   2 Dasgupta, Sanjoy, and Freund, Yoav; 
     Random projection trees and low dimensional manifolds;
     https://cseweb.ucsd.edu/~dasgupta/papers/rptree-stoc.pdf 

Provided by suh ( Stanford University's Herzenberg Lab)
License: BSD 3 clause
*/

#ifndef SUH_RPTREE_H
#define SUH_RPTREE_H
#include "suh.h"

namespace suh {
    class RpTreeNode;
    class KnnDescent;
    class RpForest;
    class RpForester;

    class RpTree {
        friend KnnDescent;
        friend RpForester;
        friend RpForest;

    public:
        ~RpTree() {
            //test_serialization();
        }
    private:
        RpTree() {}
        RpTree(std::vector<std::vector<double> > hyperplanes_in, std::vector<double> offsets_in,
               std::vector<std::vector<int>> children_in, std::vector<std::vector<int>> indices_in) : hyperplanes(
                hyperplanes_in), offsets(offsets_in), children(children_in), indices(indices_in) {};
        std::vector<std::vector<double> > hyperplanes;
        std::vector<double> offsets;
        std::vector<std::vector<int>> children;
        std::vector<std::vector<int>> indices;


        operator std::string () const;
        std::string to_string() const;
        void serialize(const char *file_name="/Users/swmeehan/flatTree.bin");
        static RpTree deserialize(const char *file_name="/Users/swmeehan/flatTree.bin");
        static RpTree *deserialize(std::istream  &in);
        void in(std::istream &in);

        inline bool operator !=(const RpTree &that){
            bool good=*this==that;
            return !good;
        }
        inline bool operator ==(const RpTree &that){
            return count_unequal(hyperplanes, that.hyperplanes)==0
                   && count_unequal(offsets, that.offsets)==0
                   && children == that.children
                   && indices == that.indices;
        }

        bool test_serialization(bool do_output=false);

        friend std::ostream &operator << (std::ostream &os, const RpTree &f);
        friend std::istream &operator >> (std::istream &is, RpTree &f);
    };


    //Random projection forest
    class RpForest{
        friend RpForester;
        RpForest(){}
        std::vector<RpTree *>grove;
        void test_serialization();
        std::string to_string();
        static RpForest *deserialize(std::string in);
        bool operator ==(const RpForest &other) const;
        long duration;
    public:
        ~RpForest();

        inline const std::vector<RpTree *>&get_grove() const {
            return grove;
        }
        inline long get_duration() const {
            return duration;
        }
    };

    using RpForestPtr=std::shared_ptr<RpForest>;

    class RpForester {
    friend KnnDescent;
        static RpForestPtr make_forest(double **, int, int, int, int, long *, bool);

        static RpTreeNode *new_tree_node(double **, int, int, long *, int, bool);

        static RpTreeNode *new_angular_tree_node(const double **, const int, const int, std::vector<int> &, long *, int);

        static RpTreeNode *new_euclidean_tree_node(const double **, const int,const  int, std::vector<int> &, long * const, const int);

        static std::tuple<std::vector<int>, std::vector<int>, std::vector<double>, double>
        angular_random_projection_split(const double **,const  int, const int, std::vector<int> &, long *);

        static std::tuple<std::vector<int>, std::vector<int>, std::vector<double>, double>
        euclidean_random_projection_split(const double **, const int, const int, std::vector<int> &, long *);

        static RpTree *new_tree(RpTreeNode *tree_node, const int leaf_size);

        static std::pair<int, int>
        recursive_flatten(const RpTreeNode *, std::vector<std::vector<double> > &, std::vector<double> &,
                          std::vector<std::vector<int> > &, std::vector<std::vector<int> > &, int, int);

        static std::vector<std::vector<int>> leaf_array(const std::vector<RpTree *> &rp_forest);

        static std::vector<int>
        search_flat_tree(const double **const, const int, const int, const std::vector<std::vector<double> > &,
                         const std::vector<double> &, const std::vector<std::vector<int>> &,
                         const std::vector<std::vector<int>> &, long *const);

        static int
        select_side(const std::vector<double> &, const double, const double **const, const int, const int, long *const);

    };

    struct RpTreeNode {
        RpTreeNode(
                const std::vector<int> indices_in,
                const bool is_leaf_in,
                const std::vector<double> hyperplane_in,
                const double offset_in,
                const RpTreeNode *left_child_in,
                const RpTreeNode *right_child_in)
                : indices(std::move(indices_in)),
                  is_leaf(is_leaf_in),
                  hyperplane(std::move(hyperplane_in)),
                  offset(offset_in),
                  left_child(left_child_in),
                  right_child(right_child_in) {
        };

        ~RpTreeNode();

        const std::vector<int> indices;
        const bool is_leaf;
        const std::vector<double> hyperplane;
        const double offset;
        const RpTreeNode *left_child;
        const RpTreeNode *right_child;

        inline int num_nodes() const {
            return num_nodes_wrapper(this);
        }

        inline int num_nodes_wrapper(const RpTreeNode *tree_node) const {
            if (tree_node == nullptr || tree_node->is_leaf) {
                return 1;
            } else {
                return 1 + num_nodes_wrapper(tree_node->left_child) + num_nodes_wrapper(tree_node->right_child);
            }
        }

        inline int num_leaves() const {
            return num_leaves_wrapper(*this);
        }

        inline int num_leaves_wrapper(const RpTreeNode &tree_node) const {
            if (tree_node.is_leaf) {
                return 1;
            } else {
                return num_leaves_wrapper(*tree_node.left_child) + num_leaves_wrapper(*tree_node.right_child);
            }
        }

    };
}

#endif //CONVERSIONS_RP_TREE_H
