//
// Created by Jonathan Ebrahimian on 7/27/20.
//
/*
   AUTHORSHIP
   Primary Developer:     Jonathan Ebrahimian <jebrahimian@mail.smu.edu>
   Math Lead:             Connor Meehan <connor.gw.meehan@gmail.com>
   Secondary Developers:  Connor Meehan 
                          Stephen Meehan <swmeehan@stanford.edu>
   Bioinformatics Lead:   Wayne Moore <wmoore@stanford.edu>

   ALGHORITHMS
   1 Dong, Wei, Charikar, Moses, Li, Kai; 
     Efficient K-Nearest Neighbor Graph Construction for Generic Similarity Measures; 
     https://www.cs.princeton.edu/cass/papers/www11.pdf 
   2 Dasgupta, Sanjoy, and Freund, Yoav; 
     Random projection trees and low dimensional manifolds;
     https://cseweb.ucsd.edu/~dasgupta/papers/rptree-stoc.pdf 

   Provided by suh ( Stanford University's Herzenberg Lab)
   License: BSD 3 clause
*/

#ifndef SUH_RANDOMPROJECTIONTREENODE_H
#define SUH_RANDOMPROJECTIONTREENODE_H
#include <string>
#include <utility>
#include <vector>
namespace suh {
    struct RandomProjectionTreeNode {
        RandomProjectionTreeNode(
                std::vector<int> indices_in, bool is_leaf_in, std::vector<double> hyperplane_in,
                double offset_in, RandomProjectionTreeNode *left_child_in,
                RandomProjectionTreeNode *right_child_in)
                : indices(std::move(indices_in)),
                  is_leaf(is_leaf_in),
                  hyperplane(std::move(hyperplane_in)),
                  offset(offset_in),
                  left_child(left_child_in),
                  right_child(right_child_in) {
        };

        ~RandomProjectionTreeNode();

        std::vector<int> indices;
        bool is_leaf;
        std::vector<double> hyperplane;
        double offset;
        RandomProjectionTreeNode *left_child;
        RandomProjectionTreeNode *right_child;

        int num_nodes() {
            return num_nodes_wrapper(this);
        }

        int num_nodes_wrapper(RandomProjectionTreeNode *tree) {
            if (tree == nullptr || tree->is_leaf) {
                return 1;
            } else {
                return 1 + num_nodes_wrapper(tree->left_child) + num_nodes_wrapper(tree->right_child);
            }
        }

        int num_leaves() {
            return num_leaves_wrapper(*this);
        }

        int num_leaves_wrapper(RandomProjectionTreeNode &tree) {
            if (tree.is_leaf) {
                return 1;
            } else {
                return num_leaves_wrapper(*tree.left_child) + num_leaves_wrapper(*tree.right_child);
            }
        }

    };
}
#endif //SUH_RANDOMPROJECTIONTREENODE_H
