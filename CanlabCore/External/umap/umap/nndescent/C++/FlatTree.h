//
// Created by Jonathan Ebrahimian on 7/28/20.
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

#ifndef SUH_FLATTREE_H
#define SUH_FLATTREE_H
#include <vector>
#include "suh.h"
namespace suh {
    class FlatTree {
    public:
        FlatTree() {}
        FlatTree(std::vector<std::vector<double> > hyperplanes_in, std::vector<double> offsets_in,
                 std::vector<std::vector<int>> children_in, std::vector<std::vector<int>> indices_in) : hyperplanes(
                hyperplanes_in), offsets(offsets_in), children(children_in), indices(indices_in) {};
        std::vector<std::vector<double> > hyperplanes;
        std::vector<double> offsets;
        std::vector<std::vector<int>> children;
        std::vector<std::vector<int>> indices;


        operator std::string () const;
        std::string to_string() const;
        void serialize(const char *file_name="/Users/swmeehan/flatTree.bin");
        static FlatTree deserialize(const char *file_name="/Users/swmeehan/flatTree.bin");
        static FlatTree *deserialize(std::istream  &in);
        void in(std::istream &in);

        inline bool operator !=(const FlatTree &that){
            bool good=*this==that;
            return !good;
        }
        inline bool operator ==(const FlatTree &that){
            return count_unequal(hyperplanes, that.hyperplanes)==0
                   && count_unequal(offsets, that.offsets)==0
                   && children == that.children
                   && indices == that.indices;
        }

        bool test_serialization(bool do_output=false);
        ~FlatTree() {
            //test_serialization();
        }
        friend std::ostream &operator << (std::ostream &os, const FlatTree &f);
        friend std::istream &operator >> (std::istream &is, FlatTree &f);
    };
}
#endif //SUH_FLATTREE_H
