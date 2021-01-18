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


#include "RpTree.h"
#include "KnnGraph.h"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <algorithm> //find_if

using namespace suh;
RpTreeNode::~RpTreeNode() {
    if (this->left_child) {
        this->left_child = nullptr;
        delete this->left_child;
    }
    if (this->right_child) {
        this->right_child = nullptr;
        delete this->right_child;
    }
}

RpForestPtr RpForester::make_forest(double ** data, int data_size, int data_size_width, int n_neighbors, int n_trees, long * rng_state, bool angular){
    int leaf_size;
    auto startLeaf = std::chrono::high_resolution_clock::now();
    if(10 > n_neighbors){
        leaf_size = 10;
    }else{
        leaf_size = n_neighbors;
    }
    if (suh::debug_timing) {
        std::cout << "Building forest n_trees=" << n_trees << ", angular=" << angular << std::endl;
    }

    RpForest *rp_forest=new RpForest();
    try{
        for(int x = 0; x < n_trees; x++){
            rp_forest->grove.push_back(
                    new_tree(new_tree_node(data, data_size, data_size_width, rng_state, leaf_size, angular), leaf_size));
        }
    }catch(int x){
        std::cout << "Random Projection forest initialisation failed due to recursion"
                     "limit being reached. Something is a little strange with your "
                     "data, and this may take longer than normal to compute." << x << std::endl;
    }

    auto stopLeaf = std::chrono::high_resolution_clock::now();
    auto durationLeaf = std::chrono::duration_cast<std::chrono::microseconds>(stopLeaf - startLeaf);
    rp_forest->duration=durationLeaf.count();

    if (suh::debug_timing) {
        std::cout << "Forest computation cost " << durationLeaf.count() << " microseconds " << std::endl;
    }

    return RpForestPtr(rp_forest);
}

RpTreeNode * RpForester::new_tree_node(double ** data, int data_size, int data_size_width, long * rng_state, int leaf_size, bool angular){
    bool is_sparse = false; //TODO figure out what scipy.sparse.isspmatrix_csr is
    std::vector<int> indices(data_size);
    for(int x = 0; x < data_size;x++){
        indices[x] = x;
    }

    if(is_sparse){
        std::cout << "TODO: look at rp_tree code. Not implementing in MVP" << std::endl;
        return nullptr;
    }else{
        if(angular){
            return new_angular_tree_node((const double  **)data, data_size, data_size_width, indices, rng_state, leaf_size);//FIXME
        }else{
            return new_euclidean_tree_node((const double **)data, data_size, data_size_width, indices, rng_state, leaf_size);
        }
    }


}

/* TODO: Check if we could be doing more by reference and less by value to speed up execution
 *
 */
RpTreeNode * RpForester::new_angular_tree_node(const double ** data,const  int data_size,const  int data_size_width, std::vector<int>& indices, long * rng_state, int leaf_size){

    if(indices.size() > leaf_size){

        std::tuple<std::vector<int>,std::vector<int>,std::vector<double>,int> results = angular_random_projection_split(data, data_size, data_size_width, indices, rng_state);
        std::vector<int> left_indices = std::get<0>(results);
        std::vector<int> right_indices = std::get<1>(results);
        std::vector<double> hyperplane = std::get<2>(results);
        int offset = std::get<3>(results);

        RpTreeNode * left_node = new_angular_tree_node(data, data_size, data_size_width, left_indices, rng_state,
                                                       leaf_size);
        RpTreeNode * right_node = new_angular_tree_node(data, data_size, data_size_width, right_indices, rng_state,
                                                        leaf_size);

        RpTreeNode * node = new RpTreeNode({}, false, hyperplane, offset, left_node, right_node);//TODO check if this works
        return node;
    }else{
        RpTreeNode * node = new RpTreeNode(indices, true, {}, {}, nullptr, nullptr);//TODO check if this works
        return node;
    }
}

RpTreeNode * RpForester::new_euclidean_tree_node(const double ** data, const int data_size, const int data_size_width, std::vector<int>& indices, long * const rng_state , const int leaf_size){
    if(indices.size() > leaf_size){
        std::tuple<std::vector<int>,std::vector<int>,std::vector<double>, double> results = euclidean_random_projection_split(data, data_size, data_size_width, indices, rng_state);
        std::vector<int> left_indices = std::get<0>(results);
        std::vector<int> right_indices = std::get<1>(results);
        std::vector<double> hyperplane = std::get<2>(results);
        double offset = std::get<3>(results);

        RpTreeNode * left_node = new_euclidean_tree_node(data, data_size, data_size_width, left_indices, rng_state,
                                                         leaf_size);
        RpTreeNode * right_node = new_euclidean_tree_node(data, data_size, data_size_width, right_indices, rng_state,
                                                          leaf_size);

        RpTreeNode * node = new RpTreeNode({}, false, hyperplane, offset, left_node, right_node);//TODO check if this works
        return node;
    }else{
        RpTreeNode * node = new RpTreeNode(indices, true, {}, {}, nullptr, nullptr);//TODO check if this works
        return node;
    }

}
std::tuple<std::vector<int>,std::vector<int>,std::vector< double>, double> RpForester::angular_random_projection_split(
        const double ** data, const int data_size, const int data_size_width, std::vector<int>& indices, long * rng_state){

    int dim = data_size;
    int left_index = KnnGraph::tau_rand_int(rng_state) % int(indices.size());
    if(left_index < 0){
        left_index = left_index + indices.size();
    }
    int right_index = KnnGraph::tau_rand_int(rng_state) % int(indices.size());
    if(right_index < 0){
        right_index = right_index + indices.size();
    }
    if(left_index == right_index){
        right_index += 1;
    }
    right_index = right_index % int(indices.size());
    if(right_index < 0){
        right_index = right_index + indices.size();
    }




    int left = indices[left_index];
    int right = indices[right_index];
    double left_norm = KnnGraph::norm(data[left], data_size_width);
    double right_norm = KnnGraph::norm(data[right], data_size_width);

    if(abs(left_norm) < 1e-8){
        left_norm = 1;
    }

    if(abs(right_norm) < 1e-8){
        right_norm = 1;
    }


    std::vector<double> hyperplane_vector(dim);

    for(int d = 0; d < dim; d++){
        hyperplane_vector[d] = ((data[left][d]/left_norm) - (data[right][d] /right_norm));
    }

    double hyperplane_norm = KnnGraph::norm(hyperplane_vector);
    if(abs(hyperplane_norm) < 1e-8){
        hyperplane_norm = 1.0;
    }

    for(int d = 0; d < dim; d++){
        hyperplane_vector[d] = hyperplane_vector[d] / hyperplane_norm;
    }

    int n_left = 0;
    int n_right = 0;

    std::vector<int> side(indices.size(),0);

//    for(int x = 0; x < indices.size(); x++){
//        side.push_back(0);
//    }
    for(int i = 0; i < indices.size(); i++){
        double margin = 0.0;
        for(int d = 0; d < dim; d++){
            margin += hyperplane_vector[d] * data[indices[i]][d];
        }

        if(abs(margin) < 1e-8){
            //std::cout << "we are in" << std::endl;
            side[i] = abs(KnnGraph::tau_rand_int(rng_state)) % 2;
            if(side[i] < 0){
                side[i] = side[i] + 2;
            }
            if(side[i] == 0){
                n_left += 1;
            }else{
                n_right += 1;
            }
        }else if(margin > 0){
            side[i] = 0;
            n_left += 1;
        }else{
            side[i] = 1;
            n_right += 1;
        }
    }
    std::vector<int> indices_left (n_left,0);
    std::vector<int> indices_right(n_right,0);
    n_left = 0;
    n_right = 0;

    for(int i = 0; i < side.size();i++){
        if(side[i] == 0){
            indices_left[n_left] = indices[i];
            n_left += 1;
        }else{
            indices_right[n_right] = indices[i];
            n_right += 1;
        }
    }
    std::tuple<std::vector<int>,std::vector<int>,std::vector<double>,int> results (indices_left,indices_right,hyperplane_vector,0);
    return results;
}
std::tuple<std::vector<int>,std::vector<int>,std::vector< double>, double> RpForester::euclidean_random_projection_split(const double ** data, const int data_size, const int data_size_width, std::vector<int>& indices, long * rng_state){
    int dim = data_size_width;

    int left_index = KnnGraph::tau_rand_int(rng_state) % int(indices.size());
    if(left_index<0){
        left_index += int(indices.size());
    }

    int right_index = KnnGraph::tau_rand_int(rng_state) % int(indices.size());
    if(right_index < 0){
        right_index += int(indices.size());
    }

    if(left_index == right_index){
        right_index += 1;
    }
    right_index = right_index % int(indices.size());

    int left = indices[left_index];
    int right = indices[right_index];

    double left_norm = KnnGraph::norm(data[left], data_size_width);
    double right_norm = KnnGraph::norm(data[right], data_size_width);

    if(abs(left_norm) < 1e-8){
        left_norm = 1.0;
    }
    if(abs(right_norm) < 1e-8){
        right_norm = 1.0;
    }
    std::vector<double> hyperplane_vector (dim);
    //long double * hyperplane_vector = new long double [dim];
    double hyperplane_offset = 0.0;

    for(int d = 0; d < dim; d++){
        hyperplane_vector[d] = ((data[left][d])) - (data[right][d]);
        hyperplane_offset -= (hyperplane_vector[d] * (data[left][d] + data[right][d]) / 2.0);
    }

    int n_left = 0;
    int n_right = 0;
    double margin;
    std::vector<int> side(indices.size(),-999999);
    for(int i = 0; i < indices.size(); i++){
        margin = hyperplane_offset;
        for(int d = 0; d < dim; d++){
            margin += hyperplane_vector[d] * data[indices[i]][d];
        }

        if(abs(margin) < 1e-8){
            side[i] = abs(KnnGraph::tau_rand_int(rng_state)) % 2;
            if(side[i] == 0){
                n_left += 1;
            }else{
                n_right += 1;
            }
        }else if(margin > 0){
            side[i] = 0;
            n_left += 1;
        }else{
            side[i] = 1;
            n_right += 1;
        }
    }


    std::vector<int> indices_left (n_left,0);
    std::vector<int> indices_right(n_right,0);

    n_left = 0;
    n_right = 0;

    for(int i = 0; i < side.size();i++){
        if(side[i] == 0){
            indices_left[n_left] = indices[i];
            n_left += 1;
        }else{
            indices_right[n_right] = indices[i];
            n_right += 1;
        }
    }


    std::tuple<std::vector<int>,std::vector<int>,std::vector<double>, double> results (indices_left,indices_right,hyperplane_vector,hyperplane_offset);
    return results;
}

RpTree * RpForester::new_tree(RpTreeNode * tree_node, const int leaf_size){

    const int n_nodes = tree_node->num_nodes();
    const int n_leaves = tree_node->num_leaves();//TODO: couldnt I just do n_nodes - 1?

    std::vector<std::vector< double > > hyperplanes;
    std::vector< double > temp;
    if(false){
        //TODO: this is the sparse case that executes if hyperplane is more than 1 dimensional
    }else{
        for(int x = 0; x < n_nodes; x++){
            for(int y = 0; y < tree_node->hyperplane.size(); y++){
                temp.push_back(0);
            }
            hyperplanes.push_back(temp);
            temp.clear();
        }
    }

    std::vector<double> offsets(n_nodes,0);

    std::vector<int> tempInt;
    std::vector<std::vector<int>> children;
    for(int x = 0; x < n_nodes; x++){
        for(int y = 0; y < 2; y++){
            tempInt.push_back(-1);
        }
        children.push_back(tempInt);
        tempInt.clear();
    }

    std::vector<std::vector<int>> indices;
    for(int x = 0; x < n_leaves; x++){
        for(int y = 0; y < leaf_size; y++){
            tempInt.push_back(-1);
        }
        indices.push_back(tempInt);
        tempInt.clear();
    }
    recursive_flatten(tree_node, hyperplanes, offsets, children, indices, 0, 0);
    delete tree_node;
    return new RpTree(hyperplanes, offsets, children, indices);
}

std::pair<int,int> RpForester::recursive_flatten(const RpTreeNode * tree, std::vector<std::vector<double> > & hyperplanes, std::vector<double>& offsets, std::vector<std::vector<int> > & children, std::vector<std::vector<int> > & indices, int node_num, int leaf_num){
    if(tree->is_leaf){
        children[node_num][0] = -1 * leaf_num;
        for(int x = 0; x < tree->indices.size();x++){
            indices[leaf_num][x] = tree->indices[x];
        }
        leaf_num += 1;
        return std::make_pair(node_num, leaf_num);
    }else{
        if(false){
            //TODO sparce case
        }else{
            hyperplanes[node_num] = tree->hyperplane;
        }
        offsets[node_num] = tree->offset;
        children[node_num][0] = node_num + 1;
        int old_node_num = node_num;
        std::pair<int,int> results = recursive_flatten(tree->left_child,hyperplanes,offsets,children,indices,node_num+1,leaf_num);
        int node_num = results.first;
        int leaf_num = results.second;
        children[old_node_num][1] = node_num + 1;
        std::pair<int,int> results2 = recursive_flatten(tree->right_child,hyperplanes,offsets,children,indices,node_num+1,leaf_num);
        int node_num2 = results2.first;
        int leaf_num2 = results2.second;
        return std::make_pair(node_num2, leaf_num2);
    }
}

RpTree::operator std::string() const {
    return to_string();
}
std::string RpTree::to_string() const {
    std::strstream s;
    s << *this;
    return s.str();
}

std::ostream &suh::operator << (std::ostream &os, const RpTree &f){
#define O1(o) \
       os<<#o<<","<<f.o.size() << "," << f.o

#define O2(o) \
       os<<#o<<","<<f.o[0].size() << "," << f.o

    os.precision(15);
    O2(hyperplanes);
    O1(offsets);
    O2(children);
    O2(indices);
    return os;
}

std::istream &suh::operator >> (std::istream &is, RpTree &f){
    f.in(is);
 return is;
}
 RpTree *RpTree::deserialize(std::istream &s){
#define GRAB(T, V)\
    std::getline(s,line);\
     suh::Deserialize<T> V(line); \
if (V.columns_==0)       \
    return nullptr;     \
if (suh::strcmpi(#V, V.name.c_str()) !=0)       { \
    std::cerr <<"deserialize expected "<< #V << " and got " <<V.name<<"??"<<std::endl;          \
    throw std::invalid_argument("Disorder in FlatTree serialization"); }


     std::string line;
     GRAB(double, hyperplanes);
    GRAB(double, offsets);
    GRAB(int, children);
     GRAB(int, indices);
     RpTree *ft=new RpTree(hyperplanes.v2, offsets.v2[0], children.v2, indices.v2);
     return ft;
}

void RpTree::in(std::istream &s){
#define GRAB2(T, V)\
    std::getline(s,line);\
     suh::Deserialize<T> V(line); \
if (V.columns_==0)       \
    return ;     \
if (suh::strcmpi(#V, V.name.c_str()) !=0)      { \
    std::cerr <<"in() expected "<< #V << " and got " <<V.name<<"??"<<std::endl;          \
    throw std::invalid_argument("Disorder in FlatTree serialization"); }\
 this->V=V.v2

    std::string line;
    //std::istream &ii=std::cin;
    //std::getline(ii, line);
    GRAB2(double, hyperplanes);
    GRAB2(double, offsets)[0];
    GRAB2(int, children);
    GRAB2(int, indices);
}


std::vector<int> &&indexesOf(std::vector<int>&in, const int find) {
    std::vector<int> found;

    std::vector<int>::iterator it = in.begin();
    while ((it = std::find_if(it, in.end(), [=](int x) { return x == find; })) != in.end()) {
        found.push_back(std::distance(in.begin(), it));
        it++;
    }
    return std::move(found);
}
bool RpTree::test_serialization(bool do_output) {

 #define OUT(V) \
    std::cout << #V << " "\
    << suh::count_unequal(V, p->V) << " dissimilarities, "\
    << suh::count_unequal(V, p->V, 0) \
    << " inequalities." << std::endl

#define OUT_INT(V) \
    std::cout << #V << " "\
    << suh::count_unequal(V, p->V) << " inequalities." << std::endl

    bool wholeTruthAndNothingButTheTruth=false;
    std::string in=(std::string)*this;
    try {
        std::stringstream s(in);
        RpTree *p=deserialize(s);
        wholeTruthAndNothingButTheTruth = *this == *p;
        if (!do_output) {
            std::cout << "good=" << wholeTruthAndNothingButTheTruth << std::endl;
            return wholeTruthAndNothingButTheTruth;
        }
        try {
            RpTree *p2 = deserialize(s);
        } catch (std::invalid_argument argument) {
            std::cerr << "...Happened at string end:  " << in << std::endl;
            //std::vector<int>found=indexesOf(p->indices, 15);
            //std::cout<<found<<std::endl;
            return false;
        }
        std::cout << "serialization on whole is " << wholeTruthAndNothingButTheTruth << std::endl;
        OUT(hyperplanes);
        OUT(offsets);
        OUT_INT(indices);
        OUT_INT(children);/*
    std::stringstream s2(in);
    FlatTree *ft=new FlatTree;

        s2 >> *ft;
        wholeTruthAndNothingButTheTruth= *this == *ft;*/
    } catch (std::invalid_argument argument) {
        std::cerr << "...ouch" << in << std::endl;
        std::stringstream s(in);
        RpTree *p = deserialize(s);
        return false;
    }
    return true;
}

RpForest::~RpForest(){
    if (suh::debug_ctor_dtor)
        std::cout << "Destructing forest of " << grove.size() << " trees" << std::endl;
    for (auto tree:grove){
        delete tree;
    }
}
using namespace std::chrono;

void RpForest::test_serialization(){
    std::string saved=to_string();
    auto startLeaf = high_resolution_clock::now();
    RpForest *newTree=deserialize(saved);
    auto stopLeaf = high_resolution_clock::now();
    auto durationLeaf = duration_cast<microseconds>(stopLeaf - startLeaf);
    std::cout << "Forest deserialization cost " << durationLeaf.count() << " microseconds " << std::endl;

    std::cout << "Forest serialization ";
    if (*this==*newTree)
        std::cout << " works" << std::endl;
    else
        std::cout << " fails" << std::endl;
}


std::string RpForest::to_string(){
    std::stringstream out;
    for (int i=0;i<grove.size();i++){
        out << *grove[i];
    }
    return out.str();
}

bool RpForest::operator ==(const RpForest &that) const{
    for (int i=0;i<grove.size();i++){
        if (*grove[i] != *that.grove[i]) return false;
    }
    return true;
}
RpForest *RpForest::deserialize(std::string in) {
    RpForest *rf=new RpForest;
    std::stringstream ss(in);
    try{
        for (RpTree *ft = RpTree::deserialize(ss);
             ft != nullptr;
             ft = RpTree::deserialize(ss)) {
            rf->grove.push_back(ft);
        }
    } catch (std::invalid_argument argument) {
        if (rf->grove.empty()) {
            delete rf;
            rf = nullptr;
        }
    }
    return rf;
}


std::vector<std::vector<int> > RpForester::leaf_array(const std::vector<RpTree *> & rp_forest){
    std::vector<std::vector<int> > leaf_array, foo;
    std::vector<int> temp;
    const int N=rp_forest.size();
    if(N> 0){
        for(int x = 0; x < N;x++){
            append<int>(leaf_array, rp_forest[x]->indices);
        }
    }else{
        temp.push_back(-1);
        leaf_array.push_back(temp);
    }

    return leaf_array;
}

std::vector<int> RpForester::search_flat_tree(
        const double** const self_data, const int pos, const int columns,
        const std::vector<std::vector<double>> & hyperplanes,
        const std::vector<double> & offsets,
        const std::vector<std::vector<int>> & children,
        const std::vector<std::vector<int>> & indices,
        long *const rng_state) {
    int node = 0;
    while(children[node][0] > 0){
        int side = select_side(hyperplanes[node], offsets[node], self_data, pos, columns, rng_state);
        if(side == 0){
            node = children[node][0];
        }else{
            node = children[node][1];
        }
    }
    return indices[0-children[node][0]];
}

int RpForester::select_side(const std::vector<double> & hyperplane, const double offset, const double ** const query_points, const int pos, const int data_size_width, long* const rng_state) {
    double margin = offset;
    for(int d = 0; d < data_size_width; d++){
        margin += hyperplane[d] * query_points[pos][d];
    }
    if(abs(margin) < 1e-8){
        int side = abs(KnnGraph::tau_rand_int(rng_state)) % 2;
        if(side < 0){
            side = side + 2;
        }
        if(side == 0){
            return 0;
        } else{
            return 1;
        }
    }else if(margin > 0){
        return 0;
    }else{
        return 1;
    }
}
