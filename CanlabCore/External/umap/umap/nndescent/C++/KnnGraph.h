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

Provided by suh (   Stanford University's Herzenberg Lab )
License: BSD 3 clause
*/

/* Nov 6 2020
 * Jonathan's C++ translation from Leland McInnes Python code gave this
 * class the same name used in Python:  utils
 * Stephen judged utils far too vague and over-general since the class's every method addresses a
 * specific 3D matrix of doubles that serves as a L nearest neighbor graph.  Thus utils became
 * KnnGraph.  Words like heap became nn_graph.  This is more meaningful to the code reader
 * and follows Google C++ naming for classes and variables
 * https://google.github.io/styleguide/cppguide.html#Naming
 *
 * At the same time Stephen committed this name change to git he move the class into the suh
 * namespace and started migrating the C++ towards
 * a) the Google C++ style guide https://google.github.io/styleguide/cppguide.html
 * b) Effective C++ idioms from Scott Myers https://www.aristeia.com/books.html
*
* The biggest visible effect of migration a) is  Google's naming conventions
* https://google.github.io/styleguide/cppguide.html#Namng
* One deviation is snake_use as does C++ std (e.g. std::string::push_back())
*  for class data names and class instance & static function names
*  instead of CamelCase or camelCase
*
* The biggest  impact of  migration b) is adding the keyword const everywhere.
*       This could help with compiler optimization for builtins that were previously non const references
*/

#ifndef SUH_NNGRAPH_H
#define SUH_NNGRAPH_H
#include "distances.h"
#include "RpTree.h"

namespace suh {
    class KnnDescent;
    class RpForester;
    class KnnGraph;
    using KnnGraphPtr = std::shared_ptr<suh::KnnGraph>;
    class KnnGraph {
        friend RpForester;
        friend KnnDescent;

        double ***  nn_graph_;

    public:
        int id=-1; // mainly for debugging purposes
        RpForestPtr rp_forest_;

        DistancePtr dist;
        const int rows_, columns_;
        KnnGraph(const KnnGraph &other);
        KnnGraph(KnnGraph &&other); // move rvalue
        inline operator double *** (){return nn_graph_;}
        inline  double *** graph(){return nn_graph_;}
        KnnGraph(double ***prior, const int rows, const int columns, DistancePtr dist, const bool copy);
        KnnGraph(const int rows, const int columns, const DistancePtr dist=std::shared_ptr<DistanceMetric>(), const bool fill=true);
        KnnGraph &operator = (const KnnGraph &);
        KnnGraph &operator = (KnnGraph &&);
        ~KnnGraph();

        static KnnGraph debug_move_copy(double ***prior, const int rows, const int columns);
        static KnnGraph debug_move_copy(const KnnGraph &in);

        virtual operator std::string() const {
            std::strstream s;
            s << "KnnGraph: 3 X " << rows_ << " X " << columns_;
            return s.str();
        }

        suh::MatrixPtr indices() const;
        suh::MatrixPtr distances() const;
        void copy_indices_column_wise(double *to) const;
        void copy_indices_MatLab(double *to) const;
        void copy_distances_column_wise(double *to) const;

    private:

        static void init(double ***const nn_graph, const int rows, const int columns);

        static int push(double ***const nn_graph, const int row, const double dist, const int index, const int flag,
                        const int columns);

        static int
        unchecked_push(double ***const nn_graph, const int row, const double dist, const int index, const int flag,
                       const int columns);

        static void rejection_sample(const int columns, const int rows, long *const rng_state, int *const indices);

        static inline int tau_rand_int(long *const rng_state) {
            rng_state[0] = (((rng_state[0] & 4294967294) << 12) & 0xFFFFFFFF) ^ (
                    (((rng_state[0] << 13) & 0xFFFFFFFF) ^ rng_state[0]) >> 19
            );

            rng_state[1] = (((rng_state[1] & 4294967288) << 4) & 0xFFFFFFFF) ^ (
                    (((rng_state[1] << 2) & 0xFFFFFFFF) ^ rng_state[1]) >> 25
            );

            rng_state[2] = (((rng_state[2] & 4294967280) << 17) & 0xFFFFFFFF) ^ (
                    (((rng_state[2] << 3) & 0xFFFFFFFF) ^ rng_state[2]) >> 11
            );
            return rng_state[0] ^ rng_state[1] ^ rng_state[2];
        }

        static inline int tau_rand_int(long *const rng_state, const int modulo) {
            int i = tau_rand_int(rng_state) % modulo;
            if (i < 0) {
                i += modulo;
            }
            return i;
        }

        static void new_build_candidates(double ***const nn_graph, const int, const int, const int, long *const rng_state,
                                         double ***const, double ***const, const double rho);

        static inline long double TauRand(long *const rng_state) {
            int integer = tau_rand_int(rng_state);
            long double tempDouble = integer;
            return abs(tempDouble / 0x7FFFFFFF);
        }

        static void sort(double ***const nn_graph, const int rows, const int columns);

        static int smallest_flagged(double ***const, const int, const int);

        static inline void sift_down(double *const nn_graph1, double *const nn_graph2, const int columns) {
            int left_child;
            int right_child;
            int swap;
            int elt = 0;
            double temp;
            double temp2;
            while (elt * 2 + 1 < columns) {
                left_child = elt * 2 + 1;
                right_child = left_child + 1;
                swap = elt;

                if (nn_graph1[swap] < nn_graph1[left_child]) {
                    swap = left_child;
                }

                if (right_child < columns && nn_graph1[swap] < nn_graph1[right_child]) {
                    swap = right_child;
                }

                if (swap == elt) {
                    break;
                } else {
                    temp = nn_graph1[swap];
                    temp2 = nn_graph1[elt];
                    nn_graph1[elt] = temp;
                    nn_graph1[swap] = temp2;

                    temp = nn_graph2[swap];
                    temp2 = nn_graph2[elt];
                    nn_graph2[elt] = temp;
                    nn_graph2[swap] = temp2;

                    elt = swap;
                }
            }
        }

        static inline void slice(double *const in, const int columns, double *const out) {
            for (int x = 0; x < columns; x++) {
                out[x] = in[x];
            }
        }

        static inline double norm(const double *const vec, const int columns) {
            double result = 0;
            for (const double *ptr = vec; ptr < vec + columns; ptr++) {
                result += pow(*ptr, 2);//TODO: make it a simple multiplication
            }
            return sqrt(result);
        }

        static inline double norm(const std::vector<double> &vec) {
            double result = 0;
            for (int i = 0; i < vec.size(); i++) {
                result += pow(vec[i], 2);//TODO: think about a bit shift
            }
            return sqrt(result);
        }


    };
}


#endif //SUH_NNGRAPH_H
