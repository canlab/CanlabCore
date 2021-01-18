//
// namespace suh contains C++ from Stanford University's Herzenberg Lab
//
// Created by Stephen Meehan on 10/29/20.
//
/*
AUTHOR
   Stephen Meehan <swmeehan@stanford.edu>

Provided by suh ( Stanford University's Herzenberg Lab)
License: BSD 3 clause
*/

/* Items in suh namespace are  migrating  towards
* a) the Google C++ style guide https://google.github.io/styleguide/cppguide.html
* b) Effective C++ idioms from Scott Meyers https://www.aristeia.com/books.html
*
* The biggest visible effect of migration a) is Google's naming conventions
* https://google.github.io/styleguide/cppguide.html#Namng
* One deviation is snake_use as does C++ std (e.g. std::string::push_back())
*  for class data names and class instance & static function names
*  instead of CamelCase or camelCase
*
* The biggest  impact of  migration b) is adding the keyword const everywhere.
*       This could help with compiler optimization for builtins that were previously non const references
*/

#ifndef SUH_H
#define SUH_H

#include <string>
#include <iostream>
#include <strstream>
#include <sstream>
#include <fstream>
#include <regex>
#include <memory>
#include <assert.h>
#include <vector>
#include <functional>
#include <cmath>
using std::pow;

namespace suh {
    static const bool debug_ctor_dtor = false, debug_file_io=true, debug_timing=false;

    using FncProgress=std::function<bool(const int iter, const int n_iters)>;

    inline double secs(long microseconds){
        return ((double)microseconds/1000000.0);
    }
    template<typename ... Args>
    std::string string_format( const std::string& format, Args ... args ){
        size_t size = snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
        if( size <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
        std::unique_ptr<char[]> buf( new char[ size ] );
        snprintf( buf.get(), size, format.c_str(), args ... );
        return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
    }

    inline std::string str_secs(long microseconds){
        return suh::string_format("%4.3f", ((double)microseconds/1000000.0));
    }

    template <typename T>std::string bank_num(T num) {
        return suh::string_format("%.2f", (double) num);
    }


    template<typename T>
    void append(std::vector<T> &a, const std::vector<T> &b) {
        a.insert(std::end(a), std::begin(b), std::end(b));
    }

    template<typename T>
    void append(std::vector<std::vector<T>> &a, const std::vector<std::vector<T>> &b) {
        a.insert(std::end(a), std::begin(b), std::end(b));
    }

    template<typename T>
    std::ostream &print( const std::vector<T> &in,  std::ostream &o = std::cout, const char *delim = ",") {
        const int N=in.size();
        for (int i=0;i<N;i++)
            if (i<N-1)
                o<<in[i]<<delim;
            else
                o<<in[i]<<std::endl;
        return o;
    }

    template <typename T> std::ostream& operator << ( std::ostream &o, const std::vector<T> &in){
        return print(in,  o);
    }

    template<typename T>
    std::ostream &print( const std::vector<std::vector<T>> &in, const bool flat=true,  std::ostream &o = std::cout, const char *delim = ",") {
        const int rows=in.size();
        if (flat)
            for (int row = 0; row < rows; row++) {
                const int columns = in[row].size();
                for (int col = 0; col < columns; col++)
                    o << in[row][col] << delim;
                if (row == rows - 1)
                    o << std::endl;
            }
        else
            for (int row=9; row < rows; row++)
                print(in[row], o, delim);
        return o;
    }

    template <typename T> std::ostream& operator << ( std::ostream &o, const std::vector<std::vector<T>> &in){
        return print(in, true, o);
    }

    template<typename T>
    double mean(const T *const nums, const int size) {
        double sum = 0;
        for (int i = 0; i < size; i++) {
            sum += nums[i];
        }
        return sum / (double) size;
    }

    inline int strcmpi(const char *l, const char *r) {
        int lenDiff = strlen(l) - strlen(r);
        if (lenDiff != 0) {
            return lenDiff;
        }
        while (*r && *l) {
            const int dif = tolower(*l) - tolower(*r);
            if (dif != 0) {
                return dif;
            }
            r++;
            l++;
        }
        return 0;
    }

    inline int  strcmpi(std::string l, std::string r) {
        int lenDiff = l.size() - r.size();
        if (lenDiff != 0) {
            return lenDiff;
        }
        for (int i = 0; i < l.size(); i++) {
            if (tolower(l[i]) != tolower(r[i]) ){
                return l[i]-r[i];
            }
        }
        return 0;
    }

    template<typename T>
    T **new_matrix(const size_t rows, const size_t columns) {
        if (rows == 0)
            throw std::invalid_argument("number of rows is 0");
        if (columns == 0)
            throw std::invalid_argument("number of columns is 0");
        T **ptr = nullptr;
        T *pool = nullptr;
        try {
            ptr = new T *[rows];  // allocate pointers (can throw here)
            pool = new T[rows * columns];  // allocate pool (can throw here)

            // now point the row pointers to the appropriate positions in
            // the memory pool
            for (unsigned i = 0; i < rows; ++i, pool += columns)
                ptr[i] = pool;

            // Done.
            return ptr;
        }
        catch (std::bad_alloc &ex) {
            delete[] ptr; // either this is nullptr or it was allocated
            throw ex;  // memory allocation error
        }
    }

    template<typename T>
    T **new_matrix(const T *pool, const size_t rows, const size_t cols) {
        if (rows == 0)
            throw std::invalid_argument("number of rows is 0");
        if (cols == 0)
            throw std::invalid_argument("number of columns is 0");
        T **ptr = nullptr;
        try {
            ptr = new T *[rows];  // allocate pointers (can throw here)
            // now point the row pointers to the appropriate positions in
            // the memory pool
            for (unsigned i = 0; i < rows; ++i, pool += cols)
                ptr[i] = (T *) pool;
            return ptr;
        }
        catch (std::bad_alloc &ex) {
            delete[] ptr; // either this is nullptr or it was allocated
            throw ex;  // memory allocation error
        }
    }

    // k prefix conforms to Google C++ style guide
    //kCopyColumnWise serves MatLab matrices
    enum class Transfer {
        kMove, kCopy, kCopyColumnWise
    };


    template<typename T>
    void copy( T **to, const T **from, const size_t rows, const size_t columns) {
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < columns; col++) {
                to[row][col] = from[row][col];
            }
        }
    }
    // by default copy MatLab matrices (column wise order)
    template<typename T>
    T **new_matrix(const T *pool, const size_t rows, const size_t columns, const Transfer transfer) {
        T **ptr = nullptr;
        if (transfer == Transfer::kMove) {
            ptr = new_matrix<T>(pool, rows, columns);
            for (unsigned i = 0; i < rows; ++i, pool += columns)
                ptr[i] = (T *) pool;
        } else {
            ptr = new_matrix<T>(rows, columns);
            if (transfer == Transfer::kCopyColumnWise) {
                for (int col = 0; col < columns; col++) {
                    for (int row = 0; row < rows; row++) {
                        ptr[row][col] = pool[row + rows * col];
                    }
                }
            } else { // Transfer::kCopy
                unsigned int N = rows * columns;
                T *const to = ptr[0], *end = to + N;
                for (int i = 0; i < N; i++) {
                    to[i] = pool[i];
                }
            }
        }
        return (T **) ptr;
    }


    template<typename T>
    T **new_matrix(const size_t rows, const size_t columns, const T default_value) {
        T **ptr = new_matrix<T>(rows, columns);
        std::fill_n(ptr[0], rows * columns, default_value);
        return ptr;
    }

    template<typename T>
    void delete_matrix(T **arr) {
        delete[] arr[0];  // remove the pool
        delete[] arr;     // remove the pointers
    }

    template<typename T>
    T ***new_matrices(const size_t n_matrices, const size_t rows, const size_t columns, const T *default_values) {
        T ***ptr = new T **[n_matrices];
        for (int i = 0; i < n_matrices; i++) {
            ptr[i] = new_matrix<T>(rows, columns, default_values[i]);
        }
        return ptr;
    }

    template<typename T>
    T ***new_matrices(const size_t n_matrices, const size_t rows, const size_t columns) {
        T ***ptr = new T **[n_matrices];
        for (int i = 0; i < n_matrices; i++) {
            ptr[i] = new_matrix<T>(rows, columns);
        }
        return ptr;
    }

    template<typename T>
    T ***new_matrices(const T ***prior, const size_t n_matrices, const size_t rows, const size_t columns) {
        T ***ptr = new T **[n_matrices];
        for (int i = 0; i < n_matrices; i++) {
            ptr[i] = new_matrix<T>(prior[i][0], rows, columns, Transfer::kCopy);
        }
        return ptr;
    }

    template<typename T>
    void
    delete_matrices(T ***matrices_ptr, const size_t n_matrices) {
        if (matrices_ptr != nullptr) {
            for (int i = 0; i < n_matrices; i++) {
                delete_matrix(matrices_ptr[i]);
            }
        }
    }

    inline std::string home_file(const char *path) {
        std::string home("/Users/");
        home.append(getenv("USER"));
        home.append(path);
        return home;
    }

    inline std::string documents_file(const char *path) {
        std::string home("/Users/");
        home.append(getenv("USER"));
        home.append("/Documents/");
        home.append(path);
        return home;
    }
    inline std::string run_umap_examples_file(const char *path) {
        std::string home("/Users/");
        home.append(getenv("USER"));
        home.append("/Documents/run_umap/examples/");
        home.append(path);
        return home;
    }

    inline std::string trim(std::string s) {
        std::regex e("^\\s+|\\s+$");   // remove leading and trailing spaces
        return regex_replace(s, e, "");
    }

    template<typename T>
    std::ostream &print(const T **const data, const int rows, const int columns, std::ostream &o = std::cout) {
        for (auto row = 0; row < rows; row++) {
            o << "#" << (row + 1) << ": ";
            for (auto col = 0; col < columns; col++) {
                if (col < columns - 1) {
                    o << data[row][col] << ",";
                } else {
                    o << data[row][col] << std::endl;
                }
            }
        }
        return o;
    }

    template<typename T>
    void print(const T *const data, const int size, std::ostream &o = std::cout) {
        for (int i = 0; i < size; i++) {
            if (i < size - 1) {
                std::cout << " " << data[i] << ',';
            } else {
                std::cout << data[i] << std::endl;
            }
        }
    }

    inline void Debug(double **nums, const int row, const int size) {
        std::cout << "#" << row << " *=" << suh::mean(nums[row], size) << ":";
        suh::print(nums[row], size);
    }

    template<typename T>
    void Debug(T *nums, const int size) {
        std::cout << " *=" << suh::mean(nums, size) << ":";
        suh::print(nums, size);
    }

    //Utility for bounds checking of Matrix
    template<class T>
    class Row {
    public:
        T *const row;
        const int columns;

        inline Row(T *const row, const int columns) : row(row), columns(columns) {
        }

        inline T operator[](int column) const {
            assert(column >= 0 && column < columns && row != nullptr);
            return row[column];
        }
    };

// Matrix is a class that contains 2D arrays allocated externally typically of builtin type
// It serves these externally builtin 2D arrays by providing a way to destruct them when the container
// loses scope.  Matrix can be a stack variable and simplify destruction.  It can also be built on the stack
// and then be switched to the heap with shared_ptr at which point the original allocation is not lost when
// the stack variable leaves scope of its closure
    template<class T>
    class Matrix {
    public:

        //Utility for bounds checking indexing of matrix otherwise use this->matrix[row][col]
        inline suh::Row<T> operator[](int row) const {
            assert(row >= 0 && row < rows_ && matrix_ != nullptr);
            return Row<T>(matrix_[row], columns_);
        }

        T **matrix_;
        size_t columns_;
        size_t rows_;

        Matrix() : matrix_(nullptr), columns_(0), rows_(0) {
        }

        Matrix(const T *pool, const size_t rows, const size_t  columns, const Transfer transfer = Transfer::kCopy)
                : matrix_(nullptr), columns_(columns), rows_(rows) {
            matrix_ = new_matrix(pool, rows, columns, transfer);
        }

        Matrix(T **const matrix, const size_t rows, const size_t columns) : matrix_(matrix), rows_(rows), columns_(columns) {
        }

        Matrix(const size_t size) : matrix_(suh::new_matrix<T>(size, size)), rows_(size), columns_(size) {
        }

        Matrix(const size_t rows, const size_t columns) : matrix_(suh::new_matrix<T>(rows, columns)), rows_(rows), columns_(columns) {
        }

        Matrix<T> &operator=(const Matrix<T> &other) = delete;

        void block_matrix_deallocation() {
            shared = this;
        }

        void clear(){
            if (delete_each_row){
                for (int row = 0; row < rows_; row++) {
                    std::memset(matrix_[row], 0, columns_*sizeof(T));
                }
            } else{
                //everybody out of the pool!
                std::memset(matrix_[0], 0, rows_*columns_*sizeof(T));
            }
        }
        virtual ~Matrix() {
            if (shared == nullptr) {
                if (matrix_ != nullptr) {

                    if (delete_each_row) {
                        for (int row = 0; row < rows_; row++) {
                            delete[] matrix_[row];
                        }
                        if (debug_ctor_dtor) std::cout << " destructing matrix & rows ";
                    } else { // delete the pool of rows X columns of T
                        delete[] matrix_[0];
                        if (debug_ctor_dtor) std::cout << " destructing matrix & pool ";
                    }
                    delete matrix_;
                    if (debug_ctor_dtor) std::cout << std::string(*this) << std::endl;
                } else {
                    if (debug_ctor_dtor)
                        std::cout << "destructing matrix with nullptr " << rows_ << " X " << columns_ << std::endl;
                }
            }
        }

        virtual operator std::string() const {
            std::strstream s;
            s << "matrix: " << rows_ << " X " << columns_;
            return s.str();
        }

        void copy_column_wise(T *to) const {
            for (int row = 0; row < rows_; row++) {
                for (int col = 0; col < columns_; col++) {
                    to[row + rows_ * col] = matrix_[row][col];
                }
            }
        }

        void copy_column_wise(T *to, const T add) const {
            for (int row = 0; row < rows_; row++) {
                for (int col = 0; col < columns_; col++) {
                    to[row + rows_ * col] = matrix_[row][col] + add;
                }
            }
        }

    private:
        Matrix<T> *shared = nullptr;
        std::shared_ptr<Matrix<T>> sh;

        Matrix(const Matrix<T> &other)
                : matrix_(other.matrix_), rows_(other.rows_),
                  columns_(other.columns_), delete_each_row(other.delete_each_row) {
        }

    protected:
        bool delete_each_row = false;

        std::vector<std::string> columnNames;

    public:

        std::shared_ptr<Matrix<T>> shared_ptr() {
            if (shared == nullptr) {
                shared = new Matrix<T>(*this);
                sh = std::shared_ptr<Matrix<T>>(shared);
            }
            return std::shared_ptr<Matrix<T>>(sh);
        }

        inline T *vector(int row = 0) const {
            if (rows_ != 1) {
                std::cerr << "vector() called for matrix with "
                          << rows_ << " rows???" << std::endl;
            }
            return matrix_[0];
        }

        std::ostream &print(std::ostream &o = std::cout) const {
            for (int col = 0; col < this->columns_; col++) {
                if (col < this->columns_ - 1) {
                    o << this->columnNames[col] << ',';
                } else {
                    o << this->columnNames[col] << std::endl;
                }
            }
            suh::print((const double **) this->matrix_, this->rows_, this->columns_, o);
            return o;
        }
        std::ostream &print(const int row, std::ostream &o = std::cout) const {
            suh::print(this->matrix_[row], this->columns_, o);
            return o;
        }
    };

    template<typename T>
    int count_unequal(const std::vector<T> v1, const std::vector<T> v2,
                      const double tolerance = .0002) {
        int dif = abs((long)(v1.size() - v2.size()));
        for (int row = 0; row < v1.size(); row++)
            if (row < v2.size()) {
                T value1 = v1[row], value2 = v2[row];
                T d = abs(value1 - value2);
                if (d > tolerance)
                    dif++;
            }
        return dif;
    }

    template<typename T>
    int count_unequal(const std::vector<std::vector<T>> v1, const std::vector<std::vector<T>> v2,
                      const double tolerance = .0002) {
        int dif = abs((long)(v1.size() - v2.size()));
        for (int row = 0; row < v1.size(); row++)
            if (row < v2.size()) {
                dif += abs((long)(v1[row].size() - v2[row].size()));
                for (int col = 0; col < v1[0].size(); col++) {
                    if (col<v2[row].size()) {
                        T value1 = v1[row][col], value2 = v2[row][col];
                        T d = abs(value1 - value2);
                        if (d > tolerance)
                            dif++;
                    }
                }
            }
        return dif;
    }

    template<typename T>
    int count_unequal(const T **thisPtr, const T **thatPtr, const size_t rows, const size_t columns) {
        int dif = 0;
        for (int row = 0; row < rows; row++)
            for (int col = 0; col < columns; col++)
                if (thisPtr[row][col] != thatPtr[row][col])
                    dif++;
        return dif;
    }

    template<typename T>
    int count_unequal(const T *thisPtr, const T *thatPtr, const size_t rows, const size_t columns) {
        int dif = 0;
        const size_t N = rows * columns;
        for (int i = 0; i < N; i++)
            if (thisPtr[i] != thatPtr[i])
                dif++;
        return dif;
    }

    template<typename T>
    int count_unequal(const Matrix<T> &this_matrix, const Matrix<T> &that_matrix) {
        const size_t rows = this_matrix.rows_ > that_matrix.rows_ ? that_matrix.rows_ : this_matrix.rows_;
        const size_t columns = this_matrix.columns_ > that_matrix.columns_ ? that_matrix.columns_ : this_matrix.columns_;
        return count_unequal((const T **) this_matrix.matrix_, (const T **) that_matrix.matrix_, rows, columns);
    }


    template<typename T>
    std::vector<int> get_unequal(const Matrix<T> &this_matrix, const Matrix<T> &that_matrix, const int row,
                                 const double tolerance = .0002) {
        const size_t columns = this_matrix.columns_ > that_matrix.columns_ ? that_matrix.columns_ : this_matrix.columns_;
        std::vector<int> notEqual;
        T *thisRow = this_matrix.matrix_[row], *thatRow = that_matrix.matrix_[row];
        for (int col = 0; col < columns; col++) {
            if (thisRow[col] != thatRow[col]) {
                if (tolerance == 0) {
                    notEqual.push_back(col);
                } else {
                    if (abs(thisRow[col] - thatRow[col]) > tolerance)
                        notEqual.push_back(col);
                }
            }
        }
        return notEqual;
    }

    template<typename T>
    int print_inequalities(const Matrix<T> &this_matrix, const Matrix<T> &that_matrix, const int row,
                           const bool print_if_equal = false, const double tolerance = .0002) {
        std::vector<int> bad = suh::get_unequal(this_matrix, that_matrix, row, tolerance);
        if (print_if_equal && !bad.empty()) {
            std::cout << "Row #" << row << " has " << bad.size() << " inequal values @ columns:" << std::endl;
            std::cout << "    ->";
            for (int i = 0; i < bad.size(); i++) {
                std::cout << " " << bad[i];
            }
            std::cout << std::endl;
        }
        return bad.size();
    }

    template<typename T>
    int
    print_inequalities(const Matrix<T> &this_matrix, const Matrix<T> &that_matrix, const bool print_if_equal = false,
                       const double tolerance = .0002) {
        const size_t rows = this_matrix.rows_ > that_matrix.rows_ ? that_matrix.rows_ : this_matrix.rows_;
        int badCnt = 0;
        for (int row = 0; row < rows; row++) {
            badCnt += print_inequalities(this_matrix, that_matrix, row);
        }
        if (badCnt > 0 || print_if_equal) {
            if (badCnt == 0) {
                std::cout << "PERFECT:  ";
            }
            std::cout << bank_num(100.0 - (100.0 * (double) badCnt / (double) (rows * this_matrix.columns_)))
                      << "% similar (" << badCnt << "/" << (rows * this_matrix.columns_);

            int inexact;
            if (tolerance > 0)
                inexact = count_unequal(this_matrix, that_matrix);
            else
                inexact = badCnt;
            if (badCnt != inexact) {
                std::cout << " more than .0002 different, " << bank_num(inexact) << " unequal";
            } else {
                std::cout << " are unequal";
            }
            std::cout << ")" << std::endl;
        }
        return badCnt;
    }

    using MatrixPtr = std::shared_ptr<suh::Matrix<double>>;
    using MatrixIntPtr = std::shared_ptr<suh::Matrix<int>>;

    template<class T>
    class CsvMatrix : public Matrix<T> {
    public:
        inline suh::Row<T> operator[](int row) const {
            return Row<T>(this->matrix_[row], this->columns_);
        }

        CsvMatrix(const std::string file_path, const std::string suffix = "")
                : Matrix<T>() {
            if (suh::debug_file_io) {
                std::cout << "Reading " << file_path << "...";
                std::cout.flush();
            }
            std::vector<T *> matrixVector;
            std::string s(file_path);
            s.append(suffix);
            std::ifstream csv_file(s);
            this->delete_each_row = true;
            if (!csv_file.good()) {
                std::cerr << "The file " << s << " does NOT exist !!!" << std::endl;
            }
            if (csv_file.is_open()) {
                std::string line, colName;
                while (std::getline(csv_file, line) && trim(line).empty());
                std::stringstream lineOfNames(line);
                while (std::getline(lineOfNames, colName, ','))
                    this->columnNames.push_back(colName);
                this->columns_ = this->columnNames.size();
                while (std::getline(csv_file, line) && trim(line).empty());
                while (!line.empty()) {
                    std::stringstream rowOfValues(line);
                    int colIdx = 0;
                    T *const row = new T[this->columns_];
                    std::string token;
                    while (std::getline(rowOfValues, token, ',')) {
                        std::stringstream tokenLine(token);
                        tokenLine >> row[colIdx];
                        if (colIdx++ == this->columns_)
                            break;
                        if (rowOfValues.peek() == ',')
                            rowOfValues.ignore();
                    }
                    matrixVector.push_back(row);
                    if (!std::getline(csv_file, line)) {
                        break;
                    }
                }
            }
            csv_file.close();
            if (this->columns_ > 1) {
                this->rows_ = matrixVector.size();
                this->matrix_ = new T *[this->rows_];
                for (int row = 0; row < this->rows_; ++row) {
                    this->matrix_[row] = matrixVector[row];
                }
            } else {
                this->rows_ = 1;
                this->columns_ = matrixVector.size();
                T *row = new T[this->columns_];
                this->matrix_ = new T *[this->rows_];
                this->matrix_[0] = row;
                for (int col = 0; col < this->columns_; ++col) {
                    T *p = matrixVector[col];
                    row[col] = p[0];
                    delete[] p;
                }
            }
            if (suh::debug_file_io)
                std::cout << this->rows_ << " X " << this->columns_ << " values read." <<std::endl;

        }

        T *operator[](const int row) {
            return this->matrix_[row];
        }

    };
    template <typename T> struct  Deserialize{
        std::string name;
        int columns_=0;
        std::vector<std::vector<T>>v2;
        Deserialize( const std::string &in){
            std::stringstream rowOfValues(in);
            std::vector<T>v1;
            std::getline(rowOfValues, name, ',');
            if (!name.empty()){
                std::string token;
                std::getline(rowOfValues, token, ',');
                std::stringstream c(token);
                c >> columns_;
                if (columns_>0) {
                    size_t colIdx = 0;
                    T value;
                    while (std::getline(rowOfValues, token, ',')) {
                        std::stringstream tokenLine(token);
                        tokenLine >> value;
                        v1.push_back(value);
                        if (colIdx == this->columns_ - 1) {
                            v2.push_back(v1);
                            v1.clear();
                            colIdx = 0;
                        } else {
                            colIdx++;
                        }
                        if (rowOfValues.peek() == ',')
                            rowOfValues.ignore();
                    }
                }
            }
        }
    };


    struct KnnSelfOtherUseCaseCsv {
        suh::CsvMatrix<int> indices, indptr;
        suh::CsvMatrix<double> self, other, correctKnnIndices, correctKnnDists;

        inline KnnSelfOtherUseCaseCsv(const std::string trainingSetPath, const std::string testSetPath) :
                other(trainingSetPath, ".csv"),
                self(testSetPath, ".csv"),
                indptr(trainingSetPath, ".indptr.csv"),
                indices(trainingSetPath, ".indices.csv"),
                correctKnnDists(testSetPath, ".knnDists.csv"),
                correctKnnIndices(testSetPath, ".knnIndices.csv") {
            assert(other.rows_ + 1 == indptr.columns_);
        }

    };



}
#endif //SUH_H