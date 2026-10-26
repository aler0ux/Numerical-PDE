#ifndef PDE_COMPARE_CSRMATRIX_HPP
#define PDE_COMPARE_CSRMATRIX_HPP

// A minimal Compressed Sparse Row matrix, built once from a list of
// (row, col, value) triplets and never mutated afterwards.
//
// This is the "standard simple technique" version of CSR: no fancy
// insertion, no column reordering heuristics -- just sort the triplets,
// sum duplicates (so assembly can add() the same entry many times), and
// lay the result out as three flat arrays. That is all a matvec-based
// iterative solver needs.

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace cmp {

struct Triplet {
    int row;
    int col;
    double value;
};

class CSRMatrix {
public:
    CSRMatrix() = default;

    // Build from triplets. `n` is the number of rows == number of columns
    // (every matrix in this project is square). Duplicate (row,col) pairs
    // are summed, which is exactly what finite-element/finite-difference
    // assembly needs: several triangles or stencils may contribute to the
    // same matrix entry.
    static CSRMatrix fromTriplets(int n, std::vector<Triplet> triplets) {
        CSRMatrix m;
        m.n_ = n;
        std::sort(triplets.begin(), triplets.end(), [](const Triplet& a, const Triplet& b) {
            if (a.row != b.row) return a.row < b.row;
            return a.col < b.col;
        });

        m.rowPtr_.assign(static_cast<size_t>(n) + 1, 0);
        m.colIdx_.reserve(triplets.size());
        m.values_.reserve(triplets.size());

        size_t i = 0;
        for (int r = 0; r < n; ++r) {
            m.rowPtr_[static_cast<size_t>(r)] = static_cast<int>(m.colIdx_.size());
            while (i < triplets.size() && triplets[i].row == r) {
                const int col = triplets[i].col;
                double sum = 0.0;
                while (i < triplets.size() && triplets[i].row == r && triplets[i].col == col) {
                    sum += triplets[i].value;
                    ++i;
                }
                // Skip exact zeros; harmless to keep but noisy to inspect.
                m.colIdx_.push_back(col);
                m.values_.push_back(sum);
            }
        }
        m.rowPtr_[static_cast<size_t>(n)] = static_cast<int>(m.colIdx_.size());
        return m;
    }

    int rows() const { return n_; }
    int nnz() const { return static_cast<int>(values_.size()); }

    const std::vector<int>& rowPtr() const { return rowPtr_; }
    const std::vector<int>& colIdx() const { return colIdx_; }
    const std::vector<double>& values() const { return values_; }

    // y = A * x
    void matvec(const std::vector<double>& x, std::vector<double>& y) const {
        if (static_cast<int>(x.size()) != n_) throw std::invalid_argument("CSRMatrix::matvec: size mismatch");
        y.assign(static_cast<size_t>(n_), 0.0);
        for (int r = 0; r < n_; ++r) {
            double sum = 0.0;
            for (int k = rowPtr_[static_cast<size_t>(r)]; k < rowPtr_[static_cast<size_t>(r) + 1]; ++k) {
                sum += values_[static_cast<size_t>(k)] * x[static_cast<size_t>(colIdx_[static_cast<size_t>(k)])];
            }
            y[static_cast<size_t>(r)] = sum;
        }
    }

    // Diagonal entries, used by the Jacobi preconditioner. 0 if a row has
    // no diagonal entry (should not happen for the matrices this project
    // builds, but better a zero than undefined behaviour).
    std::vector<double> diagonal() const {
        std::vector<double> d(static_cast<size_t>(n_), 0.0);
        for (int r = 0; r < n_; ++r) {
            for (int k = rowPtr_[static_cast<size_t>(r)]; k < rowPtr_[static_cast<size_t>(r) + 1]; ++k) {
                if (colIdx_[static_cast<size_t>(k)] == r) {
                    d[static_cast<size_t>(r)] = values_[static_cast<size_t>(k)];
                    break;
                }
            }
        }
        return d;
    }

private:
    int n_ = 0;
    std::vector<int> rowPtr_;
    std::vector<int> colIdx_;
    std::vector<double> values_;
};

}  // namespace cmp

#endif  // PDE_COMPARE_CSRMATRIX_HPP
