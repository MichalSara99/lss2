#include "lss_matrix_2d.hpp"
#include <typeinfo>
#include <vector>

namespace lss_containers
{

matrix_2d::matrix_2d(std::size_t rows, std::size_t columns, double value)
    : rows_{rows}, cols_{columns}, data_(value, rows * columns)
{
}

matrix_2d::matrix_2d(std::size_t rows, std::size_t columns) : matrix_2d(rows, columns, 0.0)
{
}

matrix_2d::~matrix_2d()
{
}

void matrix_2d::from_data(std::valarray<double> &&data)
{
    LSS_ASSERT(data.size() == total_size(), "matrix_2d: Passed data must be of the same size as total size of "
                                            "matrix.");
    data_ = std::move(data);
}

std::valarray<double> const &matrix_2d::data() const
{
    return data_;
}

std::size_t const &matrix_2d::rows() const
{
    return rows_;
}

std::size_t const &matrix_2d::columns() const
{
    return cols_;
}

std::size_t matrix_2d::total_size() const
{
    return (rows_ * cols_);
}

double &matrix_2d::operator()(std::size_t row_idx, std::size_t col_idx)
{
    LSS_ASSERT(row_idx >= 0 && row_idx < rows_, "matrix_2d: row_idx out of range.");
    LSS_ASSERT(col_idx >= 0 && col_idx < cols_, "matrix_2d: col_idx out of range.");
    return data_[col_idx + cols_ * row_idx];
}

double matrix_2d::operator()(std::size_t row_idx, std::size_t col_idx) const
{
    LSS_ASSERT(row_idx >= 0 && row_idx < rows_, "matrix_2d: row_idx out of range.");
    LSS_ASSERT(col_idx >= 0 && col_idx < cols_, "matrix_2d: col_idx out of range.");
    return data_[col_idx + cols_ * row_idx];
}

void matrix_2d::operator()(std::size_t row_idx, std::size_t col_idx, double value)
{
    LSS_ASSERT(row_idx >= 0 && row_idx < rows_, "matrix_2d: row_idx out of range.");
    LSS_ASSERT(col_idx >= 0 && col_idx < cols_, "matrix_2d: col_idx out of range.");
    data_[col_idx + cols_ * row_idx] = value;
}

std::slice_array<double> matrix_2d::row(std::size_t idx)
{
    LSS_ASSERT(idx >= 0 && idx < rows_, "matrix_2d: idx out of range.");
    return data_[std::slice(idx * cols_, cols_, 1)];
}

void matrix_2d::row(std::size_t idx, std::valarray<double> &&vals)
{
    LSS_ASSERT(idx >= 0 && idx < rows_, "matrix_2d: idx out of range.");
    LSS_ASSERT(cols_ == vals.size(), "matrix_2d: vals must have same column size");
    data_[std::slice(idx * cols_, cols_, 1)] = std::move(vals);
}

std::slice_array<double> matrix_2d::column(std::size_t idx)
{
    LSS_ASSERT(idx >= 0 && idx < cols_, "matrix_2d: idx out of range.");
    return data_[std::slice(idx, rows_, cols_)];
}

void matrix_2d::column(std::size_t idx, std::valarray<double> &&vals)
{
    LSS_ASSERT(idx >= 0 && idx < cols_, "matrix_2d: idx out of range.");
    LSS_ASSERT(rows_ == vals.size(), "matrix_2d: vals must have same row size");
    data_[std::slice(idx, rows_, cols_)] = std::move(vals);
}

} // namespace lss_containers
