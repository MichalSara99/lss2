#include "lss_matrix_3d.hpp"

namespace lss_containers
{

matrix_3d::matrix_3d(std::size_t rows, std::size_t columns, std::size_t layers, double value)
    : rows_{rows}, cols_{columns}, lays_{layers}, data_(value, rows * columns * layers)
{
}

matrix_3d::matrix_3d(std::size_t rows, std::size_t columns, std::size_t layers) : matrix_3d(rows, columns, layers, 0.0)
{
}

matrix_3d::~matrix_3d()
{
}

void matrix_3d::from_data(std::valarray<double> &&data)
{
    LSS_ASSERT(data.size() == (total_size()), "matrix_3d: Passed data must be of the same size as total size of "
                                              "matrix.");
    data_ = std::move(data);
}

std::valarray<double> const &matrix_3d::data() const
{
    return data_;
}

std::size_t const &matrix_3d::rows() const
{
    return rows_;
}

std::size_t const &matrix_3d::columns() const
{
    return cols_;
}

std::size_t const &matrix_3d::layers() const
{
    return lays_;
}

std::size_t matrix_3d::total_size() const
{
    return (rows_ * cols_ * lays_);
}

double &matrix_3d::operator()(std::size_t row_idx, std::size_t col_idx, std::size_t lay_idx)
{
    LSS_ASSERT(row_idx >= 0 && row_idx < rows_, "matrix_3d: row_idx out of range.");
    LSS_ASSERT(col_idx >= 0 && col_idx < cols_, "matrix_3d: col_idx out of range.");
    LSS_ASSERT(lay_idx >= 0 && lay_idx < lays_, "matrix_3d: lay_idx out of range.");
    return data_[col_idx + cols_ * row_idx + rows_ * cols_ * lay_idx];
}

double matrix_3d::operator()(std::size_t row_idx, std::size_t col_idx, std::size_t lay_idx) const
{
    LSS_ASSERT(row_idx >= 0 && row_idx < rows_, "matrix_3d: row_idx out of range.");
    LSS_ASSERT(col_idx >= 0 && col_idx < cols_, "matrix_3d: col_idx out of range.");
    LSS_ASSERT(lay_idx >= 0 && lay_idx < lays_, "matrix_3d: lay_idx out of range.");
    return data_[col_idx + cols_ * row_idx + rows_ * cols_ * lay_idx];
}

void matrix_3d::operator()(std::size_t row_idx, std::size_t col_idx, std::size_t lay_idx, double value)
{
    LSS_ASSERT(row_idx >= 0 && row_idx < rows_, "matrix_3d: row_idx out of range.");
    LSS_ASSERT(col_idx >= 0 && col_idx < cols_, "matrix_3d: col_idx out of range.");
    LSS_ASSERT(lay_idx >= 0 && lay_idx < lays_, "matrix_3d: lay_idx out of range.");
    data_[col_idx + cols_ * row_idx + rows_ * cols_ * lay_idx] = value;
}

std::slice_array<double> matrix_3d::row(std::size_t row_idx, std::size_t lay_idx)
{
    LSS_ASSERT(row_idx >= 0 && row_idx < rows_, "matrix_3d: row_idx out of range.");
    LSS_ASSERT(lay_idx >= 0 && lay_idx < lays_, "matrix_3d: lay_idx out of range.");
    return data_[std::slice(cols_ * row_idx + rows_ * cols_ * lay_idx, cols_, 1)];
}

void matrix_3d::row(std::size_t row_idx, std::size_t lay_idx, std::valarray<double> &&vals)
{
    LSS_ASSERT(row_idx >= 0 && row_idx < rows_, "matrix_3d: row_idx out of range.");
    LSS_ASSERT(lay_idx >= 0 && lay_idx < lays_, "matrix_3d: lay_idx out of range.");
    LSS_ASSERT(cols_ == vals.size(), "matrix_3d: vals must have same row size");
    data_[std::slice(cols_ * row_idx + rows_ * cols_ * lay_idx, cols_, 1)] = std::move(vals);
}

std::slice_array<double> matrix_3d::column(std::size_t col_idx, std::size_t lay_idx)
{
    LSS_ASSERT(col_idx >= 0 && col_idx < cols_, "matrix_3d: col_idx out of range.");
    LSS_ASSERT(lay_idx >= 0 && lay_idx < lays_, "matrix_3d: lay_idx out of range.");
    return data_[std::slice(col_idx + lay_idx * cols_ * rows_, rows_, cols_)];
}

void matrix_3d::column(std::size_t col_idx, std::size_t lay_idx, std::valarray<double> &&vals)
{
    LSS_ASSERT(col_idx >= 0 && col_idx < cols_, "matrix_3d: col_idx out of range.");
    LSS_ASSERT(lay_idx >= 0 && lay_idx < lays_, "matrix_3d: lay_idx out of range.");
    LSS_ASSERT(rows_ == vals.size(), "matrix_3d: vals must have same column size");
    data_[std::slice(col_idx + lay_idx * cols_ * rows_, rows_, cols_)] = std::move(vals);
}

std::slice_array<double> matrix_3d::layer(std::size_t row_idx, std::size_t col_idx)
{
    LSS_ASSERT(row_idx >= 0 && row_idx < rows_, "matrix_3d: row_idx out of range.");
    LSS_ASSERT(col_idx >= 0 && col_idx < cols_, "matrix_3d: col_idx out of range.");
    return data_[std::slice(col_idx + row_idx * cols_, lays_, cols_ * rows_)];
}

void matrix_3d::layer(std::size_t row_idx, std::size_t col_idx, std::valarray<double> &&vals)
{
    LSS_ASSERT(row_idx >= 0 && row_idx < rows_, "matrix_3d: row_idx out of range.");
    LSS_ASSERT(col_idx >= 0 && col_idx < cols_, "matrix_3d: col_idx out of range.");
    LSS_ASSERT(lays_ == vals.size(), "matrix_3d: vals must have same layer size");
    data_[std::slice(col_idx + row_idx * cols_, lays_, cols_ * rows_)] = std::move(vals);
}

} // namespace lss_containers
