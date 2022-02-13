#include "lss_flat_matrix.hpp"

#include <algorithm>
#include <tuple>

namespace lss_containers
{

flat_matrix::flat_matrix()
{
}

flat_matrix::~flat_matrix()
{
}

flat_matrix::flat_matrix(std::size_t nrows, std::size_t ncols) : nrows_{nrows}, ncols_{ncols}
{
    column_cnt_.resize(nrows);
    diagonal_.resize(nrows);
}

flat_matrix::flat_matrix(flat_matrix const &copy)
    : ncols_{copy.ncols_}, nrows_{copy.nrows_}, container_{copy.container_},
      column_cnt_{copy.column_cnt_}, diagonal_{copy.diagonal_}
{
}

flat_matrix::flat_matrix(flat_matrix &&other) noexcept
    : ncols_{std::move(other.ncols_)}, nrows_{std::move(other.nrows_)}, container_{std::move(other.container_)},
      column_cnt_{std::move(other.column_cnt_)}, diagonal_{std::move(other.diagonal_)}
{
}

flat_matrix &flat_matrix::operator=(flat_matrix const &copy)
{
    if (this != &copy)
    {
        ncols_ = copy.ncols_;
        nrows_ = copy.nrows_;
        container_ = copy.container_;
        column_cnt_ = copy.column_cnt_;
        diagonal_ = copy.diagonal_;
    }
    return *this;
}

flat_matrix &flat_matrix::operator=(flat_matrix &&other) noexcept
{
    if (this != &other)
    {
        ncols_ = std::move(other.ncols_);
        nrows_ = std::move(other.nrows_);
        container_ = std::move(other.container_);
        column_cnt_ = std::move(other.column_cnt_);
        diagonal_ = std::move(other.diagonal_);
    }
    return *this;
}

std::size_t flat_matrix::rows() const
{
    return nrows_;
}
std::size_t flat_matrix::columns() const
{
    return ncols_;
}
std::size_t flat_matrix::size() const
{
    return container_.size();
}
void flat_matrix::clear()
{
    container_.clear();
}
double flat_matrix::diagonal_at_row(std::size_t row_idx) const
{
    return diagonal_[row_idx];
}
std::size_t flat_matrix::non_zero_column_size(std::size_t row_idx) const
{
    return column_cnt_[row_idx];
}

void flat_matrix::emplace_back(std::size_t row_idx, std::size_t col_idx, double value)
{
    LSS_ASSERT(row_idx < nrows_, " rowIdx is outside <0," << nrows_ << ")");
    LSS_ASSERT(col_idx < ncols_, " colIdx is outside <0," << ncols_ << ")");
    container_.emplace_back(std::make_tuple(row_idx, col_idx, value));
    column_cnt_[row_idx]++;
    if (row_idx == col_idx)
        diagonal_[row_idx] = value;
}

void flat_matrix::emplace_back(std::tuple<std::size_t, std::size_t, double> tuple)
{
    LSS_ASSERT(std::get<0>(tuple) < nrows_, " rowIdx is outside <0," << nrows_ << ")");
    LSS_ASSERT(std::get<1>(tuple) < ncols_, " colIdx is outside <0," << ncols_ << ")");
    container_.emplace_back(std::move(tuple));
    column_cnt_[std::get<0>(tuple)]++;
    if (std::get<1>(tuple) == std::get<0>(tuple))
        diagonal_[std::get<1>(tuple)] = std::get<2>(tuple);
}

std::tuple<std::size_t, std::size_t, double> const &flat_matrix::at(std::size_t idx) const
{
    return container_.at(idx);
}

void flat_matrix::sort(flat_matrix_sort_enum sort)
{
    if (sort == flat_matrix_sort_enum::RowMajor)
    {
        std::sort(container_.begin(), container_.end(),
                  [this](std::tuple<std::size_t, std::size_t, double> const &lhs,
                         std::tuple<std::size_t, std::size_t, double> const &rhs) {
                      std::size_t const flat_idx_lhs = std::get<1>(lhs) + nrows_ * std::get<0>(lhs);
                      std::size_t const flat_idx_rhs = std::get<1>(rhs) + nrows_ * std::get<0>(rhs);
                      return (flat_idx_lhs < flat_idx_rhs);
                  });
    }
    else
    {
        std::sort(container_.begin(), container_.end(),
                  [this](std::tuple<std::size_t, std::size_t, double> const &lhs,
                         std::tuple<std::size_t, std::size_t, double> const &rhs) {
                      std::size_t const flat_idx_lhs = std::get<0>(lhs) + ncols_ * std::get<1>(lhs);
                      std::size_t const flat_idx_rhs = std::get<0>(rhs) + ncols_ * std::get<1>(rhs);
                      return (flat_idx_lhs < flat_idx_rhs);
                  });
    }
}

} // namespace lss_containers
