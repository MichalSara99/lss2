/**

    @file      lss_flat_matrix.hpp
    @brief     Represents flat matrix mainly for CUDA solvers
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once
#if !defined(_LSS_FLAT_MATRIX_HPP_)
#define _LSS_FLAT_MATRIX_HPP_

#include "../common/lss_enumerations.hpp"
#include "../common/lss_macros.hpp"
#include "../common/lss_utility.hpp"
#include <vector>

namespace lss_containers
{

using lss_enumerations::flat_matrix_sort_enum;

/**
    @struct flat_matrix
    @brief  Flat matrix object
**/
struct flat_matrix
{
  private:
    std::vector<std::tuple<std::size_t, std::size_t, double>> container_;
    std::vector<std::size_t> column_cnt_;
    std::vector<double> diagonal_;
    std::size_t ncols_, nrows_;

    explicit flat_matrix();

  public:
    /**
        @brief flat_matrix object constructor
        @param nrows number of rows
        @param ncols number of columns
    **/
    explicit flat_matrix(std::size_t nrows, std::size_t ncols);

    virtual ~flat_matrix();

    flat_matrix(flat_matrix const &copy);
    flat_matrix(flat_matrix &&other) noexcept;

    flat_matrix &operator=(flat_matrix const &copy);
    flat_matrix &operator=(flat_matrix &&other) noexcept;

    /**
        @brief  Number of rows.
        @retval number of rows
    **/
    std::size_t rows() const;

    /**
        @brief  Number of columns.
        @retval number of columns
    **/
    std::size_t columns() const;

    /**
        @brief  Total size.
        @retval total size
    **/
    std::size_t size() const;

    /**
        @brief Clear the matrix
    **/
    void clear();

    /**
        @brief Return diagonal element at row index.
        @param  row_idx
        @retval
    **/
    double diagonal_at_row(std::size_t row_idx) const;

    /**
        @brief Size of the non-zero container at row row_idx
        @param  row_idx
        @retval
    **/
    std::size_t non_zero_column_size(std::size_t row_idx) const;

    /**
        @brief Place value at position (row_idx,col_idx).
        @param row_idx
        @param col_idx
        @param value
    **/
    void emplace_back(std::size_t row_idx, std::size_t col_idx, double value);

    /**
        @brief Append a tuple.
        @param tuple
    **/
    void emplace_back(std::tuple<std::size_t, std::size_t, double> tuple);

    /**
        @brief Sort flat matrix.
        @param sort
    **/
    void sort(flat_matrix_sort_enum sort);

    /**
        @brief  Get a tuple from flat matrix at idx
        @param  idx
        @retval
    **/
    std::tuple<std::size_t, std::size_t, double> const &at(std::size_t idx) const;
};

} // namespace lss_containers

#endif ///_LSS_FLAT_MATRIX_HPP_
