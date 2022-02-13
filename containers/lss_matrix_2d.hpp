/**

    @file      lss_matrix_2d.hpp
    @brief     Represents 2D slicable matrix
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once
#if !defined(_LSS_MATRIX_2D_HPP_)
#define _LSS_MATRIX_2D_HPP_

#include "../common/lss_enumerations.hpp"
#include "../common/lss_macros.hpp"
#include "../common/lss_utility.hpp"
#include <typeinfo>
#include <valarray>

namespace lss_containers
{

using lss_utility::container_t;
using lss_utility::sptr_t;

/**

    @class   matrix_2d
    @brief
    @details ~

**/
class matrix_2d
{
  private:
    std::size_t rows_;
    std::size_t cols_;
    std::valarray<double> data_;

    explicit matrix_2d();

  public:
    explicit matrix_2d(std::size_t rows, std::size_t columns);
    explicit matrix_2d(std::size_t rows, std::size_t columns, double value);

    ~matrix_2d();

    matrix_2d(matrix_2d const &copy) = default;
    matrix_2d(matrix_2d &&other) noexcept = default;

    matrix_2d &operator=(matrix_2d const &copy) = default;
    matrix_2d &operator=(matrix_2d &&other) noexcept = default;

    /**
        @brief Populate matrix_2d from a passed data
        @param  data
    **/
    LSS_API void from_data(std::valarray<double> &&data);

    /**
        @brief Data from matrix_2d in row-wise fashion
        @retval data as flat vector row-wise
    **/
    LSS_API std::valarray<double> const &data() const;

    /**
        @brief Number of rows.
        @retval
    **/
    LSS_API std::size_t const &rows() const;

    /**
        @brief Number of columns.
        @retval
    **/
    LSS_API std::size_t const &columns() const;

    /**
        @brief Total size.
        @retval
    **/
    LSS_API std::size_t total_size() const;

    /**
        @brief  Assign a value to mat(i,j) = value
        @param  row_idx
        @param  col_idx
        @retval value from mat(i,j)
    **/
    LSS_API double &operator()(std::size_t row_idx, std::size_t col_idx);

    /**
        @brief  Assign a value to mat(i,j) = value
        @param  row_idx
        @param  col_idx
        @retval value from mat(i,j)
    **/
    LSS_API double operator()(std::size_t row_idx, std::size_t col_idx) const;

    /**
        @brief place value at position (row_idx,col_idx)
        @param row_idx
        @param col_idx
        @param value
    **/
    LSS_API void operator()(std::size_t row_idx, std::size_t col_idx, double value);

    /**
        @brief Return a row from matrix_2d at idx
        @param  idx
        @retval
    **/
    LSS_API std::slice_array<double> row(std::size_t idx);

    /**
        @brief Place a row into matrix_2d at idx
        @param  idx
        @param  vals
    **/
    LSS_API void row(std::size_t idx, std::valarray<double> &&vals);

    /**
        @brief Return a column from matrix_2d at idx
        @param  idx
        @retval
    **/
    LSS_API std::slice_array<double> column(std::size_t idx);

    /**
        @brief Place a column into matrix_2d at idx
        @param  idx
        @param  vals
    **/
    LSS_API void column(std::size_t idx, std::valarray<double> &&vals);
};

using matrix_2d_ptr = sptr_t<matrix_2d>;

} // namespace lss_containers

#endif ///_LSS_MATRIX_2D_HPP_
