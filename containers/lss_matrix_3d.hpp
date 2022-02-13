/**

    @file      lss_matrix_3d.hpp
    @brief     Represents 3D matrix
    @details   ~
    @author    Michal Sara
    @date      4.01.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_MATRIX_3D_HPP_)
#define _LSS_MATRIX_3D_HPP_

#include "../common/lss_enumerations.hpp"
#include "../common/lss_macros.hpp"
#include "../common/lss_utility.hpp"
#include <typeinfo>
#include <valarray>

namespace lss_containers
{

using lss_utility::sptr_t;

/**

    @class   matrix_3d
    @brief
    @details ~

**/
class matrix_3d
{
  private:
    std::size_t rows_, cols_, lays_;
    std::valarray<double> data_;

    explicit matrix_3d();

  public:
    explicit matrix_3d(std::size_t rows, std::size_t columns, std::size_t layers);
    explicit matrix_3d(std::size_t rows, std::size_t columns, std::size_t layers, double value);

    ~matrix_3d();

    matrix_3d(matrix_3d const &copy) = default;
    matrix_3d(matrix_3d &&other) noexcept = default;

    matrix_3d &operator=(matrix_3d const &copy) = default;
    matrix_3d &operator=(matrix_3d &&other) noexcept = default;

    /**
        @brief Populate matrix_3d from a passed data
        @param  data
    **/
    LSS_API void from_data(std::valarray<double> &&data);

    /**
        @brief Data from matrix_3d in row-wise fashion
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
        @brief Number of layers.
        @retval
    **/
    LSS_API std::size_t const &layers() const;

    /**
        @brief Total size.
        @retval
    **/
    LSS_API std::size_t total_size() const;

    /**
        @brief  Assign a value to mat(i,j,k) = value
        @param  row_idx
        @param  col_idx
        @param  lay_idx
        @retval value from mat(i,j,k)
    **/
    LSS_API double &operator()(std::size_t row_idx, std::size_t col_idx, std::size_t lay_idx);

    /**
        @brief  Assign a value to mat(i,j,k) = value
        @param  row_idx
        @param  col_idx
        @param  lay_idx
        @retval value from mat(i,j,k)
    **/
    LSS_API double operator()(std::size_t row_idx, std::size_t col_idx, std::size_t lay_idx) const;

    /**
        @brief place value at position (row_idx,col_idx,lay_idx)
        @param row_idx
        @param col_idx
        @param lay_idx
        @param value
    **/
    LSS_API void operator()(std::size_t row_idx, std::size_t col_idx, std::size_t lay_idx, double value);

    /**
        @brief Return a row from matrix_3d at (row_idx,lay_idx)
        @param  row_idx
        @param  lay_idx
        @retval
    **/
    LSS_API std::slice_array<double> row(std::size_t row_idx, std::size_t lay_idx);

    /**
        @brief Place a row into matrix_3d at (row_idx,lay_idx)
        @param  row_idx
        @param  lay_idx
        @param  vals
    **/
    LSS_API void row(std::size_t row_idx, std::size_t lay_idx, std::valarray<double> &&vals);

    /**
        @brief Return a column from matrix_3d at (col_idx,lay_idx)
        @param  col_idx
        @param  lay_idx
        @retval
    **/
    LSS_API std::slice_array<double> column(std::size_t col_idx, std::size_t lay_idx);

    /**
        @brief Place a column into matrix_3d at (col_idx,lay_idx)
        @param  col_idx
        @param  lay_idx
        @param  vals
    **/
    LSS_API void column(std::size_t col_idx, std::size_t lay_idx, std::valarray<double> &&vals);

    /**
        @brief Return a layer from matrix_3d at (row_idx,col_idx)
        @param  row_idx
        @param  col_idx
        @retval
    **/
    LSS_API std::slice_array<double> layer(std::size_t row_idx, std::size_t col_idx);

    /**
        @brief Place a layer into matrix_3d at (row_idx,col_idx)
        @param  row_idx
        @param  col_idx
        @param  vals
    **/
    LSS_API void layer(std::size_t row_idx, std::size_t col_idx, std::valarray<double> &&vals);
};

using matrix_3d_ptr = sptr_t<matrix_3d>;

} // namespace lss_containers

#endif ///_LSS_MATRIX_3D_HPP_
