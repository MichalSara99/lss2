/**

    @file      lss_discretization.hpp
    @brief     Represents discretizations
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once
#if !defined(_LSS_DISCRETIZATION_HPP_)
#define _LSS_DISCRETIZATION_HPP_

#include <functional>

#include "../common/lss_enumerations.hpp"
#include "../common/lss_utility.hpp"
#include "../containers/lss_matrix_2d.hpp"
#include "../containers/lss_matrix_3d.hpp"
#include "lss_grid.hpp"
#include "lss_grid_config.hpp"

using lss_containers::matrix_2d;
using lss_containers::matrix_3d;
using lss_grids::grid_config_1d_ptr;
using lss_grids::grid_config_2d_ptr;
using lss_grids::grid_config_3d_ptr;
using lss_utility::container_t;

/**
    1D discretization structure
 */
struct discretization_1d
{
  public:
    /**
     * Discretize 1D space
     *
     * \param grid_config - 1D grid config object
     * \param output - 1D container for output
     */
    static void of_space(grid_config_1d_ptr const &grid_config, container_t &output); //!<

    /**
     * Discretize function F(x) where x = first dim variable
     *
     * \param grid_config - 1D grid config object
     * \param fun - function F(x)
     * \param output - 1D container for output
     */
    static void of_function(grid_config_1d_ptr const &grid_config, std::function<double(double)> const &fun,
                            container_t &output);

    /**
     * Discretize function F(t,x) where t = time, x = first dim variable
     *
     * \param grid_config - 1D grid config object
     * \param time - time valraible t
     * \param fun - function F(t,x)
     * \param output = 1D container for output
     */
    static void of_function(grid_config_1d_ptr const &grid_config, double const &time,
                            std::function<double(double, double)> const &fun, container_t &output);
};

/**
    2D discretization structure
 */
struct discretization_2d
{
  public:
    /**
     * Discretize function F(x,y) where x=first dim variable,
     *  y = second dim variable
     *
     * \param grid_config - 2D grid config object
     * \param fun - function F(x,y)
     * \param output - 2D matrix for output
     */
    static void of_function(grid_config_2d_ptr const &grid_config, std::function<double(double, double)> const &fun,
                            matrix_2d &output);

    /**
     * Discretize function F(t,x,y) where t=time, x=first dim variable,
     *  y = second dim variable
     *
     * \param grid_config - 2D grid config object
     * \param time  - time valraible t
     * \param fun - function F(t,x,y)
     * \param output - 2D matrix for output
     */
    static void of_function(grid_config_2d_ptr const &grid_config, double const &time,
                            std::function<double(double, double, double)> const &fun, matrix_2d &output);
    /**
     * Discretize function F(t,x,y) where t=time, x=first dim variable,
     *  y = second dim variable
     *
     * \param grid_config - 2D grid config object
     * \param time - time valraible t
     * \param fun - function F(t,x,y)
     * \param rows - number of rows of the output
     * \param cols - number of columns of the output
     * \param output - 1D container for output (row-wise)
     */

    static void of_function(grid_config_2d_ptr const &grid_config, double const &time,
                            std::function<double(double, double, double)> const &fun, std::size_t const &rows,
                            std::size_t const &cols, container_t &output);
};

/**
    3D discretization structure
 */
struct discretization_3d
{
  public:
    /**
     * Discretize function F(x,y,z) where x=first dim variable,
     *  y = second dim variable, z = third dim variable
     *
     * \param grid_config - 3D grid config object
     * \param fun - function F(x,y,z)
     * \param output - 3D matrix for output
     */
    static void of_function(grid_config_3d_ptr const &grid_config,
                            std::function<double(double, double, double)> const &fun, matrix_3d &output);

    /**
     * Discretize function F(t,x,y,z) where t=time, x=first dim variable,
     *  y = second dim variable, z = third dim variable
     *
     * \param grid_config - 3D grid config object
     * \param time  - time valraible t
     * \param fun - function F(t,x,y,z)
     * \param output - 3D matrix for output
     */
    static void of_function(grid_config_3d_ptr const &grid_config, double const &time,
                            std::function<double(double, double, double, double)> const &fun, matrix_3d &output);
};

#endif ///_LSS_DISCRETIZATION_HPP_
