/**

    @file      lss_grid.hpp
    @brief     Represents Grid functions
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once

#if !defined(_LSS_GRID_HPP_)
#define _LSS_GRID_HPP_

#include <functional>

#include "../common/lss_enumerations.hpp"
#include "../common/lss_utility.hpp"
#include "lss_grid_config.hpp"
#include "lss_grid_transform_config.hpp"

namespace lss_grids
{

/**
    1D uniform_grid structure
 */

struct grid_1d
{

    static double value(grid_config_1d_ptr const &grid_config, std::size_t const &idx);

    static double step(grid_config_1d_ptr const &grid_config);

    static std::size_t index_of(grid_config_1d_ptr const &grid_config, double zeta);

    static double transformed_value(grid_transform_config_1d_ptr const &grid_trans_config, double zeta);
};

/**
    2D grid structure
 */

struct grid_2d
{
    static double value_1(grid_config_2d_ptr const &grid_config, std::size_t const &idx);

    static double value_2(grid_config_2d_ptr const &grid_config, std::size_t const &idx);

    static double step_1(grid_config_2d_ptr const &grid_config);

    static double step_2(grid_config_2d_ptr const &grid_config);

    static std::size_t index_of_1(grid_config_2d_ptr const &grid_config, double zeta);

    static std::size_t index_of_2(grid_config_2d_ptr const &grid_config, double eta);

    static double transformed_value_1(grid_transform_config_2d_ptr const &grid_trans_config, double zeta);

    static double transformed_value_2(grid_transform_config_2d_ptr const &grid_trans_config, double eta);
};

/**
    3D grid structure
 */

struct grid_3d
{
    static double value_1(grid_config_3d_ptr const &grid_config, std::size_t const &idx);

    static double value_2(grid_config_3d_ptr const &grid_config, std::size_t const &idx);

    static double value_3(grid_config_3d_ptr const &grid_config, std::size_t const &idx);

    static double step_1(grid_config_3d_ptr const &grid_config);

    static double step_2(grid_config_3d_ptr const &grid_config);

    static double step_3(grid_config_3d_ptr const &grid_config);

    static std::size_t index_of_1(grid_config_3d_ptr const &grid_config, double zeta);

    static std::size_t index_of_2(grid_config_3d_ptr const &grid_config, double eta);

    static std::size_t index_of_3(grid_config_3d_ptr const &grid_config, double ny);

    static double transformed_value_1(grid_transform_config_3d_ptr const &grid_trans_config, double zeta);

    static double transformed_value_2(grid_transform_config_3d_ptr const &grid_trans_config, double eta);

    static double transformed_value_3(grid_transform_config_3d_ptr const &grid_trans_config, double ny);
};

using grid_1d_ptr = sptr_t<grid_1d>;
using grid_2d_ptr = sptr_t<grid_2d>;
using grid_3d_ptr = sptr_t<grid_3d>;

} // namespace lss_grids
#endif ///_LSS_GRID_HPP_
