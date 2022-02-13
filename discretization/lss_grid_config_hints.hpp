/**

    @file      lss_grid_config_hints.hpp
    @brief     Hints structures for the grid configuration
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once
#if !defined(_LSS_GRID_CONFIG_HINTS_HPP_)
#define _LSS_GRID_CONFIG_HINTS_HPP_

#include "../common/lss_enumerations.hpp"
#include "../common/lss_macros.hpp"
#include "../common/lss_utility.hpp"

namespace lss_grids
{
using lss_enumerations::dimension_enum;
using lss_enumerations::grid_enum;
using lss_utility::sptr_t;

/**
    @struct grid_config_hints_1d
    @brief  1D grid_config_hints structure
**/
struct grid_config_hints_1d
{
  private:
    double accumulation_point_;
    double alpha_scale_;
    grid_enum grid_;

  public:
    /**
        @brief grid_config_hints_1d object constructor
        @param accumulation_point  - point around which the spacing gets finer
        @param alpha_scale - the higher the value the finer the spacing gets
        @param grid_type - type of the grid
    **/
    explicit grid_config_hints_1d(double accumulation_point = double(0.0), double alpha_scale = double(3.0),
                                  grid_enum grid_type = grid_enum::Uniform);

    LSS_API double accumulation_point() const;

    LSS_API double alpha_scale() const;

    LSS_API grid_enum grid() const;
};

/**
    @struct grid_config_hints_2d
    @brief  2D grid_config_hints structure
**/
struct grid_config_hints_2d
{
  private:
    double accumulation_point_;
    double alpha_scale_;
    double beta_scale_;
    grid_enum grid_;

  public:
    /**
    @brief grid_config_hints_2d object constructor
    @param accumulation_point  - point around which the spacing gets finer
    @param alpha_scale - the higher the value the finer the spacing gets for first dim
    @param alpha_scale - the higher the value the finer the spacing gets for second dim
    @param grid_type - type of the grid
    **/
    explicit grid_config_hints_2d(double accumulation_point = double(0.0), double alpha_scale = double(3.0),
                                  double beta_scale = double(50.0), grid_enum grid_type = grid_enum::Uniform);

    LSS_API double accumulation_point() const;

    LSS_API double alpha_scale() const;

    LSS_API double beta_scale() const;

    LSS_API grid_enum grid() const;
};

/**
    @struct grid_config_hints_3d
    @brief  3D grid_config_hints structure
**/
struct grid_config_hints_3d
{
  private:
    double accumulation_point_;
    double alpha_scale_;
    double beta_scale_1_;
    double beta_scale_2_;
    grid_enum grid_;

  public:
    explicit grid_config_hints_3d(double accumulation_point = double(0.0), double alpha_scale = double(3.0),
                                  double beta_scale_1 = double(50.0), double beta_scale_2 = double(50.0),
                                  grid_enum grid_type = grid_enum::Uniform);

    LSS_API double accumulation_point() const;

    LSS_API double alpha_scale() const;

    LSS_API std::pair<double, double> beta_scales() const;

    LSS_API grid_enum grid() const;
};

using grid_config_hints_1d_ptr = sptr_t<grid_config_hints_1d>;
using grid_config_hints_2d_ptr = sptr_t<grid_config_hints_2d>;
using grid_config_hints_3d_ptr = sptr_t<grid_config_hints_3d>;

///////////////////////////////////////////////////////////////////////////////////////

} // namespace lss_grids
#endif ///_LSS_GRID_CONFIG_HINTS_HPP_
