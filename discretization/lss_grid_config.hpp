/**

    @file      lss_grid_config.hpp
    @brief     Grid configuration objects
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once

#if !defined(_LSS_GRID_CONFIG_HPP_)
#define _LSS_GRID_CONFIG_HPP_

#include "../common/lss_utility.hpp"
#include "../pde_solvers/lss_pde_discretization_config.hpp"
#include "lss_discretization_config.hpp"

namespace lss_grids
{
using lss_discretization::discretization_config_1d_ptr;
using lss_pde_solvers::pde_discretization_config_2d_ptr;
using lss_pde_solvers::pde_discretization_config_3d_ptr;
using lss_utility::sptr_t;

/**
    @struct grid_config_1d
    @brief  1D Grid configuration
**/
struct grid_config_1d
{
  private:
    double step_;

  public:
    explicit grid_config_1d(discretization_config_1d_ptr const &discretization_config);

    double step() const;

    std::size_t index_of(double zeta);

    double value_for(std::size_t idx);
};

using grid_config_1d_ptr = sptr_t<grid_config_1d>;

/**
    @struct grid_config_2d
    @brief  2D Grid configuration
**/
struct grid_config_2d
{
  private:
    double step_1_;
    double step_2_;
    grid_config_1d_ptr grid_1_;
    grid_config_1d_ptr grid_2_;

  public:
    explicit grid_config_2d(pde_discretization_config_2d_ptr const &discretization_config);

    grid_config_1d_ptr const &grid_1() const;
    grid_config_1d_ptr const &grid_2() const;

    double step_1() const;
    double step_2() const;

    std::size_t index_of_1(double zeta);
    std::size_t index_of_2(double eta);

    double value_for_1(std::size_t idx);
    double value_for_2(std::size_t idx);
};

using grid_config_2d_ptr = sptr_t<grid_config_2d>;

/**
    @struct grid_config_3d
    @brief  3D Grid configuration
**/
struct grid_config_3d
{
  private:
    double step_1_;
    double step_2_;
    double step_3_;
    grid_config_1d_ptr grid_1_;
    grid_config_1d_ptr grid_2_;
    grid_config_1d_ptr grid_3_;
    grid_config_2d_ptr grid_12_;
    grid_config_2d_ptr grid_21_;
    grid_config_2d_ptr grid_13_;
    grid_config_2d_ptr grid_31_;
    grid_config_2d_ptr grid_23_;
    grid_config_2d_ptr grid_32_;

  public:
    explicit grid_config_3d(pde_discretization_config_3d_ptr const &discretization_config);

    grid_config_1d_ptr const &grid_1() const;
    grid_config_1d_ptr const &grid_2() const;
    grid_config_1d_ptr const &grid_3() const;

    grid_config_2d_ptr const &grid_12() const;
    grid_config_2d_ptr const &grid_21() const;
    grid_config_2d_ptr const &grid_13() const;
    grid_config_2d_ptr const &grid_31() const;
    grid_config_2d_ptr const &grid_23() const;
    grid_config_2d_ptr const &grid_32() const;

    double step_1() const;
    double step_2() const;
    double step_3() const;

    std::size_t index_of_1(double zeta);
    std::size_t index_of_2(double eta);
    std::size_t index_of_3(double ny);

    double value_for_1(std::size_t idx);
    double value_for_2(std::size_t idx);
    double value_for_3(std::size_t idx);
};

using grid_config_3d_ptr = sptr_t<grid_config_3d>;

///////////////////////////////////////////////////////////////////////////////////////

} // namespace lss_grids
#endif ///_LSS_GRID_CONFIG_HPP_
