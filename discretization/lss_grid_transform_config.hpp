/**

    @file      lss_grid_transform_config.hpp
    @brief     Transformed grid configurations
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once

#if !defined(_LSS_GRID_TRANSFORM_CONFIG_HPP_)
#define _LSS_GRID_TRANSFORM_CONFIG_HPP_

#include "../common/lss_utility.hpp"
#include "../pde_solvers/lss_pde_discretization_config.hpp"
#include "lss_discretization_config.hpp"
#include "lss_grid_config_hints.hpp"
#include <cmath>
#include <functional>

namespace lss_grids
{
using lss_discretization::discretization_config_1d_ptr;
using lss_pde_solvers::pde_discretization_config_2d_ptr;
using lss_pde_solvers::pde_discretization_config_3d_ptr;
using lss_utility::sptr_t;

/**
    1D grid_transform_config structure
 */
struct grid_transform_config_1d
{
  private:
    double init_;
    double alpha_;
    double c_[2];
    std::function<double(double)> a_der_;
    std::function<double(double)> b_der_;

    void initialize(double low, double high, grid_config_hints_1d_ptr const &grid_hints);

  public:
    explicit grid_transform_config_1d(discretization_config_1d_ptr const &discretization_config,
                                      grid_config_hints_1d_ptr const &grid_hints);

    std::function<double(double)> const &a_derivative() const;
    std::function<double(double)> const &b_derivative() const;

    double value_for(double zeta);
};

/**
    2D grid_config structure
 */
struct grid_transform_config_2d
{
  private:
    double init_1_;
    double init_2_;
    double alpha_;
    double beta_;
    double c_[2];
    double d_;

    std::function<double(double)> a_der_;
    std::function<double(double)> b_der_;
    std::function<double(double)> c_der_;
    std::function<double(double)> d_der_;

    void initialize(pde_discretization_config_2d_ptr const &discretization_config,
                    grid_config_hints_2d_ptr const &grid_hints);

  public:
    explicit grid_transform_config_2d(pde_discretization_config_2d_ptr const &discretization_config,
                                      grid_config_hints_2d_ptr const &grid_hints);

    std::function<double(double)> const &a_derivative() const;
    std::function<double(double)> const &b_derivative() const;
    std::function<double(double)> const &c_derivative() const;
    std::function<double(double)> const &d_derivative() const;

    double value_for_1(double zeta);
    double value_for_2(double eta);
};

/**
    3D grid_config structure
 */
struct grid_transform_config_3d
{
  private:
    double init_1_;
    double init_2_;
    double init_3_;
    double alpha_;
    double beta_1_, beta_2_;
    double c_[2];
    double e_;
    double d_;

    std::function<double(double)> a_1_der_;
    std::function<double(double)> a_2_der_;
    std::function<double(double)> b_1_der_;
    std::function<double(double)> b_2_der_;
    std::function<double(double)> c_1_der_;
    std::function<double(double)> c_2_der_;

    void initialize(pde_discretization_config_3d_ptr const &discretization_config,
                    grid_config_hints_3d_ptr const &grid_hints);

  public:
    explicit grid_transform_config_3d(pde_discretization_config_3d_ptr const &discretization_config,
                                      grid_config_hints_3d_ptr const &grid_hints);

    std::function<double(double)> const &a_1_derivative() const;
    std::function<double(double)> const &a_2_derivative() const;
    std::function<double(double)> const &b_1_derivative() const;
    std::function<double(double)> const &b_2_derivative() const;
    std::function<double(double)> const &c_1_derivative() const;
    std::function<double(double)> const &c_2_derivative() const;

    double value_for_1(double zeta);
    double value_for_2(double eta);
    double value_for_3(double ny);
};

using grid_transform_config_1d_ptr = sptr_t<grid_transform_config_1d>;
using grid_transform_config_2d_ptr = sptr_t<grid_transform_config_2d>;
using grid_transform_config_3d_ptr = sptr_t<grid_transform_config_3d>;

///////////////////////////////////////////////////////////////////////////////////////

} // namespace lss_grids
#endif ///_LSS_GRID_TRANSFORM_CONFIG_HPP_
