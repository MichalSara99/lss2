/**

    @file      lss_ode_data_transform.hpp
    @brief     ODE data transformations
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once
#if !defined(_LSS_ODE_DATA_TRANSFORM_HPP_)
#define _LSS_ODE_DATA_TRANSFORM_HPP_

#include <functional>

#include "../../common/lss_utility.hpp"
#include "../../discretization/lss_grid.hpp"
#include "../../discretization/lss_grid_transform_config.hpp"
#include "../lss_ode_data_config.hpp"

namespace lss_ode_solvers
{

using lss_grids::grid_1d;
using lss_grids::grid_2d;
using lss_grids::grid_transform_config_1d_ptr;
using lss_grids::grid_transform_config_2d_ptr;
using lss_utility::sptr_t;

/**
    @struct ode_data_transform
    @brief
**/
struct ode_data_transform
{
  private:
    bool is_nonhom_data_set_{false};
    std::function<double(double)> a_coeff_{nullptr};
    std::function<double(double)> b_coeff_{nullptr};
    std::function<double(double)> nonhom_coeff_{nullptr};

    void initialize(ode_data_config_ptr const &ode_data_config,
                    grid_transform_config_1d_ptr const grid_transform_config);

    explicit ode_data_transform() = delete;

  public:
    explicit ode_data_transform(ode_data_config_ptr const &ode_data_config,
                                grid_transform_config_1d_ptr const grid_transform_config);

    ~ode_data_transform();
    bool const &is_nonhom_data_set() const;

    std::function<double(double)> nonhom_function() const;

    std::function<double(double)> const &a_coefficient() const;

    std::function<double(double)> const &b_coefficient() const;
};

using ode_data_transform_ptr = sptr_t<ode_data_transform>;

} // namespace lss_ode_solvers

#endif ///_LSS_ODE_DATA_TRANSFORM_HPP_
