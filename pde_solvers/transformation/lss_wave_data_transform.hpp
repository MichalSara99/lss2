/**

    @file      lss_wave_data_transform.hpp
    @brief     Transformations of data for wave problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once

#if !defined(_LSS_WAVE_DATA_TRANSFORM_HPP_)
#define _LSS_WAVE_DATA_TRANSFORM_HPP_

#include <functional>

#include "../../common/lss_utility.hpp"
#include "../../discretization/lss_grid.hpp"
#include "../../discretization/lss_grid_transform_config.hpp"
#include "../lss_wave_data_config.hpp"

namespace lss_pde_solvers
{

using lss_grids::grid_1d;
using lss_grids::grid_transform_config_1d_ptr;
using lss_utility::sptr_t;

/**
    1D wave_data_transform structure
 */
struct wave_data_transform_1d
{
  private:
    bool is_wave_source_set_{false};
    std::function<double(double, double)> a_coeff_{nullptr};
    std::function<double(double, double)> b_coeff_{nullptr};
    std::function<double(double, double)> c_coeff_{nullptr};
    std::function<double(double, double)> d_coeff_{nullptr};
    std::function<double(double)> init_first_coeff_{nullptr};
    std::function<double(double)> init_second_coeff_{nullptr};
    std::function<double(double, double)> src_coeff_{nullptr};

    void initialize(wave_data_config_1d_ptr const &wave_data_config,
                    grid_transform_config_1d_ptr const grid_transform_config);

    explicit wave_data_transform_1d() = delete;

  public:
    explicit wave_data_transform_1d(wave_data_config_1d_ptr const &wave_data_config,
                                    grid_transform_config_1d_ptr const grid_transform_config);

    ~wave_data_transform_1d();

    bool const &is_wave_source_set() const;

    std::function<double(double, double)> wave_source() const;

    std::function<double(double)> const &first_initial_condition() const;

    std::function<double(double)> const &second_initial_condition() const;

    std::function<double(double, double)> const &a_coefficient() const;

    std::function<double(double, double)> const &b_coefficient() const;

    std::function<double(double, double)> const &c_coefficient() const;

    std::function<double(double, double)> const &d_coefficient() const;
};

using wave_data_transform_1d_ptr = sptr_t<wave_data_transform_1d>;

} // namespace lss_pde_solvers

#endif ///_LSS_WAVE_DATA_TRANSFORM_HPP_
