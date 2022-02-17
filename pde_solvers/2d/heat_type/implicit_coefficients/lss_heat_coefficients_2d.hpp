/**

    @file      lss_heat_coefficients_2d.hpp
    @brief     Implicit coefficients for 2D heat problems
    @details   ~
    @author    Michal Sara
    @date      14.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_COEFFICIENTS_2D_HPP_)
#define _LSS_HEAT_COEFFICIENTS_2D_HPP_

#include "../../../../common/lss_utility.hpp"
#include "../../../../discretization/lss_discretization.hpp"
#include "../../../lss_pde_discretization_config.hpp"
#include "../../../lss_splitting_method_config.hpp"
#include "../../../transformation/lss_heat_data_transform.hpp"

namespace lss_pde_solvers
{
namespace two_dimensional
{

using lss_utility::range_ptr;
using lss_utility::sptr_t;

/**
    heat_coefficients_2d object
 */
struct heat_coefficients_2d
{
  public:
    // scheme constant coefficients:
    double alpha_, beta_, gamma_, delta_, ni_, rho_, zeta_, k_;
    std::size_t space_size_x_, space_size_y_;
    range_ptr rangex_, rangey_;
    // theta variable:
    double theta_;
    // functional coefficients:
    std::function<double(double, double, double)> M_;
    std::function<double(double, double, double)> M_tilde_;
    std::function<double(double, double, double)> P_;
    std::function<double(double, double, double)> P_tilde_;
    std::function<double(double, double, double)> Z_;
    std::function<double(double, double, double)> W_;
    std::function<double(double, double, double)> C_;
    std::function<double(double, double, double)> D_;
    std::function<double(double, double, double)> E_;
    std::function<double(double, double, double)> F_;

  private:
    void initialize(pde_discretization_config_2d_ptr const &discretization_config,
                    splitting_method_config_ptr const &splitting_config);

    void initialize_coefficients(heat_data_transform_2d_ptr const &heat_data_config);

  public:
    heat_coefficients_2d() = delete;

    heat_coefficients_2d(heat_data_transform_2d_ptr const &heat_data_config,
                         pde_discretization_config_2d_ptr const &discretization_config,
                         splitting_method_config_ptr const splitting_config, double const &theta);
};

using heat_coefficients_2d_ptr = sptr_t<heat_coefficients_2d>;

} // namespace two_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_COEFFICIENTS_2D_HPP_
