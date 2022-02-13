#if !defined(_IMPLICIT_HESTON_SOLVER_2D_HPP_)
#define _IMPLICIT_HESTON_SOLVER_2D_HPP_

#include <map>
#include <vector>

#include "../../../../../boundaries/lss_boundary.hpp"
#include "../../../../../common/lss_macros.hpp"
#include "../../../../../common/lss_utility.hpp"
#include "../../../../../pde_solvers/2d/heat_type/lss_heston_equation.hpp"
#include "../../../configs/discretization_config_2d.hpp"
#include "../../../configs/grid_config_2d.hpp"
#include "../../configs/heat_data_config_2d.hpp"
#include "../../configs/heat_implicit_solver_config.hpp"
#include "../../configs/splitting_config.hpp"

namespace lss
{

using implicit_heston_solver_2d = lss_pde_solvers::two_dimensional::implicit_solvers::heston_equation;
using implicit_heston_solver_2d_ptr =
    lss_utility::sptr_t<lss_pde_solvers::two_dimensional::implicit_solvers::heston_equation>;
using boundary_2d_pair = lss_boundary::boundary_2d_pair;
using boundary_2d_ptr = lss_boundary::boundary_2d_ptr;

class implicit_heston_solver_2d_builder
{
  private:
    heat_data_config_2d_ptr heat_data_config_;
    discretization_config_2d_ptr discretization_config_;
    boundary_2d_pair horizontal_boundary_pair_;
    boundary_2d_ptr vertical_upper_boundary_;
    grid_config_2d_ptr grid_cfg_;
    splitting_config_ptr splitting_;
    heat_implicit_solver_config_ptr solver_config_;
    std::map<std::string, double> solver_config_details_;

  public:
    LSS_API explicit implicit_heston_solver_2d_builder();

    LSS_API implicit_heston_solver_2d_builder &heat_data_config(const heat_data_config_2d_ptr &heat_data_config);

    LSS_API implicit_heston_solver_2d_builder &discretization_config(
        const lss::discretization_config_2d_ptr &discretization_config);

    LSS_API implicit_heston_solver_2d_builder &vertical_upper_boundary(const boundary_2d_ptr &boundary_ptr);

    LSS_API implicit_heston_solver_2d_builder &horizontal_boundary_pair(const boundary_2d_pair &boundary_pair);

    LSS_API implicit_heston_solver_2d_builder &splitting(const splitting_config_ptr &splitting);

    LSS_API implicit_heston_solver_2d_builder &grid_hints(const grid_config_2d_ptr &grid_hints);

    LSS_API implicit_heston_solver_2d_builder &solver_config(const heat_implicit_solver_config_ptr &solver_config);

    LSS_API implicit_heston_solver_2d_builder &solver_config_details(
        const std::map<std::string, double> &solver_config_details);

    LSS_API implicit_heston_solver_2d_ptr build();
};

} // namespace lss

#endif ///_IMPLICIT_HESTON_SOLVER_2D_HPP_
