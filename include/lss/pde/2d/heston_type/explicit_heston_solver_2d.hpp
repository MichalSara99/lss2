#if !defined(_EXPLICIT_HESTON_SOLVER_2D_HPP_)
#define _EXPLICIT_HESTON_SOLVER_2D_HPP_

#include <vector>

#include "../../../../../boundaries/lss_boundary.hpp"
#include "../../../../../common/lss_macros.hpp"
#include "../../../../../common/lss_utility.hpp"
#include "../../../../../pde_solvers/2d/heat_type/lss_heston_equation.hpp"
#include "../../../configs/discretization_config_2d.hpp"
#include "../../../configs/grid_config_2d.hpp"
#include "../../configs/heat_data_config_2d.hpp"
#include "../../configs/heat_explicit_solver_config.hpp"

namespace lss
{

using explicit_heston_solver_2d = lss_pde_solvers::two_dimensional::explicit_solvers::heston_equation;
using explicit_heston_solver_2d_ptr =
    lss_utility::sptr_t<lss_pde_solvers::two_dimensional::explicit_solvers::heston_equation>;
using boundary_2d_pair = lss_boundary::boundary_2d_pair;
using boundary_2d_ptr = lss_boundary::boundary_2d_ptr;

class explicit_heston_solver_2d_builder
{
  private:
    heat_data_config_2d_ptr heat_data_config_;
    discretization_config_2d_ptr discretization_config_;
    boundary_2d_pair horizontal_boundary_pair_;
    boundary_2d_ptr vertical_upper_boundary_;
    grid_config_2d_ptr grid_cfg_;
    heat_explicit_solver_config_ptr solver_config_;

  public:
    LSS_API explicit explicit_heston_solver_2d_builder();

    LSS_API explicit_heston_solver_2d_builder &heat_data_config(const heat_data_config_2d_ptr &heat_data_config);

    LSS_API explicit_heston_solver_2d_builder &discretization_config(
        const lss::discretization_config_2d_ptr &discretization_config);

    LSS_API explicit_heston_solver_2d_builder &vertical_upper_boundary(const boundary_2d_ptr &boundary_ptr);

    LSS_API explicit_heston_solver_2d_builder &horizontal_boundary_pair(const boundary_2d_pair &boundary_pair);

    LSS_API explicit_heston_solver_2d_builder &grid_hints(const grid_config_2d_ptr &grid_hints);

    LSS_API explicit_heston_solver_2d_builder &solver_config(const heat_explicit_solver_config_ptr &solver_config);

    LSS_API explicit_heston_solver_2d_ptr build();
};

} // namespace lss

#endif ///_EXPLICIT_HESTON_SOLVER_2D_HPP_
