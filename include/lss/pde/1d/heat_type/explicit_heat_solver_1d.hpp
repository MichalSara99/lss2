#if !defined(_EXPLICIT_HEAT_SOLVER_1D_HPP_)
#define _EXPLICIT_HEAT_SOLVER_1D_HPP_

#include <vector>

#include "../../../../../boundaries/lss_boundary.hpp"
#include "../../../../../common/lss_macros.hpp"
#include "../../../../../common/lss_utility.hpp"
#include "../../../../../pde_solvers/1d/heat_type/lss_heat_equation.hpp"
#include "../../../configs/discretization_config_1d.hpp"
#include "../../../configs/grid_config_1d.hpp"
#include "../../configs/heat_data_config_1d.hpp"
#include "../../configs/heat_explicit_solver_config.hpp"

namespace lss
{

using explicit_heat_solver_1d = lss_pde_solvers::one_dimensional::explicit_solvers::heat_equation;
using explicit_heat_solver_1d_ptr =
    lss_utility::sptr_t<lss_pde_solvers::one_dimensional::explicit_solvers::heat_equation>;
using boundary_1d_pair = lss_boundary::boundary_1d_pair;

class explicit_heat_solver_1d_builder
{
  private:
    heat_data_config_1d_ptr heat_data_config_;
    discretization_config_1d_ptr discretization_config_;
    boundary_1d_pair boundary_pair_;
    grid_config_1d_ptr grid_cfg_;
    heat_explicit_solver_config_ptr solver_config_;

  public:
    LSS_API explicit explicit_heat_solver_1d_builder();

    LSS_API explicit_heat_solver_1d_builder &heat_data_config(const heat_data_config_1d_ptr &heat_data_config);

    LSS_API explicit_heat_solver_1d_builder &discretization_config(
        const lss::discretization_config_1d_ptr &discretization_config);

    LSS_API explicit_heat_solver_1d_builder &boundary_pair(const boundary_1d_pair &boundary_pair);

    LSS_API explicit_heat_solver_1d_builder &grid_hints(const grid_config_1d_ptr &grid_hints);

    LSS_API explicit_heat_solver_1d_builder &solver_config(const heat_explicit_solver_config_ptr &solver_config);

    LSS_API explicit_heat_solver_1d_ptr build();
};

} // namespace lss

#endif ///_EXPLICIT_HEAT_SOLVER_1D_HPP_
