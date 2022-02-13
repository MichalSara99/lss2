#if !defined(_IMPLICIT_WAVE_SOLVER_1D_HPP_)
#define _IMPLICIT_WAVE_SOLVER_1D_HPP_

#include <map>
#include <vector>

#include "../../../../../boundaries/lss_boundary.hpp"
#include "../../../../../common/lss_macros.hpp"
#include "../../../../../common/lss_utility.hpp"
#include "../../../../../pde_solvers/1d/wave_type/lss_wave_equation.hpp"
#include "../../../configs/discretization_config_1d.hpp"
#include "../../../configs/grid_config_1d.hpp"
#include "../../configs/wave_data_config_1d.hpp"
#include "../../configs/wave_implicit_solver_config.hpp"

namespace lss
{

using implicit_wave_solver_1d = lss_pde_solvers::one_dimensional::implicit_solvers::wave_equation;

using implicit_wave_solver_1d_ptr =
    lss_utility::sptr_t<lss_pde_solvers::one_dimensional::implicit_solvers::wave_equation>;
using boundary_1d_pair = lss_boundary::boundary_1d_pair;

class implicit_wave_solver_1d_builder
{
  private:
    wave_data_config_1d_ptr wave_data_config_;
    discretization_config_1d_ptr discretization_config_;
    boundary_1d_pair boundary_pair_;
    grid_config_1d_ptr grid_cfg_;
    wave_implicit_solver_config_ptr solver_config_;
    std::map<std::string, double> solver_config_details_;

  public:
    LSS_API explicit implicit_wave_solver_1d_builder();

    LSS_API implicit_wave_solver_1d_builder &heat_data_config(const wave_data_config_1d_ptr &wave_data_config);

    LSS_API implicit_wave_solver_1d_builder &discretization_config(
        const discretization_config_1d_ptr &discretization_config);

    LSS_API implicit_wave_solver_1d_builder &boundary_pair(const boundary_1d_pair &boundary_pair);

    LSS_API implicit_wave_solver_1d_builder &grid_hints(const grid_config_1d_ptr &grid_hints);

    LSS_API implicit_wave_solver_1d_builder &solver_config(const wave_implicit_solver_config_ptr &solver_config);

    LSS_API implicit_wave_solver_1d_builder &solver_config_details(
        const std::map<std::string, double> &solver_config_details);

    LSS_API implicit_wave_solver_1d_ptr build();
};

} // namespace lss

#endif ///_IMPLICIT_WAVE_SOLVER_1D_HPP_
