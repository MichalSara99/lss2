#if !defined(_IMPLICIT_ODE_TP_SOLVER_HPP_)
#define _IMPLICIT_ODE_TP_SOLVER_HPP_

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../ode_solvers/second_degree/lss_ode_equation.hpp"
#include "../../configs/discretization_config_ode.hpp"
#include "../../configs/grid_config_1d.hpp"
#include "../configs/ode_coefficient_data_config.hpp"
#include "../configs/ode_data_config.hpp"
#include "../configs/ode_nonhom_data_config.hpp"
#include "../configs/ode_solver_config.hpp"
#include <map>

namespace lss
{

using implicit_ode_tp_solver = lss_ode_solvers::implicit_solvers::ode_equation;
using implicit_ode_tp_solver_ptr = lss_utility::sptr_t<lss_ode_solvers::implicit_solvers::ode_equation>;
using boundary_1d_pair = lss_boundary::boundary_1d_pair;

struct implicit_ode_tp_solver_builder
{
  private:
    ode_data_config_ptr ode_data_cfg_;
    discretization_config_ode_ptr discretization_cfg_;
    boundary_1d_pair boundary_pair_;
    grid_config_1d_ptr grid_cfg_;
    ode_implicit_solver_config_ptr ode_solver_config_;
    std::map<std::string, double> ode_solver_config_details_;

  public:
    LSS_API explicit implicit_ode_tp_solver_builder();

    LSS_API implicit_ode_tp_solver_builder &ode_data_config(const ode_data_config_ptr &ode_data_config);

    LSS_API implicit_ode_tp_solver_builder &discretization_config(
        const discretization_config_ode_ptr &discretization_config);

    LSS_API implicit_ode_tp_solver_builder &boundary_pair(const boundary_1d_pair &boundary_pair);

    LSS_API implicit_ode_tp_solver_builder &grid_hints(const grid_config_1d_ptr &grid_hints);

    LSS_API implicit_ode_tp_solver_builder &solver_config(const ode_implicit_solver_config_ptr &solver_config);

    LSS_API implicit_ode_tp_solver_builder &solver_config_details(
        const std::map<std::string, double> &solver_config_details);

    LSS_API implicit_ode_tp_solver_ptr build();
};

} // namespace lss

#endif ///_IMPLICIT_ODE_TP_SOLVER_HPP_
