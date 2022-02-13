#include "lss_ode_equation_implicit_kernel.hpp"

namespace lss_ode_solvers
{

ode_equation_implicit_kernel::ode_equation_implicit_kernel(tridiagonal_solver_ptr solver,
                                                           boundary_1d_pair const &boundary_pair,
                                                           ode_data_transform_ptr const &ode_data_config,
                                                           ode_discretization_config_ptr const &discretization_config,
                                                           ode_implicit_solver_config_ptr const &solver_config,
                                                           grid_config_1d_ptr const &grid_config)
    : solver_{solver}, boundary_pair_{boundary_pair}, ode_data_cfg_{ode_data_config},
      discretization_cfg_{discretization_config}, solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void ode_equation_implicit_kernel::operator()(container_t &solution, bool is_ode_nonhom_set,
                                              std::function<double(double)> const &ode_nonhom, double omega_value)
{
    // create a ode coefficient holder:
    auto const ode_coeff_holder = std::make_shared<ode_implicit_coefficients>(ode_data_cfg_, discretization_cfg_);
    // create and set up the solver:
    solver_->set_omega(omega_value);
    solver_->set_factorization(solver_cfg_->tridiagonal_factorization());
    auto const &solver_method_ptr = std::make_shared<ode_implicit_solver_method>(solver_, ode_coeff_holder, grid_cfg_);
    if (is_ode_nonhom_set)
    {
        solver_method_ptr->solve(boundary_pair_, ode_nonhom, solution);
    }
    else
    {
        solver_method_ptr->solve(boundary_pair_, solution);
    }
}

} // namespace lss_ode_solvers
