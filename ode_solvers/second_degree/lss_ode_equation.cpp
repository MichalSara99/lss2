#include "lss_ode_equation.hpp"

#include <functional>

#include "../../discretization/lss_discretization.hpp"
#include "../../sparse_solvers/tridiagonal/cuda_solver/lss_cuda_solver.hpp"
#include "../../sparse_solvers/tridiagonal/double_sweep_solver/lss_double_sweep_solver.hpp"
#include "../../sparse_solvers/tridiagonal/sor_solver/lss_sor_solver.hpp"
#include "../../sparse_solvers/tridiagonal/sor_solver_cuda/lss_sor_solver_cuda.hpp"
#include "../../sparse_solvers/tridiagonal/thomas_lu_solver/lss_thomas_lu_solver.hpp"

namespace lss_ode_solvers
{

using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_1d_ptr;
using lss_grids::grid_config_1d;
using lss_grids::grid_config_hints_1d_ptr;
using lss_grids::grid_transform_config_1d;
using lss_grids::grid_transform_config_1d_ptr;
using lss_transformation::boundary_transform_1d;
using lss_transformation::boundary_transform_1d_ptr;

typedef discretization_1d d_1d;

namespace implicit_solvers
{
void ode_equation::initialize(ode_data_config_ptr const &ode_data_cfg,
                              grid_config_hints_1d_ptr const &grid_config_hints, boundary_1d_pair const &boundary_pair)
{
    LSS_VERIFY(ode_data_cfg, "ode_data_config must not be null");
    LSS_VERIFY(ode_discretization_cfg_, "ode_discretization_config must not be null");
    LSS_VERIFY(std::get<0>(boundary_pair), "boundary_pair.first must not be null");
    LSS_VERIFY(std::get<1>(boundary_pair), "boundary_pair.second must not be null");
    LSS_VERIFY(ode_solver_cfg_, "ode_solver_config must not be null");
    if (!ode_solver_config_details_.empty())
    {
        auto const &it = ode_solver_config_details_.find("sor_omega");
        LSS_ASSERT(it != ode_solver_config_details_.end(), "sor_omega is not defined");
    }
    // make necessary transformations:
    // create grid_transform_config:
    grid_trans_cfg_ = std::make_shared<grid_transform_config_1d>(ode_discretization_cfg_, grid_config_hints);
    // transform original heat data:
    ode_data_trans_cfg_ = std::make_shared<ode_data_transform>(ode_data_cfg, grid_trans_cfg_);
    // transform original boundary:
    boundary_ = std::make_shared<boundary_transform_1d>(boundary_pair, grid_trans_cfg_);
}

ode_equation::ode_equation(ode_data_config_ptr const &ode_data_config,
                           ode_discretization_config_ptr const &ode_discretization_config,
                           boundary_1d_pair const &boundary_pair, grid_config_hints_1d_ptr const &grid_config_hints,
                           ode_implicit_solver_config_ptr const &ode_solver_config,
                           std::map<std::string, double> const &ode_solver_config_details)
    : ode_discretization_cfg_{ode_discretization_config}, ode_solver_cfg_{ode_solver_config},
      ode_solver_config_details_{ode_solver_config_details}
{
    initialize(ode_data_config, grid_config_hints, boundary_pair);
}

ode_equation ::~ode_equation()
{
}

void ode_equation::solve(container_t &solution)
{

    LSS_ASSERT(solution.size() > 0, "The input solution container must be initialized");
    // size of space discretization:
    const std::size_t space_size = ode_discretization_cfg_->number_of_space_points();
    // This is the proper size of the container:
    LSS_ASSERT(solution.size() == space_size, "The input solution container must have the correct size");
    // grid:
    auto const &grid_cfg = std::make_shared<grid_config_1d>(ode_discretization_cfg_);
    auto const &boundary_pair = boundary_->boundary_pair();
    const bool is_ode_nonhom_set = ode_data_trans_cfg_->is_nonhom_data_set();
    // get ode_nonhom:
    auto const &ode_nonhom = ode_data_trans_cfg_->nonhom_function();
    // prepare tridiagonal_solver pointer:
    lss_tridiagonal_solver::tridiagonal_solver_ptr solver;
    double omega_value{0.0};
    if (ode_solver_cfg_->memory_space() == memory_space_enum::Device)
    {
        if (ode_solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::CUDASolver)
        {
            solver = std::make_shared<lss_cuda_solver::cuda_solver<memory_space_enum::Device>>(space_size);
        }
        else if (ode_solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::SORSolver)
        {
            solver = std::make_shared<lss_sor_solver_cuda::sor_solver_cuda>(space_size);
            LSS_ASSERT(!ode_solver_config_details_.empty(), "ode_solver_config_details map must not be empty");
            omega_value = ode_solver_config_details_["sor_omega"];
        }
        else
        {
            throw std::exception("Not supported on Device");
        }
    }
    else if (ode_solver_cfg_->memory_space() == memory_space_enum::Host)
    {
        if (ode_solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::CUDASolver)
        {
            solver = std::make_shared<lss_cuda_solver::cuda_solver<memory_space_enum::Host>>(space_size);
        }
        else if (ode_solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::SORSolver)
        {
            solver = std::make_shared<lss_sor_solver::sor_solver>(space_size);
            LSS_ASSERT(!ode_solver_config_details_.empty(), "ode_solver_config_details map must not be empty");
            omega_value = ode_solver_config_details_["sor_omega"];
        }
        else if (ode_solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::DoubleSweepSolver)
        {
            solver = std::make_shared<lss_double_sweep_solver::double_sweep_solver>(space_size);
        }
        else if (ode_solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::ThomasLUSolver)
        {
            solver = std::make_shared<lss_thomas_lu_solver::thomas_lu_solver>(space_size);
        }
        else
        {
            throw std::exception("Not supported on Host");
        }
    }
    else
    {
        throw std::exception("Unreachable");
    }

    ode_equation_implicit_kernel kernel(solver, boundary_pair, ode_data_trans_cfg_, ode_discretization_cfg_,
                                        ode_solver_cfg_, grid_cfg);
    kernel(solution, is_ode_nonhom_set, ode_nonhom, omega_value);
}

} // namespace implicit_solvers

namespace explicit_solvers
{

} // namespace explicit_solvers

} // namespace lss_ode_solvers
