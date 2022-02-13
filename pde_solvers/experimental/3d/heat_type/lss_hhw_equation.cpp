#include "lss_hhw_equation.hpp"

#include "../../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../discretization/lss_grid.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
//#include "lss_heston_equation_explicit_kernel.hpp"
#include "lss_hhw_equation_implicit_kernel.hpp"

namespace lss_pde_solvers
{

using lss_boundary::dirichlet_boundary_3d;
using lss_boundary::neumann_boundary_3d;
using lss_grids::grid_config_3d;
using lss_grids::grid_transform_config_3d;

namespace three_dimensional
{

namespace implicit_solvers
{

void hhw_equation::initialize(heat_data_config_3d_ptr const &heat_data_cfg,
                              grid_config_hints_3d_ptr const &grid_config_hints,
                              boundary_3d_pair const &x_boundary_pair, boundary_3d_ptr const &y_upper_boundary_ptr,
                              boundary_3d_pair const &z_boundary_pair)
{
    // verify and check:
    LSS_VERIFY(heat_data_cfg, "heat_data_config must not be null");
    LSS_VERIFY(discretization_cfg_, "discretization_config must not be null");

    if (auto y_ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(y_upper_boundary_ptr))
    {
        LSS_VERIFY(y_ptr, "y_upper_boundary_ptr must be of dirichlet type only");
    }

    if (auto x_ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(std::get<0>(x_boundary_pair)))
    {
        LSS_VERIFY(x_ptr, "x_boundary_pair.first must be of dirichlet type only");
    }
    if (auto x_ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(std::get<1>(x_boundary_pair)))
    {
        LSS_VERIFY(x_ptr, "x_boundary_pair.second must be of neumann type only");
    }

    if (auto z_ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(std::get<0>(z_boundary_pair)))
    {
        LSS_VERIFY(z_ptr, "z_boundary_pair.first must be of dirichlet type only");
    }
    if (auto z_ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(std::get<1>(z_boundary_pair)))
    {
        LSS_VERIFY(z_ptr, "z_boundary_pair.second must be of neumann type only");
    }

    LSS_VERIFY(splitting_method_cfg_, "splitting_method_config must not be null");
    LSS_VERIFY(solver_cfg_, "solver_config must not be null");
    LSS_VERIFY(grid_config_hints, "grid_config_hints must not be null");
    if (!solver_config_details_.empty())
    {
        auto const &it = solver_config_details_.find("sor_omega");
        LSS_ASSERT(it != solver_config_details_.end(), "sor_omega is not defined");
    }
    // make necessary transformations:
    // create grid_transform_config:
    grid_trans_cfg_ = std::make_shared<grid_transform_config_3d>(discretization_cfg_, grid_config_hints);
    // transform original heat data:
    heat_data_trans_cfg_ = std::make_shared<heat_data_transform_3d>(heat_data_cfg, grid_trans_cfg_);
    // transform original boundary:
    hhw_boundary_ = std::make_shared<hhw_boundary_transform>(x_boundary_pair, y_upper_boundary_ptr, z_boundary_pair,
                                                             grid_trans_cfg_);
}

hhw_equation::hhw_equation(heat_data_config_3d_ptr const &heat_data_config,
                           pde_discretization_config_3d_ptr const &discretization_config,
                           boundary_3d_pair const &x_boundary_pair, boundary_3d_ptr const &y_upper_boundary_ptr,
                           boundary_3d_pair const &z_boundary_pair,
                           splitting_method_config_ptr const &splitting_method_config,
                           grid_config_hints_3d_ptr const &grid_config_hints,
                           heat_implicit_solver_config_ptr const &solver_config,
                           std::map<std::string, double> const &solver_config_details)
    : discretization_cfg_{discretization_config}, splitting_method_cfg_{splitting_method_config},
      solver_cfg_{solver_config}, solver_config_details_{solver_config_details}
{
    initialize(heat_data_config, grid_config_hints, x_boundary_pair, y_upper_boundary_ptr, z_boundary_pair);
}

hhw_equation ::~hhw_equation()
{
}

void hhw_equation::solve(container_3d<by_enum::RowPlane> &solution)
{
    LSS_ASSERT((solution.rows()) > 0 && (solution.columns() > 0) && (solution.layers() > 0),
               "The input solution container must be initialized");

    // get space ranges:
    const auto &spaces = discretization_cfg_->space_range();
    // across X:
    const auto space_x = std::get<0>(spaces);
    // across Y:
    const auto space_y = std::get<1>(spaces);
    // across Z:
    const auto space_z = std::get<2>(spaces);
    // size of spaces discretization:
    const auto &space_sizes = discretization_cfg_->number_of_space_points();
    const std::size_t space_size_x = std::get<0>(space_sizes);
    const std::size_t space_size_y = std::get<1>(space_sizes);
    const std::size_t space_size_z = std::get<2>(space_sizes);
    // This is the proper size of the container:
    LSS_ASSERT((solution.columns() == space_size_y) && (solution.rows() == space_size_x) &&
                   (solution.layers() == space_size_z),
               "The input solution container must have the correct size");
    // create grid_config:
    auto const &grid_cfg = std::make_shared<grid_config_3d>(discretization_cfg_);
    auto const &y_boundary_ptr = hhw_boundary_->y_upper_boundary();
    auto const &x_boundary_pair_ptr = hhw_boundary_->x_boundary_pair();
    auto const &z_boundary_pair_ptr = hhw_boundary_->z_boundary_pair();
    // create container to carry previous solution:
    container_3d<by_enum::RowPlane> prev_sol(space_size_x, space_size_y, space_size_z, double{});
    // create container to carry next solution:
    container_3d<by_enum::RowPlane> next_sol(space_size_x, space_size_y, space_size_z, double{});
    // discretize initial condition
    d_3d::of_function(grid_cfg, heat_data_trans_cfg_->initial_condition(), prev_sol);
    // get heat_source:
    const bool is_heat_source_set = heat_data_trans_cfg_->is_heat_source_set();
    // get heat_source:
    auto const &heat_source = heat_data_trans_cfg_->heat_source();

    if (solver_cfg_->memory_space() == memory_space_enum::Device)
    {
        if (solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::CUDASolver)
        {
            typedef hhw_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::CUDASolver>
                dev_cu_solver;

            dev_cu_solver solver(x_boundary_pair_ptr, y_boundary_ptr, z_boundary_pair_ptr, heat_data_trans_cfg_,
                                 discretization_cfg_, splitting_method_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol, next_sol, is_heat_source_set, heat_source);
            solution = prev_sol;
        }
        else if (solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::SORSolver)
        {
            typedef hhw_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::SORSolver>
                dev_sor_solver;

            LSS_ASSERT(!solver_config_details_.empty(), "solver_config_details map must not be empty");
            double omega_value = solver_config_details_["sor_omega"];
            dev_sor_solver solver(x_boundary_pair_ptr, y_boundary_ptr, z_boundary_pair_ptr, heat_data_trans_cfg_,
                                  discretization_cfg_, splitting_method_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol, next_sol, is_heat_source_set, heat_source, omega_value);
            solution = prev_sol;
        }
        else
        {
            throw std::exception("Not supported on Device");
        }
    }
    else if (solver_cfg_->memory_space() == memory_space_enum::Host)
    {
        if (solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::CUDASolver)
        {
            typedef hhw_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::CUDASolver>
                host_cu_solver;

            host_cu_solver solver(x_boundary_pair_ptr, y_boundary_ptr, z_boundary_pair_ptr, heat_data_trans_cfg_,
                                  discretization_cfg_, splitting_method_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol, next_sol, is_heat_source_set, heat_source);
            solution = next_sol;
        }
        else if (solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::SORSolver)
        {
            typedef hhw_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::SORSolver>
                host_sor_solver;

            LSS_ASSERT(!solver_config_details_.empty(), "solver_config_details map must not be empty");
            double omega_value = solver_config_details_["sor_omega"];
            host_sor_solver solver(x_boundary_pair_ptr, y_boundary_ptr, z_boundary_pair_ptr, heat_data_trans_cfg_,
                                   discretization_cfg_, splitting_method_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol, next_sol, is_heat_source_set, heat_source, omega_value);
            solution = next_sol;
        }
        else if (solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::DoubleSweepSolver)
        {
            typedef hhw_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::DoubleSweepSolver>
                host_dss_solver;
            host_dss_solver solver(x_boundary_pair_ptr, y_boundary_ptr, z_boundary_pair_ptr, heat_data_trans_cfg_,
                                   discretization_cfg_, splitting_method_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol, next_sol, is_heat_source_set, heat_source);
            solution = next_sol;
        }
        else if (solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::ThomasLUSolver)
        {
            typedef hhw_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::ThomasLUSolver>
                host_lus_solver;
            host_lus_solver solver(x_boundary_pair_ptr, y_boundary_ptr, z_boundary_pair_ptr, heat_data_trans_cfg_,
                                   discretization_cfg_, splitting_method_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol, next_sol, is_heat_source_set, heat_source);
            solution = next_sol;
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
}

// void hhw_equation::solve(container_3d<by_enum::Row> &solutions)
//{
//     throw std::exception("Not implemented.");
// }

} // namespace implicit_solvers

} // namespace three_dimensional

} // namespace lss_pde_solvers
