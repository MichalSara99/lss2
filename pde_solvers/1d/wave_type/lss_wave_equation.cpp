#include "lss_wave_equation.hpp"

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../common/lss_macros.hpp"
#include "../../../containers/lss_matrix_2d.hpp"
#include "../../../discretization/lss_discretization.hpp"
#include "../../../transformation/lss_boundary_transform.hpp"
#include "../../lss_pde_discretization_config.hpp"
#include "../../lss_wave_data_config.hpp"
#include "../../lss_wave_solver_config.hpp"
#include "../../transformation/lss_wave_data_transform.hpp"
#include "lss_wave_equation_explicit_kernel.hpp"
#include "lss_wave_equation_implicit_kernel.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_enumerations::grid_enum;
using lss_utility::valcopy;
using d_1d = discretization_1d;

namespace implicit_solvers
{

void wave_equation::initialize(wave_data_config_1d_ptr const &wave_data_cfg,
                               grid_config_hints_1d_ptr const &grid_config_hints, boundary_1d_pair const &boundary_pair)
{
    LSS_VERIFY(wave_data_cfg, "wave_data_config must not be null");
    LSS_VERIFY(discretization_cfg_, "discretization_config must not be null");
    LSS_VERIFY(std::get<0>(boundary_pair), "boundary_pair.first must not be null");
    LSS_VERIFY(std::get<1>(boundary_pair), "boundary_pair.second must not be null");
    LSS_VERIFY(solver_cfg_, "solver_config must not be null");
    if (!solver_config_details_.empty())
    {
        auto const &it = solver_config_details_.find("sor_omega");
        LSS_ASSERT(it != solver_config_details_.end(), "sor_omega is not defined");
    }
    // make necessary transformations:
    // create grid_transform_config:
    grid_trans_cfg_ = std::make_shared<grid_transform_config_1d>(discretization_cfg_, grid_config_hints);
    // transform original wave data:
    wave_data_trans_cfg_ = std::make_shared<wave_data_transform_1d>(wave_data_cfg, grid_trans_cfg_);
    // transform original boundary:
    boundary_ = std::make_shared<boundary_transform_1d>(boundary_pair, grid_trans_cfg_);
}

wave_equation::wave_equation(wave_data_config_1d_ptr const &wave_data_config,
                             pde_discretization_config_1d_ptr const &discretization_config,
                             boundary_1d_pair const &boundary_pair, grid_config_hints_1d_ptr const &grid_config_hints,
                             wave_implicit_solver_config_ptr const &solver_config,
                             std::map<std::string, double> const &solver_config_details)
    : discretization_cfg_{discretization_config}, solver_cfg_{solver_config}, solver_config_details_{
                                                                                  solver_config_details}
{
    initialize(wave_data_config, grid_config_hints, boundary_pair);
}

wave_equation::~wave_equation()
{
}

void wave_equation::solve(container_t &solution)
{
    LSS_ASSERT(solution.size() > 0, "The input solution container must be initialized");
    // size of space discretization:
    const std::size_t space_size = discretization_cfg_->number_of_space_points();
    // This is the proper size of the container:
    LSS_ASSERT(solution.size() == space_size, "The input solution container must have the correct size");
    // grid:
    auto const &grid_cfg = std::make_shared<grid_config_1d>(discretization_cfg_);
    auto const &boundary_pair = boundary_->boundary_pair();
    // create container to carry previous solution:
    container_t prev_sol_0(double{}, space_size);
    // discretize initial condition
    d_1d::of_function(grid_cfg, wave_data_trans_cfg_->first_initial_condition(), prev_sol_0);
    // create container to carry second initial condition:
    container_t prev_sol_1(double{}, space_size);
    // discretize first order initial condition:
    d_1d::of_function(grid_cfg, wave_data_trans_cfg_->second_initial_condition(), prev_sol_1);
    // create container to carry new solution:
    container_t next_sol(double{}, space_size);
    // get heat_source:
    const bool is_wave_source_set = wave_data_trans_cfg_->is_wave_source_set();
    // get heat_source:
    auto const &wave_source = wave_data_trans_cfg_->wave_source();

    if (solver_cfg_->memory_space() == memory_space_enum::Device)
    {
        if (solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::CUDASolver)
        {
            typedef wave_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::CUDASolver>
                dev_cu_solver;

            dev_cu_solver solver(boundary_pair, wave_data_trans_cfg_, discretization_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol_0, prev_sol_1, next_sol, is_wave_source_set, wave_source);
            valcopy(solution, prev_sol_1);
        }
        else if (solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::SORSolver)
        {
            typedef wave_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::SORSolver>
                dev_sor_solver;

            LSS_ASSERT(!solver_config_details_.empty(), "solver_config_details map must not be empty");
            auto const omega_value = solver_config_details_["sor_omega"];
            dev_sor_solver solver(boundary_pair, wave_data_trans_cfg_, discretization_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol_0, prev_sol_1, next_sol, is_wave_source_set, wave_source, omega_value);
            valcopy(solution, prev_sol_1);
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
            typedef wave_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::CUDASolver>
                host_cu_solver;
            host_cu_solver solver(boundary_pair, wave_data_trans_cfg_, discretization_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol_0, prev_sol_1, next_sol, is_wave_source_set, wave_source);
            valcopy(solution, prev_sol_1);
        }
        else if (solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::SORSolver)
        {
            typedef wave_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::SORSolver>
                host_sor_solver;

            LSS_ASSERT(!solver_config_details_.empty(), "solver_config_details map must not be empty");
            auto const omega_value = solver_config_details_["sor_omega"];
            host_sor_solver solver(boundary_pair, wave_data_trans_cfg_, discretization_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol_0, prev_sol_1, next_sol, is_wave_source_set, wave_source, omega_value);
            valcopy(solution, prev_sol_1);
        }
        else if (solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::DoubleSweepSolver)
        {
            typedef wave_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::DoubleSweepSolver>
                host_dss_solver;
            host_dss_solver solver(boundary_pair, wave_data_trans_cfg_, discretization_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol_0, prev_sol_1, next_sol, is_wave_source_set, wave_source);
            valcopy(solution, prev_sol_1);
        }
        else if (solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::ThomasLUSolver)
        {
            typedef wave_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::ThomasLUSolver>
                host_lus_solver;
            host_lus_solver solver(boundary_pair, wave_data_trans_cfg_, discretization_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol_0, prev_sol_1, next_sol, is_wave_source_set, wave_source);
            valcopy(solution, prev_sol_1);
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

void wave_equation::solve(matrix_2d &solutions)
{
    LSS_ASSERT(((solutions.columns() > 0) && (solutions.rows() > 0)),
               "The input solution 2D container must be initialized");
    // size of space discretization:
    const std::size_t space_size = discretization_cfg_->number_of_space_points();
    // size of time discretization:
    const std::size_t time_size = discretization_cfg_->number_of_time_points();
    // This is the proper size of the container:
    LSS_ASSERT((solutions.rows() == time_size) && (solutions.columns() == space_size),
               "The input solution 2D container must have the correct size");
    // grid:
    auto const &grid_cfg = std::make_shared<grid_config_1d>(discretization_cfg_);
    auto const &boundary_pair = boundary_->boundary_pair();
    // create container to carry previous solution:
    container_t prev_sol_0(double{}, space_size);
    // discretize initial condition
    d_1d::of_function(grid_cfg, wave_data_trans_cfg_->first_initial_condition(), prev_sol_0);
    // create container to carry second initial condition:
    container_t prev_sol_1(double{}, space_size);
    // discretize first order initial condition:
    d_1d::of_function(grid_cfg, wave_data_trans_cfg_->second_initial_condition(), prev_sol_1);
    // create container to carry new solution:
    container_t next_sol(double{}, space_size);
    // get heat_source:
    const bool is_wave_source_set = wave_data_trans_cfg_->is_wave_source_set();
    // get heat_source:
    auto const &wave_source = wave_data_trans_cfg_->wave_source();

    if (solver_cfg_->memory_space() == memory_space_enum::Device)
    {
        if (solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::CUDASolver)
        {
            typedef wave_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::CUDASolver>
                dev_cu_solver;

            dev_cu_solver solver(boundary_pair, wave_data_trans_cfg_, discretization_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol_0, prev_sol_1, next_sol, is_wave_source_set, wave_source, solutions);
        }
        else if (solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::SORSolver)
        {
            typedef wave_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::SORSolver>
                dev_sor_solver;
            LSS_ASSERT(!solver_config_details_.empty(), "solver_config_details map must not be empty");
            auto const omega_value = solver_config_details_["sor_omega"];
            dev_sor_solver solver(boundary_pair, wave_data_trans_cfg_, discretization_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol_0, prev_sol_1, next_sol, is_wave_source_set, wave_source, omega_value, solutions);
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
            typedef wave_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::CUDASolver>
                host_cu_solver;
            host_cu_solver solver(boundary_pair, wave_data_trans_cfg_, discretization_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol_0, prev_sol_1, next_sol, is_wave_source_set, wave_source, solutions);
        }
        else if (solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::SORSolver)
        {
            typedef wave_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::SORSolver>
                host_sor_solver;

            LSS_ASSERT(!solver_config_details_.empty(), "solver_config_details map must not be empty");
            auto const omega_value = solver_config_details_["sor_omega"];
            host_sor_solver solver(boundary_pair, wave_data_trans_cfg_, discretization_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol_0, prev_sol_1, next_sol, is_wave_source_set, wave_source, omega_value, solutions);
        }
        else if (solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::DoubleSweepSolver)
        {
            typedef wave_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::DoubleSweepSolver>
                host_dss_solver;
            host_dss_solver solver(boundary_pair, wave_data_trans_cfg_, discretization_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol_0, prev_sol_1, next_sol, is_wave_source_set, wave_source, solutions);
        }
        else if (solver_cfg_->tridiagonal_method() == tridiagonal_method_enum::ThomasLUSolver)
        {
            typedef wave_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::ThomasLUSolver>
                host_lus_solver;
            host_lus_solver solver(boundary_pair, wave_data_trans_cfg_, discretization_cfg_, solver_cfg_, grid_cfg);
            solver(prev_sol_0, prev_sol_1, next_sol, is_wave_source_set, wave_source, solutions);
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

} // namespace implicit_solvers

namespace explicit_solvers
{
using d_1d = discretization_1d;

void wave_equation::initialize(wave_data_config_1d_ptr const &wave_data_cfg,
                               grid_config_hints_1d_ptr const &grid_config_hints, boundary_1d_pair const &boundary_pair)
{
    LSS_VERIFY(wave_data_cfg, "wave_data_config must not be null");
    LSS_VERIFY(discretization_cfg_, "discretization_config must not be null");
    LSS_VERIFY(std::get<0>(boundary_pair), "boundary_pair.first must not be null");
    LSS_VERIFY(std::get<1>(boundary_pair), "boundary_pair.second must not be null");
    LSS_VERIFY(solver_cfg_, "solver_config must not be null");

    // make necessary transformations:
    // create grid_transform_config:
    grid_trans_cfg_ = std::make_shared<grid_transform_config_1d>(discretization_cfg_, grid_config_hints);
    // transform original wave data:
    wave_data_trans_cfg_ = std::make_shared<wave_data_transform_1d>(wave_data_cfg, grid_trans_cfg_);
    // transform original boundary:
    boundary_ = std::make_shared<boundary_transform_1d>(boundary_pair, grid_trans_cfg_);
}

wave_equation::wave_equation(wave_data_config_1d_ptr const &wave_data_config,
                             pde_discretization_config_1d_ptr const &discretization_config,
                             boundary_1d_pair const &boundary_pair, grid_config_hints_1d_ptr const &grid_config_hints,
                             wave_explicit_solver_config_ptr const &solver_config)
    : discretization_cfg_{discretization_config}, solver_cfg_{solver_config}
{
    initialize(wave_data_config, grid_config_hints, boundary_pair);
}

wave_equation::~wave_equation()
{
}

void wave_equation::solve(container_t &solution)
{
    LSS_ASSERT(solution.size() > 0, "The input solution container must be initialized");
    // size of space discretization:
    const std::size_t space_size = discretization_cfg_->number_of_space_points();
    // This is the proper size of the container:
    LSS_ASSERT((solution.size() == space_size), "The input solution container must have the correct size");
    // grid:
    auto const &grid_cfg = std::make_shared<grid_config_1d>(discretization_cfg_);
    auto const &boundary_pair = boundary_->boundary_pair();
    // create container to carry previous solution:
    container_t prev_sol_0(double{}, space_size);
    // discretize initial condition
    d_1d::of_function(grid_cfg, wave_data_trans_cfg_->first_initial_condition(), prev_sol_0);
    // create container to carry second initial condition:
    container_t prev_sol_1(double{}, space_size);
    // discretize first order initial condition:
    d_1d::of_function(grid_cfg, wave_data_trans_cfg_->second_initial_condition(), prev_sol_1);
    // get wave_source:
    auto const &wave_source = wave_data_trans_cfg_->wave_source();
    // is wave_source set:
    const bool is_wave_source_set = wave_data_trans_cfg_->is_wave_source_set();
    // create container to carry new solution:
    container_t next_sol(double{}, space_size);

    if (solver_cfg_->memory_space() == memory_space_enum::Device)
    {
        typedef wave_equation_explicit_kernel<memory_space_enum::Device> device_solver;
        device_solver solver(boundary_pair, wave_data_trans_cfg_, discretization_cfg_, solver_cfg_, grid_cfg);
        solver(prev_sol_0, prev_sol_1, next_sol, is_wave_source_set, wave_source);
        valcopy(solution, next_sol);
    }
    else if (solver_cfg_->memory_space() == memory_space_enum::Host)
    {
        typedef wave_equation_explicit_kernel<memory_space_enum::Host> host_solver;
        host_solver solver(boundary_pair, wave_data_trans_cfg_, discretization_cfg_, solver_cfg_, grid_cfg);
        solver(prev_sol_0, prev_sol_1, next_sol, is_wave_source_set, wave_source);
        valcopy(solution, prev_sol_1);
    }
    else
    {
        throw std::exception("Unreachable");
    }
}

void wave_equation::solve(matrix_2d &solutions)
{
    LSS_ASSERT(((solutions.columns() > 0) && (solutions.rows() > 0)),
               "The input solution 2D container must be initialized");
    // size of space discretization:
    const std::size_t space_size = discretization_cfg_->number_of_space_points();
    // size of time discretization:
    const std::size_t time_size = discretization_cfg_->number_of_time_points();
    // This is the proper size of the container:
    LSS_ASSERT((solutions.rows() == time_size) && (solutions.columns() == space_size),
               "The input solution 2D container must have the correct size");
    // grid:
    auto const &grid_cfg = std::make_shared<grid_config_1d>(discretization_cfg_);
    auto const &boundary_pair = boundary_->boundary_pair();
    // create container to carry previous solution:
    container_t prev_sol_0(double{}, space_size);
    // discretize initial condition
    d_1d::of_function(grid_cfg, wave_data_trans_cfg_->first_initial_condition(), prev_sol_0);
    // create container to carry second initial condition:
    container_t prev_sol_1(double{}, space_size);
    // discretize first order initial condition:
    d_1d::of_function(grid_cfg, wave_data_trans_cfg_->second_initial_condition(), prev_sol_1);
    // get wave_source:
    auto const &wave_source = wave_data_trans_cfg_->wave_source();
    // is wave_source set:
    const bool is_wave_source_set = wave_data_trans_cfg_->is_wave_source_set();
    // create container to carry new solution:
    container_t next_sol(double{}, space_size);

    if (solver_cfg_->memory_space() == memory_space_enum::Device)
    {
        typedef wave_equation_explicit_kernel<memory_space_enum::Device> device_solver;
        device_solver solver(boundary_pair, wave_data_trans_cfg_, discretization_cfg_, solver_cfg_, grid_cfg);
        solver(prev_sol_0, prev_sol_1, next_sol, is_wave_source_set, wave_source, solutions);
    }
    else if (solver_cfg_->memory_space() == memory_space_enum::Host)
    {
        typedef wave_equation_explicit_kernel<memory_space_enum::Host> host_solver;
        host_solver solver(boundary_pair, wave_data_trans_cfg_, discretization_cfg_, solver_cfg_, grid_cfg);
        solver(prev_sol_0, prev_sol_1, next_sol, is_wave_source_set, wave_source, solutions);
    }
    else
    {
        throw std::exception("Unreachable");
    }
}

} // namespace explicit_solvers

} // namespace one_dimensional

} // namespace lss_pde_solvers
