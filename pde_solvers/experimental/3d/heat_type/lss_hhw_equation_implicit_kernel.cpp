#include "lss_hhw_equation_implicit_kernel.hpp"

#include "../../../../sparse_solvers/tridiagonal/cuda_solver/lss_cuda_solver.hpp"
#include "../../../../sparse_solvers/tridiagonal/double_sweep_solver/lss_double_sweep_solver.hpp"
#include "../../../../sparse_solvers/tridiagonal/sor_solver/lss_sor_solver.hpp"
#include "../../../../sparse_solvers/tridiagonal/sor_solver_cuda/lss_sor_solver_cuda.hpp"
#include "../../../../sparse_solvers/tridiagonal/thomas_lu_solver/lss_thomas_lu_solver.hpp"
#include "boundary_solver/lss_hhw_explicit_boundary_solver.hpp"
#include "implicit_coefficients/lss_heat_coefficients_3d.hpp"
//#include "splitting_method/lss_heat_craig_sneyd_method_3d.hpp"
#include "solver_method/splitting/lss_heat_douglas_rachford_method_3d.hpp"
//#include "splitting_method/lss_heat_hundsdorfer_verwer_method_3d.hpp"
//#include "splitting_method/lss_heat_modified_craig_sneyd_method_3d.hpp"
#include "solver_method/splitting/lss_heat_splitting_method_3d.hpp"
#include "time_loop/lss_implicit_time_loop_3d.hpp"

namespace lss_pde_solvers
{

using lss_cuda_solver::cuda_solver;
using lss_double_sweep_solver::double_sweep_solver;
using lss_enumerations::implicit_pde_schemes_enum;
using lss_enumerations::traverse_direction_enum;
using lss_sor_solver::sor_solver;
using lss_sor_solver_cuda::sor_solver_cuda;
using lss_thomas_lu_solver::thomas_lu_solver;

namespace three_dimensional
{

hhw_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::CUDASolver>::
    hhw_equation_implicit_kernel(boundary_3d_pair const &x_boundary_pair, boundary_3d_ptr const &y_upper_boundary_ptr,
                                 boundary_3d_pair const &z_boundary_pair,
                                 heat_data_transform_3d_ptr const &heat_data_config,
                                 pde_discretization_config_3d_ptr const &discretization_config,
                                 splitting_method_config_ptr const &splitting_config,
                                 heat_implicit_solver_config_ptr const &solver_config,
                                 grid_config_3d_ptr const &grid_config)
    : boundary_pair_x_{x_boundary_pair}, boundary_y_{y_upper_boundary_ptr}, boundary_pair_z_{z_boundary_pair},
      heat_data_cfg_{heat_data_config}, discretization_cfg_{discretization_config}, splitting_cfg_{splitting_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void hhw_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::CUDASolver>::operator()(
    matrix_3d &prev_solution, matrix_3d &next_solution, bool is_heat_sourse_set,
    std::function<double(double, double, double, double)> const &heat_source)
{
    // get time range:
    auto const &time = discretization_cfg_->time_range();
    // time step:
    const double k = discretization_cfg_->time_step();
    // size of spaces discretization:
    const auto &space_sizes = discretization_cfg_->number_of_space_points();
    const std::size_t space_size_x = std::get<0>(space_sizes);
    const std::size_t space_size_y = std::get<1>(space_sizes);
    const std::size_t space_size_z = std::get<2>(space_sizes);
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();

    // create a Heston coefficient holder:
    auto const heston_coeff_holder = std::make_shared<heat_coefficients_3d>(
        heat_data_cfg_, discretization_cfg_, splitting_cfg_, solver_cfg_->implicit_pde_scheme_value());
    heat_splitting_method_3d_ptr splitting_ptr;
    auto solver_y1 = std::make_shared<cuda_solver<memory_space_enum::Device>>(space_size_x);
    solver_y1->set_factorization(solver_cfg_->tridiagonal_factorization());
    auto solver_y2 = std::make_shared<cuda_solver<memory_space_enum::Device>>(space_size_y);
    solver_y2->set_factorization(solver_cfg_->tridiagonal_factorization());
    auto solver_u = std::make_shared<cuda_solver<memory_space_enum::Device>>(space_size_z);
    solver_u->set_factorization(solver_cfg_->tridiagonal_factorization());
    // splitting method:
    if (splitting_cfg_->splitting_method() == splitting_method_enum::DouglasRachford)
    {
        // create and set up the main solvers:
        splitting_ptr = std::make_shared<heat_douglas_rachford_method_3d>(
            solver_y1, solver_y2, solver_u, heston_coeff_holder, grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::CraigSneyd)
    {
        // splitting_ptr = std::make_shared<heat_craig_sneyd_method_3d>(solver_y1,solver_y2, solver_u,
        // heston_coeff_holder, grid_cfg_,is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::ModifiedCraigSneyd)
    {
        // splitting_ptr = std::make_shared<heat_modified_craig_sneyd_method_3d>(solver_y1,solver_y2, solver_u,
        // heston_coeff_holder,grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::HundsdorferVerwer)
    {
        // splitting_ptr = std::make_shared<heat_hundsdorfer_verwer_method_3d>(solver_y1,solver_y2, solver_u,
        // heston_coeff_holder,grid_cfg_, is_heat_sourse_set);
    }
    else
    {
        throw std::exception("Unreachable");
    }
    // create and set up lower volatility boundary solver:
    auto boundary_solver = std::make_shared<hhw_explicit_boundary_solver>(heston_coeff_holder, grid_cfg_);
    // auto boundary_solver = std::make_shared<hhw_implicit_boundary_solver>(solver_u, heston_coeff_holder, grid_cfg_);

    if (is_heat_sourse_set)
    {
        // auto scheme_function =
        //    implicit_heat_scheme<fp_type, container,
        //    allocator>::get(solver_cfg_->implicit_pde_scheme(), false);
        //// create a container to carry discretized source heat
        // container_t source_curr(space_size, NaN<fp_type>());
        // container_t source_next(space_size, NaN<fp_type>());
        // loop::run(solver, scheme_function, boundary_pair_, fun_triplet_, space,
        // time, last_time_idx, steps,
        //          traverse_dir, prev_solution, next_solution, rhs, heat_source,
        //          source_curr, source_next);
    }
    else
    {
        implicit_time_loop_3d::run(splitting_ptr, boundary_solver, boundary_pair_x_, boundary_y_, boundary_pair_z_,
                                   grid_cfg_, time, last_time_idx, k, traverse_dir, prev_solution, next_solution);
    }
}

hhw_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::SORSolver>::
    hhw_equation_implicit_kernel(boundary_3d_pair const &x_boundary_pair, boundary_3d_ptr const &y_upper_boundary_ptr,
                                 boundary_3d_pair const &z_boundary_pair,
                                 heat_data_transform_3d_ptr const &heat_data_config,
                                 pde_discretization_config_3d_ptr const &discretization_config,
                                 splitting_method_config_ptr const &splitting_config,
                                 heat_implicit_solver_config_ptr const &solver_config,
                                 grid_config_3d_ptr const &grid_config)
    : boundary_pair_x_{x_boundary_pair}, boundary_y_{y_upper_boundary_ptr}, boundary_pair_z_{z_boundary_pair},
      heat_data_cfg_{heat_data_config}, discretization_cfg_{discretization_config}, splitting_cfg_{splitting_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void hhw_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::SORSolver>::operator()(
    matrix_3d &prev_solution, matrix_3d &next_solution, bool is_heat_sourse_set,
    std::function<double(double, double, double, double)> const &heat_source, double omega_value)
{
    // get time range:
    auto const &time = discretization_cfg_->time_range();
    // time step:
    const double k = discretization_cfg_->time_step();
    // size of spaces discretization:
    const auto &space_sizes = discretization_cfg_->number_of_space_points();
    const std::size_t space_size_x = std::get<0>(space_sizes);
    const std::size_t space_size_y = std::get<1>(space_sizes);
    const std::size_t space_size_z = std::get<2>(space_sizes);
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();

    // create a Heston coefficient holder:
    auto const heston_coeff_holder = std::make_shared<heat_coefficients_3d>(
        heat_data_cfg_, discretization_cfg_, splitting_cfg_, solver_cfg_->implicit_pde_scheme_value());
    heat_splitting_method_3d_ptr splitting_ptr;
    auto solver_y1 = std::make_shared<sor_solver_cuda>(space_size_x);
    solver_y1->set_omega(omega_value);
    auto solver_y2 = std::make_shared<sor_solver_cuda>(space_size_y);
    solver_y2->set_omega(omega_value);
    auto solver_u = std::make_shared<sor_solver_cuda>(space_size_z);
    solver_u->set_omega(omega_value);
    // splitting method:
    if (splitting_cfg_->splitting_method() == splitting_method_enum::DouglasRachford)
    {
        // create and set up the main solvers:
        splitting_ptr = std::make_shared<heat_douglas_rachford_method_3d>(
            solver_y1, solver_y2, solver_u, heston_coeff_holder, grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::CraigSneyd)
    {
        // splitting_ptr = std::make_shared<heat_craig_sneyd_method_3d>(solver_y1,solver_y2, solver_u,
        // heston_coeff_holder, grid_cfg_,is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::ModifiedCraigSneyd)
    {
        // splitting_ptr = std::make_shared<heat_modified_craig_sneyd_method_3d>(solver_y1,solver_y2, solver_u,
        // heston_coeff_holder,grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::HundsdorferVerwer)
    {
        // splitting_ptr = std::make_shared<heat_hundsdorfer_verwer_method_3d>(solver_y1,solver_y2, solver_u,
        // heston_coeff_holder, grid_cfg_, is_heat_sourse_set);
    }
    else
    {
        throw std::exception("Unreachable");
    }
    // create and set up lower volatility boundary solver:
    auto boundary_solver = std::make_shared<hhw_explicit_boundary_solver>(heston_coeff_holder, grid_cfg_);
    // auto boundary_solver = std::make_shared<hhw_implicit_boundary_solver>(solver_u, heston_coeff_holder, grid_cfg_);

    if (is_heat_sourse_set)
    {
        // auto scheme_function =
        //    implicit_heat_scheme<fp_type, container,
        //    allocator>::get(solver_cfg_->implicit_pde_scheme(), false);
        //// create a container to carry discretized source heat
        // container_t source_curr(space_size, NaN<fp_type>());
        // container_t source_next(space_size, NaN<fp_type>());
        // loop::run(solver, scheme_function, boundary_pair_, fun_triplet_,
        // space, time, last_time_idx, steps,
        //          traverse_dir, prev_solution, next_solution, rhs,
        //          heat_source, source_curr, source_next);
    }
    else
    {
        implicit_time_loop_3d::run(splitting_ptr, boundary_solver, boundary_pair_x_, boundary_y_, boundary_pair_z_,
                                   grid_cfg_, time, last_time_idx, k, traverse_dir, prev_solution, next_solution);
    }
}

hhw_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::CUDASolver>::
    hhw_equation_implicit_kernel(boundary_3d_pair const &x_boundary_pair, boundary_3d_ptr const &y_upper_boundary_ptr,
                                 boundary_3d_pair const &z_boundary_pair,
                                 heat_data_transform_3d_ptr const &heat_data_config,
                                 pde_discretization_config_3d_ptr const &discretization_config,
                                 splitting_method_config_ptr const &splitting_config,
                                 heat_implicit_solver_config_ptr const &solver_config,
                                 grid_config_3d_ptr const &grid_config)
    : boundary_pair_x_{x_boundary_pair}, boundary_y_{y_upper_boundary_ptr}, boundary_pair_z_{z_boundary_pair},
      heat_data_cfg_{heat_data_config}, discretization_cfg_{discretization_config}, splitting_cfg_{splitting_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void hhw_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::CUDASolver>::operator()(
    matrix_3d &prev_solution, matrix_3d &next_solution, bool is_heat_sourse_set,
    std::function<double(double, double, double, double)> const &heat_source)
{
    // get time range:
    auto const &time = discretization_cfg_->time_range();
    // time step:
    const double k = discretization_cfg_->time_step();
    // size of spaces discretization:
    const auto &space_sizes = discretization_cfg_->number_of_space_points();
    const std::size_t space_size_x = std::get<0>(space_sizes);
    const std::size_t space_size_y = std::get<1>(space_sizes);
    const std::size_t space_size_z = std::get<2>(space_sizes);
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();
    // create a Heston coefficient holder:
    auto const heston_coeff_holder = std::make_shared<heat_coefficients_3d>(
        heat_data_cfg_, discretization_cfg_, splitting_cfg_, solver_cfg_->implicit_pde_scheme_value());
    heat_splitting_method_3d_ptr splitting_ptr;
    auto solver_y1 = std::make_shared<cuda_solver<memory_space_enum::Host>>(space_size_x);
    solver_y1->set_factorization(solver_cfg_->tridiagonal_factorization());
    auto solver_y2 = std::make_shared<cuda_solver<memory_space_enum::Host>>(space_size_y);
    solver_y2->set_factorization(solver_cfg_->tridiagonal_factorization());
    auto solver_u = std::make_shared<cuda_solver<memory_space_enum::Host>>(space_size_z);
    solver_u->set_factorization(solver_cfg_->tridiagonal_factorization());
    // splitting method:
    if (splitting_cfg_->splitting_method() == splitting_method_enum::DouglasRachford)
    {
        // create and set up the main solvers:
        splitting_ptr = std::make_shared<heat_douglas_rachford_method_3d>(
            solver_y1, solver_y2, solver_u, heston_coeff_holder, grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::CraigSneyd)
    {
        // splitting_ptr = std::make_shared<heat_craig_sneyd_method_3d>(solver_y1,solver_y2, solver_u,
        // heston_coeff_holder, grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::ModifiedCraigSneyd)
    {
        // splitting_ptr = std::make_shared<heat_modified_craig_sneyd_method_3d>(solver_y1,solver_y2, solver_u,
        // heston_coeff_holder,grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::HundsdorferVerwer)
    {
        // splitting_ptr = std::make_shared<heat_hundsdorfer_verwer_method_3d>(solver_y1,solver_y2, solver_u,
        // heston_coeff_holder,grid_cfg_, is_heat_sourse_set);
    }
    else
    {
        throw std::exception("Unreachable");
    }
    // create and set up lower volatility boundary solver:
    auto boundary_solver = std::make_shared<hhw_explicit_boundary_solver>(heston_coeff_holder, grid_cfg_);
    // auto boundary_solver = std::make_shared<hhw_implicit_boundary_solver>(solver_u, heston_coeff_holder, grid_cfg_);

    if (is_heat_sourse_set)
    {
        // auto scheme_function =
        //    implicit_heat_scheme<fp_type, container,
        //    allocator>::get(solver_cfg_->implicit_pde_scheme(), false);
        //// create a container to carry discretized source heat
        // container_t source_curr(space_size, NaN<fp_type>());
        // container_t source_next(space_size, NaN<fp_type>());
        // loop::run(solver, scheme_function, boundary_pair_, fun_triplet_,
        // space, time, last_time_idx, steps,
        //          traverse_dir, prev_solution, next_solution, rhs,
        //          heat_source, source_curr, source_next);
    }
    else
    {
        implicit_time_loop_3d::run(splitting_ptr, boundary_solver, boundary_pair_x_, boundary_y_, boundary_pair_z_,
                                   grid_cfg_, time, last_time_idx, k, traverse_dir, prev_solution, next_solution);
    }
}

hhw_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::SORSolver>::hhw_equation_implicit_kernel(
    boundary_3d_pair const &x_boundary_pair, boundary_3d_ptr const &y_upper_boundary_ptr,
    boundary_3d_pair const &z_boundary_pair, heat_data_transform_3d_ptr const &heat_data_config,
    pde_discretization_config_3d_ptr const &discretization_config, splitting_method_config_ptr const &splitting_config,
    heat_implicit_solver_config_ptr const &solver_config, grid_config_3d_ptr const &grid_config)
    : boundary_pair_x_{x_boundary_pair}, boundary_y_{y_upper_boundary_ptr}, boundary_pair_z_{z_boundary_pair},
      heat_data_cfg_{heat_data_config}, discretization_cfg_{discretization_config}, splitting_cfg_{splitting_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void hhw_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::SORSolver>::operator()(
    matrix_3d &prev_solution, matrix_3d &next_solution, bool is_heat_sourse_set,
    std::function<double(double, double, double, double)> const &heat_source, double omega_value)
{
    // get time range:
    auto const &time = discretization_cfg_->time_range();
    // time step:
    const double k = discretization_cfg_->time_step();
    // size of spaces discretization:
    const auto &space_sizes = discretization_cfg_->number_of_space_points();
    const std::size_t space_size_x = std::get<0>(space_sizes);
    const std::size_t space_size_y = std::get<1>(space_sizes);
    const std::size_t space_size_z = std::get<2>(space_sizes);
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();

    // create a Heston coefficient holder:
    auto const heston_coeff_holder = std::make_shared<heat_coefficients_3d>(
        heat_data_cfg_, discretization_cfg_, splitting_cfg_, solver_cfg_->implicit_pde_scheme_value());
    heat_splitting_method_3d_ptr splitting_ptr;
    auto solver_y1 = std::make_shared<sor_solver>(space_size_x);
    solver_y1->set_omega(omega_value);
    auto solver_y2 = std::make_shared<sor_solver>(space_size_y);
    solver_y2->set_omega(omega_value);
    auto solver_u = std::make_shared<sor_solver>(space_size_z);
    solver_u->set_omega(omega_value);
    // splitting method:
    if (splitting_cfg_->splitting_method() == splitting_method_enum::DouglasRachford)
    {
        // create and set up the main solvers:
        splitting_ptr = std::make_shared<heat_douglas_rachford_method_3d>(
            solver_y1, solver_y2, solver_u, heston_coeff_holder, grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::CraigSneyd)
    {
        // splitting_ptr = std::make_shared<heat_craig_sneyd_method_3d>(solver_y1,solver_y2, solver_u,
        // heston_coeff_holder, grid_cfg_,is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::ModifiedCraigSneyd)
    {
        // splitting_ptr = std::make_shared<heat_modified_craig_sneyd_method_3d>(solver_y1,solver_y2, solver_u,
        // heston_coeff_holder,grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::HundsdorferVerwer)
    {
        // splitting_ptr = std::make_shared<heat_hundsdorfer_verwer_method_3d>(solver_y1,solver_y2, solver_u,
        // heston_coeff_holder,grid_cfg_, is_heat_sourse_set);
    }
    else
    {
        throw std::exception("Unreachable");
    }
    // create and set up lower volatility boundary solver:
    auto boundary_solver = std::make_shared<hhw_explicit_boundary_solver>(heston_coeff_holder, grid_cfg_);
    // auto boundary_solver = std::make_shared<hhw_implicit_boundary_solver>(solver_u, heston_coeff_holder, grid_cfg_);

    if (is_heat_sourse_set)
    {
        // auto scheme_function =
        //    implicit_heat_scheme<fp_type, container,
        //    allocator>::get(solver_cfg_->implicit_pde_scheme(), false);
        //// create a container to carry discretized source heat
        // container_t source_curr(space_size, NaN<fp_type>());
        // container_t source_next(space_size, NaN<fp_type>());
        // loop::run(solver, scheme_function, boundary_pair_, fun_triplet_,
        // space, time, last_time_idx, steps,
        //          traverse_dir, prev_solution, next_solution, rhs,
        //          heat_source, source_curr, source_next);
    }
    else
    {
        implicit_time_loop_3d::run(splitting_ptr, boundary_solver, boundary_pair_x_, boundary_y_, boundary_pair_z_,
                                   grid_cfg_, time, last_time_idx, k, traverse_dir, prev_solution, next_solution);
    }
}

hhw_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::DoubleSweepSolver>::
    hhw_equation_implicit_kernel(boundary_3d_pair const &x_boundary_pair, boundary_3d_ptr const &y_upper_boundary_ptr,
                                 boundary_3d_pair const &z_boundary_pair,
                                 heat_data_transform_3d_ptr const &heat_data_config,
                                 pde_discretization_config_3d_ptr const &discretization_config,
                                 splitting_method_config_ptr const &splitting_config,
                                 heat_implicit_solver_config_ptr const &solver_config,
                                 grid_config_3d_ptr const &grid_config)
    : boundary_pair_x_{x_boundary_pair}, boundary_y_{y_upper_boundary_ptr}, boundary_pair_z_{z_boundary_pair},
      heat_data_cfg_{heat_data_config}, discretization_cfg_{discretization_config}, splitting_cfg_{splitting_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void hhw_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::DoubleSweepSolver>::operator()(
    matrix_3d &prev_solution, matrix_3d &next_solution, bool is_heat_sourse_set,
    std::function<double(double, double, double, double)> const &heat_source)
{
    // get time range:
    auto const &time = discretization_cfg_->time_range();
    // time step:
    const double k = discretization_cfg_->time_step();
    // size of spaces discretization:
    const auto &space_sizes = discretization_cfg_->number_of_space_points();
    const std::size_t space_size_x = std::get<0>(space_sizes);
    const std::size_t space_size_y = std::get<1>(space_sizes);
    const std::size_t space_size_z = std::get<2>(space_sizes);
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();
    // create a Heston coefficient holder:
    auto const heston_coeff_holder = std::make_shared<heat_coefficients_3d>(
        heat_data_cfg_, discretization_cfg_, splitting_cfg_, solver_cfg_->implicit_pde_scheme_value());
    heat_splitting_method_3d_ptr splitting_ptr;
    // create and set up the main solvers:
    auto solver_y1 = std::make_shared<double_sweep_solver>(space_size_x);
    auto solver_y2 = std::make_shared<double_sweep_solver>(space_size_y);
    auto solver_u = std::make_shared<double_sweep_solver>(space_size_z);
    // splitting method:
    if (splitting_cfg_->splitting_method() == splitting_method_enum::DouglasRachford)
    {
        splitting_ptr = std::make_shared<heat_douglas_rachford_method_3d>(
            solver_y1, solver_y2, solver_u, heston_coeff_holder, grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::CraigSneyd)
    {
        // splitting_ptr = std::make_shared<heat_craig_sneyd_method_3d>(
        //     solver_y1, solver_y2, solver_u, heston_coeff_holder, grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::ModifiedCraigSneyd)
    {
        // splitting_ptr = std::make_shared<heat_modified_craig_sneyd_method_3d>(
        //     solver_y1, solver_y2, solver_u, heston_coeff_holder, grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::HundsdorferVerwer)
    {
        // splitting_ptr = std::make_shared<heat_hundsdorfer_verwer_method_3d>(
        //     solver_y1, solver_y2, solver_u, heston_coeff_holder, grid_cfg_, is_heat_sourse_set);
    }
    else
    {
        throw std::exception("Unreachable");
    }
    // create and set up lower volatility boundary solver:
    auto boundary_solver = std::make_shared<hhw_explicit_boundary_solver>(heston_coeff_holder, grid_cfg_);
    // auto boundary_solver = std::make_shared<hhw_implicit_boundary_solver>(solver_u, heston_coeff_holder, grid_cfg_);

    if (is_heat_sourse_set)
    {
        // auto scheme_function =
        //    implicit_heat_scheme<fp_type, container,
        //    allocator>::get(solver_cfg_->implicit_pde_scheme(), false);
        //// create a container to carry discretized source heat
        // container_t source_curr(space_size, NaN<fp_type>());
        // container_t source_next(space_size, NaN<fp_type>());
        // loop::run(solver, scheme_function, boundary_pair_, fun_triplet_,
        // space, time, last_time_idx, steps,
        //          traverse_dir, prev_solution, next_solution, rhs,
        //          heat_source, source_curr, source_next);
    }
    else
    {
        implicit_time_loop_3d::run(splitting_ptr, boundary_solver, boundary_pair_x_, boundary_y_, boundary_pair_z_,
                                   grid_cfg_, time, last_time_idx, k, traverse_dir, prev_solution, next_solution);
    }
}

hhw_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::ThomasLUSolver>::
    hhw_equation_implicit_kernel(boundary_3d_pair const &x_boundary_pair, boundary_3d_ptr const &y_upper_boundary_ptr,
                                 boundary_3d_pair const &z_boundary_pair,
                                 heat_data_transform_3d_ptr const &heat_data_config,
                                 pde_discretization_config_3d_ptr const &discretization_config,
                                 splitting_method_config_ptr const &splitting_config,
                                 heat_implicit_solver_config_ptr const &solver_config,
                                 grid_config_3d_ptr const &grid_config)
    : boundary_pair_x_{x_boundary_pair}, boundary_y_{y_upper_boundary_ptr}, boundary_pair_z_{z_boundary_pair},
      heat_data_cfg_{heat_data_config}, discretization_cfg_{discretization_config}, splitting_cfg_{splitting_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void hhw_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::ThomasLUSolver>::operator()(
    matrix_3d &prev_solution, matrix_3d &next_solution, bool is_heat_sourse_set,
    std::function<double(double, double, double, double)> const &heat_source)
{
    // get time range:
    auto const &time = discretization_cfg_->time_range();
    // time step:
    const double k = discretization_cfg_->time_step();
    // size of spaces discretization:
    const auto &space_sizes = discretization_cfg_->number_of_space_points();
    const std::size_t space_size_x = std::get<0>(space_sizes);
    const std::size_t space_size_y = std::get<1>(space_sizes);
    const std::size_t space_size_z = std::get<2>(space_sizes);
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();

    // create a Heston coefficient holder:
    auto const heston_coeff_holder = std::make_shared<heat_coefficients_3d>(
        heat_data_cfg_, discretization_cfg_, splitting_cfg_, solver_cfg_->implicit_pde_scheme_value());
    heat_splitting_method_3d_ptr splitting_ptr;
    // create and set up the main solvers:
    auto solver_y1 = std::make_shared<thomas_lu_solver>(space_size_x);
    auto solver_y2 = std::make_shared<thomas_lu_solver>(space_size_y);
    auto solver_u = std::make_shared<thomas_lu_solver>(space_size_z);
    // splitting method:
    if (splitting_cfg_->splitting_method() == splitting_method_enum::DouglasRachford)
    {
        splitting_ptr = std::make_shared<heat_douglas_rachford_method_3d>(
            solver_y1, solver_y2, solver_u, heston_coeff_holder, grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::CraigSneyd)
    {
        // splitting_ptr = std::make_shared<heat_craig_sneyd_method_3d>(
        //     solver_y1, solver_y2, solver_u, heston_coeff_holder, grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::ModifiedCraigSneyd)
    {
        // splitting_ptr = std::make_shared<heat_modified_craig_sneyd_method_3d>(
        //     solver_y1, solver_y2, solver_u, heston_coeff_holder, grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::HundsdorferVerwer)
    {
        // splitting_ptr = std::make_shared<heat_hundsdorfer_verwer_method_3d>(
        //     solver_y1, solver_y2, solver_u, heston_coeff_holder, grid_cfg_, is_heat_sourse_set);
    }
    else
    {
        throw std::exception("Unreachable");
    }
    // create and set up lower volatility boundary solver:
    auto boundary_solver = std::make_shared<hhw_explicit_boundary_solver>(heston_coeff_holder, grid_cfg_);
    // auto boundary_solver = std::make_shared<hhw_implicit_boundary_solver>(solver_u, heston_coeff_holder, grid_cfg_);

    if (is_heat_sourse_set)
    {
        // auto scheme_function =
        //    implicit_heat_scheme<fp_type, container,
        //    allocator>::get(solver_cfg_->implicit_pde_scheme(), false);
        //// create a container to carry discretized source heat
        // container_t source_curr(space_size, NaN<fp_type>());
        // container_t source_next(space_size, NaN<fp_type>());
        // loop::run(solver, scheme_function, boundary_pair_, fun_triplet_,
        // space, time, last_time_idx, steps,
        //          traverse_dir, prev_solution, next_solution, rhs,
        //          heat_source, source_curr, source_next);
    }
    else
    {
        implicit_time_loop_3d::run(splitting_ptr, boundary_solver, boundary_pair_x_, boundary_y_, boundary_pair_z_,
                                   grid_cfg_, time, last_time_idx, k, traverse_dir, prev_solution, next_solution);
    }
}

} // namespace three_dimensional

} // namespace lss_pde_solvers
