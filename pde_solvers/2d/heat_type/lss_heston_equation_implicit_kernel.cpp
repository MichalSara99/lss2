#include "lss_heston_equation_implicit_kernel.hpp"

#include "../../../sparse_solvers/tridiagonal/cuda_solver/lss_cuda_solver.hpp"
#include "../../../sparse_solvers/tridiagonal/double_sweep_solver/lss_double_sweep_solver.hpp"
#include "../../../sparse_solvers/tridiagonal/sor_solver/lss_sor_solver.hpp"
#include "../../../sparse_solvers/tridiagonal/sor_solver_cuda/lss_sor_solver_cuda.hpp"
#include "../../../sparse_solvers/tridiagonal/thomas_lu_solver/lss_thomas_lu_solver.hpp"
#include "boundary_solver/lss_heston_boundary_solver.hpp"
#include "implicit_coefficients/lss_heat_coefficients_2d.hpp"
#include "solver_method/splitting/lss_heat_craig_sneyd_method.hpp"
#include "solver_method/splitting/lss_heat_douglas_rachford_method.hpp"
#include "solver_method/splitting/lss_heat_hundsdorfer_verwer_method.hpp"
#include "solver_method/splitting/lss_heat_modified_craig_sneyd_method.hpp"
#include "solver_method/splitting/lss_heat_splitting_method.hpp"
#include "time_loop/lss_implicit_time_loop_2d.hpp"

namespace lss_pde_solvers
{

using lss_cuda_solver::cuda_solver;
using lss_double_sweep_solver::double_sweep_solver;
using lss_enumerations::traverse_direction_enum;
using lss_sor_solver::sor_solver;
using lss_sor_solver_cuda::sor_solver_cuda;
using lss_thomas_lu_solver::thomas_lu_solver;

namespace two_dimensional
{

heston_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::CUDASolver>::
    heston_equation_implicit_kernel(boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                    boundary_2d_pair const &horizontal_boundary_pair,
                                    heat_data_transform_2d_ptr const &heat_data_config,
                                    pde_discretization_config_2d_ptr const &discretization_config,
                                    splitting_method_config_ptr const &splitting_config,
                                    heat_implicit_solver_config_ptr const &solver_config,
                                    grid_config_2d_ptr const &grid_config)
    : boundary_ver_{vertical_upper_boundary_ptr}, boundary_pair_hor_{horizontal_boundary_pair},
      heat_data_cfg_{heat_data_config}, discretization_cfg_{discretization_config}, splitting_cfg_{splitting_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void heston_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::CUDASolver>::operator()(
    matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
    std::function<double(double, double, double)> const &heat_source)
{
    // get time range:
    auto const &time = discretization_cfg_->time_range();
    // time step:
    const double k = discretization_cfg_->time_step();
    // size of spaces discretization:
    const auto &space_sizes = discretization_cfg_->number_of_space_points();
    const std::size_t space_size_x = std::get<0>(space_sizes);
    const std::size_t space_size_y = std::get<1>(space_sizes);
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();
    // create a Heston coefficient holder:
    auto const heston_coeff_holder = std::make_shared<heat_coefficients_2d>(
        heat_data_cfg_, discretization_cfg_, splitting_cfg_, solver_cfg_->implicit_pde_scheme_value());
    heat_splitting_method_ptr splitting_ptr;
    auto solver_y = std::make_shared<cuda_solver<memory_space_enum::Device>>(space_size_x);
    solver_y->set_factorization(solver_cfg_->tridiagonal_factorization());
    auto solver_u = std::make_shared<cuda_solver<memory_space_enum::Device>>(space_size_y);
    solver_u->set_factorization(solver_cfg_->tridiagonal_factorization());
    // splitting method:
    if (splitting_cfg_->splitting_method() == splitting_method_enum::DouglasRachford)
    {
        // create and set up the main solvers:
        splitting_ptr = std::make_shared<heat_douglas_rachford_method>(solver_y, solver_u, heston_coeff_holder,
                                                                       grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::CraigSneyd)
    {
        splitting_ptr = std::make_shared<heat_craig_sneyd_method>(solver_y, solver_u, heston_coeff_holder, grid_cfg_,
                                                                  is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::ModifiedCraigSneyd)
    {
        splitting_ptr = std::make_shared<heat_modified_craig_sneyd_method>(solver_y, solver_u, heston_coeff_holder,
                                                                           grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::HundsdorferVerwer)
    {
        splitting_ptr = std::make_shared<heat_hundsdorfer_verwer_method>(solver_y, solver_u, heston_coeff_holder,
                                                                         grid_cfg_, is_heat_sourse_set);
    }
    else
    {
        throw std::exception("Unreachable");
    }
    // create and set up lower volatility boundary solver:
    auto boundary_solver = std::make_shared<heston_boundary_solver>(heston_coeff_holder, grid_cfg_);

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
        implicit_time_loop_2d::run(splitting_ptr, boundary_solver, boundary_pair_hor_, boundary_ver_, grid_cfg_, time,
                                   last_time_idx, k, traverse_dir, prev_solution, next_solution);
    }
}

heston_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::SORSolver>::
    heston_equation_implicit_kernel(boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                    boundary_2d_pair const &horizontal_boundary_pair,
                                    heat_data_transform_2d_ptr const &heat_data_config,
                                    pde_discretization_config_2d_ptr const &discretization_config,
                                    splitting_method_config_ptr const &splitting_config,
                                    heat_implicit_solver_config_ptr const &solver_config,
                                    grid_config_2d_ptr const &grid_config)
    : boundary_ver_{vertical_upper_boundary_ptr}, boundary_pair_hor_{horizontal_boundary_pair},
      heat_data_cfg_{heat_data_config}, discretization_cfg_{discretization_config}, splitting_cfg_{splitting_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void heston_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::SORSolver>::operator()(
    matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
    std::function<double(double, double, double)> const &heat_source, double omega_value)
{
    // get time range:
    auto const &time = discretization_cfg_->time_range();
    // time step:
    const double k = discretization_cfg_->time_step();
    // size of spaces discretization:
    const auto &space_sizes = discretization_cfg_->number_of_space_points();
    const std::size_t space_size_x = std::get<0>(space_sizes);
    const std::size_t space_size_y = std::get<1>(space_sizes);
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();

    // create a Heston coefficient holder:
    auto const heston_coeff_holder = std::make_shared<heat_coefficients_2d>(
        heat_data_cfg_, discretization_cfg_, splitting_cfg_, solver_cfg_->implicit_pde_scheme_value());
    heat_splitting_method_ptr splitting_ptr;
    auto solver_y = std::make_shared<sor_solver_cuda>(space_size_x);
    solver_y->set_omega(omega_value);
    auto solver_u = std::make_shared<sor_solver_cuda>(space_size_y);
    solver_u->set_omega(omega_value);
    // splitting method:
    if (splitting_cfg_->splitting_method() == splitting_method_enum::DouglasRachford)
    {
        // create and set up the main solvers:
        splitting_ptr = std::make_shared<heat_douglas_rachford_method>(solver_y, solver_u, heston_coeff_holder,
                                                                       grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::CraigSneyd)
    {
        splitting_ptr = std::make_shared<heat_craig_sneyd_method>(solver_y, solver_u, heston_coeff_holder, grid_cfg_,
                                                                  is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::ModifiedCraigSneyd)
    {
        splitting_ptr = std::make_shared<heat_modified_craig_sneyd_method>(solver_y, solver_u, heston_coeff_holder,
                                                                           grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::HundsdorferVerwer)
    {
        splitting_ptr = std::make_shared<heat_hundsdorfer_verwer_method>(solver_y, solver_u, heston_coeff_holder,
                                                                         grid_cfg_, is_heat_sourse_set);
    }
    else
    {
        throw std::exception("Unreachable");
    }
    // create and set up lower volatility boundary solver:
    auto boundary_solver = std::make_shared<heston_boundary_solver>(heston_coeff_holder, grid_cfg_);

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
        implicit_time_loop_2d::run(splitting_ptr, boundary_solver, boundary_pair_hor_, boundary_ver_, grid_cfg_, time,
                                   last_time_idx, k, traverse_dir, prev_solution, next_solution);
    }
}

heston_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::CUDASolver>::
    heston_equation_implicit_kernel(boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                    boundary_2d_pair const &horizontal_boundary_pair,
                                    heat_data_transform_2d_ptr const &heat_data_config,
                                    pde_discretization_config_2d_ptr const &discretization_config,
                                    splitting_method_config_ptr const &splitting_config,
                                    heat_implicit_solver_config_ptr const &solver_config,
                                    grid_config_2d_ptr const &grid_config)
    : boundary_ver_{vertical_upper_boundary_ptr}, boundary_pair_hor_{horizontal_boundary_pair},
      heat_data_cfg_{heat_data_config}, discretization_cfg_{discretization_config}, splitting_cfg_{splitting_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void heston_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::CUDASolver>::operator()(
    matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
    std::function<double(double, double, double)> const &heat_source)
{
    // get time range:
    auto const &time = discretization_cfg_->time_range();
    // time step:
    const double k = discretization_cfg_->time_step();
    // size of spaces discretization:
    const auto &space_sizes = discretization_cfg_->number_of_space_points();
    const std::size_t space_size_x = std::get<0>(space_sizes);
    const std::size_t space_size_y = std::get<1>(space_sizes);
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();
    // create a Heston coefficient holder:
    auto const heston_coeff_holder = std::make_shared<heat_coefficients_2d>(
        heat_data_cfg_, discretization_cfg_, splitting_cfg_, solver_cfg_->implicit_pde_scheme_value());
    heat_splitting_method_ptr splitting_ptr;
    auto solver_y = std::make_shared<cuda_solver<memory_space_enum::Host>>(space_size_x);
    solver_y->set_factorization(solver_cfg_->tridiagonal_factorization());
    auto solver_u = std::make_shared<cuda_solver<memory_space_enum::Host>>(space_size_y);
    solver_u->set_factorization(solver_cfg_->tridiagonal_factorization());
    // splitting method:
    if (splitting_cfg_->splitting_method() == splitting_method_enum::DouglasRachford)
    {
        // create and set up the main solvers:
        splitting_ptr = std::make_shared<heat_douglas_rachford_method>(solver_y, solver_u, heston_coeff_holder,
                                                                       grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::CraigSneyd)
    {
        splitting_ptr = std::make_shared<heat_craig_sneyd_method>(solver_y, solver_u, heston_coeff_holder, grid_cfg_,
                                                                  is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::ModifiedCraigSneyd)
    {
        splitting_ptr = std::make_shared<heat_modified_craig_sneyd_method>(solver_y, solver_u, heston_coeff_holder,
                                                                           grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::HundsdorferVerwer)
    {
        splitting_ptr = std::make_shared<heat_hundsdorfer_verwer_method>(solver_y, solver_u, heston_coeff_holder,
                                                                         grid_cfg_, is_heat_sourse_set);
    }
    else
    {
        throw std::exception("Unreachable");
    }
    // create and set up lower volatility boundary solver:
    auto boundary_solver = std::make_shared<heston_boundary_solver>(heston_coeff_holder, grid_cfg_);

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
        implicit_time_loop_2d::run(splitting_ptr, boundary_solver, boundary_pair_hor_, boundary_ver_, grid_cfg_, time,
                                   last_time_idx, k, traverse_dir, prev_solution, next_solution);
    }
}

heston_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::SORSolver>::
    heston_equation_implicit_kernel(boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                    boundary_2d_pair const &horizontal_boundary_pair,
                                    heat_data_transform_2d_ptr const &heat_data_config,
                                    pde_discretization_config_2d_ptr const &discretization_config,
                                    splitting_method_config_ptr const &splitting_config,
                                    heat_implicit_solver_config_ptr const &solver_config,
                                    grid_config_2d_ptr const &grid_config)
    : boundary_ver_{vertical_upper_boundary_ptr}, boundary_pair_hor_{horizontal_boundary_pair},
      heat_data_cfg_{heat_data_config}, discretization_cfg_{discretization_config}, splitting_cfg_{splitting_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void heston_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::SORSolver>::operator()(
    matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
    std::function<double(double, double, double)> const &heat_source, double omega_value)
{
    // get time range:
    auto const &time = discretization_cfg_->time_range();
    // time step:
    const double k = discretization_cfg_->time_step();
    // size of spaces discretization:
    const auto &space_sizes = discretization_cfg_->number_of_space_points();
    const std::size_t space_size_x = std::get<0>(space_sizes);
    const std::size_t space_size_y = std::get<1>(space_sizes);
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();
    // create a Heston coefficient holder:
    auto const heston_coeff_holder = std::make_shared<heat_coefficients_2d>(
        heat_data_cfg_, discretization_cfg_, splitting_cfg_, solver_cfg_->implicit_pde_scheme_value());
    heat_splitting_method_ptr splitting_ptr;
    auto solver_y = std::make_shared<sor_solver>(space_size_x);
    solver_y->set_omega(omega_value);
    auto solver_u = std::make_shared<sor_solver>(space_size_y);
    solver_u->set_omega(omega_value);
    // splitting method:
    if (splitting_cfg_->splitting_method() == splitting_method_enum::DouglasRachford)
    {
        // create and set up the main solvers:
        splitting_ptr = std::make_shared<heat_douglas_rachford_method>(solver_y, solver_u, heston_coeff_holder,
                                                                       grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::CraigSneyd)
    {
        splitting_ptr = std::make_shared<heat_craig_sneyd_method>(solver_y, solver_u, heston_coeff_holder, grid_cfg_,
                                                                  is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::ModifiedCraigSneyd)
    {
        splitting_ptr = std::make_shared<heat_modified_craig_sneyd_method>(solver_y, solver_u, heston_coeff_holder,
                                                                           grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::HundsdorferVerwer)
    {
        splitting_ptr = std::make_shared<heat_hundsdorfer_verwer_method>(solver_y, solver_u, heston_coeff_holder,
                                                                         grid_cfg_, is_heat_sourse_set);
    }
    else
    {
        throw std::exception("Unreachable");
    }
    // create and set up lower volatility boundary solver:
    auto boundary_solver = std::make_shared<heston_boundary_solver>(heston_coeff_holder, grid_cfg_);

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
        implicit_time_loop_2d::run(splitting_ptr, boundary_solver, boundary_pair_hor_, boundary_ver_, grid_cfg_, time,
                                   last_time_idx, k, traverse_dir, prev_solution, next_solution);
    }
}

heston_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::DoubleSweepSolver>::
    heston_equation_implicit_kernel(boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                    boundary_2d_pair const &horizontal_boundary_pair,
                                    heat_data_transform_2d_ptr const &heat_data_config,
                                    pde_discretization_config_2d_ptr const &discretization_config,
                                    splitting_method_config_ptr const &splitting_config,
                                    heat_implicit_solver_config_ptr const &solver_config,
                                    grid_config_2d_ptr const &grid_config)
    : boundary_ver_{vertical_upper_boundary_ptr}, boundary_pair_hor_{horizontal_boundary_pair},
      heat_data_cfg_{heat_data_config}, discretization_cfg_{discretization_config}, splitting_cfg_{splitting_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void heston_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::DoubleSweepSolver>::operator()(
    matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
    std::function<double(double, double, double)> const &heat_source)
{
    // get time range:
    auto const &time = discretization_cfg_->time_range();
    // time step:
    const double k = discretization_cfg_->time_step();
    // size of spaces discretization:
    const auto &space_sizes = discretization_cfg_->number_of_space_points();
    const std::size_t space_size_x = std::get<0>(space_sizes);
    const std::size_t space_size_y = std::get<1>(space_sizes);
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();
    // create a Heston coefficient holder:
    auto const heston_coeff_holder = std::make_shared<heat_coefficients_2d>(
        heat_data_cfg_, discretization_cfg_, splitting_cfg_, solver_cfg_->implicit_pde_scheme_value());
    heat_splitting_method_ptr splitting_ptr;
    // create and set up the main solvers:
    auto solver_y = std::make_shared<double_sweep_solver>(space_size_x);
    auto solver_u = std::make_shared<double_sweep_solver>(space_size_y);
    // splitting method:
    if (splitting_cfg_->splitting_method() == splitting_method_enum::DouglasRachford)
    {
        splitting_ptr = std::make_shared<heat_douglas_rachford_method>(solver_y, solver_u, heston_coeff_holder,
                                                                       grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::CraigSneyd)
    {
        splitting_ptr = std::make_shared<heat_craig_sneyd_method>(solver_y, solver_u, heston_coeff_holder, grid_cfg_,
                                                                  is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::ModifiedCraigSneyd)
    {
        splitting_ptr = std::make_shared<heat_modified_craig_sneyd_method>(solver_y, solver_u, heston_coeff_holder,
                                                                           grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::HundsdorferVerwer)
    {
        splitting_ptr = std::make_shared<heat_hundsdorfer_verwer_method>(solver_y, solver_u, heston_coeff_holder,
                                                                         grid_cfg_, is_heat_sourse_set);
    }
    else
    {
        throw std::exception("Unreachable");
    }
    // create and set up lower volatility boundary solver:
    auto boundary_solver = std::make_shared<heston_boundary_solver>(heston_coeff_holder, grid_cfg_);

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
        implicit_time_loop_2d::run(splitting_ptr, boundary_solver, boundary_pair_hor_, boundary_ver_, grid_cfg_, time,
                                   last_time_idx, k, traverse_dir, prev_solution, next_solution);
    }
}

heston_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::ThomasLUSolver>::
    heston_equation_implicit_kernel(boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                    boundary_2d_pair const &horizontal_boundary_pair,
                                    heat_data_transform_2d_ptr const &heat_data_config,
                                    pde_discretization_config_2d_ptr const &discretization_config,
                                    splitting_method_config_ptr const &splitting_config,
                                    heat_implicit_solver_config_ptr const &solver_config,
                                    grid_config_2d_ptr const &grid_config)
    : boundary_ver_{vertical_upper_boundary_ptr}, boundary_pair_hor_{horizontal_boundary_pair},
      heat_data_cfg_{heat_data_config}, discretization_cfg_{discretization_config}, splitting_cfg_{splitting_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void heston_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::ThomasLUSolver>::operator()(
    matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
    std::function<double(double, double, double)> const &heat_source)
{
    // get time range:
    auto const &time = discretization_cfg_->time_range();
    // time step:
    const double k = discretization_cfg_->time_step();
    // size of spaces discretization:
    const auto &space_sizes = discretization_cfg_->number_of_space_points();
    const std::size_t space_size_x = std::get<0>(space_sizes);
    const std::size_t space_size_y = std::get<1>(space_sizes);
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();
    // create a Heston coefficient holder:
    auto const heston_coeff_holder = std::make_shared<heat_coefficients_2d>(
        heat_data_cfg_, discretization_cfg_, splitting_cfg_, solver_cfg_->implicit_pde_scheme_value());
    heat_splitting_method_ptr splitting_ptr;
    // create and set up the main solvers:
    auto solver_y = std::make_shared<thomas_lu_solver>(space_size_x);
    auto solver_u = std::make_shared<thomas_lu_solver>(space_size_y);
    // splitting method:
    if (splitting_cfg_->splitting_method() == splitting_method_enum::DouglasRachford)
    {
        splitting_ptr = std::make_shared<heat_douglas_rachford_method>(solver_y, solver_u, heston_coeff_holder,
                                                                       grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::CraigSneyd)
    {
        splitting_ptr = std::make_shared<heat_craig_sneyd_method>(solver_y, solver_u, heston_coeff_holder, grid_cfg_,
                                                                  is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::ModifiedCraigSneyd)
    {
        splitting_ptr = std::make_shared<heat_modified_craig_sneyd_method>(solver_y, solver_u, heston_coeff_holder,
                                                                           grid_cfg_, is_heat_sourse_set);
    }
    else if (splitting_cfg_->splitting_method() == splitting_method_enum::HundsdorferVerwer)
    {
        splitting_ptr = std::make_shared<heat_hundsdorfer_verwer_method>(solver_y, solver_u, heston_coeff_holder,
                                                                         grid_cfg_, is_heat_sourse_set);
    }
    else
    {
        throw std::exception("Unreachable");
    }
    // create and set up lower volatility boundary solver:
    auto boundary_solver = std::make_shared<heston_boundary_solver>(heston_coeff_holder, grid_cfg_);

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
        implicit_time_loop_2d::run(splitting_ptr, boundary_solver, boundary_pair_hor_, boundary_ver_, grid_cfg_, time,
                                   last_time_idx, k, traverse_dir, prev_solution, next_solution);
    }
}

} // namespace two_dimensional

} // namespace lss_pde_solvers
