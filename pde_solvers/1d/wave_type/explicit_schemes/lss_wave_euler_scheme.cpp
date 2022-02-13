#include "lss_wave_euler_scheme.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

bool wave_euler_scheme::is_stable()
{
    auto const &b = euler_coeffs_->b_;
    const double k = euler_coeffs_->k_;
    const double h = grid_1d::step(grid_cfg_);
    const double ratio = h / k;
    const std::size_t space_size = discretization_cfg_->number_of_space_points();
    auto const ftime = discretization_cfg_->time_range()->upper();
    double x{}, t{k};

    while (t <= ftime)
    {
        for (std::size_t i = 0; i < space_size; ++i)
        {
            x = grid_1d::value(grid_cfg_, i);
            if (b(t, x) >= ratio)
                return false;
        }
        t += k;
    }
    return true;
}

void wave_euler_scheme::initialize()
{
    LSS_ASSERT(is_stable() == true, "The chosen scheme is not stable");
}

wave_euler_scheme::wave_euler_scheme(wave_explicit_coefficients_ptr const &coefficients,
                                     boundary_1d_pair const &boundary_pair,
                                     pde_discretization_config_1d_ptr const &discretization_config,
                                     grid_config_1d_ptr const &grid_config)
    : euler_coeffs_{coefficients}, boundary_pair_{boundary_pair},
      discretization_cfg_{discretization_config}, grid_cfg_{grid_config}
{
    initialize();
}

wave_euler_scheme::~wave_euler_scheme()
{
}

void wave_euler_scheme::operator()(container_t &prev_solution_0, container_t &prev_solution_1,
                                   container_t &next_solution, bool is_wave_sourse_set,
                                   std::function<double(double, double)> const &wave_source,
                                   traverse_direction_enum traverse_dir)
{
    auto const &timer = discretization_cfg_->time_range();
    const double k = discretization_cfg_->time_step();
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    auto const &solver_method_ptr =
        std::make_shared<wave_euler_solver_method>(euler_coeffs_, grid_cfg_, is_wave_sourse_set);
    if (is_wave_sourse_set)
    {
        wave_explicit_time_loop::run(solver_method_ptr, boundary_pair_, timer, last_time_idx, k, traverse_dir,
                                     prev_solution_0, prev_solution_1, wave_source, next_solution);
    }
    else
    {
        wave_explicit_time_loop::run(solver_method_ptr, boundary_pair_, timer, last_time_idx, k, traverse_dir,
                                     prev_solution_0, prev_solution_1, next_solution);
    }
}

void wave_euler_scheme::operator()(container_t &prev_solution_0, container_t &prev_solution_1,
                                   container_t &next_solution, bool is_wave_sourse_set,
                                   std::function<double(double, double)> const &wave_source,
                                   traverse_direction_enum traverse_dir, matrix_2d &solutions)
{
    auto const &timer = discretization_cfg_->time_range();
    const double k = discretization_cfg_->time_step();
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    auto const &solver_method_ptr =
        std::make_shared<wave_euler_solver_method>(euler_coeffs_, grid_cfg_, is_wave_sourse_set);
    if (is_wave_sourse_set)
    {
        wave_explicit_time_loop::run_with_stepping(solver_method_ptr, boundary_pair_, timer, last_time_idx, k,
                                                   traverse_dir, prev_solution_0, prev_solution_1, wave_source,
                                                   next_solution, solutions);
    }
    else
    {
        wave_explicit_time_loop::run_with_stepping(solver_method_ptr, boundary_pair_, timer, last_time_idx, k,
                                                   traverse_dir, prev_solution_0, prev_solution_1, next_solution,
                                                   solutions);
    }
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
