#include "lss_explicit_time_loop_2d.hpp"

namespace lss_pde_solvers
{

namespace two_dimensional
{

void explicit_time_loop_2d::run(heat_explicit_solver_method_2d_ptr const &solver_ptr,
                                heston_boundary_solver_ptr const &boundary_solver_ptr,
                                boundary_2d_pair const &horizontal_boundary_pair,
                                boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                grid_config_2d_ptr const &grid_config, range_ptr const &time_range,
                                std::size_t const &last_time_idx, double const time_step,
                                traverse_direction_enum const &traverse_dir, matrix_2d &prev_solution,
                                matrix_2d &next_solution)
{
    const double start_time = time_range->lower();
    const double end_time = time_range->upper();
    const double k = time_step;

    double time{start_time + k};
    std::size_t time_idx{};
    if (traverse_dir == traverse_direction_enum::Forward)
    {
        time_idx = 1;
        while (time_idx <= last_time_idx)
        {
            solver_ptr->solve(prev_solution, time, next_solution);
            boundary_solver_ptr->solve(prev_solution, horizontal_boundary_pair, vertical_upper_boundary_ptr, time,
                                       next_solution);
            boundary_solver_ptr->solve(prev_solution, horizontal_boundary_pair, time, next_solution);
            prev_solution = next_solution;
            time += k;
            time_idx++;
        }
    }
    else if (traverse_dir == traverse_direction_enum::Backward)
    {
        time = end_time - k;
        time_idx = last_time_idx;
        do
        {
            time_idx--;
            solver_ptr->solve(prev_solution, time, next_solution);
            boundary_solver_ptr->solve(prev_solution, horizontal_boundary_pair, vertical_upper_boundary_ptr, time,
                                       next_solution);
            boundary_solver_ptr->solve(prev_solution, horizontal_boundary_pair, time, next_solution);
            prev_solution = next_solution;

            time -= k;
        } while (time_idx > 0);
    }
    else
    {
        throw std::exception("Unreachable");
    }
}

} // namespace two_dimensional
} // namespace lss_pde_solvers
