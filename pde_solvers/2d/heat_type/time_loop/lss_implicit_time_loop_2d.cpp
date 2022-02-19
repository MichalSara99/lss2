#include "lss_implicit_time_loop_2d.hpp"

#include "../../../../boundaries/lss_dirichlet_boundary.hpp"

namespace lss_pde_solvers
{

using lss_boundary::dirichlet_boundary_2d;

namespace two_dimensional
{

boundary_2d_pair implicit_boundary::get_vertical(grid_config_1d_ptr const &grid_config_x,
                                                 matrix_2d const &next_solution)
{
    auto const lci = next_solution.columns() - 1;
    auto lower = [=](double t, double x) -> double {
        const std::size_t i = grid_config_x->index_of(x);
        return next_solution(i, 0);
    };
    auto upper = [=](double t, double x) -> double {
        const std::size_t i = grid_config_x->index_of(x);
        return next_solution(i, lci);
    };
    auto vertical_low = std::make_shared<dirichlet_boundary_2d>(lower);
    auto vertical_high = std::make_shared<dirichlet_boundary_2d>(upper);
    return std::make_pair(vertical_low, vertical_high);
}

boundary_2d_pair implicit_boundary::get_intermed_horizontal(grid_config_1d_ptr const &grid_config_y,
                                                            matrix_2d const &next_solution)
{
    const std::size_t lri = next_solution.rows() - 1;
    auto lower = [=](double t, double y) -> double {
        const std::size_t j = grid_config_y->index_of(y);
        return next_solution(0, j);
    };
    auto upper = [=](double t, double y) -> double {
        const std::size_t j = grid_config_y->index_of(y);
        return next_solution(lri, j);
    };
    auto horizontal_low = std::make_shared<dirichlet_boundary_2d>(lower);
    auto horizontal_high = std::make_shared<dirichlet_boundary_2d>(upper);
    return std::make_pair(horizontal_low, horizontal_high);
}

void implicit_time_loop_2d::run(heat_splitting_method_ptr const &solver_ptr,
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
    boundary_2d_pair ver_boundary_pair;
    boundary_2d_pair hor_inter_boundary_pair;

    double time{};
    std::size_t time_idx{};

    if (traverse_dir == traverse_direction_enum::Forward)
    {
        time = start_time + k;
        time_idx = 1;
        while (time_idx <= last_time_idx)
        {
            boundary_solver_ptr->solve(prev_solution, horizontal_boundary_pair, vertical_upper_boundary_ptr, time,
                                       next_solution);
            ver_boundary_pair = implicit_boundary::get_vertical(grid_config->grid_1(), next_solution);
            hor_inter_boundary_pair = implicit_boundary::get_intermed_horizontal(grid_config->grid_2(), prev_solution);
            solver_ptr->solve(prev_solution, hor_inter_boundary_pair, ver_boundary_pair, time, next_solution);
            boundary_solver_ptr->solve(prev_solution, horizontal_boundary_pair, time, next_solution);

            prev_solution = next_solution;
            time += k;
            time_idx++;
        }
    }
    else
    {
        time = end_time - k;
        time_idx = last_time_idx;
        do
        {
            time_idx--;
            boundary_solver_ptr->solve(prev_solution, horizontal_boundary_pair, vertical_upper_boundary_ptr, time,
                                       next_solution);
            ver_boundary_pair = implicit_boundary::get_vertical(grid_config->grid_1(), next_solution);
            hor_inter_boundary_pair = implicit_boundary::get_intermed_horizontal(grid_config->grid_2(), prev_solution);
            solver_ptr->solve(prev_solution, hor_inter_boundary_pair, ver_boundary_pair, time, next_solution);
            boundary_solver_ptr->solve(prev_solution, horizontal_boundary_pair, time, next_solution);

            prev_solution = next_solution;
            time -= k;
        } while (time_idx > 0);
    }
}

void implicit_time_loop_2d::run_with_stepping(heat_splitting_method_ptr const &solver_ptr,
                                              heston_boundary_solver_ptr const &boundary_solver_ptr,
                                              boundary_2d_pair const &horizontal_boundary_pair,
                                              boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                              grid_config_2d_ptr const &grid_config, range_ptr const &time_range,
                                              std::size_t const &last_time_idx, double const time_step,
                                              traverse_direction_enum const &traverse_dir, matrix_2d &prev_solution,
                                              matrix_2d &next_solution, matrix_3d &solutions)
{
    const double start_time = time_range->lower();
    const double end_time = time_range->upper();
    const double k = time_step;
    boundary_2d_pair ver_boundary_pair;
    boundary_2d_pair hor_inter_boundary_pair;

    double time{};
    std::size_t time_idx{};

    if (traverse_dir == traverse_direction_enum::Forward)
    {
        // store the initial solution:
        solutions.layer_plane(0) = prev_solution.data();
        time = start_time + k;
        time_idx = 1;
        while (time_idx <= last_time_idx)
        {
            boundary_solver_ptr->solve(prev_solution, horizontal_boundary_pair, vertical_upper_boundary_ptr, time,
                                       next_solution);
            ver_boundary_pair = implicit_boundary::get_vertical(grid_config->grid_1(), next_solution);
            hor_inter_boundary_pair = implicit_boundary::get_intermed_horizontal(grid_config->grid_2(), prev_solution);
            solver_ptr->solve(prev_solution, hor_inter_boundary_pair, ver_boundary_pair, time, next_solution);
            boundary_solver_ptr->solve(prev_solution, horizontal_boundary_pair, time, next_solution);

            solutions.layer_plane(time_idx) = next_solution.data();
            prev_solution = next_solution;
            time += k;
            time_idx++;
        }
    }
    else
    {
        solutions.layer_plane(last_time_idx) = prev_solution.data();
        time = end_time - k;
        time_idx = last_time_idx;
        do
        {
            time_idx--;
            boundary_solver_ptr->solve(prev_solution, horizontal_boundary_pair, vertical_upper_boundary_ptr, time,
                                       next_solution);
            ver_boundary_pair = implicit_boundary::get_vertical(grid_config->grid_1(), next_solution);
            hor_inter_boundary_pair = implicit_boundary::get_intermed_horizontal(grid_config->grid_2(), prev_solution);
            solver_ptr->solve(prev_solution, hor_inter_boundary_pair, ver_boundary_pair, time, next_solution);
            boundary_solver_ptr->solve(prev_solution, horizontal_boundary_pair, time, next_solution);

            solutions.layer_plane(time_idx) = next_solution.data();
            prev_solution = next_solution;
            time -= k;
        } while (time_idx > 0);
    }
}

} // namespace two_dimensional
} // namespace lss_pde_solvers
