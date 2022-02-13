#include "lss_implicit_time_loop.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

void implicit_time_loop::run(implicit_solver_method_ptr const &solver_ptr, boundary_1d_pair const &boundary_pair,
                             range_ptr const &time_range, std::size_t const &last_time_idx, double const time_step,
                             traverse_direction_enum const &traverse_dir, container_t &solution)
{
    const double start_time = time_range->lower();
    const double end_time = time_range->upper();
    const double k = time_step;
    // container for next solution:
    container_t next_solution(double{}, solution.size());

    double time{start_time + k};
    std::size_t time_idx{};
    if (traverse_dir == traverse_direction_enum::Forward)
    {
        time_idx = 1;
        while (time_idx <= last_time_idx)
        {
            solver_ptr->solve(solution, boundary_pair, time, next_solution);
            solution = next_solution;
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
            solver_ptr->solve(solution, boundary_pair, time, next_solution);
            solution = next_solution;
            time -= k;
        } while (time_idx > 0);
    }
}

void implicit_time_loop::run(implicit_solver_method_ptr const &solver_ptr, boundary_1d_pair const &boundary_pair,
                             range_ptr const &time_range, std::size_t const &last_time_idx, double const time_step,
                             traverse_direction_enum const &traverse_dir,
                             std::function<double(double, double)> const &heat_source, container_t &solution)
{
    const double start_time = time_range->lower();
    const double end_time = time_range->upper();
    const double k = time_step;
    // container for next solution:
    container_t next_solution(double{}, solution.size());

    double time{start_time + k};
    double next_time{time + k};
    std::size_t time_idx{};
    if (traverse_dir == traverse_direction_enum::Forward)
    {
        time_idx = 1;
        while (time_idx <= last_time_idx)
        {
            solver_ptr->solve(solution, boundary_pair, time, next_time, heat_source, next_solution);
            solution = next_solution;
            time += k;
            next_time += k;
            time_idx++;
        }
    }
    else
    {
        time = end_time - k;
        next_time = time - k;
        time_idx = last_time_idx;
        do
        {
            time_idx--;
            solver_ptr->solve(solution, boundary_pair, time, next_time, heat_source, next_solution);
            solution = next_solution;
            time -= k;
            next_time -= k;
        } while (time_idx > 0);
    }
}

void implicit_time_loop::run_with_stepping(implicit_solver_method_ptr const &solver_ptr,
                                           boundary_1d_pair const &boundary_pair, range_ptr const &time_range,
                                           std::size_t const &last_time_idx, double const time_step,
                                           traverse_direction_enum const &traverse_dir, container_t &solution,
                                           matrix_2d &solutions)
{
    const double start_time = time_range->lower();
    const double end_time = time_range->upper();
    const double k = time_step;
    // container for next solution:
    container_t next_solution(double{}, solution.size());

    double time{start_time + k};
    std::size_t time_idx{};
    if (traverse_dir == traverse_direction_enum::Forward)
    {
        // store the initial solution:
        solutions.row(0) = solution;
        time_idx = 1;
        while (time_idx <= last_time_idx)
        {
            solver_ptr->solve(solution, boundary_pair, time, next_solution);
            solutions.row(time_idx) = next_solution;
            solution = next_solution;
            time += k;
            time_idx++;
        }
    }
    else
    {
        // store the initial solution:
        solutions.row(last_time_idx) = solution;
        time = end_time - k;
        time_idx = last_time_idx;
        do
        {
            time_idx--;
            solver_ptr->solve(solution, boundary_pair, time, next_solution);
            solutions.row(time_idx) = next_solution;
            solution = next_solution;
            time -= k;
        } while (time_idx > 0);
    }
}

void implicit_time_loop::run_with_stepping(implicit_solver_method_ptr const &solver_ptr,
                                           boundary_1d_pair const &boundary_pair, range_ptr const &time_range,
                                           std::size_t const &last_time_idx, double const time_step,
                                           traverse_direction_enum const &traverse_dir,
                                           std::function<double(double, double)> const &heat_source,
                                           container_t &solution, matrix_2d &solutions)
{
    const double start_time = time_range->lower();
    const double end_time = time_range->upper();
    const double k = time_step;
    // container for next solution:
    container_t next_solution(double{}, solution.size());

    double time{start_time + k};
    double next_time{time + k};
    std::size_t time_idx{};
    if (traverse_dir == traverse_direction_enum::Forward)
    {
        // store the initial solution:
        solutions.row(0) = solution;
        time_idx = 1;
        while (time_idx <= last_time_idx)
        {
            solver_ptr->solve(solution, boundary_pair, time, next_time, heat_source, next_solution);
            solutions.row(time_idx) = next_solution;
            solution = next_solution;
            time += k;
            next_time += k;
            time_idx++;
        }
    }
    else
    {
        // store the initial solution:
        solutions.row(last_time_idx) = solution;
        time = end_time - k;
        next_time = time - k;
        time_idx = last_time_idx;
        do
        {
            time_idx--;
            solver_ptr->solve(solution, boundary_pair, time, next_time, heat_source, next_solution);
            solutions.row(time_idx) = next_solution;
            solution = next_solution;
            time -= k;
            next_time -= k;
        } while (time_idx > 0);
    }
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
