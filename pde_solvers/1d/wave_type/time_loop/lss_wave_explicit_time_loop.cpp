#include "lss_wave_explicit_time_loop.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

void wave_explicit_time_loop::run(wave_explicit_solver_method_ptr const &solver_ptr,
                                  boundary_1d_pair const &boundary_pair, range_ptr const &time_range,
                                  std::size_t const &last_time_idx, double const time_step,
                                  traverse_direction_enum const &traverse_dir, container_t &prev_solution_0,
                                  container_t &prev_solution_1, container_t &next_solution)
{
    const double start_time = time_range->lower();
    const double end_time = time_range->upper();
    const double k = time_step;
    double time{start_time};
    double next_time{time + k};
    std::size_t time_idx{};

    if (traverse_dir == traverse_direction_enum::Forward)
    {
        // solve for initial time step:
        solver_ptr->solve_initial(prev_solution_0, prev_solution_1, boundary_pair, time, next_time, next_solution);
        prev_solution_1 = next_solution;
        time_idx = 1;

        // solve for rest of time steps:
        time += k;
        next_time += k;
        time_idx++;
        while (time_idx <= last_time_idx)
        {
            solver_ptr->solve(prev_solution_0, prev_solution_1, boundary_pair, time, next_time, next_solution);
            prev_solution_0 = prev_solution_1;
            prev_solution_1 = next_solution;
            time += k;
            next_time += k;
            time_idx++;
        }
    }
    else if (traverse_dir == traverse_direction_enum::Backward)
    {
        time_idx = last_time_idx;
        time = end_time;
        next_time = time - k;
        // solve for initial time step:
        solver_ptr->solve_terminal(prev_solution_0, prev_solution_1, boundary_pair, time, next_time, next_solution);
        prev_solution_1 = next_solution;
        time_idx--;

        // solve for rest of time steps:
        time -= k;
        next_time -= k;
        do
        {
            time_idx--;
            solver_ptr->solve(prev_solution_0, prev_solution_1, boundary_pair, time, next_time, next_solution);
            prev_solution_0 = prev_solution_1;
            prev_solution_1 = next_solution;
            time -= k;
            next_time -= k;
        } while (time_idx > 0);
    }
    else
    {
        throw std::exception("Unreachable");
    }
}

void wave_explicit_time_loop::run(wave_explicit_solver_method_ptr const &solver_ptr,
                                  boundary_1d_pair const &boundary_pair, range_ptr const &time_range,
                                  std::size_t const &last_time_idx, double const time_step,
                                  traverse_direction_enum const &traverse_dir, container_t &prev_solution_0,
                                  container_t &prev_solution_1,
                                  std::function<double(double, double)> const &wave_source, container_t &next_solution)
{
    const double start_time = time_range->lower();
    const double end_time = time_range->upper();
    const double k = time_step;
    double time{start_time};
    double next_time{time + k};
    std::size_t time_idx{};

    if (traverse_dir == traverse_direction_enum::Forward)
    {
        // solve for initial time step:
        solver_ptr->solve_initial(prev_solution_0, prev_solution_1, boundary_pair, time, next_time, wave_source,
                                  next_solution);
        prev_solution_1 = next_solution;
        time_idx = 1;

        // solve for rest of time steps:
        time += k;
        next_time += k;
        time_idx++;
        while (time_idx <= last_time_idx)
        {
            solver_ptr->solve(prev_solution_0, prev_solution_1, boundary_pair, time, next_time, wave_source,
                              next_solution);
            prev_solution_0 = prev_solution_1;
            prev_solution_1 = next_solution;
            time += k;
            next_time += k;
            time_idx++;
        }
    }
    else if (traverse_dir == traverse_direction_enum::Backward)
    {
        time_idx = last_time_idx;
        time = end_time;
        next_time = time - k;
        // solve for initial time step:
        solver_ptr->solve_terminal(prev_solution_0, prev_solution_1, boundary_pair, time, next_time, wave_source,
                                   next_solution);
        prev_solution_1 = next_solution;
        time_idx--;

        // solve for rest of time steps:
        time -= k;
        next_time -= k;
        do
        {
            time_idx--;
            solver_ptr->solve(prev_solution_0, prev_solution_1, boundary_pair, time, next_time, wave_source,
                              next_solution);
            prev_solution_0 = prev_solution_1;
            prev_solution_1 = next_solution;
            time -= k;
            next_time -= k;
        } while (time_idx > 0);
    }
    else
    {
        throw std::exception("Unreachable");
    }
}

void wave_explicit_time_loop::run_with_stepping(wave_explicit_solver_method_ptr const &solver_ptr,
                                                boundary_1d_pair const &boundary_pair, range_ptr const &time_range,
                                                std::size_t const &last_time_idx, double const time_step,
                                                traverse_direction_enum const &traverse_dir,
                                                container_t &prev_solution_0, container_t &prev_solution_1,
                                                container_t &next_solution, matrix_2d &solutions)
{
    const double start_time = time_range->lower();
    const double end_time = time_range->upper();
    const double k = time_step;
    double time{start_time};
    double next_time{time + k};
    std::size_t time_idx{};

    if (traverse_dir == traverse_direction_enum::Forward)
    {
        // store the initial solution:
        solutions.row(0) = prev_solution_0;
        // solve for initial time step:
        solver_ptr->solve_initial(prev_solution_0, prev_solution_1, boundary_pair, time, next_time, next_solution);
        prev_solution_1 = next_solution;
        time_idx = 1;
        solutions.row(time_idx) = next_solution;

        // solve for rest of time steps:
        time += k;
        next_time += k;
        time_idx++;
        while (time_idx <= last_time_idx)
        {
            solver_ptr->solve(prev_solution_0, prev_solution_1, boundary_pair, time, next_time, next_solution);
            prev_solution_0 = prev_solution_1;
            prev_solution_1 = next_solution;
            solutions.row(time_idx) = next_solution;
            time += k;
            next_time += k;
            time_idx++;
        }
    }
    else if (traverse_dir == traverse_direction_enum::Backward)
    {
        time_idx = last_time_idx;
        // store the terminal solution:
        solutions.row(last_time_idx) = prev_solution_0;
        time = end_time;
        next_time = time - k;
        // solve for terminal time step:
        solver_ptr->solve_terminal(prev_solution_0, prev_solution_1, boundary_pair, time, next_time, next_solution);
        prev_solution_1 = next_solution;
        time_idx--;
        solutions.row(time_idx) = next_solution;

        // solve for rest of time steps:
        time -= k;
        next_time -= k;
        do
        {
            time_idx--;
            solver_ptr->solve(prev_solution_0, prev_solution_1, boundary_pair, time, next_time, next_solution);
            prev_solution_0 = prev_solution_1;
            prev_solution_1 = next_solution;
            solutions.row(time_idx) = next_solution;
            time -= k;
            next_time -= k;
        } while (time_idx > 0);
    }
    else
    {
        throw std::exception("Unreachable");
    }
}

void wave_explicit_time_loop::run_with_stepping(wave_explicit_solver_method_ptr const &solver_ptr,
                                                boundary_1d_pair const &boundary_pair, range_ptr const &time_range,
                                                std::size_t const &last_time_idx, double const time_step,
                                                traverse_direction_enum const &traverse_dir,
                                                container_t &prev_solution_0, container_t &prev_solution_1,
                                                std::function<double(double, double)> const &wave_source,
                                                container_t &next_solution, matrix_2d &solutions)
{
    const double start_time = time_range->lower();
    const double end_time = time_range->upper();
    const double k = time_step;
    double time{start_time};
    double next_time{time + k};
    std::size_t time_idx{};

    if (traverse_dir == traverse_direction_enum::Forward)
    {
        // store the initial solution:
        solutions.row(0) = prev_solution_0;
        // solve for initial time step:
        solver_ptr->solve_initial(prev_solution_0, prev_solution_1, boundary_pair, time, next_time, wave_source,
                                  next_solution);
        prev_solution_1 = next_solution;
        time_idx = 1;
        solutions.row(time_idx) = next_solution;

        // solve for rest of time steps:
        time += k;
        next_time += k;
        time_idx++;
        while (time_idx <= last_time_idx)
        {
            solver_ptr->solve(prev_solution_0, prev_solution_1, boundary_pair, time, next_time, wave_source,
                              next_solution);
            prev_solution_0 = prev_solution_1;
            prev_solution_1 = next_solution;
            solutions.row(time_idx) = next_solution;
            time += k;
            next_time += k;
            time_idx++;
        }
    }
    else if (traverse_dir == traverse_direction_enum::Backward)
    {
        time_idx = last_time_idx;
        // store the terminal solution:
        solutions.row(last_time_idx) = prev_solution_0;
        time = end_time;
        next_time = time - k;
        // solve for terminal time step:
        solver_ptr->solve_terminal(prev_solution_0, prev_solution_1, boundary_pair, time, next_time, wave_source,
                                   next_solution);
        prev_solution_1 = next_solution;
        time_idx--;
        solutions.row(time_idx) = next_solution;

        // solve for rest of time steps:
        time -= k;
        next_time -= k;
        do
        {
            time_idx--;
            solver_ptr->solve(prev_solution_0, prev_solution_1, boundary_pair, time, next_time, wave_source,
                              next_solution);
            prev_solution_0 = prev_solution_1;
            prev_solution_1 = next_solution;
            solutions.row(time_idx) = next_solution;
            time -= k;
            next_time -= k;
        } while (time_idx > 0);
    }
    else
    {
        throw std::exception("Unreachable");
    }
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
