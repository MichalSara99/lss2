#include "lss_heston_euler_solver_method.hpp"

#include "../../../../discretization/lss_discretization.hpp"
#include "../../../../discretization/lss_grid.hpp"

namespace lss_pde_solvers
{

using lss_grids::grid_2d;
using d_2d = discretization_2d<std::vector, std::allocator<double>>;

namespace two_dimensional
{

void explicit_heston_scheme::rhs(heston_euler_coefficients_ptr const &cfs, grid_config_2d_ptr const &grid_config,
                                 container_2d<by_enum::Row> const &input, double const &time,
                                 container_2d<by_enum::Row> &solution)
{
    auto const one = 1.0;
    auto const &M = cfs->M_;
    auto const &M_tilde = cfs->M_tilde_;
    auto const &P = cfs->P_;
    auto const &P_tilde = cfs->P_tilde_;
    auto const &Z = cfs->Z_;
    auto const &W = cfs->W_;
    auto const &C = cfs->C_;
    auto const rows = input.rows() - 1;
    auto const cols = input.columns() - 1;

    double x{}, y{}, val{};
    for (std::size_t r = 1; r < rows; ++r)
    {
        x = grid_2d::value_1(grid_config, r);
        for (std::size_t c = 1; c < cols; ++c)
        {
            y = grid_2d::value_2(grid_config, c);
            val = C(time, x, y) * input(r - 1, c - 1) + M(time, x, y) * input(r - 1, c) -
                  C(time, x, y) * input(r - 1, c + 1) + M_tilde(time, x, y) * input(r, c - 1) +
                  (one - Z(time, x, y) - W(time, x, y)) * input(r, c) + P_tilde(time, x, y) * input(r, c + 1) -
                  C(time, x, y) * input(r + 1, c - 1) + P(time, x, y) * input(r + 1, c) +
                  C(time, x, y) * input(r + 1, c + 1);
            solution(r, c, val);
        }
    }
}

void explicit_heston_scheme::rhs_source(heston_euler_coefficients_ptr const &cfs, grid_config_2d_ptr const &grid_config,
                                        container_2d<by_enum::Row> const &input, double const &time,
                                        container_2d<by_enum::Row> const &inhom_input,
                                        container_2d<by_enum::Row> &solution)
{
    auto const one = 1.0;
    auto const &M = cfs->M_;
    auto const &M_tilde = cfs->M_tilde_;
    auto const &P = cfs->P_;
    auto const &P_tilde = cfs->P_tilde_;
    auto const &Z = cfs->Z_;
    auto const &W = cfs->W_;
    auto const &C = cfs->C_;
    auto const rho = cfs->rho_;
    auto const rows = input.rows() - 1; // without boundary
    auto const cols = input.columns() - 1;

    double x{}, y{}, val{};
    for (std::size_t r = 1; r < rows; ++r)
    {
        x = grid_2d::value_1(grid_config, r);
        for (std::size_t c = 1; c < cols; ++c)
        {
            y = grid_2d::value_2(grid_config, c);
            val = C(time, x, y) * input(r - 1, c - 1) + M(time, x, y) * input(r - 1, c) -
                  C(time, x, y) * input(r - 1, c + 1) + M_tilde(time, x, y) * input(r, c - 1) +
                  (one - Z(time, x, y) - W(time, x, y)) * input(r, c) + P_tilde(time, x, y) * input(r, c + 1) -
                  C(time, x, y) * input(r + 1, c - 1) + P(time, x, y) * input(r + 1, c) +
                  C(time, x, y) * input(r + 1, c + 1) + rho * inhom_input(r, c);
            solution(r, c, val);
        }
    }
}

void heston_euler_solver_method::initialize(bool is_heat_source_set)
{
    if (is_heat_source_set)
    {
        source_ =
            std::make_shared<container_2d<by_enum::Row>>(coefficients_->space_size_x_, coefficients_->space_size_y_);
    }
}

heston_euler_solver_method::heston_euler_solver_method(heston_euler_coefficients_ptr const &coefficients,
                                                       grid_config_2d_ptr const &grid_config, bool is_heat_source_set)
    : coefficients_{coefficients}, heston_explicit_solver_method(grid_config)
{
    initialize(is_heat_source_set);
}

heston_euler_solver_method::~heston_euler_solver_method()
{
}

void heston_euler_solver_method::solve(container_2d<by_enum::Row> &prev_solution, double const &time,
                                       container_2d<by_enum::Row> &solution)
{
    explicit_heston_scheme::rhs(coefficients_, grid_cfg_, prev_solution, time, solution);
}

void heston_euler_solver_method::solve(container_2d<by_enum::Row> &prev_solution, double const &time,
                                       std::function<double(double, double, double)> const &heat_source,
                                       container_2d<by_enum::Row> &solution)
{
    d_2d::of_function(grid_cfg_, time, heat_source, *source_);
    explicit_heston_scheme::rhs_source(coefficients_, grid_cfg_, prev_solution, time, *source_, solution);
}

} // namespace two_dimensional

} // namespace lss_pde_solvers
