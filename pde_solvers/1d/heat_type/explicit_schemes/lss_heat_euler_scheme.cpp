#include "lss_heat_euler_scheme.hpp"

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../../boundaries/lss_robin_boundary.hpp"
#include "../../../../discretization/lss_grid.hpp"
#include "../explicit_coefficients/lss_heat_euler_coefficients.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_grids::grid_1d;

bool heat_euler_scheme::is_stable(heat_coefficients_ptr const &coefficients)
{
    const double zero = 0.0;
    const double one = 1.0;
    const double two = 2.0;
    auto const &A = coefficients->A_;
    auto const &B = coefficients->B_;
    auto const &D = coefficients->D_;
    const double k = coefficients->k_;
    const double lambda = coefficients->lambda_;
    const double gamma = coefficients->gamma_;
    const double delta = coefficients->delta_;
    auto const &a = [=](double t, double x) { return ((A(t, x) + D(t, x)) / (two * lambda)); };
    auto const &b = [=](double t, double x) { return ((D(t, x) - A(t, x)) / (two * gamma)); };
    auto const &c = [=](double t, double x) { return ((lambda * a(t, x) - B(t, x)) / delta); };
    const std::size_t space_size = discretization_cfg_->number_of_space_points();
    auto const &ftime = discretization_cfg_->time_range()->upper();
    double x{}, t{k};
    while (t <= ftime)
    {
        for (std::size_t i = 0; i < space_size; ++i)
        {
            x = grid_1d::value(grid_cfg_, i);
            if (c(t, x) > zero)
                return false;
            if ((two * lambda * a(t, x) - k * c(t, x)) > one)
                return false;
            if (((gamma * std::abs(b(t, x))) * (gamma * std::abs(b(t, x)))) > (two * lambda * a(t, x)))
                return false;
        }
        t += k;
    }
    return true;
}

void heat_euler_scheme::initialize(heat_coefficients_ptr const &coefficients)
{
    LSS_ASSERT(is_stable(coefficients) == true, "The chosen scheme is not stable");
    euler_coeffs_ = std::make_shared<heat_euler_coefficients>(coefficients);
}

heat_euler_scheme::heat_euler_scheme(heat_coefficients_ptr const &coefficients, boundary_1d_pair const &boundary_pair,
                                     pde_discretization_config_1d_ptr const &discretization_config,
                                     grid_config_1d_ptr const &grid_config)
    : boundary_pair_{boundary_pair}, discretization_cfg_{discretization_config}, grid_cfg_{grid_config}
{
    initialize(coefficients);
}

heat_euler_scheme::~heat_euler_scheme()
{
}

void heat_euler_scheme::operator()(container_t &solution, bool is_heat_sourse_set,
                                   std::function<double(double, double)> const &heat_source,
                                   traverse_direction_enum traverse_dir)
{
    auto const &timer = discretization_cfg_->time_range();
    const double k = discretization_cfg_->time_step();
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    auto const &solver_method_ptr =
        std::make_shared<heat_euler_solver_method>(euler_coeffs_, grid_cfg_, is_heat_sourse_set);
    if (is_heat_sourse_set)
    {
        explicit_time_loop::run(solver_method_ptr, boundary_pair_, timer, last_time_idx, k, traverse_dir, heat_source,
                                solution);
    }
    else
    {
        explicit_time_loop::run(solver_method_ptr, boundary_pair_, timer, last_time_idx, k, traverse_dir, solution);
    }
}

void heat_euler_scheme::operator()(container_t &solution, bool is_heat_sourse_set,
                                   std::function<double(double, double)> const &heat_source,
                                   traverse_direction_enum traverse_dir, matrix_2d &solutions)
{
    auto const &timer = discretization_cfg_->time_range();
    const double k = discretization_cfg_->time_step();
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    auto const &solver_method_ptr =
        std::make_shared<heat_euler_solver_method>(euler_coeffs_, grid_cfg_, is_heat_sourse_set);
    if (is_heat_sourse_set)
    {
        explicit_time_loop::run_with_stepping(solver_method_ptr, boundary_pair_, timer, last_time_idx, k, traverse_dir,
                                              heat_source, solution, solutions);
    }
    else
    {
        explicit_time_loop::run_with_stepping(solver_method_ptr, boundary_pair_, timer, last_time_idx, k, traverse_dir,
                                              solution, solutions);
    }
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
