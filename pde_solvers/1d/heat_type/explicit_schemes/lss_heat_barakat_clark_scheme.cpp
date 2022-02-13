#include "lss_heat_barakat_clark_scheme.hpp"

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../../boundaries/lss_robin_boundary.hpp"
#include "../solver_method/lss_heat_barakat_clark_solver_method.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_1d_ptr;
using lss_boundary::dirichlet_boundary_1d;
using lss_boundary::neumann_boundary_1d;
using lss_boundary::robin_boundary_1d;

void heat_barakat_clark_scheme::initialize(heat_coefficients_ptr const &coefficients)
{
    auto const &first = boundary_pair_.first;
    if (std::dynamic_pointer_cast<neumann_boundary_1d>(first))
    {
        throw std::exception("Neumann boundary type is not supported for this scheme");
    }
    if (std::dynamic_pointer_cast<robin_boundary_1d>(first))
    {
        throw std::exception("Robin boundary type is not supported for this scheme");
    }
    auto const &second = boundary_pair_.second;
    if (std::dynamic_pointer_cast<neumann_boundary_1d>(second))
    {
        throw std::exception("Neumann boundary type is not supported for this scheme");
    }
    if (std::dynamic_pointer_cast<robin_boundary_1d>(second))
    {
        throw std::exception("Robin boundary type is not supported for this scheme");
    }
    bc_coeffs_ = std::make_shared<heat_barakat_clark_coefficients>(coefficients);
}

heat_barakat_clark_scheme::heat_barakat_clark_scheme(heat_coefficients_ptr const &coefficients,
                                                     boundary_1d_pair const &boundary_pair,
                                                     pde_discretization_config_1d_ptr const &discretization_config,
                                                     grid_config_1d_ptr const &grid_config)
    : boundary_pair_{boundary_pair}, discretization_cfg_{discretization_config}, grid_cfg_{grid_config}
{
    initialize(coefficients);
}

heat_barakat_clark_scheme::~heat_barakat_clark_scheme()
{
}

void heat_barakat_clark_scheme::operator()(container_t &solution, bool is_heat_sourse_set,
                                           std::function<double(double, double)> const &heat_source,
                                           traverse_direction_enum traverse_dir)
{
    auto const &timer = discretization_cfg_->time_range();
    const double k = discretization_cfg_->time_step();
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    auto const &solver_method_ptr =
        std::make_shared<heat_barakat_clark_solver_method>(bc_coeffs_, grid_cfg_, is_heat_sourse_set);
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

void heat_barakat_clark_scheme::operator()(container_t &solution, bool is_heat_sourse_set,
                                           std::function<double(double, double)> const &heat_source,
                                           traverse_direction_enum traverse_dir, matrix_2d &solutions)
{
    auto const &timer = discretization_cfg_->time_range();
    const double k = discretization_cfg_->time_step();
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    auto const &solver_method_ptr =
        std::make_shared<heat_barakat_clark_solver_method>(bc_coeffs_, grid_cfg_, is_heat_sourse_set);
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
