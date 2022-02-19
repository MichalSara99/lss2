#include "lss_heat_euler_scheme_2d.hpp"

#include "../../../../common/lss_macros.hpp"

namespace lss_pde_solvers
{

namespace two_dimensional
{

bool heat_euler_scheme_2d::is_stable(heat_coefficients_2d_ptr const &coefficients)
{
    // TODO: this needs to be implemented !!!
    return true;
}

void heat_euler_scheme_2d::initialize(heat_coefficients_2d_ptr const &coefficients)
{
    LSS_ASSERT(is_stable(coefficients) == true, "The chosen scheme is not stable");
    euler_coeffs_ = std::make_shared<heat_euler_coefficients_2d>(coefficients);
    // TODO: if more models present then cast this to appropriate concrete object
    heston_boundary_ = std::make_shared<heston_boundary_solver>(coefficients, grid_cfg_);
}

heat_euler_scheme_2d::heat_euler_scheme_2d(heat_coefficients_2d_ptr const &coefficients,
                                           boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                           boundary_2d_pair const &horizontal_boundary_pair,
                                           pde_discretization_config_2d_ptr const &discretization_config,
                                           grid_config_2d_ptr const &grid_config)
    : boundary_ver_{vertical_upper_boundary_ptr}, boundary_pair_hor_{horizontal_boundary_pair},
      discretization_cfg_{discretization_config}, grid_cfg_{grid_config}
{
    initialize(coefficients);
}

heat_euler_scheme_2d::~heat_euler_scheme_2d()
{
}

void heat_euler_scheme_2d::operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                                      std::function<double(double, double, double)> const &heat_source,
                                      traverse_direction_enum traverse_dir)
{
    auto const timer = discretization_cfg_->time_range();
    const double k = discretization_cfg_->time_step();
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    auto const &solver_method_ptr =
        std::make_shared<heat_euler_solver_method_2d>(euler_coeffs_, grid_cfg_, is_heat_sourse_set);
    if (is_heat_sourse_set)
    {
        throw std::exception("Not yet implemented.");
        // TODO: to be implemented!!!
        //  loop::run(solver_method_ptr, boundary_pair_, timer, last_time_idx, k,
        //  traverse_dir, heat_source, solution);
    }
    else
    {
        explicit_time_loop_2d::run(solver_method_ptr, heston_boundary_, boundary_pair_hor_, boundary_ver_, grid_cfg_,
                                   timer, last_time_idx, k, traverse_dir, prev_solution, next_solution);
    }
}

void heat_euler_scheme_2d::operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                                      std::function<double(double, double, double)> const &heat_source,
                                      traverse_direction_enum traverse_dir, matrix_3d &solutions)
{
    auto const timer = discretization_cfg_->time_range();
    const double k = discretization_cfg_->time_step();
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    auto const &solver_method_ptr =
        std::make_shared<heat_euler_solver_method_2d>(euler_coeffs_, grid_cfg_, is_heat_sourse_set);
    if (is_heat_sourse_set)
    {
        throw std::exception("Not yet implemented.");
        // TODO: to be implemented!!!
        //  loop::run(solver_method_ptr, boundary_pair_, timer, last_time_idx, k,
        //  traverse_dir, heat_source, solution);
    }
    else
    {
        explicit_time_loop_2d::run_with_stepping(solver_method_ptr, heston_boundary_, boundary_pair_hor_, boundary_ver_,
                                                 grid_cfg_, timer, last_time_idx, k, traverse_dir, prev_solution,
                                                 next_solution, solutions);
    }
}

} // namespace two_dimensional
} // namespace lss_pde_solvers
