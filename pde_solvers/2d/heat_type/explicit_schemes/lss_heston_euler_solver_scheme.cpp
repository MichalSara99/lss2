#include "lss_heston_euler_solver_scheme.hpp"

#include "../../../../common/lss_macros.hpp"

namespace lss_pde_solvers
{

namespace two_dimensional
{

bool heston_euler_scheme::is_stable(heston_implicit_coefficients_ptr const &coefficients)
{
    // TODO: this needs to be implemented !!!
    return true;
}

void heston_euler_scheme::initialize(heston_implicit_coefficients_ptr const &coefficients)
{
    LSS_ASSERT(is_stable(coefficients) == true, "The chosen scheme is not stable");
    euler_coeffs_ = std::make_shared<heston_euler_coefficients>(coefficients);
    heston_boundary_ = std::make_shared<heston_explicit_boundary_solver>(coefficients, grid_cfg_);
}

heston_euler_scheme::heston_euler_scheme(heston_implicit_coefficients_ptr const &coefficients,
                                         boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                         boundary_2d_pair const &horizontal_boundary_pair,
                                         pde_discretization_config_2d_ptr const &discretization_config,
                                         grid_config_2d_ptr const &grid_config)
    : boundary_ver_{vertical_upper_boundary_ptr}, boundary_pair_hor_{horizontal_boundary_pair},
      discretization_cfg_{discretization_config}, grid_cfg_{grid_config}
{
    initialize(coefficients);
}

heston_euler_scheme::~heston_euler_scheme()
{
}

void heston_euler_scheme::operator()(container_2d<by_enum::Row> &prev_solution,
                                     container_2d<by_enum::Row> &next_solution, bool is_heat_sourse_set,
                                     std::function<double(double, double, double)> const &heat_source,
                                     traverse_direction_enum traverse_dir)
{
    auto const timer = discretization_cfg_->time_range();
    const double k = discretization_cfg_->time_step();
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    auto const &solver_method_ptr =
        std::make_shared<heston_euler_solver_method>(euler_coeffs_, grid_cfg_, is_heat_sourse_set);
    if (is_heat_sourse_set)
    {
        // TODO: to be implemented!!!
        //  loop::run(solver_method_ptr, boundary_pair_, timer, last_time_idx, k,
        //  traverse_dir, heat_source, solution);
    }
    else
    {
        heston_explicit_time_loop::run(solver_method_ptr, heston_boundary_, boundary_pair_hor_, boundary_ver_,
                                       grid_cfg_, timer, last_time_idx, k, traverse_dir, prev_solution, next_solution);
    }
}

} // namespace two_dimensional
} // namespace lss_pde_solvers
