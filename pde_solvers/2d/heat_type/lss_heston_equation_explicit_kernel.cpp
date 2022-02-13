#include "lss_heston_equation_explicit_kernel.hpp"

#include "../../../discretization/lss_grid.hpp"
#include "explicit_schemes/lss_heston_euler_cuda_solver_scheme.hpp"
#include "explicit_schemes/lss_heston_euler_solver_scheme.hpp"

namespace lss_pde_solvers
{

using lss_enumerations::explicit_pde_schemes_enum;

namespace two_dimensional
{

heston_equation_explicit_kernel<memory_space_enum::Device>::heston_equation_explicit_kernel(
    boundary_2d_ptr const &vertical_upper_boundary_ptr, boundary_2d_pair const &horizontal_boundary_pair,
    heat_data_transform_2d_ptr const &heat_data_config, pde_discretization_config_2d_ptr const &discretization_config,
    heat_explicit_solver_config_ptr const &solver_config, grid_config_2d_ptr const &grid_config)
    : boundary_ver_{vertical_upper_boundary_ptr}, boundary_pair_hor_{horizontal_boundary_pair},
      heat_data_cfg_{heat_data_config}, discretization_cfg_{discretization_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void heston_equation_explicit_kernel<memory_space_enum::Device>::operator()(
    container_2d<by_enum::Row> &prev_solution, container_2d<by_enum::Row> &next_solution, bool is_heat_sourse_set,
    std::function<double(double, double, double)> const &heat_source)
{
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();
    // create a Heston coefficient holder:
    auto const heston_coeff_holder =
        std::make_shared<heston_implicit_coefficients>(heat_data_cfg_, discretization_cfg_, nullptr, 0.0);
    // Here make a dicision which explicit scheme to launch:
    if (solver_cfg_->explicit_pde_scheme() == explicit_pde_schemes_enum::Euler)
    {
        heston_euler_cuda_scheme euler_scheme(heston_coeff_holder, boundary_ver_, boundary_pair_hor_,
                                              discretization_cfg_, grid_cfg_);
        euler_scheme(prev_solution, next_solution, is_heat_sourse_set, heat_source, traverse_dir);
    }
    else if (solver_cfg_->explicit_pde_scheme() == explicit_pde_schemes_enum::ADEBarakatClark)
    {
        throw std::exception("Not currently supported");
    }
    else if (solver_cfg_->explicit_pde_scheme() == explicit_pde_schemes_enum::ADESaulyev)
    {
        throw std::exception("Not currently supported");
    }
    else
    {
        throw std::exception("Unreachable");
    }
}

heston_equation_explicit_kernel<memory_space_enum::Host>::heston_equation_explicit_kernel(
    boundary_2d_ptr const &vertical_upper_boundary_ptr, boundary_2d_pair const &horizontal_boundary_pair,
    heat_data_transform_2d_ptr const &heat_data_config, pde_discretization_config_2d_ptr const &discretization_config,
    heat_explicit_solver_config_ptr const &solver_config, grid_config_2d_ptr const &grid_config)
    : boundary_ver_{vertical_upper_boundary_ptr}, boundary_pair_hor_{horizontal_boundary_pair},
      heat_data_cfg_{heat_data_config}, discretization_cfg_{discretization_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void heston_equation_explicit_kernel<memory_space_enum::Host>::operator()(
    container_2d<by_enum::Row> &prev_solution, container_2d<by_enum::Row> &next_solution, bool is_heat_sourse_set,
    std::function<double(double, double, double)> const &heat_source)
{
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();
    // create a Heston coefficient holder:
    auto const heston_coeff_holder =
        std::make_shared<heston_implicit_coefficients>(heat_data_cfg_, discretization_cfg_, nullptr, 0.0);
    // Here make a dicision which explicit scheme to launch:
    if (solver_cfg_->explicit_pde_scheme() == explicit_pde_schemes_enum::Euler)
    {
        heston_euler_scheme euler_scheme(heston_coeff_holder, boundary_ver_, boundary_pair_hor_, discretization_cfg_,
                                         grid_cfg_);
        euler_scheme(prev_solution, next_solution, is_heat_sourse_set, heat_source, traverse_dir);
    }
    else if (solver_cfg_->explicit_pde_scheme() == explicit_pde_schemes_enum::ADEBarakatClark)
    {
        throw std::exception("Not currently supported");
    }
    else if (solver_cfg_->explicit_pde_scheme() == explicit_pde_schemes_enum::ADESaulyev)
    {
        throw std::exception("Not currently supported");
    }
    else
    {
        throw std::exception("Unreachable");
    }
}

} // namespace two_dimensional

} // namespace lss_pde_solvers
