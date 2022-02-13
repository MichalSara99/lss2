#include "lss_heat_equation_explicit_kernel.hpp"

#include "explicit_schemes/lss_heat_barakat_clark_scheme.hpp"
#include "explicit_schemes/lss_heat_euler_cuda_scheme.hpp"
#include "explicit_schemes/lss_heat_euler_scheme.hpp"
#include "explicit_schemes/lss_heat_saulyev_scheme.hpp"
#include "implicit_coefficients/lss_heat_coefficients.hpp"

namespace lss_pde_solvers
{
namespace one_dimensional
{

heat_equation_explicit_kernel<memory_space_enum::Device>::heat_equation_explicit_kernel(
    boundary_1d_pair const &boundary_pair, heat_data_transform_1d_ptr const &heat_data_config,
    pde_discretization_config_1d_ptr const &discretization_config, heat_explicit_solver_config_ptr const &solver_config,
    grid_config_1d_ptr const &grid_config)
    : boundary_pair_{boundary_pair}, heat_data_cfg_{heat_data_config}, discretization_cfg_{discretization_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void heat_equation_explicit_kernel<memory_space_enum::Device>::operator()(
    container_t &solution, bool is_heat_sourse_set, std::function<double(double, double)> const &heat_source)
{
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();
    // create a heat coefficient holder:
    auto const heat_coeff_holder = std::make_shared<heat_coefficients>(heat_data_cfg_, discretization_cfg_, double{});
    // Here make a dicision which explicit scheme to launch:
    if (solver_cfg_->explicit_pde_scheme() == explicit_pde_schemes_enum::Euler)
    {
        heat_euler_cuda_scheme euler_scheme(heat_coeff_holder, boundary_pair_, discretization_cfg_, grid_cfg_);
        euler_scheme(solution, is_heat_sourse_set, heat_source, traverse_dir);
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

void heat_equation_explicit_kernel<memory_space_enum::Device>::operator()(
    container_t &solution, bool is_heat_sourse_set, std::function<double(double, double)> const &heat_source,
    matrix_2d &solutions) {
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();
    // create a heat coefficient holder:
    auto const heat_coeff_holder = std::make_shared<heat_coefficients>(heat_data_cfg_, discretization_cfg_, double{});
    // Here make a dicision which explicit scheme to launch:
    if (solver_cfg_->explicit_pde_scheme() == explicit_pde_schemes_enum::Euler)
    {
        heat_euler_cuda_scheme euler_scheme(heat_coeff_holder, boundary_pair_, discretization_cfg_, grid_cfg_);
        euler_scheme(solution, is_heat_sourse_set, heat_source, traverse_dir, solutions);
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

heat_equation_explicit_kernel<memory_space_enum::Host>::heat_equation_explicit_kernel(
    boundary_1d_pair const &boundary_pair, heat_data_transform_1d_ptr const &heat_data_config,
    pde_discretization_config_1d_ptr const &discretization_config, heat_explicit_solver_config_ptr const &solver_config,
    grid_config_1d_ptr const &grid_config)
    : boundary_pair_{boundary_pair}, heat_data_cfg_{heat_data_config}, discretization_cfg_{discretization_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void heat_equation_explicit_kernel<memory_space_enum::Host>::operator()(
    container_t &solution, bool is_heat_sourse_set, std::function<double(double, double)> const &heat_source)
{
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();
    // create a heat coefficient holder:
    auto const heat_coeff_holder = std::make_shared<heat_coefficients>(heat_data_cfg_, discretization_cfg_, double{});
    // Here make a dicision which explicit scheme to launch:
    if (solver_cfg_->explicit_pde_scheme() == explicit_pde_schemes_enum::Euler)
    {
        heat_euler_scheme euler_scheme(heat_coeff_holder, boundary_pair_, discretization_cfg_, grid_cfg_);
        euler_scheme(solution, is_heat_sourse_set, heat_source, traverse_dir);
    }
    else if (solver_cfg_->explicit_pde_scheme() == explicit_pde_schemes_enum::ADEBarakatClark)
    {
        heat_barakat_clark_scheme bc_scheme(heat_coeff_holder, boundary_pair_, discretization_cfg_, grid_cfg_);
        bc_scheme(solution, is_heat_sourse_set, heat_source, traverse_dir);
    }
    else if (solver_cfg_->explicit_pde_scheme() == explicit_pde_schemes_enum::ADESaulyev)
    {
        heat_saulyev_scheme s_scheme(heat_coeff_holder, boundary_pair_, discretization_cfg_, grid_cfg_);
        s_scheme(solution, is_heat_sourse_set, heat_source, traverse_dir);
    }
    else
    {
        throw std::exception("Unreachable");
    }
}

void heat_equation_explicit_kernel<memory_space_enum::Host>::operator()(
    container_t &solution, bool is_heat_sourse_set, std::function<double(double, double)> const &heat_source,
    matrix_2d &solutions) {
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();
    // create a heat coefficient holder:
    auto const heat_coeff_holder = std::make_shared<heat_coefficients>(heat_data_cfg_, discretization_cfg_, double{});
    // Here make a dicision which explicit scheme to launch:
    if (solver_cfg_->explicit_pde_scheme() == explicit_pde_schemes_enum::Euler)
    {
        heat_euler_scheme euler_scheme(heat_coeff_holder, boundary_pair_, discretization_cfg_, grid_cfg_);
        euler_scheme(solution, is_heat_sourse_set, heat_source, traverse_dir, solutions);
    }
    else if (solver_cfg_->explicit_pde_scheme() == explicit_pde_schemes_enum::ADEBarakatClark)
    {
        heat_barakat_clark_scheme bc_scheme(heat_coeff_holder, boundary_pair_, discretization_cfg_, grid_cfg_);
        bc_scheme(solution, is_heat_sourse_set, heat_source, traverse_dir, solutions);
    }
    else if (solver_cfg_->explicit_pde_scheme() == explicit_pde_schemes_enum::ADESaulyev)
    {
        heat_saulyev_scheme s_scheme(heat_coeff_holder, boundary_pair_, discretization_cfg_, grid_cfg_);
        s_scheme(solution, is_heat_sourse_set, heat_source, traverse_dir, solutions);
    }
    else
    {
        throw std::exception("Unreachable");
    }
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
