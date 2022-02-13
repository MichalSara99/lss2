#include "lss_wave_equation_explicit_kernel.hpp"

namespace lss_pde_solvers
{
namespace one_dimensional
{

using lss_enumerations::traverse_direction_enum;

wave_equation_explicit_kernel<memory_space_enum::Device>::wave_equation_explicit_kernel(
    boundary_1d_pair const &boundary_pair, wave_data_transform_1d_ptr const &wave_data_config,
    pde_discretization_config_1d_ptr const &discretization_config, wave_explicit_solver_config_ptr const &solver_config,
    grid_config_1d_ptr const &grid_config)
    : boundary_pair_{boundary_pair}, wave_data_cfg_{wave_data_config}, discretization_cfg_{discretization_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void wave_equation_explicit_kernel<memory_space_enum::Device>::operator()(
    container_t &prev_solution_0, container_t &prev_solution_1, container_t &next_solution, bool is_wave_sourse_set,
    std::function<double(double, double)> const &wave_source)
{
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();

    // create a heat coefficient holder:
    auto wave_coeff_holder = std::make_shared<wave_explicit_coefficients>(wave_data_cfg_, discretization_cfg_);
    // get the modified wave source:
    auto const mod_wave_source =
        (is_wave_sourse_set == true) ? wave_coeff_holder->modified_wave_source(wave_source) : nullptr;

    // Here we have only Euler discretization available:
    wave_euler_cuda_scheme euler_scheme(wave_coeff_holder, boundary_pair_, discretization_cfg_, grid_cfg_);
    euler_scheme(prev_solution_0, prev_solution_1, next_solution, is_wave_sourse_set, mod_wave_source, traverse_dir);
}

void wave_equation_explicit_kernel<memory_space_enum::Device>::operator()(
    container_t &prev_solution_0, container_t &prev_solution_1, container_t &next_solution, bool is_wave_sourse_set,
    std::function<double(double, double)> const &wave_source, matrix_2d &solutions)
{
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();

    // create a heat coefficient holder:
    auto const wave_coeff_holder = std::make_shared<wave_explicit_coefficients>(wave_data_cfg_, discretization_cfg_);
    // get the modified wave source:
    auto const mod_wave_source =
        (is_wave_sourse_set == true) ? wave_coeff_holder->modified_wave_source(wave_source) : nullptr;

    // Here we have only Euler discretization available:
    wave_euler_cuda_scheme euler_scheme(wave_coeff_holder, boundary_pair_, discretization_cfg_, grid_cfg_);
    euler_scheme(prev_solution_0, prev_solution_1, next_solution, is_wave_sourse_set, mod_wave_source, traverse_dir,
                 solutions);
}

wave_equation_explicit_kernel<memory_space_enum::Host>::wave_equation_explicit_kernel(
    boundary_1d_pair const &boundary_pair, wave_data_transform_1d_ptr const &wave_data_config,
    pde_discretization_config_1d_ptr const &discretization_config, wave_explicit_solver_config_ptr const &solver_config,
    grid_config_1d_ptr const &grid_config)
    : boundary_pair_{boundary_pair}, wave_data_cfg_{wave_data_config}, discretization_cfg_{discretization_config},
      solver_cfg_{solver_config}, grid_cfg_{grid_config}
{
}

void wave_equation_explicit_kernel<memory_space_enum::Host>::operator()(
    container_t &prev_solution_0, container_t &prev_solution_1, container_t &next_solution, bool is_wave_sourse_set,
    std::function<double(double, double)> const &wave_source)
{
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();

    // create a heat coefficient holder:
    auto const wave_coeff_holder = std::make_shared<wave_explicit_coefficients>(wave_data_cfg_, discretization_cfg_);

    auto const mod_wave_source =
        (is_wave_sourse_set == true) ? wave_coeff_holder->modified_wave_source(wave_source) : nullptr;

    // Here make a dicision which explicit scheme to launch:
    wave_euler_scheme euler_scheme(wave_coeff_holder, boundary_pair_, discretization_cfg_, grid_cfg_);
    euler_scheme(prev_solution_0, prev_solution_1, next_solution, is_wave_sourse_set, mod_wave_source, traverse_dir);
}

void wave_equation_explicit_kernel<memory_space_enum::Host>::operator()(
    container_t &prev_solution_0, container_t &prev_solution_1, container_t &next_solution, bool is_wave_sourse_set,
    std::function<double(double, double)> const &wave_source, matrix_2d &solutions)
{
    // save traverse_direction
    const traverse_direction_enum traverse_dir = solver_cfg_->traverse_direction();

    // create a heat coefficient holder:
    auto const wave_coeff_holder = std::make_shared<wave_explicit_coefficients>(wave_data_cfg_, discretization_cfg_);
    auto const mod_wave_source =
        (is_wave_sourse_set == true) ? wave_coeff_holder->modified_wave_source(wave_source) : nullptr;
    // Here make a dicision which explicit scheme to launch:
    wave_euler_scheme euler_scheme(wave_coeff_holder, boundary_pair_, discretization_cfg_, grid_cfg_);
    euler_scheme(prev_solution_0, prev_solution_1, next_solution, is_wave_sourse_set, mod_wave_source, traverse_dir,
                 solutions);
}

} // namespace one_dimensional
} // namespace lss_pde_solvers
