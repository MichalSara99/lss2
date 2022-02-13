/**

    @file      lss_wave_euler_cuda_solver_method.hpp
    @brief     CUDA Euler solver method for wave problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_WAVE_EULER_CUDA_SOLVER_METHOD_HPP_)
#define _LSS_WAVE_EULER_CUDA_SOLVER_METHOD_HPP_

#include <device_launch_parameters.h>
#include <thrust/device_vector.h>

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../containers/lss_matrix_2d.hpp"
#include "../../../../discretization/lss_discretization.hpp"
#include "../../../../discretization/lss_grid.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../../../lss_pde_discretization_config.hpp"
#include "../../../lss_wave_data_config.hpp"
#include "../explicit_coefficients/lss_wave_explicit_coefficients.hpp"
#include "lss_wave_euler_solver_method.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::boundary_1d_pair;
using lss_containers::matrix_2d;
using lss_utility::container_t;

__global__ extern void wave_core_kernel(double const *a_coeff, double const *b_coeff, double const *c_coeff,
                                        double const *d_coeff, double const *input_0, double const *input_1,
                                        double *solution, const std::size_t size);

__global__ extern void wave_core_kernel(double const *a_coeff, double const *b_coeff, double const *c_coeff,
                                        double const *d_coeff, double const *input_0, double const *input_1,
                                        double const *source, double *solution, const std::size_t size);

/**
 * wave_euler_cuda_kernel object
 */
class wave_euler_cuda_kernel
{

  private:
    double k_;
    thrust::device_vector<double> d_a_;
    thrust::device_vector<double> d_b_;
    thrust::device_vector<double> d_c_;
    thrust::device_vector<double> d_d_;
    container_t h_a_;
    container_t h_b_;
    container_t h_c_;
    container_t h_d_;
    grid_config_1d_ptr grid_cfg_;
    // coefficients:
    std::function<double(double, double)> a_;
    std::function<double(double, double)> b_;
    std::function<double(double, double)> c_;
    std::function<double(double, double)> d_;

    void initialize(wave_explicit_coefficients_ptr const &coefficients);

    void discretize_coefficients(double time);

  public:
    explicit wave_euler_cuda_kernel(wave_explicit_coefficients_ptr const &coefficients,
                                    grid_config_1d_ptr const &grid_config);

    ~wave_euler_cuda_kernel();

    void launch(double time, thrust::device_vector<double> const &input_0, thrust::device_vector<double> const &input_1,
                thrust::device_vector<double> &solution);

    void launch(double time, thrust::device_vector<double> const &input_0, thrust::device_vector<double> const &input_1,
                thrust::device_vector<double> const &source, thrust::device_vector<double> &solution);
};

using wave_euler_cuda_kernel_ptr = sptr_t<wave_euler_cuda_kernel>;

/**
 * explicit_wave_cuda_scheme object
 */
class explicit_wave_cuda_scheme
{

  public:
    static void rhs(wave_explicit_coefficients_ptr const &cfs, wave_euler_cuda_kernel_ptr const &kernel,
                    grid_config_1d_ptr const &grid_config, container_t const &input_0, container_t const &input_1,
                    boundary_1d_pair const &boundary_pair, double const &time, container_t &solution);

    static void rhs_source(wave_explicit_coefficients_ptr const &cfs, wave_euler_cuda_kernel_ptr const &kernel,
                           grid_config_1d_ptr const &grid_config, container_t const &input_0,
                           container_t const &input_1, container_t const &inhom_input,
                           boundary_1d_pair const &boundary_pair, double const &time, container_t &solution);
};

/**
wave_euler_cuda_solver_method object
*/
class wave_euler_cuda_solver_method : public wave_explicit_solver_method
{
  private:
    // scheme coefficients:
    wave_explicit_coefficients_ptr coefficients_;
    // cuda kernel:
    wave_euler_cuda_kernel_ptr kernel_;
    container_t source_;

    explicit wave_euler_cuda_solver_method() = delete;

    void initialize(bool is_wave_source_set);

  public:
    explicit wave_euler_cuda_solver_method(wave_explicit_coefficients_ptr const &coefficients,
                                           grid_config_1d_ptr const &grid_config, bool is_wave_source_set);

    ~wave_euler_cuda_solver_method();

    wave_euler_cuda_solver_method(wave_euler_cuda_solver_method const &) = delete;
    wave_euler_cuda_solver_method(wave_euler_cuda_solver_method &&) = delete;
    wave_euler_cuda_solver_method &operator=(wave_euler_cuda_solver_method const &) = delete;
    wave_euler_cuda_solver_method &operator=(wave_euler_cuda_solver_method &&) = delete;

    void solve_initial(container_t &prev_solution_0, container_t &prev_solution_1,
                       boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                       container_t &solution);

    void solve_initial(container_t &prev_solution_0, container_t &prev_solution_1,
                       boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                       std::function<double(double, double)> const &wave_source, container_t &solution);

    void solve_terminal(container_t &prev_solution_0, container_t &prev_solution_1,
                        boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                        container_t &solution);

    void solve_terminal(container_t &prev_solution_0, container_t &prev_solution_1,
                        boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                        std::function<double(double, double)> const &wave_source, container_t &solution);

    void solve(container_t &prev_solution_0, container_t &prev_solution_1, boundary_1d_pair const &boundary_pair,
               double const &time, double const &next_time, container_t &solution);

    void solve(container_t &prev_solution_0, container_t &prev_solution_1, boundary_1d_pair const &boundary_pair,
               double const &time, double const &next_time, std::function<double(double, double)> const &wave_source,
               container_t &solution);
};

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_WAVE_EULER_CUDA_SOLVER_METHOD_HPP_
