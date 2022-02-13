/**

    @file      lss_heat_euler_cuda_solver_method.hpp
    @brief     Euler CUDA solver method for heat problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_EULER_CUDA_SOLVER_METHOD_HPP_)
#define _LSS_HEAT_EULER_CUDA_SOLVER_METHOD_HPP_

#include <device_launch_parameters.h>
#include <thrust/device_vector.h>

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../discretization/lss_discretization.hpp"
#include "../../../../discretization/lss_grid.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../explicit_coefficients/lss_heat_euler_coefficients.hpp"
#include "lss_heat_explicit_solver_method.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_1d_ptr;
using lss_utility::container_t;
using lss_utility::sptr_t;

__global__ extern void heat_core_kernel(double const *a_coeff, double const *b_coeff, double const *d_coeff,
                                        double const *input, double *solution, const std::size_t size);

__global__ extern void heat_core_kernel(double const *a_coeff, double const *b_coeff, double const *d_coeff,
                                        double time_step, double const *input, double const *source, double *solution,
                                        const std::size_t size);

/**
 * heat_euler_cuda_kernel object
 */
class heat_euler_cuda_kernel
{

  private:
    double k_;
    thrust::device_vector<double> d_a_;
    thrust::device_vector<double> d_b_;
    thrust::device_vector<double> d_d_;
    container_t h_a_;
    container_t h_b_;
    container_t h_d_;
    grid_config_1d_ptr grid_cfg_;
    // coefficients:
    std::function<double(double, double)> a_;
    std::function<double(double, double)> b_;
    std::function<double(double, double)> d_;

    void initialize(heat_euler_coefficients_ptr const &coefficients);

    void discretize_coefficients(double time);

  public:
    explicit heat_euler_cuda_kernel(heat_euler_coefficients_ptr const &coefficients,
                                    grid_config_1d_ptr const &grid_config);

    ~heat_euler_cuda_kernel();

    void launch(double time, thrust::device_vector<double> const &input, thrust::device_vector<double> &solution);

    void launch(double time, thrust::device_vector<double> const &input, thrust::device_vector<double> const &source,
                thrust::device_vector<double> &solution);
};

using heat_euler_cuda_kernel_ptr = sptr_t<heat_euler_cuda_kernel>;

/**
 * explicit_cuda_scheme object
 */
class explicit_cuda_scheme
{
  public:
    static void rhs(heat_euler_coefficients_ptr const &cfs, heat_euler_cuda_kernel_ptr const &kernel,
                    grid_config_1d_ptr const &grid_config, container_t const &input,
                    boundary_1d_pair const &boundary_pair, double const &time, container_t &solution);

    static void rhs_source(heat_euler_coefficients_ptr const &cfs, heat_euler_cuda_kernel_ptr const &kernel,
                           grid_config_1d_ptr const &grid_config, container_t const &input,
                           container_t const &inhom_input, boundary_1d_pair const &boundary_pair, double const &time,
                           container_t &solution);
};

/**
heat_euler_cuda_solver_method object
*/
class heat_euler_cuda_solver_method : public heat_explicit_solver_method
{

  private:
    // scheme coefficients:
    heat_euler_coefficients_ptr coefficients_;
    // cuda kernel:
    heat_euler_cuda_kernel_ptr kernel_;
    // containers:
    container_t source_;

    explicit heat_euler_cuda_solver_method() = delete;

    void initialize(bool is_heat_source_set);

  public:
    explicit heat_euler_cuda_solver_method(heat_euler_coefficients_ptr const &coefficients,
                                           grid_config_1d_ptr const &grid_config, bool is_heat_source_set);

    ~heat_euler_cuda_solver_method();

    heat_euler_cuda_solver_method(heat_euler_cuda_solver_method const &) = delete;
    heat_euler_cuda_solver_method(heat_euler_cuda_solver_method &&) = delete;
    heat_euler_cuda_solver_method &operator=(heat_euler_cuda_solver_method const &) = delete;
    heat_euler_cuda_solver_method &operator=(heat_euler_cuda_solver_method &&) = delete;

    void solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair, std::size_t const &time_idx,
               double const &time, container_t &solution) override;

    void solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair, std::size_t const &time_idx,
               double const &time, double const &next_time, std::function<double(double, double)> const &heat_source,
               container_t &solution) override;
};

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_EULER_CUDA_SOLVER_METHOD_HPP_
