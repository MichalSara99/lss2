/**

    @file      lss_heat_euler_cuda_solver_method_2d.hpp
    @brief     Explicit CUDA Euler method for 2D heat problems
    @details   ~
    @author    Michal Sara
    @date      14.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_EULER_CUDA_SOLVER_METHOD_2D_HPP_)
#define _LSS_HEAT_EULER_CUDA_SOLVER_METHOD_2D_HPP_

#include <device_launch_parameters.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../containers/lss_matrix_2d.hpp"
#include "../../../../discretization/lss_discretization.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../explicit_coefficients/lss_heat_euler_coefficients_2d.hpp"
#include "lss_heat_explicit_solver_method_2d.hpp"

namespace lss_pde_solvers
{
namespace two_dimensional
{

using lss_containers::matrix_2d;
using lss_containers::matrix_2d_ptr;
using lss_grids::grid_config_2d_ptr;
using lss_utility::container_t;
using lss_utility::range;
using lss_utility::sptr_t;

using d_2d = discretization_2d;

__global__ extern void heat_core_kernel_2d(double const *m_coeff, double const *m_tilde_coeff, double const *p_coeff,
                                           double const *p_tilde_coeff, double const *c_coeff, double const *z_coeff,
                                           double const *w_coeff, double const *input, double *solution,
                                           const std::size_t size_x, const std::size_t size_y);

__global__ extern void heat_core_kernel_2d(double const *m_coeff, double const *m_tilde_coeff, double const *p_coeff,
                                           double const *p_tilde_coeff, double const *c_coeff, double const *z_coeff,
                                           double const *w_coeff, double time_step, double const *input,
                                           double const *source, double *solution, const std::size_t size_x,
                                           const std::size_t size_y);

/**
 * heat_euler_cuda_kernel_2d object
 */
class heat_euler_cuda_kernel_2d
{

  private:
    double k_;
    std::size_t size_x_, size_y_;
    // device containers:
    thrust::device_vector<double> d_m_;
    thrust::device_vector<double> d_m_tilde_;
    thrust::device_vector<double> d_p_;
    thrust::device_vector<double> d_p_tilde_;
    thrust::device_vector<double> d_c_;
    thrust::device_vector<double> d_z_;
    thrust::device_vector<double> d_w_;
    // host containers:
    container_t h_m_;
    container_t h_m_tilde_;
    container_t h_p_;
    container_t h_p_tilde_;
    container_t h_c_;
    container_t h_z_;
    container_t h_w_;
    grid_config_2d_ptr grid_cfg_;
    // coefficients:
    std::function<double(double, double, double)> m_;
    std::function<double(double, double, double)> m_tilde_;
    std::function<double(double, double, double)> p_;
    std::function<double(double, double, double)> p_tilde_;
    std::function<double(double, double, double)> c_;
    std::function<double(double, double, double)> z_;
    std::function<double(double, double, double)> w_;

    void initialize(heat_euler_coefficients_2d_ptr const &coefficients);

    void discretize_coefficients(double time);

  public:
    explicit heat_euler_cuda_kernel_2d(heat_euler_coefficients_2d_ptr const &coefficients,
                                       grid_config_2d_ptr const &grid_config);

    ~heat_euler_cuda_kernel_2d();

    void launch(double time, thrust::device_vector<double> const &input, thrust::device_vector<double> &solution);

    void launch(double time, thrust::device_vector<double> const &input, thrust::device_vector<double> const &source,
                thrust::device_vector<double> &solution);
};

using heat_euler_cuda_kernel_2d_ptr = sptr_t<heat_euler_cuda_kernel_2d>;

class explicit_euler_cuda_rhs
{

  public:
    static void rhs(heat_euler_coefficients_2d_ptr const &cfs, heat_euler_cuda_kernel_2d_ptr const &kernel,
                    grid_config_2d_ptr const &grid_config, matrix_2d const &input, double const &time,
                    matrix_2d &solution);

    static void rhs_source(heat_euler_coefficients_2d_ptr const &cfs, heat_euler_cuda_kernel_2d_ptr const &kernel,
                           grid_config_2d_ptr const &grid_config, matrix_2d const &input, double const &time,
                           matrix_2d const &inhom_input, matrix_2d &solution);
};

/**
    heat_euler_cuda_solver_method_2d object
*/

class heat_euler_cuda_solver_method_2d : public heat_explicit_solver_method_2d
{

  private:
    // scheme coefficients:
    heat_euler_coefficients_2d_ptr coefficients_;
    // cuda kernel:
    heat_euler_cuda_kernel_2d_ptr kernel_;
    matrix_2d_ptr source_;

    explicit heat_euler_cuda_solver_method_2d() = delete;

    void initialize(bool is_heat_source_set);

  public:
    explicit heat_euler_cuda_solver_method_2d(heat_euler_coefficients_2d_ptr const &coefficients,
                                              grid_config_2d_ptr const &grid_config, bool is_heat_source_set);

    ~heat_euler_cuda_solver_method_2d();

    heat_euler_cuda_solver_method_2d(heat_euler_cuda_solver_method_2d const &) = delete;
    heat_euler_cuda_solver_method_2d(heat_euler_cuda_solver_method_2d &&) = delete;
    heat_euler_cuda_solver_method_2d &operator=(heat_euler_cuda_solver_method_2d const &) = delete;
    heat_euler_cuda_solver_method_2d &operator=(heat_euler_cuda_solver_method_2d &&) = delete;

    void solve(matrix_2d &prev_solution, double const &time, matrix_2d &solution) override;

    void solve(matrix_2d &prev_solution, double const &time,
               std::function<double(double, double, double)> const &heat_source, matrix_2d &solution) override;
};

} // namespace two_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_EULER_CUDA_SOLVER_METHOD_2D_HPP_
