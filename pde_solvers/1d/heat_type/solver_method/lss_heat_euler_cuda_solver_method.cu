#include "lss_heat_euler_cuda_solver_method.hpp"

#include "../../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../../boundaries/lss_robin_boundary.hpp"
#include <thrust/host_vector.h>

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::dirichlet_boundary_1d;
using lss_boundary::neumann_boundary_1d;
using lss_boundary::robin_boundary_1d;
using d_1d = discretization_1d;
using lss_utility::copy;

__global__ void heat_core_kernel(double const *a_coeff, double const *b_coeff, double const *d_coeff, double const *input,
                      double *solution, const std::size_t size)
{
    const std::size_t tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid == 0)
        return;
    if (tid >= (size - 1))
        return;
    solution[tid] = (d_coeff[tid] * input[tid + 1]) + (b_coeff[tid] * input[tid]) + (a_coeff[tid] * input[tid - 1]);
}

__global__ void heat_core_kernel(double const *a_coeff, double const *b_coeff,
                                double const *d_coeff, double time_step,
                      double const *input, double const *source, double *solution, const std::size_t size)
{
    const std::size_t tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid == 0)
        return;
    if (tid >= (size - 1))
        return;
    solution[tid] = (d_coeff[tid] * input[tid + 1]) + (b_coeff[tid] * input[tid]) + (a_coeff[tid] * input[tid - 1]) +
                    (time_step * source[tid]);
}


void heat_euler_cuda_kernel::initialize(heat_euler_coefficients_ptr const &coefficients) {
      const std::size_t space_size = coefficients->space_size_;
      k_ = coefficients->k_;
      a_ = coefficients->A_;
      b_ = coefficients->B_;
      d_ = coefficients->D_;
      h_a_.resize(space_size);
      h_b_.resize(space_size);
      h_d_.resize(space_size);
      d_a_.resize(space_size);
      d_b_.resize(space_size);
      d_d_.resize(space_size);
}


void heat_euler_cuda_kernel::discretize_coefficients(double time)
{
    // discretize on host
    d_1d::of_function(grid_cfg_, time, a_, h_a_);
    d_1d::of_function(grid_cfg_, time, b_, h_b_);
    d_1d::of_function(grid_cfg_, time, d_, h_d_);
    // copy to host vector:
    thrust::host_vector<double> host_a(h_a_.size());
    thrust::host_vector<double> host_b(h_b_.size());
    thrust::host_vector<double> host_d(h_d_.size());
    copy(host_a, h_a_);
    copy(host_b, h_b_);
    copy(host_d, h_d_);
    // copy to device vector:
    thrust::copy(host_a.begin(), host_a.end(), d_a_.begin());
    thrust::copy(host_b.begin(), host_b.end(), d_b_.begin());
    thrust::copy(host_d.begin(), host_d.end(), d_d_.begin());
}

heat_euler_cuda_kernel::heat_euler_cuda_kernel(heat_euler_coefficients_ptr const &coefficients,
                                               grid_config_1d_ptr const &grid_config)
    : grid_cfg_{grid_config}
{
    initialize(coefficients);
}

heat_euler_cuda_kernel::~heat_euler_cuda_kernel()
{
}

void explicit_cuda_scheme::rhs(heat_euler_coefficients_ptr const &cfs, heat_euler_cuda_kernel_ptr const &kernel,
                               grid_config_1d_ptr const &grid_config,container_t const &input,
                               boundary_1d_pair const &boundary_pair, double const &time, container_t &solution) {
    const double two = 2.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &a = cfs->A_;
    auto const &b = cfs->B_;
    auto const &d = cfs->D_;
    auto const h = grid_1d::step(grid_config);
    double x{};
    // for lower boundaries first:
    x = grid_1d::value(grid_config, 0);
    if (auto const &ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(first_bnd))
    {
        solution[0] = ptr->value(time);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        solution[0] = beta * a(time, x) + b(time, x) * input[0] + (a(time, x) + d(time, x)) * input[1];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] =
            (b(time, x) + alpha * a(time, x)) * input[0] + (a(time, x) + d(time, x)) * input[1] + beta * a(time, x);
    }
    // for upper boundaries second:
    const std::size_t N = solution.size() - 1;
    x = grid_1d::value(grid_config, N);
    if (auto const &ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(second_bnd))
    {
        solution[N] = ptr->value(time);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        solution[N] = (a(time, x) + d(time, x)) * input[N - 1] + b(time, x) * input[N] - delta * d(time, x);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (a(time, x) + d(time, x)) * input[N - 1] + (b(time, x) - gamma * d(time, x)) * input[N] -
                      delta * d(time, x);
    }
    // light-weight object with cuda kernel computing the solution:
    thrust::host_vector<double> h_input(input.size());
    thrust::host_vector<double> h_solution(solution.size());
    copy(h_input, input);
    copy(h_solution, solution);
    thrust::device_vector<double> d_input(h_input);
    thrust::device_vector<double> d_solution(h_solution);
    kernel->launch(time, d_input, d_solution);
    thrust::copy(d_solution.begin(), d_solution.end(), h_solution.begin());
    copy(solution, h_solution);
}

void explicit_cuda_scheme::rhs_source(heat_euler_coefficients_ptr const &cfs, heat_euler_cuda_kernel_ptr const &kernel,
                                      grid_config_1d_ptr const &grid_config, container_t const &input,
                                      container_t const &inhom_input,boundary_1d_pair const &boundary_pair, double const &time,
                                      container_t &solution) {
    const double two = 2.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &a = cfs->A_;
    auto const &b = cfs->B_;
    auto const &d = cfs->D_;
    auto const k = cfs->k_;
    auto const h = grid_1d::step(grid_config);
    double x{};

    // for lower boundaries first:
    x = grid_1d::value(grid_config, 0);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        solution[0] =
            beta * a(time, x) + b(time, x) * input[0] + (a(time, x) + d(time, x)) * input[1] + k * inhom_input[0];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] = (b(time, x) + alpha * a(time, x)) * input[0] + (a(time, x) + d(time, x)) * input[1] +
                      beta * a(time, x) + k * inhom_input[0];
    }
    // for upper boundaries second:
    const std::size_t N = solution.size() - 1;
    x = grid_1d::value(grid_config, N);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        solution[N] =
            (a(time, x) + d(time, x)) * input[N - 1] + b(time, x) * input[N] - delta * d(time, x) + k * inhom_input[N];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (a(time, x) + d(time, x)) * input[N - 1] + (b(time, x) - gamma * d(time, x)) * input[N] -
                      delta * d(time, x) + k * inhom_input[N];
    }
    // light-weight object with cuda kernel computing the solution:
    thrust::host_vector<double> h_input(input.size());
    thrust::host_vector<double> h_inhom_input(inhom_input.size());
    thrust::host_vector<double> h_solution(solution.size());
    copy(h_input, input);
    copy(h_inhom_input, inhom_input);
    copy(h_solution, solution);
    thrust::device_vector<double> d_input(h_input);
    thrust::device_vector<double> d_inhom_input(h_inhom_input);
    thrust::device_vector<double> d_solution(h_solution);
    kernel->launch(time, d_input, d_inhom_input, d_solution);
    thrust::copy(d_solution.begin(), d_solution.end(), h_solution.begin());
    copy(solution, h_solution);
}

void heat_euler_cuda_solver_method::initialize(bool is_heat_source_set)
{
    kernel_ = std::make_shared<heat_euler_cuda_kernel>(coefficients_, grid_cfg_);
    if (is_heat_source_set)
    {
        source_.resize(coefficients_->space_size_);
    }
}

heat_euler_cuda_solver_method::heat_euler_cuda_solver_method(heat_euler_coefficients_ptr const &coefficients,
                                                             grid_config_1d_ptr const &grid_config,
                                                             bool is_heat_source_set)
    : coefficients_{coefficients},heat_explicit_solver_method(grid_config) {
  initialize(is_heat_source_set);
}

heat_euler_cuda_solver_method::~heat_euler_cuda_solver_method()
{
}

void heat_euler_cuda_solver_method::solve(container_t &prev_solution,
    boundary_1d_pair const &boundary_pair, std::size_t const &time_idx,
    double const &time, container_t &solution) {

    explicit_cuda_scheme::rhs(coefficients_, kernel_, grid_cfg_, prev_solution, boundary_pair, time, solution);
}

void heat_euler_cuda_solver_method::solve(
    container_t &prev_solution,
    boundary_1d_pair const &boundary_pair, std::size_t const &time_idx,
    double const &time, double const &next_time,
    std::function<double(double, double)> const &heat_source,
    container_t &solution) {

    d_1d::of_function(grid_cfg_, time, heat_source, source_);
    explicit_cuda_scheme::rhs_source(coefficients_, kernel_, grid_cfg_, prev_solution, source_, boundary_pair, time,
                                     solution);
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
