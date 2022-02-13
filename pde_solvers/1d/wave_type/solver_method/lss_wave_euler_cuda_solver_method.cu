#include "lss_wave_euler_cuda_solver_method.hpp"

#include <thrust/host_vector.h>
#include "../../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../../boundaries/lss_robin_boundary.hpp"
#include"../../../../common/lss_utility.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::dirichlet_boundary_1d;
using lss_boundary::neumann_boundary_1d;
using lss_boundary::robin_boundary_1d;
using lss_utility::copy;
using d_1d = discretization_1d;

__global__ void wave_core_kernel(double const *a_coeff, double const *b_coeff, double const *c_coeff,
                                 double const *d_coeff, double const *input_0, double const *input_1, double *solution,
                                 const std::size_t size)
{
    const std::size_t tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid == 0)
        return;
    if (tid >= (size - 1))
        return;

    solution[tid] = (a_coeff[tid] * input_1[tid - 1]) + (c_coeff[tid] * input_1[tid]) +
                    (b_coeff[tid] * input_1[tid + 1]) - (d_coeff[tid] * input_0[tid]);
}

__global__ void wave_core_kernel(double const *a_coeff, double const *b_coeff, double const *c_coeff,
                                 double const *d_coeff, double const *input_0, double const *input_1,
                                 double const *source, double *solution, const std::size_t size)
{
    const std::size_t tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid == 0)
        return;
    if (tid >= (size - 1))
        return;

    solution[tid] = (a_coeff[tid] * input_1[tid - 1]) + (c_coeff[tid] * input_1[tid]) +
                    (b_coeff[tid] * input_1[tid + 1]) - (d_coeff[tid] * input_0[tid]) + source[tid];
}

void wave_euler_cuda_kernel::initialize(wave_explicit_coefficients_ptr const &coefficients)
{
    const std::size_t space_size = coefficients->space_size_;
    k_ = coefficients->k_;
    a_ = coefficients->A_;
    b_ = coefficients->B_;
    c_ = coefficients->C_;
    d_ = coefficients->D_;
    h_a_.resize(space_size);
    h_b_.resize(space_size);
    h_c_.resize(space_size);
    h_d_.resize(space_size);
    d_a_.resize(space_size);
    d_b_.resize(space_size);
    d_c_.resize(space_size);
    d_d_.resize(space_size);
}

void wave_euler_cuda_kernel::discretize_coefficients(double time)
{
    // discretize on host
    d_1d::of_function(grid_cfg_, time, a_, h_a_);
    d_1d::of_function(grid_cfg_, time, b_, h_b_);
    d_1d::of_function(grid_cfg_, time, c_, h_c_);
    d_1d::of_function(grid_cfg_, time, d_, h_d_);
    // copy to host vector:
    thrust::host_vector<double> host_a(h_a_.size());
    thrust::host_vector<double> host_b(h_b_.size());
    thrust::host_vector<double> host_c(h_c_.size());
    thrust::host_vector<double> host_d(h_d_.size());
    copy(host_a, h_a_);
    copy(host_b, h_b_);
    copy(host_c, h_c_);
    copy(host_d, h_d_);
    thrust::copy(host_a.begin(), host_a.end(), d_a_.begin());
    thrust::copy(host_b.begin(), host_b.end(), d_b_.begin());
    thrust::copy(host_c.begin(), host_c.end(), d_c_.begin());
    thrust::copy(host_d.begin(), host_d.end(), d_d_.begin());
}

wave_euler_cuda_kernel::wave_euler_cuda_kernel(wave_explicit_coefficients_ptr const &coefficients,
                                               grid_config_1d_ptr const &grid_config)
    : grid_cfg_{grid_config}
{
    initialize(coefficients);
}

wave_euler_cuda_kernel::~wave_euler_cuda_kernel()
{
}

void explicit_wave_cuda_scheme::rhs(wave_explicit_coefficients_ptr const &cfs, wave_euler_cuda_kernel_ptr const &kernel,
                                    grid_config_1d_ptr const &grid_config, container_t const &input_0,
                                    container_t const &input_1, boundary_1d_pair const &boundary_pair,
                                    double const &time, container_t &solution)
{
    const double two = 2.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &a = cfs->A_;
    auto const &b = cfs->B_;
    auto const &c = cfs->C_;
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
        solution[0] = beta * a(time, x) + c(time, x) * input_1[0] + (a(time, x) + b(time, x)) * input_1[1] -
                      d(time, x) * input_0[0];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] = beta * a(time, x) + (c(time, x) + alpha * a(time, x)) * input_1[0] +
                      (a(time, x) + b(time, x)) * input_1[1] - d(time, x) * input_0[0];
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
        solution[N] = (a(time, x) + b(time, x)) * input_1[N - 1] + c(time, x) * input_1[N] - delta * b(time, x) -
                      d(time, x) * input_0[N];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (a(time, x) + b(time, x)) * input_1[N - 1] + (c(time, x) - gamma * b(time, x)) * input_1[N] -
                      delta * b(time, x) - d(time, x) * input_0[N];
    }

    thrust::host_vector<double> h_input_0(input_0.size());
    thrust::host_vector<double> h_input_1(input_1.size());
    thrust::host_vector<double> h_solution(solution.size());
    copy(h_input_0, input_0);
    copy(h_input_1, input_1);
    copy(h_solution, solution);
    thrust::device_vector<double> d_input_0(h_input_0);
    thrust::device_vector<double> d_input_1(h_input_1);
    thrust::device_vector<double> d_solution(h_solution);
    kernel->launch(time, d_input_0, d_input_1, d_solution);
    thrust::copy(d_solution.begin(), d_solution.end(), h_solution.begin());
    copy(solution, h_solution);
}

void explicit_wave_cuda_scheme::rhs_source(wave_explicit_coefficients_ptr const &cfs,
                                           wave_euler_cuda_kernel_ptr const &kernel,
                                           grid_config_1d_ptr const &grid_config, container_t const &input_0,
                                           container_t const &input_1, container_t const &inhom_input,
                                           boundary_1d_pair const &boundary_pair, double const &time,
                                           container_t &solution)
{
    const double two = 2.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &a = cfs->A_;
    auto const &b = cfs->B_;
    auto const &c = cfs->C_;
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
        solution[0] = beta * a(time, x) + c(time, x) * input_1[0] + (a(time, x) + b(time, x)) * input_1[1] -
                      d(time, x) * input_0[0] + inhom_input[0];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] = beta * a(time, x) + (c(time, x) + alpha * a(time, x)) * input_1[0] +
                      (a(time, x) + b(time, x)) * input_1[1] - d(time, x) * input_0[0] + inhom_input[0];
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
        solution[N] = (a(time, x) + b(time, x)) * input_1[N - 1] + c(time, x) * input_1[N] - delta * b(time, x) -
                      d(time, x) * input_0[N] + inhom_input[N];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (a(time, x) + b(time, x)) * input_1[N - 1] + (c(time, x) - gamma * b(time, x)) * input_1[N] -
                      delta * b(time, x) - d(time, x) * input_0[N] + inhom_input[N];
    }
    thrust::host_vector<double> h_input_0(input_0.size());
    thrust::host_vector<double> h_input_1(input_1.size());
    thrust::host_vector<double> h_inhom_input(inhom_input.size());
    thrust::host_vector<double> h_solution(solution.size());
    copy(h_input_0, input_0);
    copy(h_input_1, input_1);
    copy(h_inhom_input, inhom_input);
    copy(h_solution, solution);
    thrust::device_vector<double> d_input_0(h_input_0);
    thrust::device_vector<double> d_input_1(h_input_1);
    thrust::device_vector<double> d_source(h_inhom_input);
    thrust::device_vector<double> d_solution(h_solution);
    kernel->launch(time, d_input_0, d_input_1, d_source, d_solution);
    thrust::copy(d_solution.begin(), d_solution.end(), h_solution.begin());
    copy(solution, h_solution);
}

void wave_euler_cuda_solver_method::initialize(bool is_wave_source_set)
{
    kernel_ = std::make_shared<wave_euler_cuda_kernel>(coefficients_, grid_cfg_);
    if (is_wave_source_set)
    {
        source_.resize(coefficients_->space_size_);
    }
}

wave_euler_cuda_solver_method::wave_euler_cuda_solver_method(wave_explicit_coefficients_ptr const &coefficients,
                                                             grid_config_1d_ptr const &grid_config,
                                                             bool is_wave_source_set)
    : coefficients_{coefficients}, wave_explicit_solver_method(grid_config)
{
    initialize(is_wave_source_set);
}

wave_euler_cuda_solver_method::~wave_euler_cuda_solver_method()
{
}

void wave_euler_cuda_solver_method::solve_initial(container_t &prev_solution_0, container_t &prev_solution_1,
                                                  boundary_1d_pair const &boundary_pair, double const &time,
                                                  double const &next_time, container_t &solution)
{
    explicit_wave_scheme::rhs_initial(coefficients_, grid_cfg_, prev_solution_0, prev_solution_1, boundary_pair, time,
                                      solution);
}

void wave_euler_cuda_solver_method::solve_initial(container_t &prev_solution_0, container_t &prev_solution_1,
                                                  boundary_1d_pair const &boundary_pair, double const &time,
                                                  double const &next_time,
                                                  std::function<double(double, double)> const &wave_source,
                                                  container_t &solution)
{
    // get the right-hand side of the scheme:
    d_1d::of_function(grid_cfg_, time, wave_source, source_);
    explicit_wave_scheme::rhs_initial_source(coefficients_, grid_cfg_, prev_solution_0, prev_solution_1, source_,
                                             boundary_pair, time, solution);
}

void wave_euler_cuda_solver_method::solve_terminal(container_t &prev_solution_0, container_t &prev_solution_1,
                                                   boundary_1d_pair const &boundary_pair, double const &time,
                                                   double const &next_time, container_t &solution)
{
    // get the right-hand side of the scheme:
    explicit_wave_scheme::rhs_terminal(coefficients_, grid_cfg_, prev_solution_0, prev_solution_1, boundary_pair, time,
                                       solution);
}

void wave_euler_cuda_solver_method::solve_terminal(container_t &prev_solution_0, container_t &prev_solution_1,
                                                   boundary_1d_pair const &boundary_pair, double const &time,
                                                   double const &next_time,
                                                   std::function<double(double, double)> const &wave_source,
                                                   container_t &solution)
{
    // get the right-hand side of the scheme:
    d_1d::of_function(grid_cfg_, time, wave_source, source_);
    explicit_wave_scheme::rhs_terminal_source(coefficients_, grid_cfg_, prev_solution_0, prev_solution_1, source_,
                                              boundary_pair, time, solution);
}

void wave_euler_cuda_solver_method::solve(container_t &prev_solution_0, container_t &prev_solution_1,
                                          boundary_1d_pair const &boundary_pair, double const &time,
                                          double const &next_time, container_t &solution)
{
    // get the right-hand side of the scheme:
    explicit_wave_cuda_scheme::rhs(coefficients_, kernel_, grid_cfg_, prev_solution_0, prev_solution_1, boundary_pair,
                                   time, solution);
}

void wave_euler_cuda_solver_method::solve(container_t &prev_solution_0, container_t &prev_solution_1,
                                          boundary_1d_pair const &boundary_pair, double const &time,
                                          double const &next_time,
                                          std::function<double(double, double)> const &wave_source,
                                          container_t &solution)
{
    // get the right-hand side of the scheme:
    d_1d::of_function(grid_cfg_, time, wave_source, source_);
    explicit_wave_cuda_scheme::rhs_source(coefficients_, kernel_, grid_cfg_, prev_solution_0, prev_solution_1, source_,
                                          boundary_pair, time, solution);
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
