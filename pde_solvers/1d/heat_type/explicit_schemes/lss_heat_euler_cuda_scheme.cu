#include "../solver_method/lss_heat_euler_cuda_solver_method.hpp"
#include "lss_heat_euler_cuda_scheme.hpp"

#include <device_launch_parameters.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#define THREADS_PER_BLOCK 256

namespace lss_pde_solvers
{

namespace one_dimensional
{


 bool heat_euler_cuda_scheme::is_stable(heat_coefficients_ptr const &coefficients) {

    const double zero = 0.0;
    const double two = 2.0;
    auto const &A = coefficients->A_;
    auto const &B = coefficients->B_;
    auto const &D = coefficients->D_;
    const double k = coefficients->k_;
    const double lambda = coefficients->lambda_;
    const double gamma = coefficients->gamma_;
    const double delta = coefficients->delta_;
    auto const &a = [=](double t, double x) {
        return ((A(t, x) + D(t, x)) / (two * lambda));
    };
    auto const &b = [=](double t, double x) {
        return ((D(t, x) - A(t, x)) / (two * gamma));
    };
    auto const &c = [=](double t, double x) {
        return ((lambda * a(t, x) - B(t, x)) / delta);
    };
    const std::size_t space_size = discretization_cfg_->number_of_space_points();
    auto const &ftime = discretization_cfg_->time_range()->upper();
    double x{}, t{k};
    while (t <= ftime) {
    for (std::size_t i = 0; i < space_size; ++i) {
        x = grid_1d::value(grid_cfg_, i);
        if (c(t, x) > zero) return false;
        if ((gamma * gamma * b(t, x) * b(t, x)) > (two * lambda * a(t, x)))
        return false;
    }
    t += k;
    }
    return true;
}

heat_euler_cuda_scheme::heat_euler_cuda_scheme(heat_coefficients_ptr const &coefficients, boundary_1d_pair const &boundary_pair,
    pde_discretization_config_1d_ptr const &discretization_config,grid_config_1d_ptr const &grid_config)
    : boundary_pair_{boundary_pair},
      discretization_cfg_{discretization_config},
      grid_cfg_{grid_config} {
  initialize(coefficients);
}

heat_euler_cuda_scheme:: ~heat_euler_cuda_scheme() {}



void heat_euler_cuda_scheme::operator()(container_t &solution, bool is_heat_sourse_set,
                                    std::function<double(double, double)> const &heat_source,
                                    traverse_direction_enum traverse_dir) {
  auto const &timer = discretization_cfg_->time_range();
  const double k = discretization_cfg_->time_step();
  // last time index:
  const std::size_t last_time_idx =
      discretization_cfg_->number_of_time_points() - 1;
  auto const &solver_method_ptr =
      std::make_shared<heat_euler_cuda_solver_method>(euler_coeffs_, grid_cfg_,
                                                      is_heat_sourse_set);
  if (is_heat_sourse_set) {
    explicit_time_loop::run(solver_method_ptr, boundary_pair_, timer,
                            last_time_idx, k,traverse_dir, heat_source, solution);
  } else {
    explicit_time_loop::run(solver_method_ptr, boundary_pair_, timer,
                            last_time_idx, k,traverse_dir, solution);
  }
}

void heat_euler_cuda_scheme::operator()(container_t &solution, bool is_heat_sourse_set,
                                std::function<double(double, double)> const &heat_source,
                                traverse_direction_enum traverse_dir, matrix_2d &solutions) {
  auto const &timer = discretization_cfg_->time_range();
  const double k = discretization_cfg_->time_step();
  // last time index:
  const std::size_t last_time_idx =
      discretization_cfg_->number_of_time_points() - 1;
  auto const &solver_method_ptr =
      std::make_shared<heat_euler_cuda_solver_method>(euler_coeffs_, grid_cfg_,
                                                      is_heat_sourse_set);
  if (is_heat_sourse_set) {
    explicit_time_loop::run_with_stepping(solver_method_ptr, boundary_pair_,
                                          timer,last_time_idx, k, traverse_dir, heat_source,
                                            solution, solutions);
  } else {
    explicit_time_loop::run_with_stepping(solver_method_ptr, boundary_pair_,
                                          timer,last_time_idx, k, traverse_dir, solution,
                                            solutions);
  }
}





void heat_euler_cuda_kernel::launch(double time, thrust::device_vector<double> const &input,
                                            thrust::device_vector<double> &solution)
{
    discretize_coefficients(time);
    const unsigned int threads_per_block = THREADS_PER_BLOCK;
    const unsigned int blocks_per_grid =
        static_cast<unsigned int>(solution.size() + threads_per_block - 1) / threads_per_block;
    double *raw_a = thrust::raw_pointer_cast(d_a_.data());
    double *raw_b = thrust::raw_pointer_cast(d_b_.data());
    double *raw_d = thrust::raw_pointer_cast(d_d_.data());
    const double *raw_input = thrust::raw_pointer_cast(input.data());
    double *raw_solution = thrust::raw_pointer_cast(solution.data());
    heat_core_kernel<<<threads_per_block, blocks_per_grid>>>(raw_a, raw_b, raw_d, raw_input, 
        raw_solution, solution.size());
}




void heat_euler_cuda_kernel::launch(double time, thrust::device_vector<double> const &input,
                                            thrust::device_vector<double> const &source,
                                            thrust::device_vector<double> &solution)
{
    discretize_coefficients(time);
    const unsigned int threads_per_block = THREADS_PER_BLOCK;
    const unsigned int blocks_per_grid =
        static_cast<unsigned int>(solution.size() + threads_per_block - 1) / threads_per_block;
    double *raw_a = thrust::raw_pointer_cast(d_a_.data());
    double *raw_b = thrust::raw_pointer_cast(d_b_.data());
    double *raw_d = thrust::raw_pointer_cast(d_d_.data());
    const double *raw_source = thrust::raw_pointer_cast(source.data());
    const double *raw_input = thrust::raw_pointer_cast(input.data());
    double *raw_solution = thrust::raw_pointer_cast(solution.data());
    heat_core_kernel<<<threads_per_block, blocks_per_grid>>>(
        raw_a, raw_b, raw_d, k_, raw_input, raw_source, raw_solution, solution.size());
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
