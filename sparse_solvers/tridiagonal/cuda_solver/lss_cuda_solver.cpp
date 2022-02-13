#include "lss_cuda_solver.hpp"

#include "../../../common/lss_utility.hpp"
#include "../../../containers/lss_flat_matrix.hpp"

namespace lss_cuda_solver
{

using lss_core_cuda_solver::flat_matrix;
using lss_utility::valcopy;

template <> void cuda_solver<memory_space_enum::Host>::initialize()
{
    const double one = static_cast<double>(1.0);
    const double step = one / static_cast<double>(discretization_size_ - 1);
    cuda_boundary_ = std::make_shared<cuda_boundary>(discretization_size_, step);
}
template <>
cuda_solver<memory_space_enum::Host>::cuda_solver(std::size_t discretization_size)
    : lss_tridiagonal_solver::tridiagonal_solver(discretization_size, factorization_enum::QRMethod)
{
    initialize();
}

template <> cuda_solver<memory_space_enum::Host>::~cuda_solver()
{
}

template <>
void cuda_solver<memory_space_enum::Host>::kernel(boundary_1d_pair const &boundary, container_t &solution,
                                                  factorization_enum factorization, double time)
{
    // get proper boundaries:
    const std::size_t N = discretization_size_ - 1;
    const auto &lowest_quad = std::make_tuple(a_[0], b_[0], c_[0], f_[0]);
    const auto &lower_quad = std::make_tuple(a_[1], b_[1], c_[1], f_[1]);
    const auto &higher_quad = std::make_tuple(a_[N - 1], b_[N - 1], c_[N - 1], f_[N - 1]);
    const auto &highest_quad = std::make_tuple(a_[N], b_[N], c_[N], f_[N]);

    cuda_boundary_->set_lowest_quad(lowest_quad);
    cuda_boundary_->set_lower_quad(lower_quad);
    cuda_boundary_->set_higher_quad(higher_quad);
    cuda_boundary_->set_highest_quad(highest_quad);
    const auto &init_coeffs = cuda_boundary_->init_coefficients(boundary, time);
    const std::size_t start_idx = cuda_boundary_->start_index();
    const auto &fin_coeffs = cuda_boundary_->final_coefficients(boundary, time);
    const std::size_t end_idx = cuda_boundary_->end_index();

    const std::size_t system_size = end_idx - start_idx + 1;
    flat_matrix mat(system_size, system_size);
    container_t rhs(double{}, system_size);

    rhs[0] = std::get<2>(init_coeffs);
    mat.emplace_back(0, 0, std::get<0>(init_coeffs));
    mat.emplace_back(0, 1, std::get<1>(init_coeffs));
    for (std::size_t t = 1; t < system_size - 1; ++t)
    {
        mat.emplace_back(t, t - 1, a_[t + start_idx]);
        mat.emplace_back(t, t, b_[t + start_idx]);
        mat.emplace_back(t, t + 1, c_[t + start_idx]);
        rhs[t] = f_[t + start_idx];
    }
    mat.emplace_back(system_size - 1, system_size - 2, std::get<0>(fin_coeffs));
    mat.emplace_back(system_size - 1, system_size - 1, std::get<1>(fin_coeffs));
    rhs[system_size - 1] = std::get<2>(fin_coeffs);

    // initialise the solver:
    real_sparse_solver_cuda<memory_space_enum::Host> rss(system_size);
    rss.initialize(system_size);
    rss.set_flat_sparse_matrix(std::move(mat));
    rss.set_rhs(rhs);
    container_t sub_solution(double{}, system_size);
    rss.solve(sub_solution, factorization);

    valcopy(solution, sub_solution, start_idx);
    // std::copy(sub_solution.begin(), sub_solution.end(), std::next(solution.begin(), start_idx));
    //  fill in the boundary values:
    if (start_idx == 1)
        solution[0] = cuda_boundary_->lower_boundary(boundary, time);
    if (end_idx == N - 1)
        solution[N] = cuda_boundary_->upper_boundary(boundary, time);
}
template <>
void cuda_solver<memory_space_enum::Host>::kernel(boundary_2d_pair const &boundary, container_t &solution,
                                                  factorization_enum factorization, double time, double space_arg)
{
    // get proper boundaries:
    const std::size_t N = discretization_size_ - 1;
    const auto &lowest_quad = std::make_tuple(a_[0], b_[0], c_[0], f_[0]);
    const auto &lower_quad = std::make_tuple(a_[1], b_[1], c_[1], f_[1]);
    const auto &higher_quad = std::make_tuple(a_[N - 1], b_[N - 1], c_[N - 1], f_[N - 1]);
    const auto &highest_quad = std::make_tuple(a_[N], b_[N], c_[N], f_[N]);

    cuda_boundary_->set_lowest_quad(lowest_quad);
    cuda_boundary_->set_lower_quad(lower_quad);
    cuda_boundary_->set_higher_quad(higher_quad);
    cuda_boundary_->set_highest_quad(highest_quad);
    const auto &init_coeffs = cuda_boundary_->init_coefficients(boundary, time, space_arg);
    const std::size_t start_idx = cuda_boundary_->start_index();
    const auto &fin_coeffs = cuda_boundary_->final_coefficients(boundary, time, space_arg);
    const std::size_t end_idx = cuda_boundary_->end_index();

    const std::size_t system_size = end_idx - start_idx + 1;
    flat_matrix mat(system_size, system_size);
    container_t rhs(double{}, system_size);

    rhs[0] = std::get<2>(init_coeffs);
    mat.emplace_back(0, 0, std::get<0>(init_coeffs));
    mat.emplace_back(0, 1, std::get<1>(init_coeffs));
    for (std::size_t t = 1; t < system_size - 1; ++t)
    {
        mat.emplace_back(t, t - 1, a_[t + start_idx]);
        mat.emplace_back(t, t, b_[t + start_idx]);
        mat.emplace_back(t, t + 1, c_[t + start_idx]);
        rhs[t] = f_[t + start_idx];
    }
    mat.emplace_back(system_size - 1, system_size - 2, std::get<0>(fin_coeffs));
    mat.emplace_back(system_size - 1, system_size - 1, std::get<1>(fin_coeffs));
    rhs[system_size - 1] = std::get<2>(fin_coeffs);

    // initialise the solver:
    real_sparse_solver_cuda<memory_space_enum::Host> rss(system_size);
    rss.initialize(system_size);
    rss.set_flat_sparse_matrix(std::move(mat));
    rss.set_rhs(rhs);
    container_t sub_solution(double{}, system_size);
    rss.solve(sub_solution, factorization);

    valcopy(solution, sub_solution, start_idx);
    // std::copy(sub_solution.begin(), sub_solution.end(), std::next(solution.begin(), start_idx));
    //  fill in the boundary values:
    if (start_idx == 1)
        solution[0] = cuda_boundary_->lower_boundary(boundary, time, space_arg);
    if (end_idx == N - 1)
        solution[N] = cuda_boundary_->upper_boundary(boundary, time, space_arg);
}

template <>
void cuda_solver<memory_space_enum::Host>::kernel(boundary_3d_pair const &boundary, container_t &solution,
                                                  factorization_enum factorization, double time, double space_1_arg,
                                                  double space_2_arg)
{
    // get proper boundaries:
    const std::size_t N = discretization_size_ - 1;
    const auto &lowest_quad = std::make_tuple(a_[0], b_[0], c_[0], f_[0]);
    const auto &lower_quad = std::make_tuple(a_[1], b_[1], c_[1], f_[1]);
    const auto &higher_quad = std::make_tuple(a_[N - 1], b_[N - 1], c_[N - 1], f_[N - 1]);
    const auto &highest_quad = std::make_tuple(a_[N], b_[N], c_[N], f_[N]);

    cuda_boundary_->set_lowest_quad(lowest_quad);
    cuda_boundary_->set_lower_quad(lower_quad);
    cuda_boundary_->set_higher_quad(higher_quad);
    cuda_boundary_->set_highest_quad(highest_quad);
    const auto &init_coeffs = cuda_boundary_->init_coefficients(boundary, time, space_1_arg, space_2_arg);
    const std::size_t start_idx = cuda_boundary_->start_index();
    const auto &fin_coeffs = cuda_boundary_->final_coefficients(boundary, time, space_1_arg, space_2_arg);
    const std::size_t end_idx = cuda_boundary_->end_index();

    const std::size_t system_size = end_idx - start_idx + 1;
    flat_matrix mat(system_size, system_size);
    container_t rhs(double{}, system_size);

    rhs[0] = std::get<2>(init_coeffs);
    mat.emplace_back(0, 0, std::get<0>(init_coeffs));
    mat.emplace_back(0, 1, std::get<1>(init_coeffs));
    for (std::size_t t = 1; t < system_size - 1; ++t)
    {
        mat.emplace_back(t, t - 1, a_[t + start_idx]);
        mat.emplace_back(t, t, b_[t + start_idx]);
        mat.emplace_back(t, t + 1, c_[t + start_idx]);
        rhs[t] = f_[t + start_idx];
    }
    mat.emplace_back(system_size - 1, system_size - 2, std::get<0>(fin_coeffs));
    mat.emplace_back(system_size - 1, system_size - 1, std::get<1>(fin_coeffs));
    rhs[system_size - 1] = std::get<2>(fin_coeffs);

    // initialise the solver:
    real_sparse_solver_cuda<memory_space_enum::Host> rss(system_size);
    rss.initialize(system_size);
    rss.set_flat_sparse_matrix(std::move(mat));
    rss.set_rhs(rhs);
    container_t sub_solution(double{}, system_size);
    rss.solve(sub_solution, factorization);

    valcopy(solution, sub_solution, start_idx);
    // std::copy(sub_solution.begin(), sub_solution.end(), std::next(solution.begin(), start_idx));
    //  fill in the boundary values:
    if (start_idx == 1)
        solution[0] = cuda_boundary_->lower_boundary(boundary, time, space_1_arg, space_2_arg);
    if (end_idx == N - 1)
        solution[N] = cuda_boundary_->upper_boundary(boundary, time, space_1_arg, space_2_arg);
}

template <> void cuda_solver<memory_space_enum::Device>::initialize()
{
    const double one = static_cast<double>(1.0);
    const double step = one / static_cast<double>(discretization_size_ - 1);
    cuda_boundary_ = std::make_shared<cuda_boundary>(discretization_size_, step);
}
template <>
cuda_solver<memory_space_enum::Device>::cuda_solver(std::size_t discretization_size)
    : lss_tridiagonal_solver::tridiagonal_solver(discretization_size, factorization_enum::QRMethod)
{
    initialize();
}

template <> cuda_solver<memory_space_enum::Device>::~cuda_solver()
{
}

template <>
void cuda_solver<memory_space_enum::Device>::kernel(boundary_1d_pair const &boundary, container_t &solution,
                                                    factorization_enum factorization, double time)
{
    // get proper boundaries:
    const std::size_t N = discretization_size_ - 1;
    const auto &lowest_quad = std::make_tuple(a_[0], b_[0], c_[0], f_[0]);
    const auto &lower_quad = std::make_tuple(a_[1], b_[1], c_[1], f_[1]);
    const auto &higher_quad = std::make_tuple(a_[N - 1], b_[N - 1], c_[N - 1], f_[N - 1]);
    const auto &highest_quad = std::make_tuple(a_[N], b_[N], c_[N], f_[N]);

    cuda_boundary_->set_lowest_quad(lowest_quad);
    cuda_boundary_->set_lower_quad(lower_quad);
    cuda_boundary_->set_higher_quad(higher_quad);
    cuda_boundary_->set_highest_quad(highest_quad);
    const auto &init_coeffs = cuda_boundary_->init_coefficients(boundary, time);
    const std::size_t start_idx = cuda_boundary_->start_index();
    const auto &fin_coeffs = cuda_boundary_->final_coefficients(boundary, time);
    const std::size_t end_idx = cuda_boundary_->end_index();

    const std::size_t system_size = end_idx - start_idx + 1;
    flat_matrix mat(system_size, system_size);
    container_t rhs(double{}, system_size);

    rhs[0] = std::get<2>(init_coeffs);
    mat.emplace_back(0, 0, std::get<0>(init_coeffs));
    mat.emplace_back(0, 1, std::get<1>(init_coeffs));
    for (std::size_t t = 1; t < system_size - 1; ++t)
    {
        mat.emplace_back(t, t - 1, a_[t + start_idx]);
        mat.emplace_back(t, t, b_[t + start_idx]);
        mat.emplace_back(t, t + 1, c_[t + start_idx]);
        rhs[t] = f_[t + start_idx];
    }
    mat.emplace_back(system_size - 1, system_size - 2, std::get<0>(fin_coeffs));
    mat.emplace_back(system_size - 1, system_size - 1, std::get<1>(fin_coeffs));
    rhs[system_size - 1] = std::get<2>(fin_coeffs);

    // initialise the solver:
    real_sparse_solver_cuda<memory_space_enum::Host> rss(system_size);
    rss.initialize(system_size);
    rss.set_flat_sparse_matrix(std::move(mat));
    rss.set_rhs(rhs);
    container_t sub_solution(double{}, system_size);
    rss.solve(sub_solution, factorization);

    valcopy(solution, sub_solution, start_idx);
    // std::copy(sub_solution.begin(), sub_solution.end(), std::next(solution.begin(), start_idx));
    // fill in the boundary values:
    if (start_idx == 1)
        solution[0] = cuda_boundary_->lower_boundary(boundary, time);
    if (end_idx == N - 1)
        solution[N] = cuda_boundary_->upper_boundary(boundary, time);
}
template <>
void cuda_solver<memory_space_enum::Device>::kernel(boundary_2d_pair const &boundary, container_t &solution,
                                                    factorization_enum factorization, double time, double space_arg)
{
    // get proper boundaries:
    const std::size_t N = discretization_size_ - 1;
    const auto &lowest_quad = std::make_tuple(a_[0], b_[0], c_[0], f_[0]);
    const auto &lower_quad = std::make_tuple(a_[1], b_[1], c_[1], f_[1]);
    const auto &higher_quad = std::make_tuple(a_[N - 1], b_[N - 1], c_[N - 1], f_[N - 1]);
    const auto &highest_quad = std::make_tuple(a_[N], b_[N], c_[N], f_[N]);

    cuda_boundary_->set_lowest_quad(lowest_quad);
    cuda_boundary_->set_lower_quad(lower_quad);
    cuda_boundary_->set_higher_quad(higher_quad);
    cuda_boundary_->set_highest_quad(highest_quad);
    const auto &init_coeffs = cuda_boundary_->init_coefficients(boundary, time, space_arg);
    const std::size_t start_idx = cuda_boundary_->start_index();
    const auto &fin_coeffs = cuda_boundary_->final_coefficients(boundary, time, space_arg);
    const std::size_t end_idx = cuda_boundary_->end_index();

    const std::size_t system_size = end_idx - start_idx + 1;
    flat_matrix mat(system_size, system_size);
    container_t rhs(double{}, system_size);

    rhs[0] = std::get<2>(init_coeffs);
    mat.emplace_back(0, 0, std::get<0>(init_coeffs));
    mat.emplace_back(0, 1, std::get<1>(init_coeffs));
    for (std::size_t t = 1; t < system_size - 1; ++t)
    {
        mat.emplace_back(t, t - 1, a_[t + start_idx]);
        mat.emplace_back(t, t, b_[t + start_idx]);
        mat.emplace_back(t, t + 1, c_[t + start_idx]);
        rhs[t] = f_[t + start_idx];
    }
    mat.emplace_back(system_size - 1, system_size - 2, std::get<0>(fin_coeffs));
    mat.emplace_back(system_size - 1, system_size - 1, std::get<1>(fin_coeffs));
    rhs[system_size - 1] = std::get<2>(fin_coeffs);

    // initialise the solver:
    real_sparse_solver_cuda<memory_space_enum::Device> rss(system_size);
    rss.initialize(system_size);
    rss.set_flat_sparse_matrix(std::move(mat));
    rss.set_rhs(rhs);
    container_t sub_solution(double{}, system_size);
    rss.solve(sub_solution, factorization);

    valcopy(solution, sub_solution, start_idx);
    // std::copy(sub_solution.begin(), sub_solution.end(), std::next(solution.begin(), start_idx));
    //  fill in the boundary values:
    if (start_idx == 1)
        solution[0] = cuda_boundary_->lower_boundary(boundary, time, space_arg);
    if (end_idx == N - 1)
        solution[N] = cuda_boundary_->upper_boundary(boundary, time, space_arg);
}

template <>
void cuda_solver<memory_space_enum::Device>::kernel(boundary_3d_pair const &boundary, container_t &solution,
                                                    factorization_enum factorization, double time, double space_1_arg,
                                                    double space_2_arg)
{
    // get proper boundaries:
    const std::size_t N = discretization_size_ - 1;
    const auto &lowest_quad = std::make_tuple(a_[0], b_[0], c_[0], f_[0]);
    const auto &lower_quad = std::make_tuple(a_[1], b_[1], c_[1], f_[1]);
    const auto &higher_quad = std::make_tuple(a_[N - 1], b_[N - 1], c_[N - 1], f_[N - 1]);
    const auto &highest_quad = std::make_tuple(a_[N], b_[N], c_[N], f_[N]);

    cuda_boundary_->set_lowest_quad(lowest_quad);
    cuda_boundary_->set_lower_quad(lower_quad);
    cuda_boundary_->set_higher_quad(higher_quad);
    cuda_boundary_->set_highest_quad(highest_quad);
    const auto &init_coeffs = cuda_boundary_->init_coefficients(boundary, time, space_1_arg, space_2_arg);
    const std::size_t start_idx = cuda_boundary_->start_index();
    const auto &fin_coeffs = cuda_boundary_->final_coefficients(boundary, time, space_1_arg, space_2_arg);
    const std::size_t end_idx = cuda_boundary_->end_index();

    const std::size_t system_size = end_idx - start_idx + 1;
    flat_matrix mat(system_size, system_size);
    container_t rhs(double{}, system_size);

    rhs[0] = std::get<2>(init_coeffs);
    mat.emplace_back(0, 0, std::get<0>(init_coeffs));
    mat.emplace_back(0, 1, std::get<1>(init_coeffs));
    for (std::size_t t = 1; t < system_size - 1; ++t)
    {
        mat.emplace_back(t, t - 1, a_[t + start_idx]);
        mat.emplace_back(t, t, b_[t + start_idx]);
        mat.emplace_back(t, t + 1, c_[t + start_idx]);
        rhs[t] = f_[t + start_idx];
    }
    mat.emplace_back(system_size - 1, system_size - 2, std::get<0>(fin_coeffs));
    mat.emplace_back(system_size - 1, system_size - 1, std::get<1>(fin_coeffs));
    rhs[system_size - 1] = std::get<2>(fin_coeffs);

    // initialise the solver:
    real_sparse_solver_cuda<memory_space_enum::Device> rss(system_size);
    rss.initialize(system_size);
    rss.set_flat_sparse_matrix(std::move(mat));
    rss.set_rhs(rhs);
    container_t sub_solution(double{}, system_size);
    rss.solve(sub_solution, factorization);

    valcopy(solution, sub_solution, start_idx);
    // std::copy(sub_solution.begin(), sub_solution.end(), std::next(solution.begin(), start_idx));
    //  fill in the boundary values:
    if (start_idx == 1)
        solution[0] = cuda_boundary_->lower_boundary(boundary, time, space_1_arg, space_2_arg);
    if (end_idx == N - 1)
        solution[N] = cuda_boundary_->upper_boundary(boundary, time, space_1_arg, space_2_arg);
}

} // namespace lss_cuda_solver
