#include "lss_karawia_solver.hpp"

#include <tuple>
#include <type_traits>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_macros.hpp"
#include "../../../common/lss_utility.hpp"
#include "lss_karawia_boundary.hpp"

namespace lss_karawia_solver
{

karawia_solver::karawia_solver(std::size_t discretization_size) : discretization_size_{discretization_size}
{
    initialize();
}

void karawia_solver::initialize()
{
    const double one = static_cast<double>(1.0);
    const double step = one / static_cast<double>(discretization_size_ - 1);
    karawia_boundary_ = std::make_shared<karawia_solver_boundary>(discretization_size_, step);
}

karawia_solver ::~karawia_solver()
{
}

void karawia_solver::solve(boundary_1d_pair const &boundary, boundary_1d_pair const &other_boundary,
                           container_t &solution)
{
    LSS_ASSERT(solution.size() == discretization_size_, "Incorrect size of solution container");
    kernel(boundary, other_boundary, solution, double{});
}

void karawia_solver::solve(boundary_1d_pair const &boundary, boundary_1d_pair const &other_boundary,
                           container_t &solution, double at_time)
{
    LSS_ASSERT(solution.size() == discretization_size_, "Incorrect size of solution container");
    kernel(boundary, other_boundary, solution, at_time);
}

void karawia_solver::solve(boundary_2d_pair const &boundary, boundary_2d_pair const &other_boundary,
                           container_t &solution, double at_time, double space_arg)
{
    LSS_ASSERT(solution.size() == discretization_size_, "Incorrect size of solution container");
    kernel(boundary, other_boundary, solution, at_time, space_arg);
}

void karawia_solver::set_diagonals(container_t lowest_diagonal, container_t lower_diagonal, container_t diagonal,
                                   container_t upper_diagonal, container_t uppest_diagonal)
{
    LSS_ASSERT(lowest_diagonal.size() == discretization_size_, "Inncorect size for lowerDiagonal");
    LSS_ASSERT(lower_diagonal.size() == discretization_size_, "Inncorect size for lowerDiagonal");
    LSS_ASSERT(diagonal.size() == discretization_size_, "Inncorect size for diagonal");
    LSS_ASSERT(upper_diagonal.size() == discretization_size_, "Inncorect size for upperDiagonal");
    LSS_ASSERT(uppest_diagonal.size() == discretization_size_, "Inncorect size for upperDiagonal");
    a_ = std::move(lowest_diagonal);
    b_ = std::move(lower_diagonal);
    c_ = std::move(diagonal);
    d_ = std::move(upper_diagonal);
    e_ = std::move(uppest_diagonal);
}

void karawia_solver::set_rhs(container_t const &rhs)
{
    LSS_ASSERT(rhs.size() == discretization_size_, "Inncorect size for right-hand side");
    f_ = rhs;
}

void lss_karawia_solver::karawia_solver::kernel(boundary_1d_pair const &boundary,
                                                boundary_1d_pair const &other_boundary, container_t &solution,
                                                double time)
{
    // check the diagonal dominance:
    // LSS_ASSERT(is_diagonally_dominant() == true, "Tridiagonal matrix must be
    // diagonally dominant.");

    //// clear the working containers:
    // alpha_.clear();
    // beta_.clear();
    // gamma_.clear();
    // mu_.clear();

    // resize the working containers:
    alpha_.resize(discretization_size_);
    beta_.resize(discretization_size_);
    gamma_.resize(discretization_size_);
    mu_.resize(discretization_size_);

    // get proper boundaries:
    const std::size_t N = discretization_size_ - 1;
    const auto &lowest_sixta = std::make_tuple(a_[2], b_[2], c_[2], d_[2], e_[2], f_[2]);
    const auto &lower_sixta = std::make_tuple(a_[3], b_[3], c_[3], d_[3], e_[3], f_[3]);
    const auto &higher_sixta = std::make_tuple(a_[N - 3], b_[N - 3], c_[N - 3], d_[N - 3], e_[N - 3], f_[N - 3]);
    const auto &highest_sixta = std::make_tuple(a_[N - 2], b_[N - 2], c_[N - 2], d_[N - 2], e_[N - 2], f_[N - 2]);

    const auto &init_coeffs = karawia_boundary_->init_coefficients(boundary, other_boundary, time);
    const std::size_t start_idx = karawia_boundary_->start_index();
    const auto &fin_coeffs = karawia_boundary_->final_coefficients(boundary, other_boundary, time);
    const std::size_t end_idx = karawia_boundary_->end_index();

    // init values for the working containers:
    mu_[start_idx] = c_[start_idx];
    alpha_[start_idx] = d_[start_idx] / mu_[start_idx];
    beta_[start_idx] = e_[start_idx] / mu_[start_idx];
    f_[start_idx] = std::get<0>(init_coeffs) / mu_[start_idx];
    //
    const std::size_t next_idx = start_idx + 1;
    gamma_[next_idx] = b_[next_idx];
    mu_[next_idx] = c_[next_idx] - alpha_[start_idx] * gamma_[next_idx];
    alpha_[next_idx] = (d_[next_idx] - beta_[start_idx] * gamma_[next_idx]) / mu_[next_idx];
    beta_[next_idx] = e_[next_idx] / mu_[next_idx];
    f_[next_idx] = (std::get<1>(init_coeffs) - f_[start_idx] * gamma_[next_idx]) / mu_[next_idx];

    for (std::size_t t = next_idx + 1; t <= end_idx - 2; ++t)
    {
        gamma_[t] = b_[t] - alpha_[t - 2] * a_[t];
        mu_[t] = c_[t] - beta_[t - 2] * a_[t] - alpha_[t - 1] * gamma_[t];
        alpha_[t] = (d_[t] - beta_[t - 1] * gamma_[t]) / mu_[t];
        beta_[t] = e_[t] / mu_[t];
        f_[t] = (f_[t] - f_[t - 2] * a_[t] - f_[t - 1] * gamma_[t]) / mu_[t];
    }

    gamma_[end_idx - 1] = b_[end_idx - 1] - alpha_[end_idx - 3] * a_[end_idx - 1];
    mu_[end_idx - 1] =
        c_[end_idx - 1] - beta_[end_idx - 3] * a_[end_idx - 1] - alpha_[end_idx - 2] * gamma_[end_idx - 1];
    alpha_[end_idx - 1] = (d_[end_idx - 1] - beta_[end_idx - 2] * gamma_[end_idx - 1]) / mu_[end_idx - 1];
    //
    gamma_[end_idx] = b_[end_idx - 1] - alpha_[end_idx - 2] * a_[end_idx];
    mu_[end_idx] = c_[end_idx] - beta_[end_idx - 2] * a_[end_idx] - alpha_[end_idx - 1] * gamma_[end_idx];
    f_[end_idx - 1] =
        (std::get<0>(fin_coeffs) - f_[end_idx - 2] * a_[end_idx - 1] - f_[end_idx - 2] * gamma_[end_idx - 1]) /
        mu_[end_idx - 1];
    f_[end_idx] =
        (std::get<1>(fin_coeffs) - f_[end_idx - 1] * a_[end_idx] - f_[end_idx - 1] * gamma_[end_idx]) / mu_[end_idx];

    solution[end_idx] = f_[end_idx];
    solution[end_idx - 1] = f_[end_idx - 1] - alpha_[end_idx - 1] * solution[end_idx];
    for (std::size_t t = end_idx - 1; t-- > start_idx /*&& t >= 0*/; /*t--*/)
    {
        solution[t] = f_[t] - alpha_[t] * solution[t + 1] - beta_[t] * solution[t + 2];
    }

    // fill in the boundary values:
    solution[0] = karawia_boundary_->lower_boundary(boundary, time);
    solution[1] = karawia_boundary_->lower_boundary(other_boundary, time);
    solution[N - 1] = karawia_boundary_->upper_boundary(other_boundary, time);
    solution[N] = karawia_boundary_->upper_boundary(boundary, time);
}

void karawia_solver::kernel(boundary_2d_pair const &boundary, boundary_2d_pair const &other_boundary,
                            container_t &solution, double time, double space_args)
{
    // check the diagonal dominance:
    // LSS_ASSERT(is_diagonally_dominant() == true, "Tridiagonal matrix must be
    // diagonally dominant.");

    //// clear the working containers:
    // alpha_.clear();
    // beta_.clear();
    // gamma_.clear();
    // mu_.clear();

    // resize the working containers:
    alpha_.resize(discretization_size_);
    beta_.resize(discretization_size_);
    gamma_.resize(discretization_size_);
    mu_.resize(discretization_size_);

    // get proper boundaries:
    const std::size_t N = discretization_size_ - 1;
    const auto &lowest_sixta = std::make_tuple(a_[2], b_[2], c_[2], d_[2], e_[2], f_[2]);
    const auto &lower_sixta = std::make_tuple(a_[3], b_[3], c_[3], d_[3], e_[3], f_[3]);
    const auto &higher_sixta = std::make_tuple(a_[N - 3], b_[N - 3], c_[N - 3], d_[N - 3], e_[N - 3], f_[N - 3]);
    const auto &highest_sixta = std::make_tuple(a_[N - 2], b_[N - 2], c_[N - 2], d_[N - 2], e_[N - 2], f_[N - 2]);

    const auto &init_coeffs = karawia_boundary_->init_coefficients(boundary, other_boundary, time, space_args);
    const std::size_t start_idx = karawia_boundary_->start_index();
    const auto &fin_coeffs = karawia_boundary_->final_coefficients(boundary, other_boundary, time, space_args);
    const std::size_t end_idx = karawia_boundary_->end_index();

    // init values for the working containers:
    mu_[start_idx] = c_[start_idx];
    alpha_[start_idx] = d_[start_idx] / mu_[start_idx];
    beta_[start_idx] = e_[start_idx] / mu_[start_idx];
    f_[start_idx] = std::get<0>(init_coeffs) / mu_[start_idx];
    //
    const std::size_t next_idx = start_idx + 1;
    gamma_[next_idx] = b_[next_idx];
    mu_[next_idx] = c_[next_idx] - alpha_[start_idx] * gamma_[next_idx];
    alpha_[next_idx] = (d_[next_idx] - beta_[start_idx] * gamma_[next_idx]) / mu_[next_idx];
    beta_[next_idx] = e_[next_idx] / mu_[next_idx];
    f_[next_idx] = (std::get<1>(init_coeffs) - f_[start_idx] * gamma_[next_idx]) / mu_[next_idx];

    for (std::size_t t = next_idx + 1; t <= end_idx - 2; ++t)
    {
        gamma_[t] = b_[t] - alpha_[t - 2] * a_[t];
        mu_[t] = c_[t] - beta_[t - 2] * a_[t] - alpha_[t - 1] * gamma_[t];
        alpha_[t] = (d_[t] - beta_[t - 1] * gamma_[t]) / mu_[t];
        beta_[t] = e_[t] / mu_[t];
        f_[t] = (f_[t] - f_[t - 2] * a_[t] - f_[t - 1] * gamma_[t]) / mu_[t];
    }

    gamma_[end_idx - 1] = b_[end_idx - 1] - alpha_[end_idx - 3] * a_[end_idx - 1];
    mu_[end_idx - 1] =
        c_[end_idx - 1] - beta_[end_idx - 3] * a_[end_idx - 1] - alpha_[end_idx - 2] * gamma_[end_idx - 1];
    alpha_[end_idx - 1] = (d_[end_idx - 1] - beta_[end_idx - 2] * gamma_[end_idx - 1]) / mu_[end_idx - 1];
    //
    gamma_[end_idx] = b_[end_idx - 1] - alpha_[end_idx - 2] * a_[end_idx];
    mu_[end_idx] = c_[end_idx] - beta_[end_idx - 2] * a_[end_idx] - alpha_[end_idx - 1] * gamma_[end_idx];
    f_[end_idx - 1] =
        (std::get<0>(fin_coeffs) - f_[end_idx - 2] * a_[end_idx - 1] - f_[end_idx - 2] * gamma_[end_idx - 1]) /
        mu_[end_idx - 1];
    f_[end_idx] =
        (std::get<1>(fin_coeffs) - f_[end_idx - 1] * a_[end_idx] - f_[end_idx - 1] * gamma_[end_idx]) / mu_[end_idx];

    solution[end_idx] = f_[end_idx];
    solution[end_idx - 1] = f_[end_idx - 1] - alpha_[end_idx - 1] * solution[end_idx];
    for (std::size_t t = end_idx - 1; t-- > start_idx;)
    {
        solution[t] = f_[t] - alpha_[t] * solution[t + 1] - beta_[t] * solution[t + 2];
    }

    // fill in the boundary values:
    solution[0] = karawia_boundary_->lower_boundary(boundary, time, space_args);
    solution[1] = karawia_boundary_->lower_boundary(other_boundary, time, space_args);
    solution[N - 1] = karawia_boundary_->upper_boundary(other_boundary, time, space_args);
    solution[N] = karawia_boundary_->upper_boundary(boundary, time, space_args);
}

} // namespace lss_karawia_solver
