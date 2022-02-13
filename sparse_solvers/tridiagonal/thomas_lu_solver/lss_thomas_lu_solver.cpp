#include "lss_thomas_lu_solver.hpp"

#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_macros.hpp"

namespace lss_thomas_lu_solver
{

using lss_utility::valcopy;

thomas_lu_solver::thomas_lu_solver(std::size_t discretization_size)
    : lss_tridiagonal_solver::tridiagonal_solver(discretization_size, factorization_enum::None)
{
    initialize();
}

thomas_lu_solver::~thomas_lu_solver()
{
}

void thomas_lu_solver::initialize()
{
    const double one = static_cast<double>(1.0);
    const double step = one / static_cast<double>(discretization_size_ - 1);
    tlu_boundary_ = std::make_shared<thomas_lu_solver_boundary>(discretization_size_, step);
}

bool lss_thomas_lu_solver::thomas_lu_solver::is_diagonally_dominant() const
{
    // if (std::abs(b_[0]) < std::abs(c_[0]))
    //    return false;
    // if (std::abs(b_[system_size_ - 1]) < std::abs(a_[system_size_ - 1]))
    //    return false;

    // for (std::size_t t = 0; t < system_size_ - 1; ++t)
    //    if (std::abs(b_[t]) < (std::abs(a_[t]) + std::abs(c_[t])))
    //        return false;
    return true;
}

void lss_thomas_lu_solver::thomas_lu_solver::kernel(boundary_1d_pair const &boundary, container_t &solution,
                                                    factorization_enum factorization, double time)
{
    // check the diagonal dominance:
    LSS_ASSERT(is_diagonally_dominant() == true, "Tridiagonal matrix must be diagonally dominant.");

    //// clear the working containers:
    // beta_.clear();
    // gamma_.clear();

    // resize the working containers:
    beta_.resize(discretization_size_);
    gamma_.resize(discretization_size_);

    // get proper boundaries:
    const std::size_t N = discretization_size_ - 1;
    const auto &lowest_quad = std::make_tuple(a_[0], b_[0], c_[0], f_[0]);
    const auto &lower_quad = std::make_tuple(a_[1], b_[1], c_[1], f_[1]);
    const auto &higher_quad = std::make_tuple(a_[N - 1], b_[N - 1], c_[N - 1], f_[N - 1]);
    const auto &highest_quad = std::make_tuple(a_[N], b_[N], c_[N], f_[N]);
    tlu_boundary_->set_lowest_quad(lowest_quad);
    tlu_boundary_->set_lower_quad(lower_quad);
    tlu_boundary_->set_higher_quad(higher_quad);
    tlu_boundary_->set_highest_quad(highest_quad);

    const auto &init_coeffs = tlu_boundary_->init_coefficients(boundary, time);
    const std::size_t start_idx = tlu_boundary_->start_index();
    const auto &fin_coeffs = tlu_boundary_->final_coefficients(boundary, time);
    const std::size_t end_idx = tlu_boundary_->end_index();
    const double a = std::get<0>(fin_coeffs);
    const double b = std::get<1>(fin_coeffs);
    const double r = std::get<2>(fin_coeffs);

    // init values for the working containers:
    beta_[start_idx] = std::get<0>(init_coeffs);
    gamma_[start_idx] = std::get<1>(init_coeffs);

    for (std::size_t t = start_idx + 1; t < end_idx; ++t)
    {
        beta_[t] = b_[t] - (a_[t] * gamma_[t - 1]);
        gamma_[t] = c_[t] / beta_[t];
    }
    beta_[end_idx] = b - (a * gamma_[end_idx - 1]);

    solution[start_idx] = std::get<3>(init_coeffs);
    for (std::size_t t = start_idx + 1; t < end_idx; ++t)
    {
        solution[t] = (f_[t] - (a_[t] * solution[t - 1])) / beta_[t];
    }
    solution[end_idx] = (r - (a * solution[end_idx - 1])) / beta_[end_idx];

    f_[end_idx] = solution[end_idx];
    for (std::size_t t = end_idx; t-- > start_idx;)
    {
        f_[t] = solution[t] - (gamma_[t] * f_[t + 1]);
    }
    // first copy to solution container:
    valcopy(solution, f_);
    // std::copy(f_.begin(), f_.end(), solution.begin());
    // fill in the boundary values:
    if (start_idx == 1)
        solution[0] = tlu_boundary_->lower_boundary(boundary, time);
    if (end_idx == N - 1)
        solution[N] = tlu_boundary_->upper_boundary(boundary, time);
}

void lss_thomas_lu_solver::thomas_lu_solver::kernel(boundary_2d_pair const &boundary, container_t &solution,
                                                    factorization_enum factorization, double time, double space_arg)
{
    // check the diagonal dominance:
    LSS_ASSERT(is_diagonally_dominant() == true, "Tridiagonal matrix must be diagonally dominant.");

    //// clear the working containers:
    // beta_.clear();
    // gamma_.clear();

    // resize the working containers:
    beta_.resize(discretization_size_);
    gamma_.resize(discretization_size_);

    // get proper boundaries:
    const std::size_t N = discretization_size_ - 1;
    const auto &lowest_quad = std::make_tuple(a_[0], b_[0], c_[0], f_[0]);
    const auto &lower_quad = std::make_tuple(a_[1], b_[1], c_[1], f_[1]);
    const auto &higher_quad = std::make_tuple(a_[N - 1], b_[N - 1], c_[N - 1], f_[N - 1]);
    const auto &highest_quad = std::make_tuple(a_[N], b_[N], c_[N], f_[N]);
    tlu_boundary_->set_lowest_quad(lowest_quad);
    tlu_boundary_->set_lower_quad(lower_quad);
    tlu_boundary_->set_higher_quad(higher_quad);
    tlu_boundary_->set_highest_quad(highest_quad);

    const auto &init_coeffs = tlu_boundary_->init_coefficients(boundary, time, space_arg);
    const std::size_t start_idx = tlu_boundary_->start_index();
    const auto &fin_coeffs = tlu_boundary_->final_coefficients(boundary, time, space_arg);
    const std::size_t end_idx = tlu_boundary_->end_index();
    const double a = std::get<0>(fin_coeffs);
    const double b = std::get<1>(fin_coeffs);
    const double r = std::get<2>(fin_coeffs);

    // init values for the working containers:
    beta_[start_idx] = std::get<0>(init_coeffs);
    gamma_[start_idx] = std::get<1>(init_coeffs);

    for (std::size_t t = start_idx + 1; t < end_idx; ++t)
    {
        beta_[t] = b_[t] - (a_[t] * gamma_[t - 1]);
        gamma_[t] = c_[t] / beta_[t];
    }
    beta_[end_idx] = b - (a * gamma_[end_idx - 1]);

    solution[start_idx] = std::get<3>(init_coeffs);
    for (std::size_t t = start_idx + 1; t < end_idx; ++t)
    {
        solution[t] = (f_[t] - (a_[t] * solution[t - 1])) / beta_[t];
    }
    solution[end_idx] = (r - (a * solution[end_idx - 1])) / beta_[end_idx];

    f_[end_idx] = solution[end_idx];
    for (std::size_t t = end_idx; t-- > start_idx;)
    {
        f_[t] = solution[t] - (gamma_[t] * f_[t + 1]);
    }
    // first copy to solution container:
    valcopy(solution, f_);
    // std::copy(f_.begin(), f_.end(), solution.begin());
    // fill in the boundary values:
    if (start_idx == 1)
        solution[0] = tlu_boundary_->lower_boundary(boundary, time, space_arg);
    if (end_idx == N - 1)
        solution[N] = tlu_boundary_->upper_boundary(boundary, time, space_arg);
}

void lss_thomas_lu_solver::thomas_lu_solver::kernel(boundary_3d_pair const &boundary, container_t &solution,
                                                    factorization_enum factorization, double time, double space_1_arg,
                                                    double space_2_arg)
{
    // check the diagonal dominance:
    LSS_ASSERT(is_diagonally_dominant() == true, "Tridiagonal matrix must be diagonally dominant.");

    //// clear the working containers:
    // beta_.clear();
    // gamma_.clear();

    // resize the working containers:
    beta_.resize(discretization_size_);
    gamma_.resize(discretization_size_);

    // get proper boundaries:
    const std::size_t N = discretization_size_ - 1;
    const auto &lowest_quad = std::make_tuple(a_[0], b_[0], c_[0], f_[0]);
    const auto &lower_quad = std::make_tuple(a_[1], b_[1], c_[1], f_[1]);
    const auto &higher_quad = std::make_tuple(a_[N - 1], b_[N - 1], c_[N - 1], f_[N - 1]);
    const auto &highest_quad = std::make_tuple(a_[N], b_[N], c_[N], f_[N]);
    tlu_boundary_->set_lowest_quad(lowest_quad);
    tlu_boundary_->set_lower_quad(lower_quad);
    tlu_boundary_->set_higher_quad(higher_quad);
    tlu_boundary_->set_highest_quad(highest_quad);

    const auto &init_coeffs = tlu_boundary_->init_coefficients(boundary, time, space_1_arg, space_2_arg);
    const std::size_t start_idx = tlu_boundary_->start_index();
    const auto &fin_coeffs = tlu_boundary_->final_coefficients(boundary, time, space_1_arg, space_2_arg);
    const std::size_t end_idx = tlu_boundary_->end_index();
    const double a = std::get<0>(fin_coeffs);
    const double b = std::get<1>(fin_coeffs);
    const double r = std::get<2>(fin_coeffs);

    // init values for the working containers:
    beta_[start_idx] = std::get<0>(init_coeffs);
    gamma_[start_idx] = std::get<1>(init_coeffs);

    for (std::size_t t = start_idx + 1; t < end_idx; ++t)
    {
        beta_[t] = b_[t] - (a_[t] * gamma_[t - 1]);
        gamma_[t] = c_[t] / beta_[t];
    }
    beta_[end_idx] = b - (a * gamma_[end_idx - 1]);

    solution[start_idx] = std::get<3>(init_coeffs);
    for (std::size_t t = start_idx + 1; t < end_idx; ++t)
    {
        solution[t] = (f_[t] - (a_[t] * solution[t - 1])) / beta_[t];
    }
    solution[end_idx] = (r - (a * solution[end_idx - 1])) / beta_[end_idx];

    f_[end_idx] = solution[end_idx];
    for (std::size_t t = end_idx; t-- > start_idx;)
    {
        f_[t] = solution[t] - (gamma_[t] * f_[t + 1]);
    }
    // first copy to solution container:
    valcopy(solution, f_);
    // std::copy(f_.begin(), f_.end(), solution.begin());
    // fill in the boundary values:
    if (start_idx == 1)
        solution[0] = tlu_boundary_->lower_boundary(boundary, time, space_1_arg, space_2_arg);
    if (end_idx == N - 1)
        solution[N] = tlu_boundary_->upper_boundary(boundary, time, space_1_arg, space_2_arg);
}

} // namespace lss_thomas_lu_solver
