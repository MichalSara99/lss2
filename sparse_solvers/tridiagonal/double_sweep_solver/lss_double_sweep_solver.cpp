#include "lss_double_sweep_solver.hpp"

#include <tuple>
#include <type_traits>
#include <vector>

namespace lss_double_sweep_solver
{

using lss_utility::valcopy;

double_sweep_solver::double_sweep_solver(std::size_t discretization_size)
    : lss_tridiagonal_solver::tridiagonal_solver(discretization_size, factorization_enum::None)
{
    initialize();
}

double_sweep_solver::~double_sweep_solver()
{
}

void double_sweep_solver::initialize()
{
    const double one = static_cast<double>(1.0);
    const double step = one / static_cast<double>(discretization_size_ - 1);
    dss_boundary_ = std::make_shared<double_sweep_boundary>(discretization_size_, step);
}

void double_sweep_solver::kernel(boundary_1d_pair const &boundary, container_t &solution,
                                 factorization_enum factorization, double time)
{
    //// clear coefficients:
    // K_.clear();
    // L_.clear();
    //  resize coefficients:
    K_.resize(discretization_size_);
    L_.resize(discretization_size_);
    // get proper boundaries:
    const std::size_t N = discretization_size_ - 1;
    const auto &low_quad = std::make_tuple(a_[0], b_[0], c_[0], f_[0]);
    dss_boundary_->set_low_quad(low_quad);
    // init coefficients:
    const auto pair = dss_boundary_->coefficients(boundary, time);
    const std::size_t start_index = dss_boundary_->start_index();
    const std::size_t end_index = dss_boundary_->end_index(boundary);

    L_[0] = std::get<1>(pair);
    K_[0] = std::get<0>(pair);

    double tmp{};
    double mone = static_cast<double>(-1.0);
    for (std::size_t t = 1; t <= end_index; ++t)
    {
        tmp = b_[t] + (a_[t] * L_[t - 1]);
        L_[t] = mone * c_[t] / tmp;
        K_[t] = (f_[t] - (a_[t] * K_[t - 1])) / tmp;
    }

    f_[N] = dss_boundary_->upper_boundary(boundary, K_[N - 1], K_[N], L_[N - 1], L_[N], time);

    for (std::size_t t = N; t-- > start_index;)
    {
        f_[t] = (L_[t] * f_[t + 1]) + K_[t];
    }
    if (start_index == 1)
        f_[0] = dss_boundary_->lower_boundary(boundary, time);

    valcopy(solution, f_);
    // std::copy(f_.begin(), f_.end(), solution.begin());
}

void double_sweep_solver::kernel(boundary_2d_pair const &boundary, container_t &solution,
                                 factorization_enum factorization, double time, double space_arg)
{
    //// clear coefficients:
    // K_.clear();
    // L_.clear();
    //  resize coefficients:
    K_.resize(discretization_size_);
    L_.resize(discretization_size_);
    // get proper boundaries:
    const std::size_t N = discretization_size_ - 1;
    const auto &low_quad = std::make_tuple(a_[0], b_[0], c_[0], f_[0]);
    dss_boundary_->set_low_quad(low_quad);
    // init coefficients:
    const auto pair = dss_boundary_->coefficients(boundary, time, space_arg);
    const std::size_t start_index = dss_boundary_->start_index();
    const std::size_t end_index = dss_boundary_->end_index(boundary);

    L_[0] = std::get<1>(pair);
    K_[0] = std::get<0>(pair);

    double tmp{};
    double mone = static_cast<double>(-1.0);
    for (std::size_t t = 1; t <= end_index; ++t)
    {
        tmp = b_[t] + (a_[t] * L_[t - 1]);
        L_[t] = mone * c_[t] / tmp;
        K_[t] = (f_[t] - (a_[t] * K_[t - 1])) / tmp;
    }

    f_[N] = dss_boundary_->upper_boundary(boundary, K_[N - 1], K_[N], L_[N - 1], L_[N], time, space_arg);

    for (std::size_t t = N; t-- > start_index;)
    {
        f_[t] = (L_[t] * f_[t + 1]) + K_[t];
    }
    if (start_index == 1)
        f_[0] = dss_boundary_->lower_boundary(boundary, time, space_arg);

    valcopy(solution, f_);
    // std::copy(f_.begin(), f_.end(), solution.begin());
}

void double_sweep_solver::kernel(boundary_3d_pair const &boundary, container_t &solution,
                                 factorization_enum factorization, double time, double space_1_arg, double space_2_arg)
{
    //// clear coefficients:
    // K_.clear();
    // L_.clear();
    //  resize coefficients:
    K_.resize(discretization_size_);
    L_.resize(discretization_size_);
    // get proper boundaries:
    const std::size_t N = discretization_size_ - 1;
    const auto &low_quad = std::make_tuple(a_[0], b_[0], c_[0], f_[0]);
    dss_boundary_->set_low_quad(low_quad);
    // init coefficients:
    const auto pair = dss_boundary_->coefficients(boundary, time, space_1_arg, space_2_arg);
    const std::size_t start_index = dss_boundary_->start_index();
    const std::size_t end_index = dss_boundary_->end_index(boundary);

    L_[0] = std::get<1>(pair);
    K_[0] = std::get<0>(pair);

    double tmp{};
    double mone = static_cast<double>(-1.0);
    for (std::size_t t = 1; t <= end_index; ++t)
    {
        tmp = b_[t] + (a_[t] * L_[t - 1]);
        L_[t] = mone * c_[t] / tmp;
        K_[t] = (f_[t] - (a_[t] * K_[t - 1])) / tmp;
    }

    f_[N] = dss_boundary_->upper_boundary(boundary, K_[N - 1], K_[N], L_[N - 1], L_[N], time, space_1_arg, space_2_arg);

    for (std::size_t t = N; t-- > start_index;)
    {
        f_[t] = (L_[t] * f_[t + 1]) + K_[t];
    }
    if (start_index == 1)
        f_[0] = dss_boundary_->lower_boundary(boundary, time, space_1_arg, space_2_arg);

    valcopy(solution, f_);
    // std::copy(f_.begin(), f_.end(), solution.begin());
}

// double precision specialization:

} // namespace lss_double_sweep_solver
