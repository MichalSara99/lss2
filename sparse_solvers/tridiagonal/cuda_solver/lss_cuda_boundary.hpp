/**

    @file      lss_cuda_boundary.hpp
    @brief     Boundary for tridiagonal CUDA solver
    @details   ~
    @author    Michal Sara
    @date      12.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_CUDA_BOUNDARY_HPP_)
#define _LSS_CUDA_BOUNDARY_HPP_

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../common/lss_utility.hpp"

namespace lss_cuda_solver
{

using quad_t = std::tuple<double, double, double, double>;
using triplet_t = std::tuple<double, double, double>;

using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_3d_pair;
using lss_utility::sptr_t;

class cuda_boundary
{
  private:
    quad_t lowest_quad_;
    quad_t lower_quad_;
    quad_t higher_quad_;
    quad_t highest_quad_;
    double space_step_;
    double b_init_, c_init_, f_init_;
    double a_end_, b_end_, f_end_;
    std::size_t start_index_, end_index_;
    std::size_t discretization_size_;

    explicit cuda_boundary() = delete;

    void initialise(boundary_1d_pair const &boundary, double time);

    void initialise(boundary_2d_pair const &boundary, double time, double space_arg);

    void initialise(boundary_3d_pair const &boundary, double time, double space_1_arg, double space_2_arg);

    void finalise(boundary_1d_pair const &boundary, double time);

    void finalise(boundary_2d_pair const &boundary, double time, double space_arg);

    void finalise(boundary_3d_pair const &boundary, double time, double space_1_arg, double space_2_arg);

  public:
    explicit cuda_boundary(const std::size_t discretization_size, const double &space_step);

    ~cuda_boundary();

    void set_lowest_quad(const quad_t &lowest_quad);

    void set_lower_quad(const quad_t &lower_quad);

    void set_highest_quad(const quad_t &highest_quad);

    void set_higher_quad(const quad_t &higher_quad);

    const triplet_t init_coefficients(boundary_1d_pair const &boundary, double time);

    const triplet_t init_coefficients(boundary_2d_pair const &boundary, double time, double space_arg);

    const triplet_t init_coefficients(boundary_3d_pair const &boundary, double time, double space_1_arg,
                                      double space_2_arg);

    const triplet_t final_coefficients(boundary_1d_pair const &boundary, double time);

    const triplet_t final_coefficients(boundary_2d_pair const &boundary, double time, double space_arg);

    const triplet_t final_coefficients(boundary_3d_pair const &boundary, double time, double space_1_arg,
                                       double space_2_arg);

    std::size_t start_index() const;

    std::size_t end_index() const;

    const double upper_boundary(boundary_1d_pair const &boundary, double time);

    const double upper_boundary(boundary_2d_pair const &boundary, double time, double space_arg);

    const double upper_boundary(boundary_3d_pair const &boundary, double time, double space_1_arg, double space_2_arg);

    const double lower_boundary(boundary_1d_pair const &boundary, double time);

    const double lower_boundary(boundary_2d_pair const &boundary, double time, double space_arg);

    const double lower_boundary(boundary_3d_pair const &boundary, double time, double space_1_arg, double space_2_arg);
};

using cuda_boundary_ptr = sptr_t<cuda_boundary>;

} // namespace lss_cuda_solver

#endif ///_LSS_CUDA_BOUNDARY_HPP_
