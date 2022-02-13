/**

    @file      lss_double_sweep_boundary.hpp
    @brief     Boundary for tridiagonal Double Sweep Solver
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_DOUBLE_SWEEP_BOUNDARY_HPP_)
#define _LSS_DOUBLE_SWEEP_BOUNDARY_HPP_

#pragma warning(disable : 4244)

#include <type_traits>
#include <vector>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../common/lss_utility.hpp"

namespace lss_double_sweep_solver
{
using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_3d_pair;
using lss_utility::sptr_t;
using quad_t = std::tuple<double, double, double, double>;

class double_sweep_boundary
{
  private:
    quad_t low_quad_;
    double space_step_, l_, k_, upper_;
    std::size_t start_index_;
    std::size_t discretization_size_;

    explicit double_sweep_boundary() = delete;

    void initialise(boundary_1d_pair const &boundary, double time);

    void initialise(boundary_2d_pair const &boundary, double time, double space_arg);

    void initialise(boundary_3d_pair const &boundary, double time, double space_1_arg, double space_2_arg);

    void finalise(boundary_1d_pair const &boundary, const double &k_nm1, const double &k_n, const double &l_nm1,
                  const double &l_n, double time);

    void finalise(boundary_2d_pair const &boundary, const double &k_nm1, const double &k_n, const double &l_nm1,
                  const double &l_n, double time, double space_arg);

    void finalise(boundary_3d_pair const &boundary, const double &k_nm1, const double &k_n, const double &l_nm1,
                  const double &l_n, double time, double space_1_arg, double space_2_arg);

  public:
    explicit double_sweep_boundary(const std::size_t &discretization_size, const double &space_step);

    ~double_sweep_boundary();

    void set_low_quad(const quad_t &quad);

    std::size_t start_index() const;

    std::size_t end_index(boundary_1d_pair const &boundary) const;

    std::size_t end_index(boundary_2d_pair const &boundary) const;

    std::size_t end_index(boundary_3d_pair const &boundary) const;

    const std::pair<double, double> coefficients(boundary_1d_pair const &boundary, double time);

    const std::pair<double, double> coefficients(boundary_2d_pair const &boundary, double time, double space_arg);

    const std::pair<double, double> coefficients(boundary_3d_pair const &boundary, double time, double space_1_arg,
                                                 double space_2_arg);

    const double upper_boundary(boundary_1d_pair const &boundary, const double &k_nm1, const double &k_n,
                                const double &l_nm1, const double &l_n, double time);

    const double upper_boundary(boundary_2d_pair const &boundary, const double &k_nm1, const double &k_n,
                                const double &l_nm1, const double &l_n, double time, double space_arg);

    const double upper_boundary(boundary_3d_pair const &boundary, const double &k_nm1, const double &k_n,
                                const double &l_nm1, const double &l_n, double time, double space_1_arg,
                                double space_2_arg);

    const double lower_boundary(boundary_1d_pair const &boundary, double time);

    const double lower_boundary(boundary_2d_pair const &boundary, double time, double space_arg);

    const double lower_boundary(boundary_3d_pair const &boundary, double time, double space_1_arg, double space_2_arg);
};

using double_sweep_boundary_ptr = sptr_t<double_sweep_boundary>;

} // namespace lss_double_sweep_solver

#endif ///_LSS_DOUBLE_SWEEP_BOUNDARY_HPP_
