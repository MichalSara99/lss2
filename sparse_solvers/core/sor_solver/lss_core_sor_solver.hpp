/**

    @file      lss_core_sor_solver.hpp
    @brief     Core SOR solver
    @details   ~
    @author    Michal Sara
    @date      12.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_CORE_SOR_SOLVER_HPP_)
#define _LSS_CORE_SOR_SOLVER_HPP_

#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_macros.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../containers/lss_flat_matrix.hpp"

namespace lss_core_sor_solver
{

using lss_containers::flat_matrix;
using lss_utility::container_t;
using lss_utility::sptr_t;
using lss_utility::uptr_t;

class core_sor_solver
{

  private:
    container_t b_;
    uptr_t<flat_matrix> matrix_data_ptr_;
    std::size_t system_size_;
    double omega_;

    explicit core_sor_solver();

    void kernel(container_t &solution);
    bool is_diagonally_dominant();

  public:
    explicit core_sor_solver(std::size_t system_size);

    virtual ~core_sor_solver();

    core_sor_solver(core_sor_solver const &) = delete;
    core_sor_solver(core_sor_solver &&) = delete;
    core_sor_solver &operator=(core_sor_solver const &) = delete;
    core_sor_solver &operator=(core_sor_solver &&) = delete;

    void set_flat_sparse_matrix(flat_matrix flat_matrix);

    void set_rhs(container_t const &rhs);

    void set_omega(double value);

    void solve(container_t &solution);

    container_t const solve();
};

using core_sor_solver_ptr = sptr_t<core_sor_solver>;

} // namespace lss_core_sor_solver

#endif ///_LSS_CORE_SOR_SOLVER_HPP_
