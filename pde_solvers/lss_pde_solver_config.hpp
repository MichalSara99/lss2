/**

    @file      lss_pde_solver_config.hpp
    @brief     PDE solver configuration objects
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_PDE_SOLVER_CONFIG_HPP_)
#define _LSS_PDE_SOLVER_CONFIG_HPP_

#include <map>

#include "../common/lss_enumerations.hpp"
#include "../common/lss_macros.hpp"
#include "../common/lss_utility.hpp"

namespace lss_pde_solvers
{

using lss_enumerations::dimension_enum;
using lss_enumerations::explicit_pde_schemes_enum;
using lss_enumerations::factorization_enum;
using lss_enumerations::implicit_pde_schemes_enum;
using lss_enumerations::memory_space_enum;
using lss_enumerations::traverse_direction_enum;
using lss_enumerations::tridiagonal_method_enum;
using lss_utility::sptr_t;

/**
     pde_implicit_solver_config base class
 */
struct pde_implicit_solver_config
{
  private:
    memory_space_enum memory_space_;
    traverse_direction_enum traverse_direction_;
    tridiagonal_method_enum tridiagonal_method_;
    factorization_enum tridiagonal_factorization_;

    explicit pde_implicit_solver_config() = delete;

    void initialize();

  public:
    explicit pde_implicit_solver_config(memory_space_enum const &memory_space,
                                        traverse_direction_enum const &traverse_direction,
                                        tridiagonal_method_enum const &tridiagonal_method,
                                        factorization_enum const &tridiagonal_factorization);
    virtual ~pde_implicit_solver_config();

    /**
        @brief Memory space enum value
        @retval
    **/
    LSS_API memory_space_enum memory_space() const;

    /**
        @brief Traverse direction enum value
        @retval
    **/
    LSS_API traverse_direction_enum traverse_direction() const;

    /**
        @brief Tridiagonal method enum value.
        @retval
    **/
    LSS_API tridiagonal_method_enum tridiagonal_method() const;

    /**
        @brief Tridiagonal factorization enum value
        @retval
    **/
    LSS_API factorization_enum tridiagonal_factorization() const;
};

/**
     pde_explicit_solver_config base class
 */
struct pde_explicit_solver_config
{
  private:
    memory_space_enum memory_space_;
    traverse_direction_enum traverse_direction_;

    explicit pde_explicit_solver_config() = delete;

  public:
    explicit pde_explicit_solver_config(memory_space_enum const &memory_space,
                                        traverse_direction_enum const &traverse_direction);
    ~pde_explicit_solver_config();

    /**
        @brief Memory space enum value
        @retval
    **/
    LSS_API memory_space_enum memory_space() const;

    /**
        @brief Traverse direction enum value
        @retval
    **/
    LSS_API traverse_direction_enum traverse_direction() const;
};

using pde_implicit_solver_config_ptr = sptr_t<pde_implicit_solver_config>;

using pde_explicit_solver_config_ptr = sptr_t<pde_explicit_solver_config>;

} // namespace lss_pde_solvers

#endif ///_LSS_PDE_SOLVER_CONFIG_HPP_
