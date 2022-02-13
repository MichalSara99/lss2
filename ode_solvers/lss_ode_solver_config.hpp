/**

    @file      lss_ode_solver_config.hpp
    @brief      Represents ODE solver configurations
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once
#if !defined(_LSS_ODE_SOLVER_CONFIG_HPP_)
#define _LSS_ODE_SOLVER_CONFIG_HPP_

#include "../common/lss_enumerations.hpp"
#include "../common/lss_utility.hpp"

namespace lss_ode_solvers
{

using lss_enumerations::factorization_enum;
using lss_enumerations::memory_space_enum;
using lss_enumerations::tridiagonal_method_enum;
using lss_utility::sptr_t;

/**
    @struct ode_implicit_solver_config
    @brief
**/
struct ode_implicit_solver_config
{
  private:
    memory_space_enum memory_space_;
    tridiagonal_method_enum tridiagonal_method_;
    factorization_enum tridiagonal_factorization_;

    explicit ode_implicit_solver_config() = delete;

    void initialize();

  public:
    explicit ode_implicit_solver_config(memory_space_enum const &memory_space,
                                        tridiagonal_method_enum const &tridiagonal_method,
                                        factorization_enum const &tridiagonal_factorization);
    ~ode_implicit_solver_config();

    /**
        @brief Memory space enum value
        @retval
    **/
    LSS_API memory_space_enum memory_space() const;

    /**
        @brief Tridiagonal method enum value
        @retval
    **/
    LSS_API tridiagonal_method_enum tridiagonal_method() const;

    /**
        @brief Tridiagonal factorization enum value
        @retval
    **/
    LSS_API factorization_enum tridiagonal_factorization() const;
};

using ode_implicit_solver_config_ptr = sptr_t<ode_implicit_solver_config>;

// =================================================
// ===== some default implicit solver configs ======
// =================================================
namespace default_ode_solver_configs
{

static auto dev_cusolver_qr_solver_config_ptr = std::make_shared<ode_implicit_solver_config>(
    memory_space_enum::Device, tridiagonal_method_enum::CUDASolver, factorization_enum::QRMethod);

static auto dev_cusolver_lu_solver_config_ptr = std::make_shared<ode_implicit_solver_config>(
    memory_space_enum::Device, tridiagonal_method_enum::CUDASolver, factorization_enum::LUMethod);

static auto host_cusolver_qr_solver_config_ptr = std::make_shared<ode_implicit_solver_config>(
    memory_space_enum::Host, tridiagonal_method_enum::CUDASolver, factorization_enum::QRMethod);

static auto host_cusolver_lu_solver_config_ptr = std::make_shared<ode_implicit_solver_config>(
    memory_space_enum::Host, tridiagonal_method_enum::CUDASolver, factorization_enum::LUMethod);

static auto dev_sorsolver_solver_config_ptr = std::make_shared<ode_implicit_solver_config>(
    memory_space_enum::Device, tridiagonal_method_enum::SORSolver, factorization_enum::None);

static auto host_sorsolver_solver_config_ptr = std::make_shared<ode_implicit_solver_config>(
    memory_space_enum::Host, tridiagonal_method_enum::SORSolver, factorization_enum::None);

static auto host_dssolver_solver_config_ptr = std::make_shared<ode_implicit_solver_config>(
    memory_space_enum::Host, tridiagonal_method_enum::DoubleSweepSolver, factorization_enum::None);

static auto host_tlusolver_solver_config_ptr = std::make_shared<ode_implicit_solver_config>(
    memory_space_enum::Host, tridiagonal_method_enum::ThomasLUSolver, factorization_enum::None);

} // namespace default_ode_solver_configs
} // namespace lss_ode_solvers

#endif ///_LSS_ODE_SOLVER_CONFIG_HPP_
