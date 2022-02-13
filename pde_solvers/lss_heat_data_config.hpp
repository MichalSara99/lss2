/**

    @file      lss_heat_data_config.hpp
    @brief     Data Configuration object for heat problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_DATA_CONFIG_HPP_)
#define _LSS_HEAT_DATA_CONFIG_HPP_

#include <functional>

#include "../common/lss_macros.hpp"
#include "../common/lss_utility.hpp"

namespace lss_pde_solvers
{

using lss_utility::sptr_t;

/**
    1D heat_initial_data_config structure
 */
struct heat_coefficient_data_config_1d
{
  private:
    std::function<double(double, double)> a_coeff_;
    std::function<double(double, double)> b_coeff_;
    std::function<double(double, double)> c_coeff_;

    explicit heat_coefficient_data_config_1d() = delete;

    void initialize();

  public:
    explicit heat_coefficient_data_config_1d(std::function<double(double, double)> const &a_coefficient,
                                             std::function<double(double, double)> const &b_coefficient,
                                             std::function<double(double, double)> const &c_coefficient);

    /**
        @brief  Diffusion coefficient function
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double)> const &a_coefficient() const;

    /**
        @brief  Convection coefficient function
        @retval function for convection
    **/
    LSS_API std::function<double(double, double)> const &b_coefficient() const;

    /**
        @brief Coefficient function
        @retval function
    **/
    LSS_API std::function<double(double, double)> const &c_coefficient() const;
};

/**
    2D heat_initial_data_config structure
 */
struct heat_coefficient_data_config_2d
{
  private:
    std::function<double(double, double, double)> a_coeff_;
    std::function<double(double, double, double)> b_coeff_;
    std::function<double(double, double, double)> c_coeff_;
    std::function<double(double, double, double)> d_coeff_;
    std::function<double(double, double, double)> e_coeff_;
    std::function<double(double, double, double)> f_coeff_;

    explicit heat_coefficient_data_config_2d() = delete;

    void initialize();

  public:
    explicit heat_coefficient_data_config_2d(std::function<double(double, double, double)> const &a_coefficient,
                                             std::function<double(double, double, double)> const &b_coefficient,
                                             std::function<double(double, double, double)> const &c_coefficient,
                                             std::function<double(double, double, double)> const &d_coefficient,
                                             std::function<double(double, double, double)> const &e_coefficient,
                                             std::function<double(double, double, double)> const &f_coefficient);

    /**
        @brief  Diffusion coefficient function for first space variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double)> const &a_coefficient() const;

    /**
        @brief  Diffusion coefficient function for second space variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double)> const &b_coefficient() const;

    /**
        @brief  Diffusion coefficient function for mixed first and second variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double)> const &c_coefficient() const;

    /**
        @brief  Convection coefficient function for first space variable
        @retval function for convection
    **/
    LSS_API std::function<double(double, double, double)> const &d_coefficient() const;

    /**
        @brief  Convection coefficient function for second space variable
        @retval function for convection
    **/
    LSS_API std::function<double(double, double, double)> const &e_coefficient() const;

    /**
        @brief Coefficient function
        @retval function
    **/
    LSS_API std::function<double(double, double, double)> const &f_coefficient() const;
};

/**
    3D heat_initial_data_config structure
 */
struct heat_coefficient_data_config_3d
{
  private:
    std::function<double(double, double, double, double)> a_coeff_;
    std::function<double(double, double, double, double)> b_coeff_;
    std::function<double(double, double, double, double)> c_coeff_;
    std::function<double(double, double, double, double)> d_coeff_;
    std::function<double(double, double, double, double)> e_coeff_;
    std::function<double(double, double, double, double)> f_coeff_;
    std::function<double(double, double, double, double)> g_coeff_;
    std::function<double(double, double, double, double)> h_coeff_;
    std::function<double(double, double, double, double)> i_coeff_;
    std::function<double(double, double, double, double)> j_coeff_;

    explicit heat_coefficient_data_config_3d() = delete;

    void initialize();

  public:
    explicit heat_coefficient_data_config_3d(
        std::function<double(double, double, double, double)> const &a_coefficient,
        std::function<double(double, double, double, double)> const &b_coefficient,
        std::function<double(double, double, double, double)> const &c_coefficient,
        std::function<double(double, double, double, double)> const &d_coefficient,
        std::function<double(double, double, double, double)> const &e_coefficient,
        std::function<double(double, double, double, double)> const &f_coefficient,
        std::function<double(double, double, double, double)> const &g_coefficient,
        std::function<double(double, double, double, double)> const &h_coefficient,
        std::function<double(double, double, double, double)> const &i_coefficient,
        std::function<double(double, double, double, double)> const &j_coefficient);

    /**
        @brief  Diffusion coefficient function for first space variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double, double)> const &a_coefficient() const;

    /**
        @brief  Diffusion coefficient function for second space variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double, double)> const &b_coefficient() const;

    /**
        @brief  Diffusion coefficient function for third space variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double, double)> const &c_coefficient() const;

    /**
        @brief  Diffusion coefficient function for mixed first and second
    variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double, double)> const &d_coefficient() const;

    /**
        @brief  Diffusion coefficient function for mixed first and third
    variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double, double)> const &e_coefficient() const;

    /**
        @brief  Diffusion coefficient function for mixed second and third
    variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double, double)> const &f_coefficient() const;

    /**
        @brief  Convection coefficient function for first space variable
        @retval  function for convection
    **/
    LSS_API std::function<double(double, double, double, double)> const &g_coefficient() const;

    /**
        @brief  Convection coefficient function for second space variable
        @retval  function for convection
    **/
    LSS_API std::function<double(double, double, double, double)> const &h_coefficient() const;

    /**
        @brief  Convection coefficient function for third space variable
        @retval  function for convection
    **/
    LSS_API std::function<double(double, double, double, double)> const &i_coefficient() const;

    /**
        @brief  coefficient function
        @retval  function
    **/
    LSS_API std::function<double(double, double, double, double)> const &j_coefficient() const;
};

using heat_coefficient_data_config_1d_ptr = sptr_t<heat_coefficient_data_config_1d>;

using heat_coefficient_data_config_2d_ptr = sptr_t<heat_coefficient_data_config_2d>;

using heat_coefficient_data_config_3d_ptr = sptr_t<heat_coefficient_data_config_3d>;

/**
    1D heat_initial_data_config structure
 */
struct heat_initial_data_config_1d
{
  private:
    std::function<double(double)> initial_condition_;

    explicit heat_initial_data_config_1d() = delete;

  public:
    explicit heat_initial_data_config_1d(std::function<double(double)> const &initial_condition);

    /**
        @brief  Initial/terminal condition function
        @retval function
    **/
    LSS_API std::function<double(double)> const &initial_condition() const;
};

/**
    2D heat_initial_data_config structure
 */
struct heat_initial_data_config_2d
{
  private:
    std::function<double(double, double)> initial_condition_;

    explicit heat_initial_data_config_2d() = delete;

  public:
    explicit heat_initial_data_config_2d(std::function<double(double, double)> const &initial_condition);

    /**
        @brief  Initial/terminal condition function
        @retval  function
    **/
    LSS_API std::function<double(double, double)> const &initial_condition() const;
};

/**
    3D heat_initial_data_config structure
 */
struct heat_initial_data_config_3d
{
  private:
    std::function<double(double, double, double)> initial_condition_;

    explicit heat_initial_data_config_3d() = delete;

  public:
    explicit heat_initial_data_config_3d(std::function<double(double, double, double)> const &initial_condition);

    /**
        @brief  Initial/terminal condition function
        @retval  function
    **/
    LSS_API std::function<double(double, double, double)> const &initial_condition() const;
};

using heat_initial_data_config_1d_ptr = sptr_t<heat_initial_data_config_1d>;

using heat_initial_data_config_2d_ptr = sptr_t<heat_initial_data_config_2d>;

using heat_initial_data_config_3d_ptr = sptr_t<heat_initial_data_config_3d>;

/**
    1D heat_source_data_config structure
 */
struct heat_source_data_config_1d
{
  private:
    std::function<double(double, double)> heat_source_;

    explicit heat_source_data_config_1d() = delete;

  public:
    explicit heat_source_data_config_1d(std::function<double(double, double)> const &heat_source);

    /**
        @brief  Source function
        @retval   function
    **/
    LSS_API std::function<double(double, double)> const &heat_source() const;
};

/**
    2D heat_source_data_config structure
 */
struct heat_source_data_config_2d
{
  private:
    std::function<double(double, double, double)> heat_source_;

    explicit heat_source_data_config_2d() = delete;

  public:
    explicit heat_source_data_config_2d(std::function<double(double, double, double)> const &heat_source);

    /**
        @brief  Source function
        @retval   function
    **/
    LSS_API std::function<double(double, double, double)> const &heat_source() const;
};

/**
    3D heat_source_data_config structure
 */
struct heat_source_data_config_3d
{
  private:
    std::function<double(double, double, double, double)> heat_source_;

    explicit heat_source_data_config_3d() = delete;

  public:
    explicit heat_source_data_config_3d(std::function<double(double, double, double, double)> const &heat_source);

    /**
        @brief  Source function
        @retval  function
    **/
    LSS_API std::function<double(double, double, double, double)> const &heat_source() const;
};

using heat_source_data_config_1d_ptr = sptr_t<heat_source_data_config_1d>;

using heat_source_data_config_2d_ptr = sptr_t<heat_source_data_config_2d>;

using heat_source_data_config_3d_ptr = sptr_t<heat_source_data_config_3d>;

/**
    1D heat_data_config structure
 */
struct heat_data_config_1d
{
  private:
    heat_coefficient_data_config_1d_ptr coefficient_data_cfg_;
    heat_initial_data_config_1d_ptr initial_data_cfg_;
    heat_source_data_config_1d_ptr source_data_cfg_;

    void initialize();
    explicit heat_data_config_1d() = delete;

  public:
    explicit heat_data_config_1d(heat_coefficient_data_config_1d_ptr const &coefficient_data_config,
                                 heat_initial_data_config_1d_ptr const &initial_data_config,
                                 heat_source_data_config_1d_ptr const &source_data_config = nullptr);

    ~heat_data_config_1d();

    /**
        @brief Source data configuration
        @retval
    **/
    LSS_API heat_source_data_config_1d_ptr const &source_data_config() const;

    /**
        @brief Initial/terminal condition function
        @retval  function
    **/
    LSS_API std::function<double(double)> const &initial_condition() const;

    /**
        @brief  Diffusion coefficient function
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double)> const &a_coefficient() const;

    /**
        @brief  Convection coefficient function
        @retval  function for convection
    **/
    LSS_API std::function<double(double, double)> const &b_coefficient() const;

    /**
        @brief  Coefficient function
        @retval  function
    **/
    LSS_API std::function<double(double, double)> const &c_coefficient() const;
};

/**
    2D heat_data_config structure
 */
struct heat_data_config_2d
{
  private:
    heat_coefficient_data_config_2d_ptr coefficient_data_cfg_;
    heat_initial_data_config_2d_ptr initial_data_cfg_;
    heat_source_data_config_2d_ptr source_data_cfg_;

    void initialize();

    explicit heat_data_config_2d() = delete;

  public:
    explicit heat_data_config_2d(heat_coefficient_data_config_2d_ptr const &coefficient_data_config,
                                 heat_initial_data_config_2d_ptr const &initial_data_config,
                                 heat_source_data_config_2d_ptr const &source_data_config = nullptr);

    ~heat_data_config_2d();

    /**
        @brief Source data configuration
        @retval
    **/
    LSS_API heat_source_data_config_2d_ptr const &source_data_config() const;

    /**
        @brief Initial/terminal condition function
        @retval  function
    **/
    LSS_API std::function<double(double, double)> const &initial_condition() const;

    /**
        @brief  Diffusion coefficient function for first space variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double)> const &a_coefficient() const;

    /**
        @brief  Diffusion coefficient function for second space variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double)> const &b_coefficient() const;

    /**
        @brief  Diffusion coefficient function for mixed first and second space variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double)> const &c_coefficient() const;

    /**
        @brief  Convection coefficient function for first space variable
        @retval  function for convection
    **/
    LSS_API std::function<double(double, double, double)> const &d_coefficient() const;

    /**
        @brief  Convection coefficient function for second space variable
        @retval  function for convection
    **/
    LSS_API std::function<double(double, double, double)> const &e_coefficient() const;

    /**
        @brief   coefficient function
        @retval  function
    **/
    LSS_API std::function<double(double, double, double)> const &f_coefficient() const;
};

struct heat_data_config_3d
{
  private:
    heat_coefficient_data_config_3d_ptr coefficient_data_cfg_;
    heat_initial_data_config_3d_ptr initial_data_cfg_;
    heat_source_data_config_3d_ptr source_data_cfg_;

    void initialize();

    explicit heat_data_config_3d() = delete;

  public:
    explicit heat_data_config_3d(heat_coefficient_data_config_3d_ptr const &coefficient_data_config,
                                 heat_initial_data_config_3d_ptr const &initial_data_config,
                                 heat_source_data_config_3d_ptr const &source_data_config = nullptr);

    ~heat_data_config_3d();

    /**
        @brief Source data configuration
        @retval
    **/
    LSS_API heat_source_data_config_3d_ptr const &source_data_config() const;

    /**
        @brief Initial/terminal condition function
        @retval  function
    **/
    LSS_API std::function<double(double, double, double)> const &initial_condition() const;

    /**
        @brief  Diffusion coefficient function for first space variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double, double)> const &a_coefficient() const;

    /**
        @brief  Diffusion coefficient function for second space variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double, double)> const &b_coefficient() const;

    /**
        @brief  Diffusion coefficient function for third space variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double, double)> const &c_coefficient() const;

    /**
        @brief  Diffusion coefficient function for mixed first and second
    variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double, double)> const &d_coefficient() const;

    /**
        @brief  Diffusion coefficient function for mixed first and third
    variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double, double)> const &e_coefficient() const;

    /**
        @brief  Diffusion coefficient function for mixed second and third
    variable
        @retval  function for diffusion
    **/
    LSS_API std::function<double(double, double, double, double)> const &f_coefficient() const;

    /**
        @brief  Convection coefficient function for first space variable
        @retval  function for convection
    **/
    LSS_API std::function<double(double, double, double, double)> const &g_coefficient() const;

    /**
        @brief  Convection coefficient function for second space variable
        @retval  function for convection
    **/
    LSS_API std::function<double(double, double, double, double)> const &h_coefficient() const;

    /**
        @brief  Convection coefficient function for third space variable
        @retval  function for convection
    **/
    LSS_API std::function<double(double, double, double, double)> const &i_coefficient() const;

    /**
        @brief   coefficient function
        @retval  function
    **/
    LSS_API std::function<double(double, double, double, double)> const &j_coefficient() const;
};

using heat_data_config_1d_ptr = sptr_t<heat_data_config_1d>;

using heat_data_config_2d_ptr = sptr_t<heat_data_config_2d>;

using heat_data_config_3d_ptr = sptr_t<heat_data_config_3d>;

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_DATA_CONFIG_HPP_
