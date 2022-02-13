#include "lss_heat_data_config.hpp"

#include <functional>

#include "../common/lss_macros.hpp"
#include "../common/lss_utility.hpp"

namespace lss_pde_solvers
{

void heat_coefficient_data_config_1d::initialize()
{
    LSS_VERIFY(a_coeff_, "a_coefficient must not be null");
    LSS_VERIFY(b_coeff_, "b_coefficient must not be null");
    LSS_VERIFY(c_coeff_, "c_coefficient must not be null");
}

heat_coefficient_data_config_1d::heat_coefficient_data_config_1d(
    std::function<double(double, double)> const &a_coefficient,
    std::function<double(double, double)> const &b_coefficient,
    std::function<double(double, double)> const &c_coefficient)
    : a_coeff_{a_coefficient}, b_coeff_{b_coefficient}, c_coeff_{c_coefficient}
{
    initialize();
}

std::function<double(double, double)> const &heat_coefficient_data_config_1d::a_coefficient() const
{
    return a_coeff_;
}

std::function<double(double, double)> const &heat_coefficient_data_config_1d::b_coefficient() const
{
    return b_coeff_;
}

std::function<double(double, double)> const &heat_coefficient_data_config_1d::c_coefficient() const
{
    return c_coeff_;
}

void heat_coefficient_data_config_2d::initialize()
{
    LSS_VERIFY(a_coeff_, "a_coefficient must not be null");
    LSS_VERIFY(b_coeff_, "b_coefficient must not be null");
    LSS_VERIFY(c_coeff_, "c_coefficient must not be null");
    LSS_VERIFY(d_coeff_, "d_coefficient must not be null");
    LSS_VERIFY(e_coeff_, "e_coefficient must not be null");
    LSS_VERIFY(f_coeff_, "f_coefficient must not be null");
}

heat_coefficient_data_config_2d::heat_coefficient_data_config_2d(
    std::function<double(double, double, double)> const &a_coefficient,
    std::function<double(double, double, double)> const &b_coefficient,
    std::function<double(double, double, double)> const &c_coefficient,
    std::function<double(double, double, double)> const &d_coefficient,
    std::function<double(double, double, double)> const &e_coefficient,
    std::function<double(double, double, double)> const &f_coefficient)
    : a_coeff_{a_coefficient}, b_coeff_{b_coefficient}, c_coeff_{c_coefficient}, d_coeff_{d_coefficient},
      e_coeff_{e_coefficient}, f_coeff_{f_coefficient}
{
    initialize();
}

std::function<double(double, double, double)> const &heat_coefficient_data_config_2d::a_coefficient() const
{
    return a_coeff_;
}

std::function<double(double, double, double)> const &heat_coefficient_data_config_2d::b_coefficient() const
{
    return b_coeff_;
}

std::function<double(double, double, double)> const &heat_coefficient_data_config_2d::c_coefficient() const
{
    return c_coeff_;
}

std::function<double(double, double, double)> const &heat_coefficient_data_config_2d::d_coefficient() const
{
    return d_coeff_;
}

std::function<double(double, double, double)> const &heat_coefficient_data_config_2d::e_coefficient() const
{
    return e_coeff_;
}

std::function<double(double, double, double)> const &heat_coefficient_data_config_2d::f_coefficient() const
{
    return f_coeff_;
}

void heat_coefficient_data_config_3d::initialize()
{
    LSS_VERIFY(a_coeff_, "a_coefficient must not be null");
    LSS_VERIFY(b_coeff_, "b_coefficient must not be null");
    LSS_VERIFY(c_coeff_, "c_coefficient must not be null");
    LSS_VERIFY(d_coeff_, "d_coefficient must not be null");
    LSS_VERIFY(e_coeff_, "e_coefficient must not be null");
    LSS_VERIFY(f_coeff_, "f_coefficient must not be null");
    LSS_VERIFY(g_coeff_, "d_coefficient must not be null");
    LSS_VERIFY(h_coeff_, "e_coefficient must not be null");
    LSS_VERIFY(i_coeff_, "f_coefficient must not be null");
    LSS_VERIFY(j_coeff_, "f_coefficient must not be null");
}

heat_coefficient_data_config_3d::heat_coefficient_data_config_3d(
    std::function<double(double, double, double, double)> const &a_coefficient,
    std::function<double(double, double, double, double)> const &b_coefficient,
    std::function<double(double, double, double, double)> const &c_coefficient,
    std::function<double(double, double, double, double)> const &d_coefficient,
    std::function<double(double, double, double, double)> const &e_coefficient,
    std::function<double(double, double, double, double)> const &f_coefficient,
    std::function<double(double, double, double, double)> const &g_coefficient,
    std::function<double(double, double, double, double)> const &h_coefficient,
    std::function<double(double, double, double, double)> const &i_coefficient,
    std::function<double(double, double, double, double)> const &j_coefficient)
    : a_coeff_{a_coefficient}, b_coeff_{b_coefficient}, c_coeff_{c_coefficient}, d_coeff_{d_coefficient},
      e_coeff_{e_coefficient}, f_coeff_{f_coefficient}, g_coeff_{g_coefficient}, h_coeff_{h_coefficient},
      i_coeff_{i_coefficient}, j_coeff_{j_coefficient}
{
    initialize();
}

std::function<double(double, double, double, double)> const &heat_coefficient_data_config_3d::a_coefficient() const
{
    return a_coeff_;
}

std::function<double(double, double, double, double)> const &heat_coefficient_data_config_3d::b_coefficient() const
{
    return b_coeff_;
}

std::function<double(double, double, double, double)> const &heat_coefficient_data_config_3d::c_coefficient() const
{
    return c_coeff_;
}

std::function<double(double, double, double, double)> const &heat_coefficient_data_config_3d::d_coefficient() const
{
    return d_coeff_;
}

std::function<double(double, double, double, double)> const &heat_coefficient_data_config_3d::e_coefficient() const
{
    return e_coeff_;
}

std::function<double(double, double, double, double)> const &heat_coefficient_data_config_3d::f_coefficient() const
{
    return f_coeff_;
}

std::function<double(double, double, double, double)> const &heat_coefficient_data_config_3d::g_coefficient() const
{
    return g_coeff_;
}

std::function<double(double, double, double, double)> const &heat_coefficient_data_config_3d::h_coefficient() const
{
    return h_coeff_;
}

std::function<double(double, double, double, double)> const &heat_coefficient_data_config_3d::i_coefficient() const
{
    return i_coeff_;
}

std::function<double(double, double, double, double)> const &heat_coefficient_data_config_3d::j_coefficient() const
{
    return j_coeff_;
}

heat_initial_data_config_1d::heat_initial_data_config_1d(std::function<double(double)> const &initial_condition)
    : initial_condition_{initial_condition}
{
    LSS_VERIFY(initial_condition_, "initial_condition must not be null");
}

std::function<double(double)> const &heat_initial_data_config_1d::initial_condition() const
{
    return initial_condition_;
}

heat_initial_data_config_2d::heat_initial_data_config_2d(std::function<double(double, double)> const &initial_condition)
    : initial_condition_{initial_condition}
{
    LSS_VERIFY(initial_condition_, "initial_condition must not be null");
}

std::function<double(double, double)> const &heat_initial_data_config_2d::initial_condition() const
{
    return initial_condition_;
}

heat_initial_data_config_3d::heat_initial_data_config_3d(
    std::function<double(double, double, double)> const &initial_condition)
    : initial_condition_{initial_condition}
{
    LSS_VERIFY(initial_condition_, "initial_condition must not be null");
}

std::function<double(double, double, double)> const &heat_initial_data_config_3d::initial_condition() const
{
    return initial_condition_;
}

heat_source_data_config_1d::heat_source_data_config_1d(std::function<double(double, double)> const &heat_source)
    : heat_source_{heat_source}
{
    LSS_VERIFY(heat_source_, "heat_source must not be null");
}

std::function<double(double, double)> const &heat_source_data_config_1d::heat_source() const
{
    return heat_source_;
}

heat_source_data_config_2d::heat_source_data_config_2d(std::function<double(double, double, double)> const &heat_source)
    : heat_source_{heat_source}
{
    LSS_VERIFY(heat_source_, "heat_source must not be null");
}

std::function<double(double, double, double)> const &heat_source_data_config_2d::heat_source() const
{
    return heat_source_;
}

heat_source_data_config_3d::heat_source_data_config_3d(
    std::function<double(double, double, double, double)> const &heat_source)
    : heat_source_{heat_source}
{
    LSS_VERIFY(heat_source_, "heat_source must not be null");
}

std::function<double(double, double, double, double)> const &heat_source_data_config_3d::heat_source() const
{
    return heat_source_;
}

void heat_data_config_1d::initialize()
{
    LSS_VERIFY(coefficient_data_cfg_, "coefficient_data_config must not be null");
    LSS_VERIFY(initial_data_cfg_, "initial_data_config must not be null");
}

heat_data_config_1d::heat_data_config_1d(heat_coefficient_data_config_1d_ptr const &coefficient_data_config,
                                         heat_initial_data_config_1d_ptr const &initial_data_config,
                                         heat_source_data_config_1d_ptr const &source_data_config)
    : coefficient_data_cfg_{coefficient_data_config}, initial_data_cfg_{initial_data_config}, source_data_cfg_{
                                                                                                  source_data_config}
{
    initialize();
}

heat_data_config_1d ::~heat_data_config_1d()
{
}

heat_source_data_config_1d_ptr const &heat_data_config_1d::source_data_config() const
{
    return source_data_cfg_;
}

std::function<double(double)> const &heat_data_config_1d::initial_condition() const
{
    return initial_data_cfg_->initial_condition();
}

std::function<double(double, double)> const &heat_data_config_1d::a_coefficient() const
{
    return coefficient_data_cfg_->a_coefficient();
}

std::function<double(double, double)> const &heat_data_config_1d::b_coefficient() const
{
    return coefficient_data_cfg_->b_coefficient();
}

std::function<double(double, double)> const &heat_data_config_1d::c_coefficient() const
{
    return coefficient_data_cfg_->c_coefficient();
}

void heat_data_config_2d::initialize()
{
    LSS_VERIFY(coefficient_data_cfg_, "coefficient_data_config must not be null");
    LSS_VERIFY(initial_data_cfg_, "initial_data_config must not be null");
}

heat_data_config_2d::heat_data_config_2d(heat_coefficient_data_config_2d_ptr const &coefficient_data_config,
                                         heat_initial_data_config_2d_ptr const &initial_data_config,
                                         heat_source_data_config_2d_ptr const &source_data_config)
    : coefficient_data_cfg_{coefficient_data_config}, initial_data_cfg_{initial_data_config}, source_data_cfg_{
                                                                                                  source_data_config}
{
    initialize();
}

heat_data_config_2d::~heat_data_config_2d()
{
}

heat_source_data_config_2d_ptr const &heat_data_config_2d::source_data_config() const
{
    return source_data_cfg_;
}

std::function<double(double, double)> const &heat_data_config_2d::initial_condition() const
{
    return initial_data_cfg_->initial_condition();
}

std::function<double(double, double, double)> const &heat_data_config_2d::a_coefficient() const
{
    return coefficient_data_cfg_->a_coefficient();
}

std::function<double(double, double, double)> const &heat_data_config_2d::b_coefficient() const
{
    return coefficient_data_cfg_->b_coefficient();
}

std::function<double(double, double, double)> const &heat_data_config_2d::c_coefficient() const
{
    return coefficient_data_cfg_->c_coefficient();
}

std::function<double(double, double, double)> const &heat_data_config_2d::d_coefficient() const
{
    return coefficient_data_cfg_->d_coefficient();
}

std::function<double(double, double, double)> const &heat_data_config_2d::e_coefficient() const
{
    return coefficient_data_cfg_->e_coefficient();
}

std::function<double(double, double, double)> const &heat_data_config_2d::f_coefficient() const
{
    return coefficient_data_cfg_->f_coefficient();
}

void heat_data_config_3d::initialize()
{
    LSS_VERIFY(coefficient_data_cfg_, "coefficient_data_config must not be null");
    LSS_VERIFY(initial_data_cfg_, "initial_data_config must not be null");
}

heat_data_config_3d::heat_data_config_3d(heat_coefficient_data_config_3d_ptr const &coefficient_data_config,
                                         heat_initial_data_config_3d_ptr const &initial_data_config,
                                         heat_source_data_config_3d_ptr const &source_data_config)
    : coefficient_data_cfg_{coefficient_data_config}, initial_data_cfg_{initial_data_config}, source_data_cfg_{
                                                                                                  source_data_config}
{
    initialize();
}

heat_data_config_3d::~heat_data_config_3d()
{
}

heat_source_data_config_3d_ptr const &heat_data_config_3d::source_data_config() const
{
    return source_data_cfg_;
}

std::function<double(double, double, double)> const &heat_data_config_3d::initial_condition() const
{
    return initial_data_cfg_->initial_condition();
}

std::function<double(double, double, double, double)> const &heat_data_config_3d::a_coefficient() const
{
    return coefficient_data_cfg_->a_coefficient();
}

std::function<double(double, double, double, double)> const &heat_data_config_3d::b_coefficient() const
{
    return coefficient_data_cfg_->b_coefficient();
}

std::function<double(double, double, double, double)> const &heat_data_config_3d::c_coefficient() const
{
    return coefficient_data_cfg_->c_coefficient();
}

std::function<double(double, double, double, double)> const &heat_data_config_3d::d_coefficient() const
{
    return coefficient_data_cfg_->d_coefficient();
}

std::function<double(double, double, double, double)> const &heat_data_config_3d::e_coefficient() const
{
    return coefficient_data_cfg_->e_coefficient();
}

std::function<double(double, double, double, double)> const &heat_data_config_3d::f_coefficient() const
{
    return coefficient_data_cfg_->f_coefficient();
}

std::function<double(double, double, double, double)> const &heat_data_config_3d::g_coefficient() const
{
    return coefficient_data_cfg_->g_coefficient();
}

std::function<double(double, double, double, double)> const &heat_data_config_3d::h_coefficient() const
{
    return coefficient_data_cfg_->h_coefficient();
}

std::function<double(double, double, double, double)> const &heat_data_config_3d::i_coefficient() const
{
    return coefficient_data_cfg_->i_coefficient();
}

std::function<double(double, double, double, double)> const &heat_data_config_3d::j_coefficient() const
{
    return coefficient_data_cfg_->j_coefficient();
}

} // namespace lss_pde_solvers
