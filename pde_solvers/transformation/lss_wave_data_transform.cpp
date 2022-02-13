#include "lss_wave_data_transform.hpp"

namespace lss_pde_solvers
{

void wave_data_transform_1d::initialize(wave_data_config_1d_ptr const &wave_data_config,
                                        grid_transform_config_1d_ptr const grid_transform_config)
{
    auto const A = wave_data_config->a_coefficient();
    auto const B = wave_data_config->b_coefficient();
    auto const C = wave_data_config->c_coefficient();
    auto const D = wave_data_config->d_coefficient();
    auto const a = grid_transform_config->a_derivative();
    auto const b = grid_transform_config->b_derivative();
    auto const init_first = wave_data_config->first_initial_condition();
    auto const init_second = wave_data_config->second_initial_condition();
    auto const &src_cfg = wave_data_config->source_data_config();
    std::function<double(double, double)> src = nullptr;
    if (src_cfg)
    {
        src = src_cfg->wave_source();
        is_wave_source_set_ = true;
        src_coeff_ = [=](double t, double zeta) {
            auto const x = grid_1d::transformed_value(grid_transform_config, zeta);
            return src(t, x);
        };
    }

    a_coeff_ = [=](double t, double zeta) {
        auto const x = grid_1d::transformed_value(grid_transform_config, zeta);
        return A(t, x);
    };

    b_coeff_ = [=](double t, double zeta) {
        auto const x = grid_1d::transformed_value(grid_transform_config, zeta);
        return (B(t, x) / (a(zeta) * a(zeta)));
    };

    c_coeff_ = [=](double t, double zeta) {
        auto const x = grid_1d::transformed_value(grid_transform_config, zeta);
        auto const a_val = a(zeta);
        auto const first = C(t, x) / a_val;
        auto const second = (B(t, x) * b(zeta)) / (a_val * a_val * a_val);
        return (first - second);
    };

    d_coeff_ = [=](double t, double zeta) {
        auto const x = grid_1d::transformed_value(grid_transform_config, zeta);
        return D(t, x);
    };

    init_first_coeff_ = [=](double zeta) {
        auto const x = grid_1d::transformed_value(grid_transform_config, zeta);
        return init_first(x);
    };

    init_second_coeff_ = [=](double zeta) {
        auto const x = grid_1d::transformed_value(grid_transform_config, zeta);
        return init_second(x);
    };
}

wave_data_transform_1d::wave_data_transform_1d(wave_data_config_1d_ptr const &wave_data_config,
                                               grid_transform_config_1d_ptr const grid_transform_config)
{
    initialize(wave_data_config, grid_transform_config);
}

wave_data_transform_1d::~wave_data_transform_1d()
{
}

bool const &wave_data_transform_1d::is_wave_source_set() const
{
    return is_wave_source_set_;
}

std::function<double(double, double)> wave_data_transform_1d::wave_source() const
{
    return (is_wave_source_set() == true) ? src_coeff_ : nullptr;
}

std::function<double(double)> const &wave_data_transform_1d::first_initial_condition() const
{
    return init_first_coeff_;
}

std::function<double(double)> const &wave_data_transform_1d::second_initial_condition() const
{
    return init_second_coeff_;
}

std::function<double(double, double)> const &wave_data_transform_1d::a_coefficient() const
{
    return a_coeff_;
}

std::function<double(double, double)> const &wave_data_transform_1d::b_coefficient() const
{
    return b_coeff_;
}

std::function<double(double, double)> const &wave_data_transform_1d::c_coefficient() const
{
    return c_coeff_;
}

std::function<double(double, double)> const &wave_data_transform_1d::d_coefficient() const
{
    return d_coeff_;
}

} // namespace lss_pde_solvers
