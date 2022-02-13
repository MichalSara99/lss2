#include "lss_ode_data_transform.hpp"

#include "../../common/lss_macros.hpp"

namespace lss_ode_solvers
{
ode_data_transform::ode_data_transform(ode_data_config_ptr const &ode_data_config,
                                       grid_transform_config_1d_ptr const grid_transform_config)
{
    initialize(ode_data_config, grid_transform_config);
}

ode_data_transform ::~ode_data_transform()
{
}

bool const &ode_data_transform::is_nonhom_data_set() const
{
    return is_nonhom_data_set_;
}

std::function<double(double)> ode_data_transform::nonhom_function() const
{
    return (is_nonhom_data_set() == true) ? nonhom_coeff_ : nullptr;
}

std::function<double(double)> const &ode_data_transform::a_coefficient() const
{
    return a_coeff_;
}

std::function<double(double)> const &ode_data_transform::b_coefficient() const
{
    return b_coeff_;
}

void ode_data_transform::initialize(ode_data_config_ptr const &ode_data_config,
                                    grid_transform_config_1d_ptr const grid_transform_config)
{
    auto const A = ode_data_config->a_coefficient();
    auto const B = ode_data_config->b_coefficient();
    auto const a = grid_transform_config->a_derivative();
    auto const b = grid_transform_config->b_derivative();
    auto const &nonhom_cfg = ode_data_config->nonhom_data_config();
    std::function<double(double)> nonhom = nullptr;
    if (nonhom_cfg)
    {
        nonhom = nonhom_cfg->nonhom_function();
        is_nonhom_data_set_ = true;
        nonhom_coeff_ = [=](double zeta) {
            auto const x = grid_1d::transformed_value(grid_transform_config, zeta);
            auto const a_val = a(zeta);
            return (a_val * a_val * nonhom(x));
        };
    }
    a_coeff_ = [=](double zeta) {
        auto const x = grid_1d::transformed_value(grid_transform_config, zeta);
        auto const first = a(zeta) * A(x);
        auto const second = b(zeta) / a(zeta);
        return (first - second);
    };

    b_coeff_ = [=](double zeta) {
        auto const x = grid_1d::transformed_value(grid_transform_config, zeta);
        auto const a_val = a(zeta);
        return (a_val * a_val * B(x));
    };
}
} // namespace lss_ode_solvers
