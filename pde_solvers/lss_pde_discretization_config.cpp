#include "lss_pde_discretization_config.hpp"

namespace lss_pde_solvers
{

pde_discretization_config_1d::pde_discretization_config_1d(range_ptr const &space_range,
                                                           std::size_t const &number_of_space_points,
                                                           range_ptr const &time_range,
                                                           std::size_t const &number_of_time_points)
    : lss_discretization::discretization_config_1d(space_range, number_of_space_points), time_range_{time_range},
      number_of_time_points_{number_of_time_points}
{
}
pde_discretization_config_1d::~pde_discretization_config_1d()
{
}
range_ptr const &pde_discretization_config_1d::time_range() const
{
    return time_range_;
}

std::size_t pde_discretization_config_1d::number_of_time_points() const
{
    return number_of_time_points_;
}

double pde_discretization_config_1d::time_step() const
{
    return ((time_range_->spread()) / static_cast<double>(number_of_time_points_ - 1));
}

pde_discretization_config_2d::pde_discretization_config_2d(
    range_ptr const &space_range_1, range_ptr const &space_range_2, std::size_t const &number_of_space_points_1,
    std::size_t const &number_of_space_points_2, range_ptr const &time_range, std::size_t const &number_of_time_points)
    : space_range_1_{space_range_1}, space_range_2_{space_range_2}, number_of_space_points_1_{number_of_space_points_1},
      number_of_space_points_2_{number_of_space_points_2}, time_range_{time_range}, number_of_time_points_{
                                                                                        number_of_time_points}
{
}

pde_discretization_config_2d::~pde_discretization_config_2d()
{
}

sptr_t<pde_discretization_config_1d> const pde_discretization_config_2d::pde_discretization_1() const
{
    return std::make_shared<pde_discretization_config_1d>(space_range_1_, number_of_space_points_1_, time_range_,
                                                          number_of_time_points_);
}

sptr_t<pde_discretization_config_1d> const pde_discretization_config_2d::pde_discretization_2() const
{
    return std::make_shared<pde_discretization_config_1d>(space_range_2_, number_of_space_points_2_, time_range_,
                                                          number_of_time_points_);
}

std::pair<range_ptr, range_ptr> const pde_discretization_config_2d::space_range() const
{
    return std::make_pair(space_range_1_, space_range_2_);
}

range_ptr const &pde_discretization_config_2d::time_range() const
{
    return time_range_;
}

std::pair<std::size_t, std::size_t> const pde_discretization_config_2d::number_of_space_points() const
{
    return std::make_pair(number_of_space_points_1_, number_of_space_points_2_);
}

std::size_t pde_discretization_config_2d::number_of_time_points() const
{
    return number_of_time_points_;
}

std::pair<double, double> pde_discretization_config_2d::space_step() const
{
    return std::make_pair(((space_range_1_->spread()) / static_cast<double>(number_of_space_points_1_ - 1)),
                          ((space_range_2_->spread()) / static_cast<double>(number_of_space_points_2_ - 1)));
}

double pde_discretization_config_2d::time_step() const
{
    return ((time_range_->spread()) / static_cast<double>(number_of_time_points_ - 1));
}

pde_discretization_config_3d::pde_discretization_config_3d(
    range_ptr const &space_range_1, range_ptr const &space_range_2, range_ptr const &space_range_3,
    std::size_t const &number_of_space_points_1, std::size_t const &number_of_space_points_2,
    std::size_t const &number_of_space_points_3, range_ptr const &time_range, std::size_t const &number_of_time_points)
    : space_range_1_{space_range_1}, space_range_2_{space_range_2}, space_range_3_{space_range_3},
      number_of_space_points_1_{number_of_space_points_1}, number_of_space_points_2_{number_of_space_points_2},
      number_of_space_points_3_{number_of_space_points_3}, time_range_{time_range}, number_of_time_points_{
                                                                                        number_of_time_points}
{
}

pde_discretization_config_3d::~pde_discretization_config_3d()
{
}

sptr_t<pde_discretization_config_1d> const pde_discretization_config_3d::pde_discretization_1() const
{
    return std::make_shared<pde_discretization_config_1d>(space_range_1_, number_of_space_points_1_, time_range_,
                                                          number_of_time_points_);
}

sptr_t<pde_discretization_config_1d> const pde_discretization_config_3d::pde_discretization_2() const
{
    return std::make_shared<pde_discretization_config_1d>(space_range_2_, number_of_space_points_2_, time_range_,
                                                          number_of_time_points_);
}

sptr_t<pde_discretization_config_1d> const pde_discretization_config_3d::pde_discretization_3() const
{
    return std::make_shared<pde_discretization_config_1d>(space_range_3_, number_of_space_points_3_, time_range_,
                                                          number_of_time_points_);
}

sptr_t<pde_discretization_config_2d> const pde_discretization_config_3d::pde_discretization_12() const
{
    return std::make_shared<pde_discretization_config_2d>(space_range_1_, space_range_2_, number_of_space_points_1_,
                                                          number_of_space_points_2_, time_range_,
                                                          number_of_time_points_);
}

sptr_t<pde_discretization_config_2d> const pde_discretization_config_3d::pde_discretization_21() const
{
    return std::make_shared<pde_discretization_config_2d>(space_range_2_, space_range_1_, number_of_space_points_2_,
                                                          number_of_space_points_1_, time_range_,
                                                          number_of_time_points_);
}

sptr_t<pde_discretization_config_2d> const pde_discretization_config_3d::pde_discretization_13() const
{
    return std::make_shared<pde_discretization_config_2d>(space_range_1_, space_range_3_, number_of_space_points_1_,
                                                          number_of_space_points_3_, time_range_,
                                                          number_of_time_points_);
}

sptr_t<pde_discretization_config_2d> const pde_discretization_config_3d::pde_discretization_31() const
{
    return std::make_shared<pde_discretization_config_2d>(space_range_3_, space_range_1_, number_of_space_points_3_,
                                                          number_of_space_points_1_, time_range_,
                                                          number_of_time_points_);
}

sptr_t<pde_discretization_config_2d> const pde_discretization_config_3d::pde_discretization_23() const
{
    return std::make_shared<pde_discretization_config_2d>(space_range_2_, space_range_3_, number_of_space_points_2_,
                                                          number_of_space_points_3_, time_range_,
                                                          number_of_time_points_);
}

sptr_t<pde_discretization_config_2d> const pde_discretization_config_3d::pde_discretization_32() const
{
    return std::make_shared<pde_discretization_config_2d>(space_range_3_, space_range_2_, number_of_space_points_3_,
                                                          number_of_space_points_2_, time_range_,
                                                          number_of_time_points_);
}

std::tuple<range_ptr, range_ptr, range_ptr> const pde_discretization_config_3d::space_range() const
{
    return std::make_tuple(space_range_1_, space_range_2_, space_range_3_);
}

range_ptr const &pde_discretization_config_3d::time_range() const
{
    return time_range_;
}

std::tuple<std::size_t, std::size_t, std::size_t> const pde_discretization_config_3d::number_of_space_points() const
{
    return std::make_tuple(number_of_space_points_1_, number_of_space_points_2_, number_of_space_points_3_);
}

std::size_t pde_discretization_config_3d::number_of_time_points() const
{
    return number_of_time_points_;
}

std::tuple<double, double, double> pde_discretization_config_3d::space_step() const
{
    return std::make_tuple(((space_range_1_->spread()) / static_cast<double>(number_of_space_points_1_ - 1)),
                           ((space_range_2_->spread()) / static_cast<double>(number_of_space_points_2_ - 1)),
                           ((space_range_3_->spread()) / static_cast<double>(number_of_space_points_3_ - 1)));
}

double pde_discretization_config_3d::time_step() const
{
    return ((time_range_->spread()) / static_cast<double>(number_of_time_points_ - 1));
}

} // namespace lss_pde_solvers
