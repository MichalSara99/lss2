#if !defined(_ROBIN_2D_HPP_)
#define _ROBIN_2D_HPP_

#include <functional>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../boundaries/lss_robin_boundary.hpp"
#include "../../../common/lss_macros.hpp"

namespace lss
{

using robin_2d_ptr = lss_boundary::robin_boundary_2d_ptr;
using robin_2d = lss_boundary::robin_boundary_2d;

struct robin_2d_builder
{
  private:
    std::function<double(double, double)> linear_value_;
    std::function<double(double, double)> value_;

  public:
    LSS_API explicit robin_2d_builder();

    LSS_API robin_2d_builder &linear_value(const std::function<double(double, double)> &linear_value);

    LSS_API robin_2d_builder &value(const std::function<double(double, double)> &value);

    LSS_API robin_2d_ptr build();
};

} // namespace lss

#endif ///_ROBIN_2D_HPP_
