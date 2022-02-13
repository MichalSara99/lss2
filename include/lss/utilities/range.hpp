#if !defined(_RANGE_HPP_)
#define _RANGE_HPP_

#include <functional>

#include "../../../common/lss_macros.hpp"
#include "../../../common/lss_range.hpp"

namespace lss
{

using range_ptr = lss_utility::range_ptr;
using range = lss_utility::range;

struct range_builder
{
  private:
    double l_, u_;

  public:
    LSS_API explicit range_builder();

    LSS_API range_builder &lower(double value);

    LSS_API range_builder &upper(double value);

    LSS_API range_ptr build();
};

} // namespace lss

#endif ///_RANGE_HPP_
