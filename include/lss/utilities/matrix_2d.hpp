#if !defined(_MATRIX_2D_HPP_)
#define _MATRIX_2D_HPP_

#include <functional>

#include "../../../common/lss_macros.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../containers/lss_matrix_2d.hpp"

namespace lss
{

using matrix_2d = lss_containers::matrix_2d;
using matrix_2d_ptr = lss_containers::matrix_2d_ptr;

struct matrix_2d_builder
{
  private:
    std::size_t rows_;
    std::size_t columns_;
    double dvalue_;

  public:
    LSS_API explicit matrix_2d_builder();

    LSS_API matrix_2d_builder &rows(std::size_t rows);

    LSS_API matrix_2d_builder &columns(std::size_t columns);

    LSS_API matrix_2d_builder &default_value(double value);

    LSS_API matrix_2d_ptr build();
};

} // namespace lss

#endif ///_MATRIX_2D_HPP_
