#include "matrix_2d.hpp"

namespace lss
{

matrix_2d_builder::matrix_2d_builder()
{
}

matrix_2d_builder &matrix_2d_builder::rows(std::size_t rows)
{
    rows_ = rows;
    return *this;
}

matrix_2d_builder &matrix_2d_builder::columns(std::size_t columns)
{
    columns_ = columns;
    return *this;
}

matrix_2d_builder &matrix_2d_builder::default_value(double value)
{
    dvalue_ = value;
    return *this;
}

matrix_2d_ptr matrix_2d_builder::build()
{
    return std::make_shared<matrix_2d>(rows_, columns_, dvalue_);
}

} // namespace lss
