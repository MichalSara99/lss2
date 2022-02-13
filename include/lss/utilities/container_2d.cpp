#include "container_2d.hpp"

namespace lss
{

rmatrix_2d_builder::rmatrix_2d_builder()
{
}

rmatrix_2d_builder &rmatrix_2d_builder::rows(std::size_t rows)
{
    rows_ = rows;
    return *this;
}

rmatrix_2d_builder &rmatrix_2d_builder::columns(std::size_t columns)
{
    columns_ = columns;
    return *this;
}

rmatrix_2d_builder &rmatrix_2d_builder::default_value(double value)
{
    dvalue_ = value;
    return *this;
}

rmatrix_2d_ptr rmatrix_2d_builder::build()
{
    return std::make_shared<rmatrix_2d>(rows_, columns_, dvalue_);
}

cmatrix_2d_builder::cmatrix_2d_builder()
{
}

cmatrix_2d_builder &cmatrix_2d_builder::rows(std::size_t rows)
{
    rows_ = rows;
    return *this;
}

cmatrix_2d_builder &cmatrix_2d_builder::columns(std::size_t columns)
{
    columns_ = columns;
    return *this;
}

cmatrix_2d_builder &cmatrix_2d_builder::default_value(double value)
{
    dvalue_ = value;
    return *this;
}

cmatrix_2d_ptr cmatrix_2d_builder::build()
{
    return std::make_shared<cmatrix_2d>(rows_, columns_, dvalue_);
}

} // namespace lss
