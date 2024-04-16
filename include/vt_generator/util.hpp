#pragma once 

#include <vt_generator/type.hpp>

namespace vt_generator
{
inline Scalar calCurvature(const std::vector<Scalar>& xm, const std::vector<Scalar>& ym, const int index)
{
    int former = (index - 1 + xm.size())%xm.size();
    int later = (index + 1)%xm.size();

    Scalar area4 = 2*((xm[former] - xm[index])*(ym[later] - ym[index])
                    - (ym[former] - ym[index])*(xm[later] - xm[index]));

    Scalar denominator = std::sqrt(std::pow(xm[former] - xm[index], 2.0) + std::pow(ym[former] - ym[index], 2.0))
                        +std::sqrt(std::pow(xm[later]  - xm[index], 2.0) + std::pow(ym[later]  - ym[index], 2.0))
                        +std::sqrt(std::pow(xm[former] - xm[later], 2.0) + std::pow(ym[former] - ym[later], 2.0));

    return area4/denominator;
}

inline Vector calNormalVector(const std::vector<Scalar>& xm, const std::vector<Scalar>& ym, const int index)
{
    Vector2 normal;
    int former = (index - 1 + xm.size())%xm.size();
    int later = (index + 1)%xm.size();
    Scalar area4 = -2*((xm[former] - xm[index])*(ym[later] - ym[index])
                    - (ym[former] - ym[index])*(xm[later] - xm[index]));

    Scalar x_numerator1 = (ym[index] - ym[former])*(std::pow(xm[later], 2.0)
                                                  - std::pow(xm[index], 2.0)
                                                  + std::pow(ym[later], 2.0)
                                                  - std::pow(xm[index], 2.0));
    
    Scalar x_numerator2 = (ym[index] - ym[later])*(std::pow(xm[former], 2.0)
                                                 - std::pow(xm[index], 2.0)
                                                 + std::pow(ym[former], 2.0)
                                                 - std::pow(xm[index], 2.0));

    Scalar y_numerator1 = (xm[index] - xm[later])*(std::pow(xm[former], 2.0)
                                                 - std::pow(xm[index], 2.0)
                                                 + std::pow(ym[former], 2.0)
                                                 - std::pow(ym[index], 2.0));

    Scalar y_numerator2 = (xm[index] - xm[former])*(std::pow(xm[later], 2.0)
                                                  - std::pow(xm[index], 2.0)
                                                  + std::pow(ym[later], 2.0)
                                                  - std::pow(ym[index], 2.0));

    Scalar center_x = (x_numerator1 - x_numerator2)/area4;
    Scalar center_y = (y_numerator1 - y_numerator2)/area4;

    normal.x() = xm[index] - center_x;
    normal.y() = ym[index] - center_y;

    normal.normalize();

    Vector2 tangent;
    tangent.x() = normal.y();
    tangent.y() = normal.x();

    return normal.normalized();
}

}