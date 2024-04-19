#pragma once 

#include <vt_generator/type.hpp>

namespace vt_generator
{
inline Scalar distance(const std::vector<Scalar>& xm, const std::vector<Scalar>& ym, const int i1, const int i2)
{
    return std::sqrt(std::pow(xm[i1] - xm[i2], 2.0) + std::pow(ym[i1] - ym[i2], 2.0));
}

inline std::vector<int> thinout(const std::vector<Scalar>& xm, const std::vector<Scalar>& ym, const Scalar dist_thres)
{
    std::vector<int> indexes{0};
    Scalar dist = 0.0;
    Scalar total = 0.0;
    for (int i=0; i<xm.size(); ++i)
    {   
        total +=distance(xm, ym, i, (i+1)%xm.size()); 
        dist += distance(xm, ym, i, (i+1)%xm.size());
        if (dist > dist_thres)
        {
            indexes.push_back((i+1)%xm.size());
            dist = 0.0;
        };
    }
    std::cout << total << std::endl;
    return indexes;
}

inline Scalar calCurvature(const std::vector<Scalar>& xm, const std::vector<Scalar>& ym, const int index, int width = 1)
{
    int former = (index - width + xm.size())%xm.size();
    int later = (index + width)%xm.size();

    Scalar area4 = 2*((xm[index] - xm[former])*(ym[later] - ym[index])
                    - (ym[index] - ym[former])*(xm[later] - xm[index]));

    Scalar denominator = std::sqrt(std::pow(xm[former] - xm[index], 2.0) + std::pow(ym[former] - ym[index], 2.0))
                        +std::sqrt(std::pow(xm[later]  - xm[index], 2.0) + std::pow(ym[later]  - ym[index], 2.0))
                        +std::sqrt(std::pow(xm[former] - xm[later], 2.0) + std::pow(ym[former] - ym[later], 2.0));

    return area4/(denominator + 1e-8*std::log(1.0 + std::exp(-2*denominator*1e8)));
}

inline Scalar calDs(const std::vector<Scalar>& xm, const std::vector<Scalar>& ym, const int index)
{
    int later = (index + 1)%xm.size();
    Vector2 cur{xm[index], ym[index]};
    Vector2 lat{xm[later], ym[later]};

    return (lat - cur).norm();
}

inline Vector calNormalVector(const std::vector<Scalar>& xm, const std::vector<Scalar>& ym, const int index, const int width = 1)
{
    Vector2 normal;
    int former = (index - width + xm.size())%xm.size();
    int later = (index + width)%xm.size();
    Scalar area4 = 2*((xm[index] - xm[later])*(ym[index] - ym[former])
                    - (xm[index] - xm[former])*(ym[index] - xm[later]));

    Scalar x_numerator1 = (ym[index] - ym[later])* (std::pow(xm[index], 2.0)
                                                  - std::pow(xm[former], 2.0)
                                                  + std::pow(ym[index], 2.0)
                                                  - std::pow(xm[former], 2.0));
    
    Scalar x_numerator2 = (ym[index] - ym[former])*(std::pow(xm[index], 2.0)
                                                 - std::pow(xm[former], 2.0)
                                                 + std::pow(ym[index], 2.0)
                                                 - std::pow(xm[former], 2.0));

    Scalar y_numerator1 = (xm[index] - xm[later])*(std::pow(xm[index], 2.0)
                                                 - std::pow(xm[former], 2.0)
                                                 + std::pow(ym[index], 2.0)
                                                 - std::pow(ym[former], 2.0));

    Scalar y_numerator2 = (xm[index] - xm[former])*(std::pow(xm[index], 2.0)
                                                  - std::pow(xm[later], 2.0)
                                                  + std::pow(ym[index], 2.0)
                                                  - std::pow(ym[later], 2.0));

    Scalar center_x = (x_numerator1 - x_numerator2)/area4;
    Scalar center_y = (y_numerator1 - y_numerator2)/area4;

    normal.x() = xm[index] - center_x;
    normal.y() = ym[index] - center_y;

    // normal.normalize();

    // Vector2 tangent;
    // tangent.x() = normal.y();
    // tangent.y() =-normal.x();

    return normal.normalized();
}

}