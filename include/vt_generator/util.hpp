#pragma once 

#include <vt_generator/type.hpp>

namespace vt_generator
{
inline Scalar distance(const std::vector<Scalar>& xm, const std::vector<Scalar>& ym, const int i1, const int i2)
{
    return std::sqrt(std::pow(xm[i1] - xm[i2], 2.0) + std::pow(ym[i1] - ym[i2], 2.0));
}

inline Scalar outer(const Vector2& v1, const Vector2& v2)
{
    return v1.x()*v2.y() - v1.y()*v2.x();
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
                    - (xm[index] - xm[former])*(ym[index] - ym[later]));

    Scalar x_numerator1 = (ym[index] - ym[later])* (std::pow(xm[index], 2.0)
                                                  - std::pow(xm[former], 2.0)
                                                  + std::pow(ym[index], 2.0)
                                                  - std::pow(ym[former], 2.0));
    
    Scalar x_numerator2 = (ym[index] - ym[former])*(std::pow(xm[index], 2.0)
                                                  - std::pow(xm[later], 2.0)
                                                  + std::pow(ym[index], 2.0)
                                                  - std::pow(ym[later], 2.0));

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

    Vector2 direction;
    direction.x() = xm[later] - xm[index];
    direction.y() = ym[later] - ym[index];

    if (outer(direction, normal) < 0.0) normal = -normal;

    return normal.normalized();
}

inline Vector calNormalVector2(const std::vector<Scalar>& xm, const std::vector<Scalar>& ym, const int index, const int width = 1)
{
    Vector2 normal;
    int former = (index - width + xm.size())%xm.size();
    int later = (index + width)%xm.size();

    Vector2 tangent1;
    tangent1.x() = (xm[index] - xm[former]);
    tangent1.y() = (ym[index] - ym[former]);
    // tangent1.normalize();

    Vector2 tangent2;
    tangent2.x() = (xm[later] - xm[index]);
    tangent2.y() = (ym[later] - ym[index]);
    // tangent2.normalize();

    Vector2 t = (tangent1 + tangent2) / 2;
    // Vector2 t = tangent2;

    normal.x() =-t.y();
    normal.y() = t.x();
    // if (index == 24)
    // {
    //     std::cout << "former: " << xm[former] << ", " << ym[former] << std::endl;
    //     std::cout << "index: " << xm[index] << ", " << ym[index] << std::endl;
    //     std::cout << "later: " << xm[later] << ", " << ym[later] << std::endl;
    //     std::cout << "t1: " << tangent1.x() << ", " << tangent2.y() << std::endl;
    //     std::cout << "t2: " << tangent2.x() << ", " << tangent2.y() << std::endl;
    //     std::cout << "normal: " << normal.x() << ", " << normal.y() << std::endl;
    // }

    return normal.normalized();
}


inline Vector2 getClosestPointOnLine(const Vector2& normal, const Vector2& center, const std::vector<Scalar>& line_x, const std::vector<Scalar>& line_y)
{
    Scalar min_distance = 1e9;
    Vector2 min_distance_pos;
    for (int i=0; i<line_x.size(); ++i)
    {
        Vector2 sp, ep;
        sp.x() = line_x.at(i);
        sp.y() = line_y.at(i);
        ep.x() = line_x.at((i+1)%line_x.size());
        ep.x() = line_y.at((i+1)%line_y.size());

        Vector2 vs = sp - center;
        Vector2 ve = ep - center;

        // intersect
        // if (outer(normal, vs)*outer(normal, ve) >= 0) continue;
        // Scalar amp = (outer(normal, sp) - (outer(normal, center)))/outer(ep - sp, normal);
        // Vector2 cp = (ep - sp)*amp + sp;
        Vector2 cp = sp;

        if ((cp - center).norm() >= min_distance) continue;

        // update 
        min_distance = (cp - center).norm();
        min_distance_pos = cp;
    }
    
    return center + normal*min_distance;
}

}