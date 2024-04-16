#include <vt_generator/csv/csv_reader.hpp>
#include <vt_generator/util.hpp>

using namespace vt_generator;

int main()
{
    std::string csv_name;
    csv::Reader csv_reader(csv_name);
    
    std::vector<Scalar> curvature;
    for (int i=0; i<csv_reader.size(); ++i)
    {
        curvature.push_back(calCurvature(csv_reader.xm(), csv_reader.ym(), i));
    }

    Scalar m;
    Scalar Iz;
    Scalar Cf;
    Scalar Cr;
    Scalar Lf;
    Scalar Lr;  
    

    return 1;
}