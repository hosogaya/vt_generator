#include <vt_generator/rmpc/rmpc.hpp>
#include <vt_generator/util.hpp>
#include <vt_generator/csv/csv_reader_opt_curvature.hpp>
#include <vt_generator/csv/csv_writer_rmpc.hpp>
#include <vt_generator/csv/csv_writer.hpp>

using namespace vt_generator;
int main()
{
    std::string file_name = "line_modified.csv";
    csv::Reader reader(file_name);

    Scalar m = 3.95;
    Scalar Iz = 0.04712;
    Scalar Cf = 4.7180;
    Scalar Cr = 5.4562;
    Scalar Lf = 0.162;
    Scalar Lr = 0.162;

    Scalar lon_acc_min =-1.0e1;
    Scalar lon_acc_max = 1.0e1;
    Scalar lat_acc_min =-1.0e1;
    Scalar lat_acc_max = 1.0e1;
    Scalar steer_min =-M_PI*15.0/180;
    Scalar steer_max = M_PI*15.0/180;
    Scalar beta_min =-M_PI*10/180;
    Scalar beta_max = M_PI*10/180;

    std::vector<Scalar> ds;
    for (int i=0; i<reader.size(); ++i)
    {
        ds.push_back(calDs(reader.xm(), reader.ym(), i));
    }
    // make bounds
    VecBound bounds_input(2);
    bounds_input.at(0) = Bound(lon_acc_min, lon_acc_max);
    bounds_input.at(1) = Bound(steer_min, steer_max);
    
    std::vector<VecBound> bounds_constraints(reader.size());
    for (int i=0; i<reader.size(); ++i)
    {
        auto& bound = bounds_constraints.at(i);
        bound.resize(7);
        bound.at(0) = Bound(-reader.outer_w().at(i), reader.inner_w().at(i)); // normal
        bound.at(1) = Bound(-M_PI*20/180, M_PI*20/180); // xi
        bound.at(2) = Bound(1.0, 10.0); // vx
        bound.at(3) = Bound(-1.0, 1.0); // vy
        bound.at(4) = Bound(-2*M_PI, M_PI); // w(omega)
        bound.at(5) = Bound(lat_acc_min, lat_acc_max); // lateral acc
        bound.at(6) = Bound(beta_min, beta_max); // beta (slip angle)
    }

    rmpc::RMPC prob(m, Iz, Cf, Cr, Lf, Lr, reader.curvature(), ds, bounds_input, bounds_constraints, 5, 10);
    if (!prob.ok())
    {
        std::cout << "the problem not consistent" << std::endl;
        return -1;
    }

    std::vector<Scalar> init_x(5);
    init_x.at(0) = 0.0; // n
    init_x.at(1) = 0.0; // xi
    init_x.at(2) = 2.0; // vx
    init_x.at(3) = 0.0; // vy
    init_x.at(4) = 0.0; // w (omega)

    if (!prob.solve(init_x))
    {
        std::cout << "failed to solve" << std::endl;
    }

    csv::WriterRMPC writer("rmpc.csv");
    writer.writeResult(reader.outer_w(), reader.inner_w(), 
                    reader.outer_x(), reader.outer_y(), reader.inner_x(), reader.inner_x(), 
                    prob.getOptState(), prob.getOptInput(), reader.curvature());

    csv::Writer writer_tra("rmpc_tar.csv");
    std::vector<Vector2> normal;
    for (int i=0; i<reader.size(); ++i)
    {
        normal.push_back(calNormalVector2(reader.xm(), reader.ym(), i, 1));
    }
    std::vector<Scalar> xm, ym, ref_v;
    for (int i=0; i<reader.size(); ++i)
    {
        Vector2 center{reader.xm()[i], reader.ym()[i]};
        xm.push_back(center.x() + normal.at(i).x()*prob.getOptState().at(i).at(0));
        ym.push_back(center.y() + normal.at(i).y()*prob.getOptState().at(i).at(0));
        ref_v.push_back(prob.getOptState().at(i).at(2));
    }
    writer_tra.writeResult(xm, ym, reader.outer_w(), reader.inner_w(), 
                        reader.center_x(), reader.center_y(),
                        reader.outer_x(), reader.outer_y(), reader.inner_x(), reader.inner_y(), 
                        reader.curvature(), ref_v, reader.size());


    return 1;
}