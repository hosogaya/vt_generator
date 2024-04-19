#include <vt_generator/csv/csv_reader.hpp>
#include <vt_generator/util.hpp>
#include <vt_generator/vt_generator.hpp>

#include <vt_generator/csv/csv_writer.hpp>

#include <vt_generator/path_smoother/cost_curv_ad.hpp>
#include <vt_generator/path_smoother/variable.hpp>
#include <vt_generator/path_smoother/max_curvature.hpp>

using namespace vt_generator;

int main()
{
    std::vector<Scalar> xm, ym, tr_left, tr_right, outer_x, outer_y, inner_x, inner_y;
    std::vector<Scalar> x_smooth, y_smooth, curvature;
    std::vector<Vector2> normal;
    std::string csv_name = "../data/icra2023_1/icra2023_1_course_info.csv";
    const int width = 1;
    Scalar min_trun_radius = 4.0;
    csv::Reader csv_reader(csv_name);
    if (!csv_reader.isOpen()) return -1;

    // auto indexes = thinout(csv_reader.xm(), csv_reader.ym(), 0.01);
    std::vector<int> indexes(csv_reader.xm().size());
    for (int i=0; i<indexes.size(); ++i) indexes[i] = i;

    for (size_t i=0; i<indexes.size(); ++i)
    {
        xm.push_back(csv_reader.xm()            [indexes.at(i)]);
        ym.push_back(csv_reader.ym()            [indexes.at(i)]);
        tr_left.push_back(csv_reader.tr_left()  [indexes.at(i)]);
        tr_right.push_back(csv_reader.tr_right()[indexes.at(i)]);
        outer_x.push_back(csv_reader.outer_x()  [indexes.at(i)]);
        outer_y.push_back(csv_reader.outer_y()  [indexes.at(i)]);
        inner_x.push_back(csv_reader.inner_x()  [indexes.at(i)]);
        inner_y.push_back(csv_reader.inner_y()  [indexes.at(i)]);
    }

    // for (int i=0; i<xm.size(); ++i)
    // {
    //     int former = (i-1+xm.size())%xm.size();
    //     int later  = (i+1)%xm.size();
    //     Vector2 pf{xm.at(former), ym.at(former)};
    //     Vector2 p {xm.at(i),      ym.at(i)};
    //     Vector2 pl{xm.at(later),  ym.at(later)};
    //     std::cout << (pl - p).dot(p - pf)/(pl - p).norm()/(p - pf).norm();
    // }


    VecBound b_pos(indexes.size()*2); // xy
    for (size_t i=0; i<indexes.size(); ++i)
    {
        b_pos[2*i].lower_ = inner_x[i]; 
        b_pos[2*i].upper_ = outer_x[i]; 
        b_pos[2*i + 1].lower_ = inner_y[i];
        b_pos[2*i + 1].upper_ = outer_y[i];
    }
    
    VecBound b_const(indexes.size());
    for (size_t i=0; i<indexes.size(); ++i)
    {
        b_const[i].lower_ =-1.0/min_trun_radius;
        b_const[i].upper_ = 1.0/min_trun_radius;
    }

    using namespace vt_generator;
    path_smoother::decorder_.setBounds(b_pos);
    std::shared_ptr<path_smoother::Variable> var = std::make_shared<path_smoother::Variable>(indexes.size());
    std::shared_ptr<path_smoother::Cost>    cost = std::make_shared<path_smoother::Cost>(    indexes.size());
    std::shared_ptr<path_smoother::MaxCurvature> curv = std::make_shared<path_smoother::MaxCurvature>(indexes.size(), b_const);

    ifopt::Problem prob;
    prob.AddCostSet(cost);
    prob.AddVariableSet(var);
    prob.AddConstraintSet(curv);

    ifopt::IpoptSolver solver;
    solver.SetOption("print_level", 5);
    solver.SetOption("max_cpu_time", 1.0e20);
    solver.SetOption("max_iter", 10000);
    solver.SetOption("tol", 1.0e-4);
    solver.SetOption("acceptable_iter", 100);
    solver.SetOption("acceptable_tol", 1.0e-2);
    solver.SetOption("acceptable_constr_viol_tol", 0.001);
    solver.Solve(prob);

    // if (solver.GetReturnStatus() == 0 || solver.GetReturnStatus() == 1)
    {
        x_smooth.resize(xm.size());
        y_smooth.resize(ym.size());
        curvature.resize(xm.size());
        csv::Writer writer("smooth_path.csv");
        auto opt_var = prob.GetOptVariables()->GetValues();
        for (size_t i=0; i<xm.size(); ++i)
        {
            x_smooth[i] = path_smoother::decorder_.decodeX(opt_var(i), i);
            y_smooth[i] = path_smoother::decorder_.decodeY(opt_var(i), i);
        }
        std::cout << std::endl;
        for (size_t i=0; i<x_smooth.size(); ++i)
        {
            curvature[i] = calCurvature(x_smooth, y_smooth, i, width);
        }
        writer.writeResult(x_smooth, y_smooth, tr_left, tr_right, xm, ym,
                            outer_x, outer_y, inner_x, inner_y, curvature, x_smooth, xm.size());
    }

    return 1;
}