#include <vt_generator/csv/csv_reader.hpp>
#include <vt_generator/util.hpp>
#include <vt_generator/vt_generator.hpp>

#include <vt_generator/csv/csv_writer.hpp>

#include <vt_generator/path_smoother/cost.hpp>
#include <vt_generator/path_smoother/variable.hpp>

using namespace vt_generator;

int main()
{
    std::vector<Scalar> xm, ym, tr_left, tr_right, outer_x, outer_y, inner_x, inner_y;
    std::vector<Scalar> curvature;
    std::vector<Vector2> normal;
    std::string csv_name = "../data/icra20getADfun23_1/icra2023_1_course_info.csv";
    const int width = 10;
    {
        csv::Reader csv_reader(csv_name);
        if (!csv_reader.isOpen()) return -1;

        auto indexes = thinout(csv_reader.xm(), csv_reader.ym(), 0.02);
        // std::vector<int> indexes(csv_reader.xm().size());
        // for (int i=0; i<indexes.size(); ++i) indexes[i] = i;

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
    }

    return 1;
    const int horizon = xm.size();

    std::cout << "horizon: " << horizon << std::endl;
    
    // std::vector<Scalar> curvature;
    // for (int i=0; i<horizon; ++i)
    // {
    //     curvature.push_back(calCurvature(xm, ym, i, width));
    // }
    for (int i=0; i<horizon; ++i)
    {
        normal.push_back(calNormalVector(xm, ym, i, width));
        if (curvature[i] > 0) normal[i] = -normal[i];
    } 
    std::vector<Scalar> ds;
    Scalar total = 0.0;
    for (int i=0; i<horizon; ++i)
    {
        ds.push_back(calDs(xm, ym, i));
    }

    for (const auto k: curvature)
    {
        if (std::isnan(k)) 
        {
            std::cout << "there is nan in curvature" << std::endl;
            return 1;
        }
    }

    for (const auto k: ds)
    {
        if (std::isnan(k)) 
        {
            std::cout << "there is nan in ds" << std::endl;
            return 1;
        }
    }

    ifopt::Problem prob;
    {
        Scalar m = 3.95;
        Scalar Iz = 0.04712;
        Scalar Cf = 4.7180;
        Scalar Cr = 5.4562;
        Scalar Lf = 0.162;
        Scalar Lr = 0.162;

        Scalar lon_acc_min =-1.5e1;
        Scalar lon_acc_max = 1.5e1;
        Scalar lat_acc_min =-1.0e1;
        Scalar lat_acc_max = 1.0e1;
        Scalar steer_min =-M_PI*20.0/180;
        Scalar steer_max = M_PI*20.0/180;
        Scalar tread = 0.2; // m
        Scalar beta_min =-M_PI*10/180;
        Scalar beta_max = M_PI*10/180;

        VecBound b_var(horizon*denormalizer_.size());
        for (int i=0; i<horizon; ++i)
        {
            b_var[denormalizer_.acc(i)] = Bound(lon_acc_min, lon_acc_max);
            b_var[denormalizer_.steer(i)] = Bound(steer_min, steer_max);
            b_var[denormalizer_.n(i)] = Bound(-(tr_right[i] - tread/2.0), tr_left[i] - tread/2.0);
            b_var[denormalizer_.xi(i)] = Bound(-M_PI*30/180, M_PI*30/180);
            b_var[denormalizer_.vx(i)] = Bound(1.0, 20.0);
            b_var[denormalizer_.vy(i)] = Bound(-10.0, 10.0);
            b_var[denormalizer_.w(i)] =  Bound(-4*M_PI, 4*M_PI);
        }

        VecBound b_costrants(horizon*constraints_size_);
        for (int i=0; i<horizon; ++i)
        {
            // state equation
            b_costrants[i*constraints_size_ + 0] = ifopt::BoundZero;
            b_costrants[i*constraints_size_ + 1] = ifopt::BoundZero;
            b_costrants[i*constraints_size_ + 2] = ifopt::BoundZero;
            b_costrants[i*constraints_size_ + 3] = ifopt::BoundZero;
            b_costrants[i*constraints_size_ + 4] = ifopt::BoundZero;
            // lateral acc
            b_costrants[i*constraints_size_ + 5] = Bound(lat_acc_min, lat_acc_max);
            // slip angle
            b_costrants[i*constraints_size_ + 6] = Bound(beta_min, beta_max);
        }

        std::shared_ptr<Variable> var = std::make_shared<Variable>(horizon, b_var);
        std::shared_ptr<KBM> constraint = std::make_shared<KBM>(m, Iz, Cf, Cr, Lf, Lr, curvature, ds, b_costrants, horizon);
        std::shared_ptr<TimeCost> cost = std::make_shared<TimeCost>(horizon, curvature, ds);

        Vector x_init(horizon*denormalizer_.size());
        for (int i=0; i<horizon; ++i)
        {
            x_init(denormalizer_.acc(i)) = 0.0;
            x_init(denormalizer_.steer(i)) = 0.0;
            x_init(denormalizer_.n(i)) = 0.0;
            x_init(denormalizer_.xi(i)) = 0.0;
            x_init(denormalizer_.vx(i)) = -1.0;
            x_init(denormalizer_.vy(i)) = 0.0;
            x_init(denormalizer_.w(i)) = 0.0;
        }

        var->SetVariables(x_init);
        prob.AddVariableSet(var);
        prob.AddConstraintSet(constraint);
        prob.AddCostSet(cost);
    }

    ifopt::IpoptSolver solver;
    solver.SetOption("print_level", 5);
    solver.SetOption("max_cpu_time", 10.0e20);
    solver.SetOption("max_iter", 10000);
    solver.SetOption("limited_memory_max_history", 6); // default 6
    solver.SetOption("tol", 1.0e-4);
    solver.SetOption("constr_viol_tol", 0.0001);
    solver.SetOption("max_soft_resto_iters", 2);
    solver.SetOption("evaluate_orig_obj_at_resto_trial", "no");
    // solver.SetOption("least_square_init_primal", "yes");
    // solver.SetOption("least_square_init_duals", "yes");

    solver.SetOption("acceptable_iter", 100);
    solver.SetOption("acceptable_tol", 1.0e-2);
    solver.SetOption("acceptable_constr_viol_tol", 0.001);
    solver.Solve(prob);

    if (solver.GetReturnStatus() == 0)
    {
        csv::Writer writer("opt_x_success.csv");
        writer.writeResult(prob.GetOptVariables()->GetValues(), horizon);   
    }
    else if (solver.GetReturnStatus() == 1)
    {
        csv::Writer writer("opt_x_acceptable.csv");
        writer.writeResult(prob.GetOptVariables()->GetValues(), horizon);
    }
    else
    {
        csv::Writer writer("opt_x_failed.csv");
        writer.writeResult(prob.GetOptVariables()->GetValues(), horizon);
    }

    std::vector<Scalar> opt_var(horizon*denormalizer_.size());
    for (int i=0; i<prob.GetOptVariables()->GetValues().size(); ++i)
    {
        opt_var[i] = prob.GetOptVariables()->GetValues()(i);
    }
    
    // opt_x, opt_y, outer_width, inner_width
    // center_x, center_y, outer_x, outer_y inner_x, inner_y
    // curvature, ref_v
    std::vector<Scalar> opt_x  (xm.size()), opt_y  (ym.size()),
                        ref_v  (xm.size());
    for (size_t i=0; i<xm.size(); ++i)
    {
        auto normal = calNormalVector(xm, ym, i, width);
        if (curvature[i] > 0) normal = -normal;
        if (normal.norm() < 0.5) std::cout << "failed to cal normal" << std::endl;

        opt_x[i]   = xm.at(i) + normal.x()*denormalizer_.denormalizeN(opt_var, i);
        opt_y[i]   = ym.at(i) + normal.y()*denormalizer_.denormalizeN(opt_var, i);
        ref_v[i]   = denormalizer_.denormalizeVx(opt_var, i);
    }
    csv::Writer path_writer("opt_trajectory.csv");
    path_writer.writeResult(opt_x,     opt_y, 
                            tr_right,  tr_left, 
                            xm,        ym, 
                            outer_x,   outer_y, 
                            inner_x,   inner_y, 
                            curvature, ref_v, 
                            horizon);

    return 1;
}