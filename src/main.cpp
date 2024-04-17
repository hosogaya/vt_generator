#include <vt_generator/csv/csv_reader.hpp>
#include <vt_generator/util.hpp>
#include <vt_generator/vt_generator.hpp>

using namespace vt_generator;

int main()
{
    std::vector<Scalar> xm, ym, tr_left, tr_right;
    {
        std::string csv_name = "../icra2023_1_course_info.csv";
        csv::Reader csv_reader(csv_name);

        auto indexes = thinout(csv_reader.xm(), csv_reader.ym(), 0.05);

        for (size_t i=0; i<indexes.size(); ++i)
        {
            xm.push_back(csv_reader.xm()[indexes[i]]);
            ym.push_back(csv_reader.ym()[indexes[i]]);
            tr_left.push_back(csv_reader.tr_left()[indexes[i]]);
            tr_right.push_back(csv_reader.tr_right()[indexes[i]]);
        }
    }

    const int horizon = xm.size();

    std::cout << "horizon: " << horizon << std::endl;
    
    std::vector<Scalar> curvature;
    for (int i=0; i<horizon; ++i)
    {
        curvature.push_back(calCurvature(xm, ym, i, 2));
    }
    std::vector<Scalar> ds;
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

        Scalar lon_acc_min =-10.0;
        Scalar lon_acc_max = 10.0;
        Scalar lat_acc_min =-1.0e4;
        Scalar lat_acc_max = 1.0e4;
        Scalar steer_min =-M_PI*45.0/180;
        Scalar steer_max = M_PI*45.0/180;
        Scalar tread = 0.2; // m
        Scalar beta_min =-M_PI*10/180;
        Scalar beta_max = M_PI*10/180;

        Scalar jark_min =-100.0;
        Scalar jark_max = 100.0;
        Scalar steer_diff_min =-M_PI*90/180;
        Scalar steer_diff_max = M_PI*90/180;

        VecBound b_var(horizon*denormalizer_.size());
        for (int i=0; i<horizon; ++i)
        {
            b_var[denormalizer_.acc(i)] = Bound(lon_acc_min, lon_acc_max);
            b_var[denormalizer_.steer(i)] = Bound(steer_min, steer_max);
            b_var[denormalizer_.n(i)] = Bound(-(tr_left[i] - tread/2.0), tr_right[i] - tread/2.0);
            b_var[denormalizer_.xi(i)] = Bound(-M_PI*90/180, M_PI*90/180);
            b_var[denormalizer_.vx(i)] = Bound(0.0, 20.0);
            b_var[denormalizer_.vy(i)] = Bound(-100.0, 100.0);
            b_var[denormalizer_.w(i)] =  Bound(-100*M_PI, 100*M_PI);
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
            // slip angle
            b_costrants[i*constraints_size_ + 5] = Bound(beta_min, beta_max);
            // lateral acc
            b_costrants[i*constraints_size_ + 6] = Bound(lat_acc_min, lat_acc_max);
            // input 
            // b_costrants[i*9 + 7] = Bound(jark_min, jark_max);
            // b_costrants[i*9 + 8] = Bound(steer_diff_min, steer_diff_max);
        }

        std::cout << "make var, const, cost" << std::endl;
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
            x_init(denormalizer_.vx(i)) = 1.0;
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
    solver.SetOption("max_iter", 3000);
    solver.SetOption("limited_memory_max_history", 4); // default 6
    solver.SetOption("tol", 1.0e-5);
    solver.Solve(prob);
    
    return 1;
}