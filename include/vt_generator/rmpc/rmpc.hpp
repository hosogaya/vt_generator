#pragma ocne

#include <vt_generator/type.hpp>

namespace vt_generator
{
namespace rmpc
{

class RMPC
{
public:
    RMPC(const Scalar m, const Scalar Iz, 
        const Scalar Cf, const Scalar Cr,
        const Scalar Lf, const Scalar Lr, 
        const std::vector<Scalar>& curvature,
        const std::vector<Scalar>& ds,
        const VecBound& bounds_input, 
        const std::vector<VecBound>& bounds_constraints,
        const int prediction_horizon, const int sampling_num);
    ~RMPC();

    bool ok() const {
        if (ds_.size() != curvature_.size()) return false;
        if (ds_.size() != bounds_constraints_.size()) return false;
        if (bounds_input_.size() != 2) return false;
        for (const auto& b: bounds_constraints_)
        {
            if (b.size() != 7)
            {
                return false;
            }
        }
        return true;
    }
    const std::vector<std::vector<Scalar>>& getOptState() const {return opt_state_;}
    const std::vector<std::pair<Scalar, Scalar>>& getOptInput() const {return opt_input_;}
    bool solve(const std::vector<Scalar>& init_state);
private:    
    std::vector<Scalar> calNextState(const std::vector<Scalar>& state, const std::pair<Scalar, Scalar>& input, const int index);
    bool evalConstratins(const std::vector<Scalar>& state, const std::pair<Scalar, Scalar>& input, const int index);
    Scalar evalCosts(const std::vector<Scalar>& state, const std::pair<Scalar, Scalar>& input, const int index);

    const VecBound bounds_input_;
    const std::vector<VecBound> bounds_constraints_;
    const int prediction_horizon_;
    const int sampling_num_;
    const std::vector<Scalar> curvature_;
    const std::vector<Scalar> ds_;

    // vehicle param
    const Scalar m_;
    const Scalar Iz_;
    const Scalar Cf_;
    const Scalar Cr_;
    const Scalar Lf_;
    const Scalar Lr_;

    // solution
    std::vector<std::vector<Scalar>> opt_state_;
    std::vector<std::pair<Scalar, Scalar>> opt_input_;
};

RMPC::RMPC(const Scalar m, const Scalar Iz, 
        const Scalar Cf, const Scalar Cr,
        const Scalar Lf, const Scalar Lr, 
        const std::vector<Scalar>& curvature,
        const std::vector<Scalar>& ds,
        const VecBound& bounds_input, 
        const std::vector<VecBound>& bounds_constraints,
        const int prediction_horizon, const int sampling_num)
: m_(m), Iz_(Iz), Cf_(Cf), Cr_(Cr), Lf_(Lf), Lr_(Lr), 
curvature_(curvature), ds_(ds), 
bounds_input_(bounds_input), bounds_constraints_(bounds_constraints), 
prediction_horizon_(prediction_horizon), sampling_num_(sampling_num)
{
    
}

RMPC::~RMPC() {}


bool RMPC::solve(const std::vector<Scalar>& init_state)
{
    if (init_state.size() != 5) 
    {
        std::cout << "init state has wrong size" << std::endl;
        return false;
    }
    std::vector<Scalar> acc_lists_(sampling_num_);
    std::vector<Scalar> steer_lists(sampling_num_);

    const auto& acc_bound = bounds_input_.at(0);
    const auto& steer_bound = bounds_input_.at(1);
    for (int i=0; i<sampling_num_; ++i)
    {
        acc_lists_.at(i) = acc_bound.lower_+(acc_bound.upper_ - acc_bound.lower_)*i/sampling_num_;
        steer_lists.at(i) = acc_bound.lower_+(acc_bound.upper_ - acc_bound.lower_)*i/sampling_num_;
    }

    std::vector<std::pair<Scalar, Scalar>> input_lists;
    for (size_t j=0; j<acc_lists_.size(); ++j)
    {
        for (size_t k=0; k<steer_lists.size(); ++k)
        {
            input_lists.emplace_back(std::make_pair(acc_lists_[j], steer_lists[k]));
        }
    }

    std::cout << "loop start" << std::endl;

    // loop for points
    std::vector<std::vector<Scalar>> opt_state{init_state};
    std::vector<std::pair<Scalar, Scalar>> opt_input;
    for (size_t p_index=0; p_index<ds_.size(); ++p_index)
    {
        Scalar opt_cost = 1e10;
        std::vector<int> local_opt_input_indexes(prediction_horizon_);  
        std::vector<int> input_indexes(prediction_horizon_);
        for (size_t i=0; i<prediction_horizon_; ++i) input_indexes.at(i) = 0;
        
        bool is_solution = false;
        std::cout << "search opt input for the point" << std::endl;
        while (input_indexes.at(0) <= input_lists.size())
        {
            std::vector<Scalar> state = opt_state.at(p_index);
            Scalar cost = 0.0;
            bool is_success = true;
            for (size_t i=0; i<prediction_horizon_; ++i) 
            {
                // constraints
                std::cout << "eval constraints" << std::endl;
                if (!evalConstratins(state, input_lists.at(input_indexes.at(i)), p_index+i))
                {
                    is_success = false;
                    break;
                }
                // evaluation
                std::cout << "eval cost" << std::endl;
                cost += evalCosts(state, input_lists.at(input_indexes.at(i)), p_index+i);
                // state equation
                std::cout << "cal next state" << std::endl;
                state = calNextState(state, input_lists.at(input_indexes.at(i)), p_index+i);
            }
            
            if (cost < opt_cost && is_success) 
            {
                local_opt_input_indexes = input_indexes;
                opt_cost = cost;
                is_solution = true;
            }
            ++input_indexes.at(input_indexes.size() - 1);
            for (int i=input_indexes.size(); i>=0; --i)
            {
                if (input_indexes.at(i) >= input_lists.size())
                {
                    input_indexes.at(i) = 0;
                    if (i>0) ++input_indexes.at(i-1);
                }
            }
        }
        if (!is_solution) 
        {
            std::cout << "There is no solution which stasfies the constraints" << std::endl;
            return false;
        }
        std::cout << "save the opt input and state" << std::endl;
        opt_input.push_back(input_lists.at(local_opt_input_indexes.at(0)));
        opt_state.push_back(calNextState(opt_state.at(p_index), input_lists.at(local_opt_input_indexes.at(0)), p_index));
    }
    // save the results
    opt_input_ = opt_input;
    opt_state_ = opt_state;

    return true;
}


std::vector<Scalar> RMPC::calNextState(const std::vector<Scalar>& state, const std::pair<Scalar, Scalar>& input, const int index)
{
    std::vector<Scalar> next(state.size());

    const auto& acc = input.first;
    const auto& steer = input.second;
    const auto& n = state.at(0);
    const auto& xi = state.at(1);
    const auto& vx = state.at(2);
    const auto& vy = state.at(3);
    const auto& w = state.at(4);
    const auto& curvature = curvature_.at(index);
    const auto& ds = ds_.at(index);
    
    auto v    = std::sqrt(std::pow(vx,2.0) + std::pow(vy, 2.0));
    auto beta = std::atan(vy/vx);
    auto dtds = (1.0 - n*curvature) / (v*std::cos(xi + beta));
    auto ff   = -Cf_*((vy - Lf_*w)/vx - steer);
    auto fr   = -Cr_*(vy - Lr_*w)/vx;

    next.at(0) = n  + (std::tan(xi + beta)*(1-n*curvature))*ds; 
    next.at(1) = xi + (w - curvature*v*std::cos(xi + beta)/ (1.0 - n*curvature))*dtds*ds;
    next.at(2) = vx + (acc   - ff*std::sin(steer)/m_ + vy*w)*dtds*ds;
    next.at(3) = vy + (fr/m_ + ff*std::cos(steer)/m_ - vx*w)*dtds*ds;
    next.at(4) = w  + (ff*Lf_*std::cos(steer) - fr*Lr_)*dtds*ds;

    return next;
}

bool RMPC::evalConstratins(const std::vector<Scalar>& state, const std::pair<Scalar, Scalar>& input, const int index)
{
    const auto& acc = input.first;
    const auto& steer = input.second;
    const auto& n = state.at(0);
    const auto& xi = state.at(1);
    const auto& vx = state.at(2);
    const auto& vy = state.at(3);
    const auto& w = state.at(4);

    auto v    = std::sqrt(std::pow(vx,2.0) + std::pow(vy, 2.0));
    auto beta = std::atan(vy/vx);
    auto ff   = -Cf_*((vy - Lf_*w)/vx - steer);
    auto fr   = -Cr_*(vy - Lr_*w)/vx;

    const auto& bounds = bounds_constraints_.at(index);

    std::vector<Scalar> constraints(bounds.size());
    constraints.at(0) = n;
    constraints.at(1) = xi;
    constraints.at(2) = vx;
    constraints.at(3) = vy;
    constraints.at(4) = w;
    constraints.at(5) = (ff + fr) /m_;
    constraints.at(6) = beta;

    for (size_t i=0; i<bounds.size(); ++i)
    {
        if (constraints.at(i) < bounds.at(i).lower_) return false;
        if (constraints.at(i) > bounds.at(i).upper_) return false;
    }

    return true;
}

Scalar RMPC::evalCosts(const std::vector<Scalar>& state, const std::pair<Scalar, Scalar>& input, const int index)
{
    const auto& acc = input.first;
    const auto& steer = input.second;
    const auto& n = state.at(0);
    const auto& xi = state.at(1);
    const auto& vx = state.at(2);
    const auto& vy = state.at(3);
    const auto& w = state.at(4);
    const auto& curvature = curvature_.at(index);
    const auto& ds = ds_.at(index);

    auto v    = std::sqrt(std::pow(vx,2.0) + std::pow(vy, 2.0));
    auto beta = std::atan(vy/vx);
    auto dtds = (1.0 - n*curvature) / (v*std::cos(xi + beta));
    auto ff   = -Cf_*((vy - Lf_*w)/vx - steer);
    auto fr   = -Cr_*(vy - Lr_*w)/vx;

    const auto& bounds = bounds_constraints_.at(index);
    std::vector<Scalar> constraints(bounds.size());
    constraints.at(0) = n;
    constraints.at(1) = xi;
    constraints.at(2) = vx;
    constraints.at(3) = vy;
    constraints.at(4) = w;
    constraints.at(5) = (ff + fr) /m_;
    constraints.at(6) = beta;

    Scalar value = dtds*ds;
    for (size_t i=0; i<constraints.size(); ++i)
    {
        auto& bound = bounds.at(i);
        auto& g = constraints.at(i);
        if (g - bound.lower_ < bound.upper_ - g)
        {
            value += (g - bound.lower_)/(bound.upper_ - bound.lower_);
        }
        else
        {
            value += (bound.upper_ - g)/(bound.upper_ - bound.lower_);
        }
    }

    return value;
}
}
}