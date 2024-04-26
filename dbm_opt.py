from casadi import Opti, sumsqr, norm_2, power, sqrt, exp, atan, sin, cos, atan2
from tabulate import tabulate
import pandas as pd
import numpy as np
import csv

def box_print(text) -> None:
    print(f"\n{tabulate(text, tablefmt='grid')}\n")
    
def distance(p1, p2):
    return sqrt(pow(p1[0]-p2[0], 2) + pow(p1[1]-p2[1], 2))

def calCurvature(p1, p2, p3):
    numerator = 2*(
        (p2[0]-p1[0])*(p3[1]-p2[1]) 
      - (p2[1]-p1[1])*(p3[0]-p2[0])
    )  # double the area of the triangle
    divisor = (
          distance(p1, p2)
        * distance(p3, p2)
        * distance(p1, p3)
    )  # perimeter of the triangle
    return numerator / divisor

class DbmPathOpt:
    def __init__(self, line_list: list, init_path: list) ->None:
        self.m = 3.95
        self.Iz = 0.04712
        self.Cr = 4.7180
        self.Cf = 5.4562
        self.Lf = 0.152
        self.Lr = 0.172
        self.line_list = line_list
        self.init_path = init_path
        
        self.min_turn_radius = 0.8
        self.n_cal_curv_point = 3
        # longitudinal acc
        self.max_lon_acc = 10.0
        self.min_lon_acc = -7
        # lateral acc
        self.min_lat_acc = -5.0
        self.max_lat_acc =  5.0
        # longitudinal velocity
        self.min_vx =1.5
        self.max_vx =7.0
        # lateral velocity
        self.min_vy = -3.0
        self.max_vy = 3.0
        
    def denomalize(self, norm_pos_list):
        denorm_pos_list = []
        for idx, line in enumerate(self.line_list):
            denorm_pos_list.append(line[0] + (line[1] - line[0])*norm_pos_list[idx])
            
        return denorm_pos_list
    
    def normalize(self, denorm_pos_list):
        norm_pos_list = []
        for idx, line in enumerate(self.line_list):
            norm_pos_list.append(
                norm_2(denorm_pos_list[idx] - line[0])
                /norm_2(line[1] - line[0])
            )
        return norm_pos_list
            
    def angle(self, v):
        return atan2(v[1], v[0])
        
    def optimize(self, optimizer: Opti):
        horizon = len(self.line_list)
        self.norm_pos_list = optimizer.variable(horizon)
        self.yaw_list = optimizer.variable(horizon)
        self.vx_list = optimizer.variable(horizon)
        self.vy_list = optimizer.variable(horizon)
        self.omega_list = optimizer.variable(horizon)
        self.acc_list = optimizer.variable(horizon)
        self.steer_list = optimizer.variable(horizon)
        self.dt_list = optimizer.variable(horizon)


        # set initial values
        init_norm_pos = self.normalize(self.init_path)
        init_yaw_list = []
        for i in range(horizon):
            init_yaw_list.append(self.angle(self.init_path[(i+1)%horizon] - self.init_path[i]))
        
        for i in range(horizon):
            if i==0: continue
            if init_yaw_list[i] - init_yaw_list[i-1] > np.pi:
                init_yaw_list[i] = init_yaw_list[i] - 2*np.pi
            elif init_yaw_list[i] - init_yaw_list[i-1] < -np.pi:
                init_yaw_list[i] = init_yaw_list[i] + 2*np.pi
                
        
        for i in range(horizon):
            optimizer.set_initial(self.norm_pos_list[i], init_norm_pos[i])
            optimizer.set_initial(self.yaw_list[i], init_yaw_list[i])
            optimizer.set_initial(self.vx_list[i], 2.0)
            optimizer.set_initial(self.vy_list[i], 0.0)
            optimizer.set_initial(self.omega_list[i], (init_yaw_list[(i+1)%horizon] - init_yaw_list[i]) / 0.05)
            optimizer.set_initial(self.acc_list[i], 0.0)
            optimizer.set_initial(self.steer_list[i], init_yaw_list[(i+1)%horizon] - init_yaw_list[i])
            optimizer.set_initial(self.dt_list[i], 0.05)
        
        self.denorm_pos_list = self.denomalize(self.norm_pos_list)
        evaluation = 0.0
        for i in range(horizon):
            slip_f = (self.vy_list[i] - self.Lf*self.omega_list[i])/self.vx_list[i] - self.steer_list[i]
            slip_r = (self.vy_list[i] - self.Lr*self.omega_list[i])/self.vx_list[i]
            ff = -self.Cf*slip_f
            fr = -self.Cr*slip_r
            
            evaluation += 2*self.dt_list[i]
            # evaluation += power(self.acc_list[(i+1)%horizon] - self.acc_list[i], 2.0) 
            evaluation += power(self.steer_list[(i+1)%horizon] - self.steer_list[i], 2.0)
            evaluation += power(self.vx_list[(i+1)%horizon] - self.vx_list[i], 2.0)
            # evaluation += power(self.vy_list[(i+1)%horizon] - self.vy_list[i], 2.0)
            
            ## state equation
            # x
            optimizer.subject_to(self.denorm_pos_list[(i+1)%horizon][0] ==
                                self.denorm_pos_list[i][0] +
                                (
                                    self.vx_list[i]*cos(self.yaw_list[i])
                                    -self.vy_list[i]*sin(self.yaw_list[i])
                                )*self.dt_list[i])
            
            # y
            optimizer.subject_to(self.denorm_pos_list[(i+1)%horizon][1] ==
                                self.denorm_pos_list[i][1] +
                                (
                                    self.vx_list[i]*sin(self.yaw_list[i])
                                    +self.vy_list[i]*cos(self.yaw_list[i])
                                )*self.dt_list[i])
            
            # yaw
            if i+1 == horizon:
                optimizer.subject_to(self.yaw_list[(i+1)%horizon] + 2*np.pi== 
                                    self.yaw_list[i] + 
                                    self.omega_list[i]*self.dt_list[i])
            else:
                optimizer.subject_to(self.yaw_list[(i+1)%horizon] ==
                                    self.yaw_list[i] + 
                                    self.omega_list[i]*self.dt_list[i])


            # vx
            optimizer.subject_to(self.vx_list[(i+1)%horizon] == 
                                 self.vx_list[i] + 
                                 (
                                     self.acc_list[i]
                                    -ff*sin(self.steer_list[i])/self.m
                                    +self.vy_list[i]*self.omega_list[i]
                                )*self.dt_list[i])
            
            # # vy
            optimizer.subject_to(self.vy_list[(i+1)%horizon] == 
                                 self.vy_list[i] + 
                                (
                                    fr/self.m
                                    +ff*cos(self.steer_list[i])/self.m
                                    -self.vx_list[i]*self.omega_list[i]
                                )*self.dt_list[i])
            # omega
            optimizer.subject_to(self.omega_list[(i+1)%horizon] == 
                                 self.omega_list[i] + 
                                (
                                    ff*self.Lf*cos(self.steer_list[i])
                                    -fr*self.Lr
                                )*self.dt_list[i]/self.Iz)
            
            # state
            optimizer.subject_to(self.norm_pos_list[i] >= 0.05)
            optimizer.subject_to(self.norm_pos_list[i] <= 0.95)
            # optimizer.subject_to(self.yaw_list[i] >-3.14)
            # optimizer.subject_to(self.yaw_list[i] < 3.14)
            optimizer.subject_to(self.vx_list[i] > self.min_vx)
            optimizer.subject_to(self.vx_list[i] < self.max_vx)
            optimizer.subject_to(self.vy_list[i] > self.min_vy)
            optimizer.subject_to(self.vy_list[i] < self.max_vy)
            optimizer.subject_to(self.omega_list[i] >-5*3.14)
            optimizer.subject_to(self.omega_list[i] < 5*3.14)
            optimizer.subject_to(self.dt_list[i] > 0)
            # input
            optimizer.subject_to(self.acc_list[i] > self.min_lon_acc)
            optimizer.subject_to(self.acc_list[i] < self.max_lon_acc)
            optimizer.subject_to(self.steer_list[i] >-0.436) # 25 deg
            optimizer.subject_to(self.steer_list[i] < 0.436) # 25 deg
            
            # slip angle
            sigmoid = 0.087266*1./1.+exp(-(self.vx_list[i] - 1.5))
            optimizer.subject_to(slip_f >-0.17453 - sigmoid) # 10 ~ 15deg
            optimizer.subject_to(slip_f < 0.17453 + sigmoid)
            optimizer.subject_to(slip_r >-0.17453 - sigmoid) # 10 ~ 15deg
            optimizer.subject_to(slip_r < 0.17453 + sigmoid)
            
            # lateral acc
            optimizer.subject_to((ff*cos(self.steer_list[i]) + fr) > self.min_lat_acc*self.m)
            optimizer.subject_to((ff*cos(self.steer_list[i]) + fr) < self.max_lat_acc*self.m)
            
            # curvature
            optimizer.subject_to(power(calCurvature(self.denorm_pos_list[(i-self.n_cal_curv_point)], 
                                                    self.denorm_pos_list[i], 
                                                    self.denorm_pos_list[(i+self.n_cal_curv_point)%len(self.denorm_pos_list)])
                                        , 2.0)
                                 < power(1.0 / self.min_turn_radius, 2.0))
            
            
        optimizer.minimize(evaluation)
        
        try:
            self.solution = optimizer.solve()
            total_time = sum(self.solution.value(self.dt_list))
            print("total time: {}".format(total_time))
        except RuntimeError:
            box_print([["Optimization failed"]])
            self.solution = optimizer.debug
            self.optimization_failed = True
            
            
    def solve(self):
        optimizer = Opti()
        optimizer.solver(
            "ipopt", 
            # casadi option
            {"print_time": True},
            # solver option
            {
                "max_iter": 5000,
                "print_level": 5, 
                "sb": "yes", 
                "tol": 1e-5,
                "print_frequency_time": 1,
                "evaluate_orig_obj_at_resto_trial": "no"
            }
        )
        
        self.optimize(optimizer=optimizer)
        
    def save_result(self, file_name):
        with open(file_name, "w", newline="") as f:
            labels = [
                "opt_x",
                "opt_y",
                "outer_width",
                "inner_width",
                "center_x",
                "center_y",
                "outer_x",
                "outer_y",
                "inner_x",
                "inner_y",
                "curvature",
                "ref_v",
                "yaw",
                "vx",
                "vy",
                "omega",
                "acc",
                "steer",
                "dt",
            ]
            writer = csv.DictWriter(f, fieldnames=labels)
            writer.writeheader()
            pos_list = self.denomalize(self.solution.value(self.norm_pos_list))
            vel_list = self.solution.value(self.vx_list)
            curvature_list = []
            for i, pos in enumerate(pos_list):
                curvature_list.append(calCurvature(pos_list[i-1], pos_list[i], pos_list[(i+1)%len(pos_list)]))
                
            yaw_list = self.solution.value(self.yaw_list)
            vy_list = self.solution.value(self.vy_list)
            omega_list = self.solution.value(self.omega_list)
            acc_list = self.solution.value(self.acc_list)
            steer_list = self.solution.value(self.steer_list)
            dt_list = self.solution.value(self.dt_list)
            for (opt_point, line, vel, center, curv, 
                 yaw, vx, vy, omega, acc, steer, dt) in zip(pos_list, self.line_list, vel_list, self.init_path, curvature_list, 
                                                            yaw_list, vel_list, vy_list,omega_list, acc_list, steer_list, dt_list):
                data = {
                    "opt_x": opt_point[0],
                    "opt_y": opt_point[1],
                    "outer_width": distance(center, line[1]),
                    "inner_width": distance(center, line[0]),
                    "center_x": center[0],
                    "center_y": center[1],
                    "outer_x": line[1][0],
                    "outer_y": line[1][1],
                    "inner_x": line[0][0],
                    "inner_y": line[0][1],
                    "curvature": curv,
                    "ref_v": vel,
                    "yaw": yaw,
                    "vx": vx,
                    "vy": vy,
                    "omega": omega,
                    "acc": acc,
                    "steer": steer,
                    "dt": dt,
                }        
                writer.writerow(data)
                
        
    
def read_csv(file_name):
    df = pd.read_csv(file_name)
    
    path = []
    line_pos = []
    for i in range(df.shape[0]):
        path.append(np.array([df["opt_x"][i], df["opt_y"][i]]))
        line_pos.append([
                        np.array([df["inner_x"][i], df["inner_y"][i]]), 
                        np.array([df["outer_x"][i], df["outer_y"][i]])
                        ])
        
    return path, line_pos


def main():
    file_name = "build/line_modified.csv"
    path, line_pos = read_csv(file_name=file_name)
    
    prob = DbmPathOpt(line_list=line_pos, init_path=path)
    prob.solve()
    prob.save_result("ego_ref_waypoint.csv")

if __name__ == "__main__":
    main()