import cv2
import csv
import fire
import os
from course_detection import detect_course, get_yaml_obj, check_paths
from casadi import Opti, sumsqr, norm_2, power, sqrt, exp
from tabulate import tabulate
from math import dist, floor
from numpy import ndarray, linalg, array, flipud, uint16
import numpy as np
from typing import Tuple, List


class PathOptimizer:
    def __init__(self, lateral_line_list, max_iter=100000, tolerance=0.1) -> None:
        self.line_list = lateral_line_list
        self.max_iter = max_iter
        self.tolerance = tolerance
        self.optimization_failed = False

    def set_params(
        self,
        max_vel,
        lon_acc_max,
        lon_acc_min,
        lat_acc_abs_max,
        use_min_turn_r,
        min_turn_r,
        n_curv_calc,
        use_dir_lock,
        use_dec,
        use_dyn_lat,
        lat_acc_slope_center,
        lat_acc_slope_inclination,
        lat_acc_max_divisor,
        use_dyn_dec,
        dec_slope_center,
        dec_slope_inclination,
        dec_max_divisor,
    ) -> None:
        self.max_vel = max_vel
        self.lon_acc_max = lon_acc_max
        self.lon_acc_min = lon_acc_min
        self.lat_acc_abs_max = lat_acc_abs_max
        self.use_min_turn_r = use_min_turn_r
        self.min_turn_r = min_turn_r
        self.n_curv_calc = n_curv_calc
        self.use_dir_lock = use_dir_lock
        self.use_dec = use_dec
        self.use_dyn_lat = use_dyn_lat
        self.lat_acc_slope_center = lat_acc_slope_center
        self.lat_acc_slope_inclination = lat_acc_slope_inclination
        self.lat_acc_max_divisor = lat_acc_max_divisor
        self.use_dyn_dec = use_dyn_dec
        self.dec_slope_center = dec_slope_center
        self.dec_slope_inclination = dec_slope_inclination
        self.dec_max_divisor = dec_max_divisor

    def set_max_iter(self, max_iter: int) -> None:
        self.max_iter = max_iter

    def set_tolerance(self, tolerance: float) -> None:
        self.tolerance = tolerance

    def get_curvature(self, point_s, point_m, point_e) -> float:
        numerator = (point_s[0] - point_m[0]) * (point_s[1] - point_e[1]) - (
            point_s[1] - point_m[1]
        ) * (
            point_s[0] - point_e[0]
        )  # double the area of the triangle
        divisor = (
            norm_2(point_s - point_m)
            * norm_2(point_m - point_e)
            * norm_2(point_s - point_e)
        )  # perimeter of the triangle
        return 2 * numerator / divisor

    def denormalize(self, norm_pos_list) -> list:
        denorm_pos_list = []
        for idx, line in enumerate(self.line_list):
            denorm_pos_list.append(
                line[0] * (1 - norm_pos_list[idx]) + line[1] * norm_pos_list[idx]
            )
        return denorm_pos_list

    def optimize_curvature(self, optimizer: Opti):
        # declare constants
        N = len(self.line_list)

        # declare variables
        self.norm_pos_list = optimizer.variable(N)  # normalized position on the line
        self.curv_list = optimizer.variable(N)  # curvature on the point

        # set initial values
        optimizer.set_initial(self.norm_pos_list, 0.5)

        evaluation = 0
        denorm_pos_list = self.denormalize(self.norm_pos_list)
        for i in range(N):
            # declare objective
            evaluation += pow(self.curv_list[i] - self.curv_list[i-1], 2)# + pow(self.curv_list[i], 2.0)

            # declare constraints
            optimizer.subject_to(
                self.curv_list[i]
                == self.get_curvature(
                    denorm_pos_list[(i + self.n_curv_calc) % N],
                    denorm_pos_list[i],
                    denorm_pos_list[i - self.n_curv_calc],
                )
            )

            optimizer.subject_to(self.norm_pos_list[i] >= 0.3)
            optimizer.subject_to(self.norm_pos_list[i] <=0.7)
            optimizer.subject_to(
                pow(self.curv_list[i], 2.0) <= pow(1 / self.min_turn_r, 2.0)
            )

        optimizer.minimize(evaluation)
        try:
            self.solution = optimizer.solve()
        except RuntimeError:
            box_print([["Optimization failed"]])
            self.solution = optimizer.debug
            self.optimization_failed = True
            
    def optimize_lap_time(self, optimizer: Opti):
        # declare constants
        N = len(self.line_list)

        # declare variables
        self.norm_pos_list = optimizer.variable(N)  # normalized position on the line
        self.vel_list = optimizer.variable(N)  # velocity on the point
        self.curv_list = optimizer.variable(N)# curvature on the point
        self.long_acc_list = optimizer.variable(N) # longitudinal acceleration on the point
        self.lat_acc_list = optimizer.variable(N) # lateral acceleration on the point

        # set initial values
        optimizer.set_initial(self.norm_pos_list, 0.5)
        optimizer.set_initial(self.vel_list, 1.0)

        lap_time = 0
        denorm_pos_list = self.denormalize(self.norm_pos_list)
        for i in range(N):
            # declare objective
            lap_time += (
                (
                    norm_2(denorm_pos_list[(i + 1) % N] - denorm_pos_list[i])
                    + norm_2(denorm_pos_list[i] - denorm_pos_list[i - 1])
                )
                / 2
            ) / self.vel_list[i]

            # declare constraints
            optimizer.subject_to(self.norm_pos_list[i] >= 0)
            optimizer.subject_to(self.norm_pos_list[i] <= 1)
            optimizer.subject_to(self.vel_list[i] > 0)
            optimizer.subject_to(self.vel_list[i] <= self.max_vel)
            optimizer.subject_to(self.long_acc_list[i] <= self.lon_acc_max)
            optimizer.subject_to(self.long_acc_list[i] >= self.lon_acc_min)
            if self.use_dyn_lat:
                divisor = 1 + (self.lat_acc_max_divisor - 1) / (
                    1
                    + exp(
                        -self.lat_acc_slope_inclination
                        * (self.vel_list[i] - self.lat_acc_slope_center)
                    )
                )
                optimizer.subject_to(
                    pow(self.lat_acc_list[i], 2.0) <= pow(self.lat_acc_abs_max / divisor, 2.0)
                )
            else:
                optimizer.subject_to(
                    pow(self.lat_acc_list[i], 2.0) <= pow(self.lat_acc_abs_max, 2.0)
                )

            if self.use_dyn_dec:
                divisor = 1 + (self.dec_max_divisor - 1) / (
                    1
                    + exp(
                        -self.dec_slope_inclination
                        * (self.vel_list[i] - self.dec_slope_center)
                    )
                )
                optimizer.subject_to(self.long_acc_list[i] >= self.lon_acc_min / divisor)
            else:
                optimizer.subject_to(self.long_acc_list[i] >= self.lon_acc_min)

            optimizer.subject_to(
                self.long_acc_list[i]
                == (pow(self.vel_list[i],2.0) - pow(self.vel_list[i - 1],2.0))
                / (2 * norm_2(denorm_pos_list[i] - denorm_pos_list[i - 1]))
            )

            if self.use_dir_lock:
                vector_1 = denorm_pos_list[i] - denorm_pos_list[i - 1]
                vector_2 = denorm_pos_list[i - 1] - denorm_pos_list[i - 2]
                dot = vector_1[0] * vector_2[0] + vector_1[1] * vector_2[1]
                optimizer.subject_to(dot > 0.0)


            optimizer.subject_to(
                self.curv_list[i]
                == self.get_curvature(
                    denorm_pos_list[(i + self.n_curv_calc + 1) % N],
                    denorm_pos_list[ i],
                    denorm_pos_list[ i - self.n_curv_calc - 1],
                )
            )
            optimizer.subject_to(
                self.lat_acc_list[i]
                == pow(self.vel_list[i], 2.0) * self.curv_list[i]
            )
            if self.use_dec:
                optimizer.subject_to(
                    self.long_acc_list[i] / (self.lon_acc_max)
                    < sqrt(1 - pow(self.lat_acc_list[i] / (self.lat_acc_abs_max), 2.0))
                )
                optimizer.subject_to(
                    self.long_acc_list[i] / (self.lon_acc_min)
                    > -sqrt(1 - pow(self.lat_acc_list[i] / (self.lat_acc_abs_max), 2.0))
                )
            else:
                optimizer.subject_to(
                      pow(self.long_acc_list[i]   / (self.lon_acc_max),     2.0)
                    + pow(self.lat_acc_list[i] / (self.lat_acc_abs_max), 2.0)
                    <= 1.0
                )
            if self.use_min_turn_r:
                optimizer.subject_to(
                    pow(self.curv_list[i],2.0) <= pow(1 / (self.min_turn_r* (1.0 - pow(self.vel_list[i]/self.max_vel, 2.0))), 2.0)
                )

        optimizer.minimize(lap_time)

        try:
            self.solution = optimizer.solve()
        except RuntimeError:
            box_print([["Optimization failed"]])
            self.solution = optimizer.debug
            self.optimization_failed = True

    def solve(self, method="lap time") -> None:
        optimizer = Opti()
        optimizer.solver(
            "ipopt",  # select solver
            {"print_time": True},  # casadi options
            {
                "max_iter": self.max_iter,
                "print_level": 5,
                "sb": "yes",
                "tol": self.tolerance,
                "print_frequency_time": 1,
            },  # solver options
        )

        if method == "lap time":
            self.optimize_lap_time(optimizer)
        elif method == "curvature":
            self.optimize_curvature(optimizer)

    def get_results(self, method="lap time") -> Tuple[List[ndarray], List[float]]:
        norm_pos_list = self.solution.value(self.norm_pos_list)
        curv_list = self.solution.value(self.curv_list)
        if method == "lap time":
            vel_list = self.solution.value(self.vel_list)
        else:
            vel_list = []
            for i in range(len(curv_list)):
                vel_list.append(0.0)
        optimized_path = []
        for idx, (norm_pos, line) in enumerate(zip(norm_pos_list, self.line_list)):
            optimized_path.append(line[0] * (1 - norm_pos) + line[1] * norm_pos)
        max_curv_list = curv_list
        # try:
        #     for curv_vals in curv_list:
        #         max_curv_list.append(max(curv_vals, key=abs))
        # except:
        #     max_curv_list.append(curv_vals)
        return optimized_path, max_curv_list, vel_list

    def print_stats(self) -> None:
        stats = self.solution.stats()
        return_status = stats["return_status"]
        iter_count = stats["iter_count"]
        t_proc = stats["t_proc_total"]
        t_wall = stats["t_wall_total"]
        print(
            f"\nstatus: {return_status},   iteration count: {iter_count},   processor time: {t_proc:.5f}s,   wall time: {t_wall:.5f}s\n"
        )


def box_print(text) -> None:
    print(f"\n{tabulate(text, tablefmt='grid')}\n")


def save_course_info(
    file_path, opt_path, lateral_line_list, curv_list, vel_list, no_margined_line_list
) -> None:
    with open(file_path, "w", newline="") as f:
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
            "true_outer_x",
            "true_outer_y",
            "true_inner_x",
            "true_inner_y",
        ]
        writer = csv.DictWriter(f, fieldnames=labels)
        writer.writeheader()
        for opt_point, lat_line, curvature, velocity, nm_line in zip(
            opt_path, lateral_line_list, curv_list, vel_list, no_margined_line_list
        ):
            outer_width = dist(opt_point, lat_line[0])
            inner_width = dist(opt_point, lat_line[1])
            center_x = (lat_line[0][0] + lat_line[1][0]) / 2
            center_y = (lat_line[0][1] + lat_line[1][1]) / 2
            data = {
                "opt_x": opt_point[0],
                "opt_y": opt_point[1],
                "outer_width": outer_width,
                "inner_width": inner_width,
                "center_x": center_x,
                "center_y": center_y,
                "outer_x": lat_line[0][0],
                "outer_y": lat_line[0][1],
                "inner_x": lat_line[1][0],
                "inner_y": lat_line[1][1],
                "curvature": curvature,
                "ref_v": velocity,
                "true_outer_x": nm_line[0][0],
                "true_outer_y": nm_line[0][1],
                "true_inner_x": nm_line[1][0],
                "true_inner_y": nm_line[1][1],
            }
            writer.writerow(data)


def draw_arrow(
    img, pt_1, pt_2, length=30, bgr=(0, 255, 0), thickness=2, tip_length=0.2
) -> ndarray:
    pt_1 = array(pt_1)
    pt_2 = array(pt_2)
    delta = (pt_2 - pt_1) / linalg.norm(pt_2 - pt_1)
    vert = flipud(delta)
    pt_2 = pt_1 + delta * length
    pt_1 += vert * length * tip_length
    pt_2 += vert * length * tip_length
    pt_1 = pt_1.astype(uint16).tolist()
    pt_2 = pt_2.astype(uint16).tolist()
    cv2.arrowedLine(
        img=img,
        pt1=pt_1,
        pt2=pt_2,
        color=bgr,
        thickness=thickness,
        tipLength=tip_length,
    )
    return img


def optimize_path(map_name: str = None) -> Tuple[List[ndarray], List[float]]:
    map_yaml_file, out_dir, config_file = check_paths(map_name)

    # load course prep parameters
    path_yaml_obj = get_yaml_obj(config_file)
    method = path_yaml_obj["path_optimization"]["method"]
    params_list = path_yaml_obj["path_optimization"]["speed_params"]
    use_dec = path_yaml_obj["path_optimization"]["use_deceleration"]
    use_min_turn_r = path_yaml_obj["path_optimization"]["use_minimum_turning_radius"]
    min_turn_r = path_yaml_obj["path_optimization"]["minimum_turning_radius"]
    use_dir_lock = path_yaml_obj["path_optimization"]["use_direction_lock"]
    n_curv_calc = path_yaml_obj["path_optimization"]["curvature_calculation_points"]
    max_iter = path_yaml_obj["path_optimization"]["max_iterations"]
    tolerance = path_yaml_obj["path_optimization"]["tolerance"]
    min_vel = path_yaml_obj["path_optimization"]["minimum_velocity"]
    use_dyn_lat = path_yaml_obj["path_optimization"]["use_dynamic_lat_acc"]
    lat_acc_slope_center = path_yaml_obj["path_optimization"]["lat_acc_slope_center"]
    lat_acc_slope_inclination = path_yaml_obj["path_optimization"][
        "lat_acc_slope_inclination"
    ]
    lat_acc_max_divisor = path_yaml_obj["path_optimization"]["lat_acc_max_divisor"]
    use_dyn_dec = path_yaml_obj["path_optimization"]["use_dynamic_dec"]
    dec_slope_center = path_yaml_obj["path_optimization"]["dec_slope_center"]
    dec_slope_inclination = path_yaml_obj["path_optimization"]["dec_slope_inclination"]
    dec_max_divisor = path_yaml_obj["path_optimization"]["dec_max_divisor"]
    bgr_list = path_yaml_obj["visualization"]["bgr_list"]
    line_thickness = path_yaml_obj["visualization"]["line_thickness"]
    m_crop_length = path_yaml_obj["course_prep"]["crop_length"]
    

    # preperate course
    lat_line_creator, line_list_no_margin = detect_course(map_name, 0.0)
    lat_line_creator, line_list_thinned = detect_course(map_name, m_crop_length)

    print("margined lat_line_list size: {}".format(len(line_list_thinned)))
    print("no marginedlat_line_list size: {}".format(len(line_list_no_margin)))
    
    if (len(line_list_no_margin) > len(line_list_thinned)):
        index_list = np.linspace(0, len(line_list_no_margin)-1, len(line_list_thinned)).astype(int)
        adjusted_list = [line for idx, line in enumerate(line_list_no_margin)  if idx in index_list]
        print(len(adjusted_list))

    # optimize path
    path_optimizer = PathOptimizer(line_list_thinned)
    path_optimizer.set_max_iter(max_iter)
    path_optimizer.set_tolerance(tolerance)
    opt_path_list: List[List[ndarray]] = []
    curv_lists: List[float] = []
    vel_lists: List[float] = []

    if not use_dec:
        for params in params_list:
            params[2] = -params[1]  # lon_acc_min = -lon_acc_max

    for max_vel, lon_acc_max, lon_acc_min, lat_acc_abs_max in params_list:
        path_optimizer.set_params(
            max_vel,
            lon_acc_max,
            lon_acc_min,
            lat_acc_abs_max,
            use_min_turn_r,
            min_turn_r,
            n_curv_calc,
            use_dir_lock,
            use_dec,
            use_dyn_lat,
            lat_acc_slope_center,
            lat_acc_slope_inclination,
            lat_acc_max_divisor,
            use_dyn_dec,
            dec_slope_center,
            dec_slope_inclination,
            dec_max_divisor,
        )
        if method == "lap time":
            box_print(
                [
                    [f"method: {method}"],
                    [
                        f"max vel: {max_vel}m/s   long acc: {lon_acc_max}m/s^2   long dec: {lon_acc_min}m/s^2   lat acc: {lat_acc_abs_max}m/s^2"
                    ],
                ]
            )
        elif method == "curvature":
            box_print(
                [
                    [f"method: {method}"],
                ]
            )
        path_optimizer.solve(method)
        path_optimizer.print_stats()
        opt_path, curv_list, vel_list = path_optimizer.get_results(method=method)
        corrected_vel_list = []
        for vel in vel_list:
            # corrected_vel_list.append(max(vel, min_vel))
            corrected_vel_list.append(vel)
                
        opt_path_list.append(opt_path)
        curv_lists.append(curv_list)
        vel_lists.append(corrected_vel_list)
        

    # save image
    img = lat_line_creator.map_img.copy()
    legend_text = f"method: {method} optimization"
    text_scale = lat_line_creator.map_width / 4000
    text_scale = max(text_scale, 0.3)
    cv2.putText(
        img,
        legend_text,
        (int(30 * text_scale), int(40 * text_scale)),
        cv2.QT_FONT_NORMAL,
        1 * text_scale,
        (0, 0, 0),
    )
    for idx, (opt_path, params) in enumerate(zip(opt_path_list, params_list)):
        bgr = bgr_list[idx % len(bgr_list)]
        img = lat_line_creator.draw_polyline(img, opt_path, bgr, line_thickness)
        img = lat_line_creator.draw_points(img, opt_path, size=1)
        if method == "lap time":
            legend_text = f"max vel: {params[0]}m/s ,  long acc: {params[1]}m/s^2 ,  long dec: {params[2]}m/s^2 ,  lat acc: {params[3]}m/s^2"
            cv2.putText(
                img,
                legend_text,
                (int(30 * text_scale), int(40 * (idx + 2) * text_scale)),
                cv2.QT_FONT_NORMAL,
                1 * text_scale,
                bgr,
            )
    pt_1 = opt_path_list[0][0]
    pt_2 = opt_path_list[0][1]
    p_pt_1 = lat_line_creator.metric2pixel_point(pt_1)
    p_pt_2 = lat_line_creator.metric2pixel_point(pt_2)
    img = draw_arrow(img, p_pt_1, p_pt_2)
    if path_optimizer.optimization_failed:
        file_path = os.path.join(out_dir, f"debug.png")
    else:
        file_path = os.path.join(out_dir, f"opt_path.png")
    cv2.imwrite(file_path, img)
    print("Saved opt paths as image.")

    # save opt paths and course info as csv
    file_path_list = []
    for opt_path, curv_list, vel_list, params in zip(
        opt_path_list, curv_lists, vel_lists, params_list
    ):
        if path_optimizer.optimization_failed:
            file_path = os.path.join(out_dir, f"debug.csv")
        else:
            if method == "lap time":
                file_path = os.path.join(
                    out_dir,
                    f"opt_path_lap_time_{params[0]}_{params[1]}_{params[2]}_{params[3]}.csv",
                )
            elif method == "curvature":
                file_path = os.path.join(out_dir, f"opt_path_curvature.csv")
        save_course_info(file_path, opt_path, line_list_thinned, curv_list, vel_list, adjusted_list)
        file_path_list.append(file_path)
    print("Saved opt paths and course info as csv files.")

    return file_path_list


def main(map_name: str = None) -> None:
    optimize_path(map_name)


if __name__ == "__main__":
    fire.Fire(main)
