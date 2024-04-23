# Vhicle trajectory generator 

# Requirements
* cpp
  * cppad
  * [ifopt](https://github.com/ethz-adrl/ifopt.git)
* python
  * casadi
  * numpy
  * opencv-python
  * opencv-contrib-python
  * fire
  * tabulate

# How to use
1. Copy your `<map_name>.png` and `<map_name>.yaml` to `map` directory
2. Plan the smooth path with a command `python3 path_optimizer.py <map_name>` (smooth path is a path which is minimized the curvature change). The output file is `data/<map_name>/opt_path_curvature.csv`
3. Calculate the normal vector of smooth path and bounds of the course. The output file is `build/line_modified.csv`
```shell
mkdir -p build
cd build 
cmake ..
make
./line_modifier
```
4. Plan optimal path considering Dynamic Bycicle model by `python3 dbm_opt.py`
The output file is `ego_ref_waypoint.csv`.

# Theory

# Assumption
* The vehicle runs the path under the Dynamic Bicycle Model.
* The position of waypoints is on the line segment which is parallel to the normal vector of the reference path (center line or pre-planned path).
* The vehicle runs the course counter-clock-wise. If the direction is clock-wise, please edit
```python
# yaw
if i+1 == horizon:
    optimizer.subject_to(self.yaw_list[(i+1)%horizon] + 2*np.pi== 
                        self.yaw_list[i] + 
                        self.omega_list[i]*self.dt_list[i])
```
to 
```python
# yaw
if i+1 == horizon:
    optimizer.subject_to(self.yaw_list[(i+1)%horizon] - 2*np.pi== 
                        self.yaw_list[i] + 
                        self.omega_list[i]*self.dt_list[i])
```

# Optmization Problem

## Smooth path planning 
This problem plans the reference path for the Lap time optimization. 
* Exploring variable
  * waypoint positions
  * curvature of the planned path
* Evaluation
  * minimize the difference between curvatures at adjacent waypoints 
* Constraints
  * turning radius
  * distance from the course outlines. 

## Lap time Optimization
* Exploring variable
  * waypoint positoins
  * velocity (longitudinal and lateral)
  * yaw angle, yaw rate (named omega in code)
  * longitudinal acceleration
  * steer angle
  * time to arrive the next waypoint
* Cost
  * minimize lap time
  * smoothness of acceleration
* Constraints
  * State equation under Dynamic Bycicle Model
  * Lateral acceration limit
  * turning radius
  * Slip angle of the front/rear tires. 