# Vhicle trajectory generator 

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

