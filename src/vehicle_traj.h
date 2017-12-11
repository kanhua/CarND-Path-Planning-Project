//
// Created by Kan-Hua Lee on 2017/11/04.
//

#ifndef PATH_PLANNING_VEHICLE_TRAJ_H
#define PATH_PLANNING_VEHICLE_TRAJ_H

#endif //PATH_PLANNING_VEHICLE_TRAJ_H
#include <vector>
#include <iostream>
#include "spline.h"
#include "util.h"

const double SPEEDLIMIT = 21.5;


struct map_data {
  std::vector<double> map_waypoints_x;
  std::vector<double> map_waypoints_y;
  std::vector<double> map_waypoints_s;
  std::vector<double> map_waypoints_dx;
  std::vector<double> map_waypoints_dy;

};

struct car_state {
    double car_x = 0;
    double car_y = 0;
    double car_s = 0;
    double car_d = 0;
    double car_yaw = 0;
    double car_speed = 0;
};

enum CarStateNum {
  left_lane = 0,
  middle_lane = 1,
  right_lane = 2,
  stopping = 4,
};

//The mapping of the state. The first is the incumbent state.
const std::vector<std::vector<CarStateNum>> state_map = {{left_lane, middle_lane, stopping},
                                                         {middle_lane, left_lane, right_lane, stopping},
                                                         {right_lane, middle_lane, stopping}};

std::vector<double> arange(double lower_bound, double higher_bound, double delta_t);

std::vector<double> global_to_car(double global_map_x, double global_map_y,
                             double global_car_x, double global_car_y, double car_yaw);

std::vector<double> car_to_global(double local_map_x, double local_map_y,
                                  double global_car_x, double global_car_y, double car_yaw);

std::vector<double> JMT(std::vector<double> start, std::vector<double> end, double T);

void fill_spline(const std::vector<double> &map_x,
                 const std::vector<double> &map_y,
                 int points_to_generate,
                 double desired_speed,
                 double acceleration,
                 double car_speed,
                 double read_in_interval,
                 std::vector<double> &traj_x,
                 std::vector<double> &traj_y);


void
gen_next_traj(const car_state &cstate,
              const std::vector<double> &previous_path_x,
              const std::vector<double> &previous_path_y,
              const std::vector<std::vector<double>> &sensor_fusion,
              const std::vector<double> &map_waypoints_x,
              const std::vector<double> &map_waypoints_y,
              const std::vector<double> &map_waypoints_dx,
              const std::vector<double> &map_waypoints_dy,
              const std::vector<double> &map_waypoints_s,
              double end_path_s,
              double end_path_d,
              std::vector<double> &next_x_vals,
              std::vector<double> &next_y_vals);

void fit_waypoint_spline(const std::vector<double> &map_waypoints_x, const std::vector<double> &map_waypoints_y,
                         const std::vector<double> &map_waypoints_s, const std::vector<double> &map_waypoints_dx,
                         const std::vector<double> &map_waypoints_dy, int start_index, int end_index,
                         tk::spline &s_x, tk::spline &s_y, tk::spline &s_dx, tk::spline &s_dy);

void fill_jmt(double start_s,
              double v_i,
              double v_F,
              double start_d,
              double end_d,
              int points_to_generate,
              std::vector<double> &traj_s,
              std::vector<double> &traj_d);

std::vector<double> getXY(double s, double d, const std::vector<double> &maps_s,
                          const std::vector<double> &maps_x, const std::vector<double> &maps_y);


// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
std::vector<double> getFrenet(double x, double y, double theta_in_rad,
                              const std::vector<double> &maps_x, const std::vector<double> &maps_y);

void map_to_car_coords_array(const car_state &cstate, std::vector<double> &next_map_x, std::vector<double> &next_map_y);

void car_to_map_cords_array(const car_state &cstate,
                            std::vector<double> &next_x_vals,
                            std::vector<double> &next_y_vals);

double
eval_cost(double car_x,
          double car_y,
          double car_theta,
          double delta_t,
          const std::vector<std::vector<double>> &sensor_fusion,
          const std::vector<double> &map_x,
          const std::vector<double> &map_y,
          const double lane_width);

std::vector<double> fill_poly_traj(std::vector<double> a_vec, std::vector<double> t_vec);


