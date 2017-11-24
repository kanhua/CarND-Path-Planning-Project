//
// Created by Kan-Hua Lee on 2017/11/04.
//

#ifndef PATH_PLANNING_VEHICLE_TRAJ_H
#define PATH_PLANNING_VEHICLE_TRAJ_H

#endif //PATH_PLANNING_VEHICLE_TRAJ_H
#include <vector>
#include <iostream>
#include "Eigen-3.3/Eigen/Dense"
#include "spline.h"


using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

constexpr double pi();

double deg2rad(double x);

double rad2deg(double x);


vector<double> map_to_car_coords(double global_map_x, double global_map_y,
                                 double global_car_x, double global_car_y, double car_yaw);


vector<double> car_to_map_coords(double local_map_x, double local_map_y,
                                 double global_car_x, double global_car_y, double car_yaw);



vector<double> JMT(vector< double> start, vector <double> end, double T);


void fill_spline_2(vector<double> &map_x, vector<double> &map_y, vector<double> &traj_x, vector<double> &traj_y,
                   int points_to_generate, double desired_speed, double acceleration, double car_speed);


void print_map(const vector<double> &map_x, const vector<double> &map_y, int number);

void gen_traj_from_spline(double car_x, double car_y, double car_s, double car_d, double car_speed, double car_yaw,
                          const vector<double> &previous_path_x, const vector<double> &previous_path_y,
                          const vector<vector<double>> &sensor_fusion, const vector<double> &map_waypoints_x,
                          const vector<double> &map_waypoints_y, const vector<double> &map_waypoints_dx,
                          const vector<double> &map_waypoints_dy, vector<double> &next_x_vals,
                          vector<double> &next_y_vals);


void gen_traj_from_jmt(double car_x, double car_y, double car_s, double car_d, double car_speed, double car_yaw,
                       const vector<double> &previous_path_x, const vector<double> &previous_path_y,
                       const vector<vector<double>> &sensor_fusion, const vector<double> &map_waypoints_x,
                       const vector<double> &map_waypoints_y, const vector<double> &map_waypoints_dx,
                       const vector<double> &map_waypoints_dy, vector<double> &next_x_vals, vector<double> &next_y_vals,
                       vector<double> &map_waypoints_s);

void fill_jmt(vector<double> &map_x, vector<double> &map_y, vector<double> &traj_x, vector<double> &traj_y,
              int points_to_generate, double desired_speed, double acceleration, double car_speed);

vector<double> getXY(double s, double d, const vector<double> &maps_s,
                     const vector<double> &maps_x, const vector<double> &maps_y);