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

double mph2mps(double x);

double mps2mph(double x);

struct map_data {
    vector<double> map_waypoints_x;
    vector<double> map_waypoints_y;
    vector<double> map_waypoints_s;
    vector<double> map_waypoints_dx;
    vector<double> map_waypoints_dy;

};

struct car_state {
    double car_x = 0;
    double car_y = 0;
    double car_s = 0;
    double car_d = 0;
    double car_yaw = 0;
    double car_speed = 0;
};


vector<double> map_to_car_coords(double global_map_x, double global_map_y,
                                 double global_car_x, double global_car_y, double car_yaw);


vector<double> car_to_map_coords(double local_map_x, double local_map_y,
                                 double global_car_x, double global_car_y, double car_yaw);

void
get_lane_cost(double car_s, const vector<vector<double>> &sensor_fusion, const int lane_width, int prev_lane_number,
              vector<double> &lane_cost);


vector<double> JMT(vector< double> start, vector <double> end, double T);


void fill_spline(vector<double> &map_x, vector<double> &map_y, vector<double> &traj_x, vector<double> &traj_y,
                 int points_to_generate, double desired_speed, double acceleration, double car_speed);


void print_map(const vector<double> &map_x, const vector<double> &map_y, int number);

void gen_traj_from_spline(car_state &cstate,
                          vector<double> &previous_path_x,
                          vector<double> &previous_path_y,
                          const vector<vector<double>> &sensor_fusion, const vector<double> &map_waypoints_x,
                          const vector<double> &map_waypoints_y, const vector<double> &map_waypoints_dx,
                          const vector<double> &map_waypoints_dy, vector<double> &next_x_vals,
                          vector<double> &next_y_vals);


void fill_jmt(vector<double> &map_x, vector<double> &map_y, vector<double> &traj_x, vector<double> &traj_y,
              int points_to_generate, double desired_speed, double acceleration, double car_speed);

vector<double> getXY(double s, double d, const vector<double> &maps_s,
                     const vector<double> &maps_x, const vector<double> &maps_y);