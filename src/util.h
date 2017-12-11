//
// Created by Kan-Hua Lee on 2017/10/21.
//

#ifndef PATH_PLANNING_UTIL_H
#define PATH_PLANNING_UTIL_H

#endif //PATH_PLANNING_UTIL_H

#include <vector>
#include <cmath>
#include "spdlog/spdlog.h"

double pi();
double deg2rad(double x);

double rad2deg(double x);

double mph2mps(double x);

double mps2mph(double x);

double distance(double x1, double y1, double x2, double y2);

double getmax(double x1, double x2);

double getmin(double x1, double x2);

std::vector<double> arange(double lower_bound, double higher_bound, double delta_t);

std::vector<double> arange_pt(double lower_bound, int num_points, double interval);

void print_map(const std::vector<double> &map_x, const std::vector<double> &map_y, int number);

void log_map(std::shared_ptr<spdlog::logger> spdlogger, const std::vector<double> &map_x,
             const std::vector<double> &map_y, int number);



