//
// Created by Kan-Hua Lee on 2017/12/11.
//
#include <iostream>
#include "spdlog/spdlog.h"
#include "util.h"

double pi() { return 3.14159265359; }

double deg2rad(double x) { return x * pi() / 180; }

double rad2deg(double x) { return x * 180 / pi(); }

double mph2mps(double x) { return x * 0.44704; }

double mps2mph(double x) { return x * 2.23694; }

double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

double getmax(double x1, double x2) {
  if (x1 > x2) {
    return x1;
  } else {
    return x2;
  }
}

double getmin(double x1, double x2) {
  if (x1 < x2) {
    return x1;
  } else return x2;
}

std::vector<double> arange(double lower_bound, double higher_bound, double delta_t) {
  std::vector<double> t;
  t.push_back(lower_bound);
  double current_t = lower_bound + delta_t;
  while (current_t < higher_bound) {
    t.push_back(current_t);
    current_t += delta_t;
  }
  t.push_back(higher_bound);
  return t;

}

std::vector<double> arange_pt(double lower_bound, int num_points, double interval) {
  std::vector<double> t;

  for (int i = 0; i <= num_points; i++) {

    t.push_back(lower_bound + i * interval);
  }
  return t;

}

void print_map(const std::vector<double> &map_x, const std::vector<double> &map_y, int number) {
  unsigned int number_to_print;
  if (number == -1) number_to_print = map_x.size();
  else {
    number_to_print = number;
  }

  std::cout << "map:value" << std::endl;
  for (int i = 0; i < number_to_print; i++) {
    std::cout << map_x[i] << "," << map_y[i] << std::endl;
  }

}

void log_map(std::shared_ptr<spdlog::logger> spdlogger,
             const std::vector<double> &map_x,
             const std::vector<double> &map_y,
             int number) {
  unsigned int number_to_print;
  if (number == -1) number_to_print = map_x.size();
  else {
    number_to_print = number;
  }

  for (int i = 0; i < number_to_print; i++) {
    spdlogger->info("{:03.2f},{:03.2f}", map_x[i], map_y[i]);
  }

}