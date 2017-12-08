//
// Created by Kan-Hua Lee on 2017/12/01.
//

//
// Created by Kan-Hua Lee on 2017/11/09.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "plstream.h"
#include "vehicle_traj.h"
using namespace std;

int main() {
  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }

  vector<vector<int>> state_map = {{0, 1, 4},
                                   {0, 1, 2, 4},
                                   {1, 2, 4}};

  cout << state_map[0][2] << endl;

  double test_angle = -pi() / 2;
  while (test_angle < pi() / 2) {
    vector<double> nc = getFrenet(963, 1131.95, test_angle += 0.1, map_waypoints_x, map_waypoints_y);
    cout << "s:" << nc[0] << endl;
    cout << "d:" << nc[1] << endl;
  }

  return 0;
}