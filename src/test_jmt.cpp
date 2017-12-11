//
// Created by Kan-Hua Lee on 2017/11/09.
//

#include <fstream>
#include <sstream>

#include "plstream.h"
//#include "spline.h"
#include "vehicle_traj.h"
#include "Eigen-3.3/Eigen/Dense"
#include "spdlog/spdlog.h"

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

void plot2d(const vector<double> &x, const vector<double> &y, const char *filename) {
  plstream *pls;
  pls = new plstream();

  const int NSIZE = x.size();

  PLFLT plt_x[NSIZE], plt_y[NSIZE];

  for (int i = 0; i < NSIZE; i++) {
    plt_x[i] = x[i];
    plt_y[i] = y[i];
  }


  // Parse and process command line arguments
  //pls->parseopts( &argc, argv, PL_PARSE_FULL );

  pls->sfnam(filename);       // file name
  pls->sdev("png");

  // Initialize plplot
  pls->init();

  // Create a labelled box to hold the plot.

  auto xmin_it = min_element(x.begin(), x.end());
  auto xmax_it = max_element(x.begin(), x.end());

  auto ymin_it = min_element(y.begin(), y.end());
  auto ymax_it = max_element(y.begin(), y.end());

  pls->col0(3);
  pls->env(*xmin_it, *xmax_it, *ymin_it, *ymax_it, 0, 0);
  pls->lab("x", "y=100 x#u2#d", "Simple PLplot demo of a 2D line plot");

  pls->col0(2);
  // Plot the data that was prepared above.
  //pls->line( x.size(), plt_x, plt_y );
  pls->poin(x.size(), plt_x, plt_y, 9);


  // In C++ we don't call plend() to close PLplot library
  // this is handled by the destructor
  delete pls;
}

void test_s_only() {
  vector<double> a = JMT({0, 10, 0}, {10, 20, 1}, 1);

  VectorXd t = VectorXd::LinSpaced(100, 0, 1);

  vector<double> new_t(t.size());
  vector<double> s_traj(t.size());

  VectorXd::Map(&new_t[0], t.size()) = t;

  s_traj = fill_poly_traj(a, new_t);
  plot2d(new_t, s_traj, "jmt_s.png");

}

void test_s_d() {
  vector<double> s_a = JMT({0, 10, 0}, {10, 20, 0}, 1);
  vector<double> d_a = JMT({0, 0, 0}, {5, 0, 0}, 1);

  VectorXd t = VectorXd::LinSpaced(100, 0, 1);

  vector<double> new_t(t.size());
  vector<double> s_traj(t.size());
  vector<double> d_traj(t.size());

  VectorXd::Map(&new_t[0], t.size()) = t;

  s_traj = fill_poly_traj(s_a, new_t);
  d_traj = fill_poly_traj(d_a, new_t);
  plot2d(new_t, s_traj, "jmt_s_1.png");
  plot2d(new_t, d_traj, "jmt_d_1.png");
  plot2d(s_traj, d_traj, "jmt_sd.png");

}

void load_map(vector<double> &map_waypoints_x,
              vector<double> &map_waypoints_y,
              vector<double> &map_waypoints_s,
              vector<double> &map_waypoints_dx,
              vector<double> &map_waypoints_dy);
void switch_lane() {
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  load_map(map_waypoints_x, map_waypoints_y, map_waypoints_s, map_waypoints_dx, map_waypoints_dy);

  vector<double> s_a = JMT({30, 10, 0}, {60, 20, 0}, 1);
  vector<double> d_a = JMT({0, 0, 0}, {5, 0, 0}, 1);

  VectorXd t = VectorXd::LinSpaced(100, 0, 1);

  vector<double> new_t(t.size());
  vector<double> s_traj(t.size());
  vector<double> d_traj(t.size());

  VectorXd::Map(&new_t[0], t.size()) = t;

  s_traj = fill_poly_traj(s_a, new_t);
  d_traj = fill_poly_traj(d_a, new_t);

  vector<double> x_traj(s_traj.size()), y_traj(s_traj.size());
  for (int i = 0; i < s_traj.size(); i++) {
    vector<double> xy = getXY(s_traj[i], d_traj[i],
                              map_waypoints_s, map_waypoints_x, map_waypoints_y);
    x_traj[i] = (xy[0]);
    y_traj[i] = (xy[1]);
  }

  plot2d(new_t, x_traj, "jmt_x_1.png");
  plot2d(new_t, y_traj, "jmt_y_1.png");
  plot2d(x_traj, y_traj, "jmt_xy.png");

}

void fit_waypoint_spline2(const vector<double> &map_waypoints_x, const vector<double> &map_waypoints_y,
                          const vector<double> &map_waypoints_s, const vector<double> &map_waypoints_dx,
                          const vector<double> &map_waypoints_dy, int start_index, int end_index,
                          tk::spline &s_x, tk::spline &s_y, tk::spline &s_dx, tk::spline &s_dy) {


  //TODO this implementation does not deal with end_index> max_index
  vector<double> x(map_waypoints_x.begin() + start_index, map_waypoints_x.begin() + end_index);
  vector<double> y(map_waypoints_y.begin() + start_index, map_waypoints_y.begin() + end_index);
  vector<double> s(map_waypoints_s.begin() + start_index, map_waypoints_s.begin() + end_index);
  vector<double> dx(map_waypoints_dx.begin() + start_index, map_waypoints_dx.begin() + end_index);
  vector<double> dy(map_waypoints_dy.begin() + start_index, map_waypoints_dy.begin() + end_index);

  s_x.set_points(s, x);
  s_y.set_points(s, y);
  s_dx.set_points(s, dx);
  s_dy.set_points(s, dy);

}

void start_up() {
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  load_map(map_waypoints_x, map_waypoints_y, map_waypoints_s, map_waypoints_dx, map_waypoints_dy);

  vector<double> way_s;
  vector<double> way_d;

  way_s.insert(way_s.begin(), map_waypoints_s.begin(), map_waypoints_s.begin() + 50);

  tk::spline sx;
  tk::spline sy;
  tk::spline sdx;
  tk::spline sdy;

  fill_poly_traj({0, 1}, {0, 1});

  fit_waypoint_spline2(map_waypoints_x,
                       map_waypoints_y,
                       map_waypoints_s,
                       map_waypoints_dx,
                       map_waypoints_dy,
                       0,
                       5,
                       sx,
                       sy,
                       sdx,
                       sdy);


  for (int i = 0; i < way_s.size(); i++) way_d.push_back(6);

  vector<double> test_s = arange(way_s[0], way_s[4], 0.1);

  vector<double> next_x_vals;
  vector<double> next_y_vals;
  for (int i = 0; i < test_s.size(); i++) {
    double &s = test_s[i];
    next_x_vals.push_back(sx(s) + 6 * sdx(s));
    next_y_vals.push_back(sy(s) + 6 * sdy(s));
  }




  int points_to_generate = 200;

  vector<double> time(next_x_vals.size());

  for (int i = 0; i < next_x_vals.size(); i++) {
    time[i] = i;
  }

  plot2d(time, next_x_vals, "jmt_st.png");

  //print_map(next_s_vals,next_d_vals,50);
  plot2d(next_x_vals, next_y_vals, "jmt_xy.png");
  vector<double> map_x(map_waypoints_x.begin(), map_waypoints_x.begin() + 10);
  vector<double> map_y(map_waypoints_y.begin(), map_waypoints_y.begin() + 10);
  vector<double> map_s(map_waypoints_s.begin(), map_waypoints_s.begin() + 10);
  vector<double> map_dx(map_waypoints_dx.begin(), map_waypoints_dx.begin() + 10);
  plot2d(map_x, map_y,
         "jmt_origxy.png");
  plot2d(map_s, map_dx, "s_dx_fig.png");

}

void load_map(vector<double> &map_waypoints_x,
              vector<double> &map_waypoints_y,
              vector<double> &map_waypoints_s,
              vector<double> &map_waypoints_dx,
              vector<double> &map_waypoints_dy) {// Load up map values for waypoint's x,y,s and d normalized normal vectorsvector<double> map_waypoints_dx;

// Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
// The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ios_base::in);

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
}

int main() {
  //start_up();

  auto console = spdlog::stdout_color_mt("console");
  auto sd_logger = spdlog::basic_logger_mt("sd_logger", "sd_log.txt");
  sd_logger->info("start of the message");

  vector<double> a = {0, 1, 2};
  vector<double> b = {0, 1, 2};

  //spdlog::set_pattern("%v");
  console->set_pattern("%v");

  log_map(sd_logger, a, b, -1);


  return 0;

}