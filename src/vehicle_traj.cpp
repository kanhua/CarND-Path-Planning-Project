//
// Created by Kan-Hua Lee on 2017/11/04.
//

#include <fstream>
#include "spline.h"
#include "vehicle_traj.h"
#include "Eigen-3.3/Eigen/Dense"
#include "util.h"
#include "jmt.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

void fit_waypoint_spline(const vector<double> &map_waypoints_x, const vector<double> &map_waypoints_y,
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

int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y) {

  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for (int i = 0; i < maps_x.size(); i++) {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x, y, map_x, map_y);
    if (dist < closestLen) {
      closestLen = dist;
      closestWaypoint = i;
    }

  }

  return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta_in_rad, const vector<double> &maps_x, const vector<double> &maps_y) {

  int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2((map_y - y), (map_x - x));

  double angle = abs(theta_in_rad - heading);

  if (angle > pi() / 4) {
    closestWaypoint++;
  }

  // Use modulus to deal with boundary
  return closestWaypoint % maps_x.size();

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta_in_rad,
                         const vector<double> &maps_x, const vector<double> &maps_y) {
  int next_wp = NextWaypoint(x, y, theta_in_rad, maps_x, maps_y);

  int prev_wp;
  prev_wp = next_wp - 1;
  if (next_wp == 0) {
    prev_wp = maps_x.size() - 1;
  }

  double n_x = maps_x[next_wp] - maps_x[prev_wp];
  double n_y = maps_y[next_wp] - maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
  double proj_x = proj_norm * n_x;
  double proj_y = proj_norm * n_y;

  double frenet_d = distance(x_x, x_y, proj_x, proj_y);

  //see if d value is positive or negative by comparing it to a center point

  double center_x = 1000 - maps_x[prev_wp];
  double center_y = 2000 - maps_y[prev_wp];
  double centerToPos = distance(center_x, center_y, x_x, x_y);
  double centerToRef = distance(center_x, center_y, proj_x, proj_y);

  if (centerToPos <= centerToRef) {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for (int i = 0; i < prev_wp; i++) {
    frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
  }

  frenet_s += distance(0, 0, proj_x, proj_y);

  return {frenet_s, frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s,
                     const vector<double> &maps_x, const vector<double> &maps_y) {
  int prev_wp = -1;

  while (s > maps_s[prev_wp + 1] && (prev_wp < (int) (maps_s.size() - 1))) {
    prev_wp++;
  }

  int wp2 = (prev_wp + 1) % maps_x.size();

  double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s - maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
  double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

  double perp_heading = heading - pi() / 2;

  double x = seg_x + d * cos(perp_heading);
  double y = seg_y + d * sin(perp_heading);

  return {x, y};

}



/// Do the spline fitting from map_x and map_y. This procedure takes into account the car_speed,
/// desired_speed and acceleration to generate the trajectory.
/// \param map_x
/// \param map_y
/// \param points_to_generate
/// \param desired_speed speed (m/s) that the car tries to approach
/// \param acceleration
/// \param car_speed
/// \param read_in_interval
/// \param traj_x array of x that will be filled into.
/// \param traj_y array of y that will be filled into.
void fill_spline(const vector<double> &map_x, const vector<double> &map_y, int points_to_generate, double desired_speed,
                 double acceleration, double car_speed, double read_in_interval, vector<double> &traj_x,
                 vector<double> &traj_y) {


  // time interval between each reading of the simulator
  //double desired_speed=20; // desired speed in m/s


  tk::spline s;
  s.set_points(map_x, map_y);

  vector<double> &new_traj_x = traj_x;
  vector<double> &new_traj_y = traj_y;

  double current_x = map_x[1];

  double dx = map_x[1] - map_x[0];
  double dy = map_y[1] - map_y[0];

  double inst_speed = 0; // instaneous speed

  if (car_speed < 0.1) {
    inst_speed = car_speed + 0.05;
  } else {
    inst_speed = sqrt(dx * dx + dy * dy) / read_in_interval;
  }

  double theta = atan2(dy, dx);

  for (int i = 0; i < points_to_generate; i++) {

    //TODO ugly code here
    if (acceleration > 0) {
      if (inst_speed < desired_speed) {
        inst_speed += acceleration * read_in_interval;
      } else {
        inst_speed = desired_speed;
      }
    } else {
      if (inst_speed > desired_speed) {
        inst_speed += acceleration * read_in_interval;
      } else {
        inst_speed = desired_speed;
      }

    }

    // check that the car does not go backwards
    assert(inst_speed > -0.5);

    current_x += inst_speed * cos(theta) * read_in_interval;
    double y = s(current_x);

    new_traj_x.push_back(current_x);
    new_traj_y.push_back(y);
  }

}



void fill_jmt(double start_s,
              double v_i,
              double v_F,
              double start_d,
              double end_d,
              int points_to_generate,
              std::vector<double> &traj_s,
              std::vector<double> &traj_d) {

  assert(points_to_generate > 0);
  const double read_in_interval = 0.02; // time interval between each reading of the simulator
  const double a_max = 8;

  //cout << "start_sd:" << start_s << "," << start_d << endl;

  double delta_t, end_s, final_v, a;


  // no need to change the lane
  delta_t = points_to_generate * read_in_interval;

  a = (v_F - v_i) / delta_t;
  if (abs(a) > a_max) {
    if (a > 0) a = a_max;
    else a = -a_max;
  }

  // calculate the final boundary conditions

  final_v = v_i + delta_t * a;
  //assert(final_v <23);

  end_s = start_s + v_i * delta_t + 0.5 * a * pow(delta_t, 2);

  cout << "calculate final car speed (v_f):" << final_v << endl;

  vector<double> s_coef, sp_coef, spp_coef;
  vector<double> d_coef, dp_coef, dpp_coef;

  full_coef_JMT({start_s, v_i, a}, {end_s, final_v, a}, delta_t, s_coef, sp_coef, spp_coef);
  full_coef_JMT({start_d, 0, 0}, {end_d, 0, 0}, delta_t, d_coef, dp_coef, dpp_coef);


  vector<double> t = arange_pt(0, points_to_generate, read_in_interval);

  traj_s = fill_poly_traj(s_coef, t);
  traj_d = fill_poly_traj(d_coef, t);

  vector<double> traj_sp = fill_poly_traj(sp_coef, t);
  vector<double> traj_spp = fill_poly_traj(spp_coef, t);

  vector<double> traj_dp = fill_poly_traj(dp_coef, t);
  vector<double> traj_dpp = fill_poly_traj(dpp_coef, t);

  vector<double> ntraj_x(traj_s.begin() + 1, traj_s.begin() + points_to_generate + 1);
  vector<double> ntraj_y(traj_d.begin() + 1, traj_d.begin() + points_to_generate + 1);

  traj_s = ntraj_x;
  traj_d = ntraj_y;


  // write the calculated trajectory to a file for debugging
  static ofstream fout("jmt_traj_log.txt");
  static unsigned int iteration = 0;

  for (int i = 0; i < points_to_generate; i++) {
    fout << iteration << "," << traj_s[i] << "," << traj_sp[i] << "," << traj_spp[i] << endl;
  }

  iteration++;

}

vector<double> global_to_car(double global_map_x, double global_map_y,
                             double global_car_x, double global_car_y, double car_yaw) {
  double local_map_x;
  double local_map_y;

  double dx = global_map_x - global_car_x;
  double dy = global_map_y - global_car_y;

  local_map_x = cos(car_yaw) * dx + sin(car_yaw) * dy;
  local_map_y = -sin(car_yaw) * dx + cos(car_yaw) * dy;

  return {local_map_x, local_map_y};
}

vector<double> car_to_global(double local_map_x, double local_map_y,
                             double global_car_x, double global_car_y, double car_yaw) {

  double global_map_x;
  double global_map_y;

  global_map_x = cos(car_yaw) * local_map_x - sin(car_yaw) * local_map_y + global_car_x;
  global_map_y = sin(car_yaw) * local_map_x + cos(car_yaw) * local_map_y + global_car_y;

  return {global_map_x, global_map_y};

}

double traj_end_eval_state(double delta_t, const vector<double> &next_x_val, const vector<double> &next_y_val,
                           const vector<vector<double>> &sensor_fusion, const vector<double> &map_x,
                           const vector<double> &map_y,
                           const car_state &cstate) {

  const double lane_width = 4;

  //int path_index = floor(delta_t / simulator_interval);
  int path_index = next_x_val.size() - 1;

  assert(path_index > 0 && path_index < next_x_val.size());
  double car_next_x = next_x_val[path_index];
  double car_next_y = next_y_val[path_index];

  double car_theta = atan2(next_y_val[path_index] - next_y_val[path_index - 1],
                           next_x_val[path_index] - next_x_val[path_index - 1]);

  return eval_cost(car_next_x, car_next_y, car_theta, delta_t, sensor_fusion, map_x, map_y, lane_width);

}

/// Predicit the state of the car after delta_t with given acceleration
/// \param current_car_state
/// \param delta_t
/// \param acc
/// \param next_car_lane
/// \param map_x
/// \param map_y
/// \param map_s
/// \return
car_state guess_next_car_state(const car_state &current_car_state,
                               const double delta_t,
                               const double acc,
                               int next_car_lane,
                               const vector<double> &map_x,
                               const vector<double> &map_y,
                               const vector<double> &map_s) {

  //check input
  assert (next_car_lane >= -1 && next_car_lane < 3);
  double lane_width = 4.0;

  const double &s = current_car_state.car_s;
  const double &a = acc;
  const double &v0 = current_car_state.car_speed;

  int current_lane = floor(current_car_state.car_d / 4.0);//TODO fix magic number here
  if (next_car_lane == -1) {
    next_car_lane = current_lane;
  }

  double final_s = s + v0 * delta_t + 0.5 * acc * pow(delta_t, 2);

  car_state next_car_state;

  next_car_state.car_s = final_s;

  next_car_state.car_d = 0.5 * lane_width + lane_width * next_car_lane;

  vector<double> nc = getXY(next_car_state.car_s, next_car_state.car_d,
                            map_s, map_x, map_y);
  next_car_state.car_x = nc[0];
  next_car_state.car_y = nc[1];

  next_car_state.car_yaw = atan2(nc[1] - current_car_state.car_y, nc[0] - current_car_state.car_x);
  next_car_state.car_speed = v0 + acc * delta_t;

  return next_car_state;

}

/// Check if the car can switch lane to target_lane_num safely
/// \param curent_car_state
/// \param target_lane_num The lane number that the car wants to switch to
/// \param sensor_fusion
/// \param lane_width
/// \return
double safe_switch(const car_state &curent_car_state, const int target_lane_num,
                   const vector<vector<double>> &sensor_fusion,
                   const double lane_width) {

  int min_back_car_dist = 100;
  for (int i = 0; i < sensor_fusion.size(); i++) {
    double neighbor_car_d = sensor_fusion[i][6];
    double neighbor_car_s = sensor_fusion[i][5];

    double neighbor_car_x = sensor_fusion[i][1];
    double neighbor_car_y = sensor_fusion[i][2];

    double vx = sensor_fusion[i][3];
    double vy = sensor_fusion[i][4];
    double neighbor_car_theta = atan2(vy, vx);
    double vs = sqrt(vy * vy + vx * vx);

    int neighbor_car_lane = floor(neighbor_car_d / lane_width);

    if (target_lane_num == neighbor_car_lane) {
      double current_car_dist = abs(curent_car_state.car_s - neighbor_car_s);
      if ((current_car_dist < min_back_car_dist)) {
        min_back_car_dist = current_car_dist;
      }
    }

  }
  if (min_back_car_dist < 10) {
    cout << "Warning: dangerous to switch" << endl;
    return 1000;
  } else {
    return 0; //add 0.1 to avoid 1/0
  }

}

/// Calculate the time required for the car to hit the car in the front,
/// assuming that both the car and the car in the front drive in constant speed.
/// \param curent_car_state The current state of the car
/// \param sensor_fusion
/// \param map_x
/// \param map_y
/// \param lane_width Width of the lane
/// \return the required time in seconds
double eval_next_collision(car_state curent_car_state, const vector<vector<double>> &sensor_fusion,
                           const vector<double> &map_x, const vector<double> &map_y, const double lane_width) {

  double cost = 0;

  // find the nearest front car
  int front_car_index = -1;
  int car_lane = floor(curent_car_state.car_d / lane_width);

  double collision_time = 1000;

  for (int i = 0; i < sensor_fusion.size(); i++) {
    double neighbor_car_d = sensor_fusion[i][6];
    double neighbor_car_s = sensor_fusion[i][5];

    double neighbor_car_x = sensor_fusion[i][1];
    double neighbor_car_y = sensor_fusion[i][2];

    double vx = sensor_fusion[i][3];
    double vy = sensor_fusion[i][4];
    double neighbor_car_theta = atan2(vy, vx);
    double vs = sqrt(vy * vy + vx * vx);

    int neighbor_car_lane = floor(neighbor_car_d / lane_width);

    if (car_lane == neighbor_car_lane) {
      double current_car_dist = neighbor_car_s - curent_car_state.car_s;
      if (current_car_dist > 0) {
        double ct = current_car_dist / (curent_car_state.car_speed);
        if (ct < collision_time) collision_time = ct;
      }
    }

  }
  return collision_time;
}

double
eval_cost(double car_x, double car_y, double car_theta, double delta_t, const vector<vector<double>> &sensor_fusion,
          const vector<double> &map_x, const vector<double> &map_y, const double lane_width) {
  vector<double> car_nc = getFrenet(car_x, car_y, car_theta, map_x, map_y);
  double car_next_s = car_nc[0];
  double car_next_d = car_nc[1];

  int car_lane = floor(car_next_d / lane_width);
  assert (car_lane >= 0 && car_lane < 4);

  double cost = 0;
  for (int i = 0; i < sensor_fusion.size(); i++) {
    double neighbor_car_d = sensor_fusion[i][6];
    double neighbor_car_s = sensor_fusion[i][5];

    double neighbor_car_x = sensor_fusion[i][1];
    double neighbor_car_y = sensor_fusion[i][2];

    double vx = sensor_fusion[i][3];
    double vy = sensor_fusion[i][4];
    double neighbor_car_theta = atan2(vy, vx);

    vector<double> nc = getFrenet(neighbor_car_x + delta_t * vx, neighbor_car_y + delta_t * vy,
                                  neighbor_car_theta, map_x, map_y);

    double neighbor_car_next_s = nc[0];
    double neighbor_car_next_d = nc[1];

    int neighbor_car_lane = floor(neighbor_car_next_d / lane_width);

    if (car_lane == neighbor_car_lane) {
      double future_car_dist = neighbor_car_next_s - car_next_s;
      double current_car_dist = neighbor_car_s - car_next_s;

      double current_cost = 1 / current_car_dist;
      if (current_cost > cost) cost = current_cost;

    }

  }

  return cost;
}

void initialize_reference_points(const car_state &cstate,
                                 const vector<double> &previous_path_x,
                                 const vector<double> &previous_path_y,
                                 double &before_next_path_start_x,
                                 double &before_next_path_start_y,
                                 double &next_path_start_x,
                                 double &next_path_start_y,
                                 double &ref_x,
                                 double &ref_y,
                                 double &ref_yaw,
                                 double &ref_speed);

void gen_next_map_waypoints_by_wpindex(const vector<double> &map_waypoints_x,
                                       const vector<double> &map_waypoints_y,
                                       const vector<double> &map_waypoints_dx,
                                       const vector<double> &map_waypoints_dy,
                                       int prev_lane_number,
                                       int next_lane_number,
                                       int num_next_index,
                                       const int lane_width,
                                       double ref_x,
                                       double ref_y,
                                       double ref_yaw,
                                       vector<double> &next_map_waypoints_x,
                                       vector<double> &next_map_waypoints_y);

void gen_traj_from_jmt(const car_state &cstate,
                       const int state_number,
                       const vector<double> &previous_path_x,
                       const vector<double> &previous_path_y,
                       const vector<vector<double>> &sensor_fusion,
                       const vector<double> &map_waypoints_x,
                       const vector<double> &map_waypoints_y,
                       const vector<double> &map_waypoints_s,
                       const vector<double> &map_waypoints_dx,
                       const vector<double> &map_waypoints_dy,
                       vector<double> &next_x_vals,
                       vector<double> &next_y_vals,
                       int total_future_points,
                       double read_in_interval) {

  assert(state_number >= 0 && state_number <= 4);
  assert(next_x_vals.size() == 0);
  assert(next_y_vals.size() == 0);

  static double s_end_path_s = 0;
  static double s_end_path_d = 0;

  double acceleration = 6;

  double desired_speed = 20;
  int next_lane_number = -1;

  const double lane_width = 4.0; // the width of the lane

  // Tell which lane that the car currently stays
  // The lane number that the car stays on. The lane next to the center line is zero.
  int prev_lane_number = floor(cstate.car_d / lane_width);

  if (state_number >= 0 && state_number < 3) {
    next_lane_number = state_number;
  } else {
    next_lane_number = prev_lane_number;
  }

  assert(next_lane_number >= 0 && next_lane_number < 3);

  int prev_points = previous_path_x.size();
  assert(previous_path_x.size() == previous_path_y.size());

  cout << "prev points left:" << previous_path_x.size() << endl;
  //cout << "start xy value:" << next_path_start_x << "," << next_path_start_y << endl;

  double px_0, px_n1, py_0, py_n1, ref_x, ref_y, ref_yaw, ref_speed;
  initialize_reference_points(cstate,
                              previous_path_x,
                              previous_path_y,
                              px_n1,
                              py_n1,
                              px_0,
                              py_0,
                              ref_x,
                              ref_y,
                              ref_yaw,
                              ref_speed);

  cout << "ref point speed:" << ref_speed << endl;

  double next_map_waypoints_s;
  double next_map_waypoints_d;

  if (previous_path_x.size() > 1) {

    next_map_waypoints_s = s_end_path_s;
    next_map_waypoints_d = s_end_path_d;
  } else {
    next_map_waypoints_s = cstate.car_s;
    next_map_waypoints_d = cstate.car_d;
  }


  // spline fitting:
  tk::spline sx, sy, sdx, sdy;
  sx.set_points(map_waypoints_s, map_waypoints_x);
  sy.set_points(map_waypoints_s, map_waypoints_y);
  sdx.set_points(map_waypoints_s, map_waypoints_dx);
  sdy.set_points(map_waypoints_s, map_waypoints_dy);


  // Find next points
  int points_to_generate = total_future_points - prev_points;

  double target_speed = 0;
  if (state_number == stopping) {
    //deacceleration
    if (ref_speed < SPEEDLIMIT) target_speed = ref_speed * 0.5;
    else {
      target_speed = 0;
    }

  } else {

    target_speed = desired_speed;
  }

  if ((prev_lane_number != next_lane_number) && (points_to_generate > 0)) {
    points_to_generate = 100;
  }

  if (points_to_generate < 0) points_to_generate = 2;

  fill_jmt(next_map_waypoints_s,
           ref_speed,
           target_speed,
           next_map_waypoints_d,
           (next_lane_number + 0.5) * lane_width,
           points_to_generate,
           next_x_vals,
           next_y_vals);
  s_end_path_s = next_x_vals[next_x_vals.size() - 1];
  s_end_path_d = next_y_vals[next_y_vals.size() - 1];
  cout << "last sd value:" << next_x_vals[next_x_vals.size() - 1] << "," << next_y_vals[next_y_vals.size() - 1]
       << endl;
  for (int i = 0; i < next_x_vals.size(); i++) {
    double s = next_x_vals[i];
    double d = next_y_vals[i];

    next_x_vals[i] = sx(s) + d * sdx(s);
    next_y_vals[i] = sy(s) + d * sdy(s);
  }
  if (prev_points > 2) {
    next_x_vals.insert(next_x_vals.begin(), previous_path_x.begin(), previous_path_x.end());
    next_y_vals.insert(next_y_vals.begin(), previous_path_y.begin(), previous_path_y.end());
  }

}

void gen_next_map_waypoints_by_interp(const vector<double> &map_waypoints_x,
                                      const vector<double> &map_waypoints_y,
                                      const vector<double> &map_waypoints_s,
                                      double ref_x,
                                      double ref_y,
                                      double ref_yaw,
                                      int next_lane_number,
                                      const int lane_width,
                                      vector<double> &next_map_waypoints_x,
                                      vector<double> &next_map_waypoints_y);
void gen_traj_from_spline(const car_state &cstate,
                          const int state_number,
                          const vector<double> &previous_path_x,
                          const vector<double> &previous_path_y,
                          const vector<vector<double>> &sensor_fusion,
                          const vector<double> &map_waypoints_x,
                          const vector<double> &map_waypoints_y,
                          const vector<double> &map_waypoints_s,
                          const vector<double> &map_waypoints_dx,
                          const vector<double> &map_waypoints_dy,
                          double endpoint_s,
                          double endpoint_d,
                          vector<double> &next_x_vals,
                          vector<double> &next_y_vals,
                          int total_future_points,
                          double read_in_interval) {

  assert(state_number >= 0 && state_number <= 4);
  assert(next_x_vals.size() == 0);
  assert(next_y_vals.size() == 0);

  double acceleration = 6;

  double before_next_path_start_x = 0;
  double before_next_path_start_y = 0;

  double next_path_start_x = 0;
  double next_path_start_y = 0;

  double ref_x = 0;
  double ref_y = 0;
  double ref_yaw = 0;
  double ref_speed = 0;

  double desired_speed = 20;
  int next_lane_number = -1;


  //Get the map coordinates of the next few points
  //Assuming staying on the second lane at the moment

  int num_next_index = 6;

  const int lane_width = 4; // the width of the lane

  // Tell which lane that the car currently stays
  // The lane number that the car stays on. The lane next to the center line is zero.
  int prev_lane_number = floor(cstate.car_d / lane_width);

  if (state_number >= 0 && state_number < 3) {
    next_lane_number = state_number;
  } else {
    next_lane_number = prev_lane_number;
  }

  assert(next_lane_number >= 0 && next_lane_number < 3);

  int prev_points = previous_path_x.size();
  assert(previous_path_x.size() == previous_path_y.size());
  initialize_reference_points(cstate,
                              previous_path_x,
                              previous_path_y,
                              before_next_path_start_x,
                              before_next_path_start_y,
                              next_path_start_x,
                              next_path_start_y,
                              ref_x,
                              ref_y,
                              ref_yaw,
                              ref_speed);




  vector<double> next_map_waypoints_x;
  vector<double> next_map_waypoints_y;

  next_map_waypoints_x.push_back(before_next_path_start_x);
  next_map_waypoints_x.push_back(next_path_start_x);

  next_map_waypoints_y.push_back(before_next_path_start_y);
  next_map_waypoints_y.push_back(next_path_start_y);

  gen_next_map_waypoints_by_interp(map_waypoints_x,
                                   map_waypoints_y,
                                   map_waypoints_s,
                                   ref_x,
                                   ref_y,
                                   ref_yaw,
                                   next_lane_number,
                                   lane_width,
                                   next_map_waypoints_x,
                                   next_map_waypoints_y);



  /*
  gen_next_map_waypoints_by_wpindex(map_waypoints_x,
                                    map_waypoints_y,
                                    map_waypoints_dx,
                                    map_waypoints_dy,
                                    prev_lane_number,
                                    next_lane_number,
                                    num_next_index,
                                    lane_width,
                                    ref_x,
                                    ref_y,
                                    ref_yaw,
                                    next_map_waypoints_x,
                                    next_map_waypoints_y);
                                    */

  map_to_car_coords_array(cstate, next_map_waypoints_x, next_map_waypoints_y);


  // Find next points
  int points_to_generate = total_future_points - prev_points;

  if (state_number == stopping) {
    //deacceleration
    fill_spline(next_map_waypoints_x,
                next_map_waypoints_y,
                points_to_generate,
                (cstate.car_speed + 1) / 2,
                -acceleration,
                cstate.car_speed,
                read_in_interval,
                next_x_vals,
                next_y_vals);
  } else {

    fill_spline(next_map_waypoints_x,
                next_map_waypoints_y,
                points_to_generate,
                desired_speed,
                acceleration,
                cstate.car_speed,
                read_in_interval,
                next_x_vals,
                next_y_vals);
  }


  // Convert the coordinates back to the map coordinates

  car_to_map_cords_array(cstate, next_x_vals, next_y_vals);

  if (prev_points > 2) {
    next_x_vals.insert(next_x_vals.begin(), previous_path_x.begin(), previous_path_x.end());
    next_y_vals.insert(next_y_vals.begin(), previous_path_y.begin(), previous_path_y.end());
  }

  vector<double>().swap(next_map_waypoints_x);
  vector<double>().swap(next_map_waypoints_y);

}
void gen_next_map_waypoints_by_interp(const vector<double> &map_waypoints_x,
                                      const vector<double> &map_waypoints_y,
                                      const vector<double> &map_waypoints_s,
                                      double ref_x,
                                      double ref_y,
                                      double ref_yaw,
                                      int next_lane_number,
                                      const int lane_width,
                                      vector<double> &next_map_waypoints_x,
                                      vector<double> &next_map_waypoints_y) {
  vector<double> sd = getFrenet(ref_x, ref_y, ref_yaw, map_waypoints_x, map_waypoints_y);

  for (int i = 0; i < 3; i++) {
    vector<double> xy = getXY(sd[0] + 30 * (i + 1), (next_lane_number + 0.5) * lane_width,
                              map_waypoints_s, map_waypoints_x, map_waypoints_y);

    next_map_waypoints_x.push_back(xy[0]);
    next_map_waypoints_y.push_back(xy[1]);

  }
}
void gen_next_map_waypoints_by_wpindex(const vector<double> &map_waypoints_x,
                                       const vector<double> &map_waypoints_y,
                                       const vector<double> &map_waypoints_dx,
                                       const vector<double> &map_waypoints_dy,
                                       int prev_lane_number,
                                       int next_lane_number,
                                       int num_next_index,
                                       const int lane_width,
                                       double ref_x,
                                       double ref_y,
                                       double ref_yaw,
                                       vector<double> &next_map_waypoints_x,
                                       vector<double> &next_map_waypoints_y) {// Construct the array for spline fitting

  //Find the closest index

  int closest_index = NextWaypoint(ref_x, ref_y,
                                   ref_yaw, map_waypoints_x, map_waypoints_y);

  // Use the next index to start if changing lane. This is to avoid the instability when switching lanes
  if (next_lane_number != prev_lane_number) {
    closest_index++;
  }



  for (int i = 0; i < num_next_index; i++) {
    int waypoints_index = (closest_index + i) % map_waypoints_x.size();

    next_map_waypoints_x.push_back(map_waypoints_x[waypoints_index] +
        ((next_lane_number + 0.5) * lane_width) * map_waypoints_dx[waypoints_index]);
    next_map_waypoints_y.push_back(map_waypoints_y[waypoints_index] +
        ((next_lane_number + 0.5) * lane_width) * map_waypoints_dy[waypoints_index]);

  }
}

/// Set the starting point of the next generated path
/// \param cstate
/// \param previous_path_x
/// \param previous_path_y
/// \param before_next_path_start_x
/// \param before_next_path_start_y
/// \param next_path_start_x
/// \param next_path_start_y
/// \param ref_x
/// \param ref_y
/// \param ref_yaw
/// \param ref_speed
void initialize_reference_points(const car_state &cstate,
                                 const vector<double> &previous_path_x,
                                 const vector<double> &previous_path_y,
                                 double &before_next_path_start_x,
                                 double &before_next_path_start_y,
                                 double &next_path_start_x,
                                 double &next_path_start_y,
                                 double &ref_x,
                                 double &ref_y,
                                 double &ref_yaw,
                                 double &ref_speed) {

  int prev_points = previous_path_x.size();
  assert(previous_path_x.size() == previous_path_y.size());

  if (previous_path_x.size() > 2) {
//Use the end points of previous_path_x
    before_next_path_start_x = previous_path_x[prev_points - 2];
    before_next_path_start_y = previous_path_y[prev_points - 2];

    next_path_start_x = previous_path_x[prev_points - 1];
    next_path_start_y = previous_path_y[prev_points - 1];

    ref_x = next_path_start_x;
    ref_y = next_path_start_y;
    ref_yaw = atan2(next_path_start_y - before_next_path_start_y, next_path_start_x - before_next_path_start_x);

    ref_speed = distance(before_next_path_start_x, before_next_path_start_y,
                         next_path_start_x, next_path_start_y) / 0.02;

  } else {

//Use the current coordinates of the car
    next_path_start_x = cstate.car_x;
    next_path_start_y = cstate.car_y;

    before_next_path_start_x = cstate.car_x - cos(cstate.car_yaw);
    before_next_path_start_y = cstate.car_y - sin(cstate.car_yaw);

    ref_x = cstate.car_x;
    ref_y = cstate.car_y;
    ref_yaw = cstate.car_yaw;

    ref_speed = cstate.car_speed;

  }
}

vector<double> eval_state_cost(const car_state &cstate,
                               const vector<vector<double>> &sensor_fusion,
                               const vector<double> &map_waypoints_x,
                               const vector<double> &map_waypoints_y,
                               const vector<double> &map_waypoints_s,
                               const int lane_width,
                               const vector<CarStateNum> &state_to_try,
                               const double delta_t,
                               const double acceleration);
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
              std::vector<double> &next_y_vals) {

  const int lane_width = 4;
  int prev_lane_number = floor(cstate.car_d / lane_width);
  //cout << "lane num:" << prev_lane_number << endl;

  assert(prev_lane_number < 3 && prev_lane_number >= 0);

  //Determine which state that the car is by calculating the lane that the car is staying.
  //Therefore, strictly speaking, "stopping" is not a state.
  int car_lane = floor(cstate.car_d / lane_width);
  assert (car_lane >= 0 && car_lane <= 2);

  const vector<CarStateNum> state_to_try = state_map[car_lane];
  int num_states = state_to_try.size();

  const double delta_t = 0.2;
  const double acceleration = 6;

  vector<double> state_cost = eval_state_cost(cstate,
                                              sensor_fusion,
                                              map_waypoints_x,
                                              map_waypoints_y,
                                              map_waypoints_s,
                                              lane_width,
                                              state_to_try,
                                              delta_t,
                                              acceleration);

  //print out state cost:
  for (int i = 0; i < num_states; i++) {
    cout << "state" << state_to_try[i] << ":" << state_cost[i] << endl;
  }

  auto result_state = min_element(state_cost.begin(), state_cost.end());
  int min_state_index = result_state - state_cost.begin();


  //TODO magic number 50 and 0.02
  cout << "recommended next state:" << state_to_try[min_state_index] << endl;

  /*
  if (previous_path_x.size() < 50) {
    //Only genereate the trajectory points when not enough trajectory points in the buffer
    gen_traj_from_jmt(cstate,
                      state_to_try[min_state_index],
                      previous_path_x,
                      previous_path_y,
                      sensor_fusion,
                      map_waypoints_x,
                      map_waypoints_y,
                      map_waypoints_s,
                      map_waypoints_dx,
                      map_waypoints_dy,
                      next_x_vals,
                      next_y_vals,
                      50,
                      0.02);
  } else {
    next_x_vals = previous_path_x;
    next_y_vals = previous_path_y;
  }
*/

  gen_traj_from_spline(cstate,
                       state_to_try[min_state_index],
                       previous_path_x,
                       previous_path_y,
                       sensor_fusion,
                       map_waypoints_x,
                       map_waypoints_y,
                       map_waypoints_s,
                       map_waypoints_dx,
                       map_waypoints_dy,
                       end_path_s,
                       end_path_d,
                       next_x_vals,
                       next_y_vals,
                       50,
                       0.02);

}
/// Evaluate the cost of the next possible states
/// \param cstate
/// \param sensor_fusion
/// \param map_waypoints_x
/// \param map_waypoints_y
/// \param map_waypoints_s
/// \param lane_width
/// \param state_to_try
/// \param delta_t
/// \param acceleration
/// \return a vector of cost of each state
vector<double> eval_state_cost(const car_state &cstate,
                               const vector<vector<double>> &sensor_fusion,
                               const vector<double> &map_waypoints_x,
                               const vector<double> &map_waypoints_y,
                               const vector<double> &map_waypoints_s,
                               const int lane_width,
                               const vector<CarStateNum> &state_to_try,
                               const double delta_t,
                               const double acceleration) {
  car_state next_car_state;
  int num_states = state_to_try.size();
  vector<double> state_cost(num_states);

  for (int i = 0; i < num_states; i++) {
    state_cost[i] = 0;
  }

  double next_coll_time = 100;
  for (int i = 0; i < num_states; i++) {
    if (i == 0) {
      //evaluate the incumbent state
      next_car_state =
          guess_next_car_state(cstate,
                               delta_t,
                               acceleration,
                               state_to_try[i],
                               map_waypoints_x,
                               map_waypoints_y,
                               map_waypoints_s);

    } else if (state_to_try[i] >= 0 && state_to_try[i] <= 2) {
      //evaluate switching lanes
      next_car_state =
          guess_next_car_state(cstate,
                               delta_t,
                               acceleration,
                               state_to_try[i],
                               map_waypoints_x,
                               map_waypoints_y,
                               map_waypoints_s);

//Add the cost of switching lanes because of the car in the back of the neighbor lane
      state_cost[i] += safe_switch(cstate, state_to_try[i], sensor_fusion, lane_width);

    } else if (state_to_try[i] == stopping) {
      next_car_state =
          guess_next_car_state(cstate,
                               delta_t,
                               -acceleration,
                               state_to_try[0],
                               map_waypoints_x,
                               map_waypoints_y,
                               map_waypoints_s);

    }

    next_coll_time = eval_next_collision(next_car_state, sensor_fusion, map_waypoints_x, map_waypoints_y,
                                         lane_width);

//cout << "time to collide of state" << state_to_try[i] << ":" << next_coll_time << endl;


    double coll_inv = 1 / next_coll_time;
    //if (i != 0 || state_to_try[i] != stopping) {
    //  coll_inv += 0.2;
    //}

    state_cost[i] += (coll_inv + 0.05 * (20 - next_car_state.car_speed));

  }
  return state_cost;
}

void car_to_map_cords_array(const car_state &cstate, vector<double> &next_x_vals, vector<double> &next_y_vals) {
  for (int i = 0; i < next_x_vals.size(); i++) {
    vector<double> nc = car_to_global(next_x_vals[i], next_y_vals[i],
                                      cstate.car_x, cstate.car_y, cstate.car_yaw);

    next_x_vals[i] = nc[0];
    next_y_vals[i] = nc[1];
  }
}

void map_to_car_coords_array(const car_state &cstate, vector<double> &next_map_x, vector<double> &next_map_y) {
  for (int i = 0; i < next_map_x.size(); i++) {
    //convert the map points to car coordinates
    vector<double> nc = global_to_car(next_map_x[i], next_map_y[i],
                                      cstate.car_x, cstate.car_y, cstate.car_yaw);

    next_map_x[i] = nc[0];
    next_map_y[i] = nc[1];

  }
}


