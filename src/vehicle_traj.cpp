//
// Created by Kan-Hua Lee on 2017/11/04.
//


#include "vehicle_traj.h"


void
get_lane_cost(double car_s, const vector<vector<double>> &sensor_fusion, const int lane_width, int prev_lane_number,
              vector<double> &lane_cost);

void map_to_car_coords_array(const car_state &cstate, vector<double> &next_map_x, vector<double> &next_map_y);

void car_to_map_cords_array(const car_state &cstate, vector<double> &next_x_vals, vector<double> &next_y_vals);

constexpr double pi() { return 3.14159265359; }

double deg2rad(double x) { return x * pi() / 180; }

double rad2deg(double x) { return x * 180 / pi(); }

double mph2mps(double x) { return x * 0.44704; }

double mps2mph(double x) { return x * 2.23694; }

double distance(double x1, double y1, double x2, double y2) {
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
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

    return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta,
                         const vector<double> &maps_x, const vector<double> &maps_y) {
    int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

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


vector<double> JMT(vector<double> start, vector<double> end, double T) {
    /*
    Calculate the Jerk Minimizing Trajectory that connects the initial state
    to the final state in time T.

    INPUTS

    start - the vehicles start location given as a length three array
        corresponding to initial values of [s, s_dot, s_double_dot]

    end   - the desired end state for vehicle. Like "start" this is a
        length three array.

    T     - The duration, in seconds, over which this maneuver should occur.

    OUTPUT
    an array of length 6, each value corresponding to a coefficent in the polynomial
    s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

    EXAMPLE

    > JMT( [0, 10, 0], [10, 10, 0], 1)
    [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
    */

    MatrixXd A = MatrixXd(3, 3);
    A << T * T * T, T * T * T * T, T * T * T * T * T,
            3 * T * T, 4 * T * T * T, 5 * T * T * T * T,
            6 * T, 12 * T * T, 20 * T * T * T;

    MatrixXd B = MatrixXd(3, 1);
    B << end[0] - (start[0] + start[1] * T + .5 * start[2] * T * T),
            end[1] - (start[1] + start[2] * T),
            end[2] - start[2];

    MatrixXd Ai = A.inverse();

    MatrixXd C = Ai * B;

    vector<double> result = {start[0], start[1], .5 * start[2]};
    for (int i = 0; i < C.size(); i++) {
        result.push_back(C.data()[i]);
    }

    return result;

}

void fill_spline(vector<double> &map_x, vector<double> &map_y, vector<double> &traj_x, vector<double> &traj_y,
                 int points_to_generate, double desired_speed, double acceleration, double car_speed) {


    const double read_in_interval = 0.02; // time interval between each reading of the simulator
    //double desired_speed=20; // desired speed in m/s


    tk::spline s;
    s.set_points(map_x, map_y);


    vector<double> new_traj_x;
    vector<double> new_traj_y;

    double current_x = map_x[1];

    double dx = map_x[1] - map_x[0];
    double dy = map_y[1] - map_y[0];

    double inst_speed = 0; // instaneous speed

    if (car_speed < 0.1) {
        inst_speed = car_speed;
    } else {
        inst_speed = sqrt(dx * dx + dy * dy) / read_in_interval;
    }


    double theta = atan2(dy, dx);

    for (int i = 0; i < points_to_generate; i++) {

        if (inst_speed < desired_speed) {
            inst_speed += acceleration * read_in_interval;
        } else {
            inst_speed = desired_speed;
        }

        current_x += inst_speed * cos(theta) * read_in_interval;
        double y = s(current_x);

        new_traj_x.push_back(current_x);
        new_traj_y.push_back(y);
    }

    traj_x = new_traj_x;
    traj_y = new_traj_y;

}


vector<double> fill_poly_traj(vector<double> a_vec, vector<double> t_vec) {
    vector<double> s_traj(t_vec.size());

    for (int i = 0; i < t_vec.size(); i++) {
        double s = 0;
        for (int j = 0; j < a_vec.size(); j++) {
            s += a_vec[j] * pow(t_vec[i], j);
        }
        s_traj[i] = s;

    }
    return s_traj;
}

vector<double> arange(double lower_bound, double higher_bound, double delta_t) {
    vector<double> t;
    t.push_back(lower_bound);
    double current_t = lower_bound + delta_t;
    while (current_t < higher_bound) {
        t.push_back(current_t);
        current_t += delta_t;
    }

    return t;

}


void fill_jmt(vector<double> &map_x, vector<double> &map_y, vector<double> &traj_x, vector<double> &traj_y,
              int points_to_generate, double desired_speed, double acceleration, double car_speed) {

    // x-coordinate of the target point
    double max_x_shift = 30;

    const double read_in_interval = 0.02; // time interval between each reading of the simulator
    //double desired_speed=20; // desired speed in m/s

    //For now, use the first section to acc. to the full desired speed
    vector<double> new_traj_x;
    vector<double> new_traj_y;

    double start_x = map_x[1];
    double start_y = map_y[1];


    double next_x = map_x[2];
    double next_y = map_y[2];

    for (int i = 2; i < map_x.size(); i++) {

        if (i == 2) {

            double ds = sqrt(pow(next_x - start_x, 2) + pow(next_y - start_y, 2));

            double t_required = ds / ((car_speed + desired_speed) / 2);
            vector<double> x_coef = JMT({start_x, desired_speed, 0}, {next_x, desired_speed, 0}, t_required);
            vector<double> y_coef = JMT({start_y, 0, 0}, {next_y, 0, 0}, t_required);

            vector<double> t = arange(0, t_required, read_in_interval);

            vector<double> x = fill_poly_traj(x_coef, t);
            vector<double> y = fill_poly_traj(y_coef, t);

            new_traj_x.insert(new_traj_x.end(), x.begin(), x.end());
            new_traj_y.insert(new_traj_y.end(), y.begin(), y.end());
        }

    }

    traj_x.clear();
    traj_y.clear();

    for (int i = 0; i < points_to_generate; i++) {
        traj_x.push_back(new_traj_x[i]);
        traj_y.push_back(new_traj_y[i]);

    }

    print_map(traj_x, traj_y, 10);


}


vector<double> map_to_car_coords(double global_map_x, double global_map_y,
                                 double global_car_x, double global_car_y, double car_yaw) {
    double local_map_x;
    double local_map_y;

    double dx = global_map_x - global_car_x;
    double dy = global_map_y - global_car_y;

    local_map_x = cos(car_yaw) * dx + sin(car_yaw) * dy;
    local_map_y = -sin(car_yaw) * dx + cos(car_yaw) * dy;

    return {local_map_x, local_map_y};
}

vector<double> car_to_map_coords(double local_map_x, double local_map_y,
                                 double global_car_x, double global_car_y, double car_yaw) {

    double global_map_x;
    double global_map_y;

    global_map_x = cos(car_yaw) * local_map_x - sin(car_yaw) * local_map_y + global_car_x;
    global_map_y = sin(car_yaw) * local_map_x + cos(car_yaw) * local_map_y + global_car_y;

    return {global_map_x, global_map_y};


}

void print_map(const vector<double> &map_x, const vector<double> &map_y, int number) {
    int number_to_print;
    if (number == -1) {
        number_to_print = map_x.size();
    } else {
        number_to_print = number;
    }

    cout << "map:value" << endl;
    for (int i = 0; i < number_to_print; i++) {
        cout << map_x[i] << "," << map_y[i] << endl;
    }


}


double eval_state(double delta_t,
                  const vector<double> &next_x_val,
                  const vector<double> &next_y_val,
                  const vector<vector<double>> &sensor_fusion, const vector<double> &map_x,
                  const vector<double> &map_y) {
    const double lane_width = 4;
    const double simulator_interval = 0.02;
    vector<double> lane_cost(3);

    int path_index = floor(simulator_interval * delta_t);
    assert(path_index > 0);
    double car_next_x = next_x_val[path_index];
    double car_next_y = next_y_val[path_index];

    double car_theta = atan2(next_x_val[path_index] - next_x_val[path_index - 1],
                             next_y_val[path_index] - next_y_val[path_index - 1]);

    vector<double> car_nc = getFrenet(car_next_x, car_next_y, car_theta, map_x, map_y);
    double car_next_s = car_nc[0];
    double car_next_d = car_nc[1];

    int car_lane = floor(car_next_d / lane_width);


    for (int i = 0; i < sensor_fusion.size(); i++) {
        double neighbor_car_d = sensor_fusion[i][6];
        double neighbor_car_s = sensor_fusion[i][5];

        double neighbor_car_x = sensor_fusion[i][1];
        double neighbor_car_y = sensor_fusion[i][2];

        double vx = sensor_fusion[i][3];
        double vy = sensor_fusion[i][4];
        double neighbor_car_theta = atan2(vx, vy);

        vector<double> nc = getFrenet(neighbor_car_x + delta_t * vx, neighbor_car_y + delta_t * vy,
                                      neighbor_car_theta, map_x, map_y);

        double neighbor_car_next_s = nc[0];
        double neighbor_car_next_d = nc[1];


        int on_lane = floor(neighbor_car_next_d / lane_width);

        double car_dist = neighbor_car_next_s - car_next_s;

        double cost = 1 / car_dist;
        if (lane_cost[on_lane] < cost && cost > 0) {
            //only consider the cars in front
            lane_cost[on_lane] = cost;
        }

    }

    return lane_cost[car_lane];

}


void gen_traj_from_spline_x(car_state &cstate,
                            int lane_number,
                            vector<double> &previous_path_x,
                            vector<double> &previous_path_y,
                            const vector<vector<double>> &sensor_fusion, const vector<double> &map_waypoints_x,
                            const vector<double> &map_waypoints_y, const vector<double> &map_waypoints_dx,
                            const vector<double> &map_waypoints_dy, vector<double> &next_x_vals,
                            vector<double> &next_y_vals) {
    next_x_vals.clear();
    next_y_vals.clear();

    double acceleration = 6;

    int total_future_points = 50;

    double next_path_start_x = 0;
    double next_path_start_y = 0;

    double next_path_start2_x = 0;
    double next_path_start2_y = 0;

    cout << "previous path size: " << previous_path_x.size() << endl;

    double ref_x = 0;
    double ref_y = 0;
    double ref_yaw = 0;

    double desired_speed = 20;


    const int lane_width = 4; // the width of the lane

    // Tell which lane that the car currently stays
    // The lane number that the car stays on. The lane next to the center line is zero.
    int prev_lane_number = floor(cstate.car_d / lane_width);

    assert(lane_number >= -1 && lane_number < 3);

    if (lane_number == -1) {
        lane_number = prev_lane_number;
    }

    int prev_points = previous_path_x.size();
    assert(previous_path_x.size() == previous_path_y.size());

    // Set the starting point of the next generated path
    if (previous_path_x.size() > 2) {
        //Use the end points of previous_path_x
        next_path_start_x = previous_path_x[prev_points - 2];
        next_path_start_y = previous_path_y[prev_points - 2];

        next_path_start2_x = previous_path_x[prev_points - 1];
        next_path_start2_y = previous_path_y[prev_points - 1];

        ref_x = next_path_start2_x;
        ref_y = next_path_start2_y;
        ref_yaw = atan2(next_path_start2_y - next_path_start_y, next_path_start2_x - next_path_start_x);


    } else {

        //Use the current coordinates of the car
        next_path_start2_x = cstate.car_x;
        next_path_start2_y = cstate.car_y;

        next_path_start_x = cstate.car_x - cos(cstate.car_yaw);
        next_path_start_y = cstate.car_y - sin(cstate.car_yaw);

        ref_x = cstate.car_x;
        ref_y = cstate.car_y;
        ref_yaw = cstate.car_yaw;

    }


    //Find the closest index

    int closest_index = NextWaypoint(ref_x, ref_y,
                                     ref_yaw, map_waypoints_x, map_waypoints_y);

    // Use the next index to start if changing lane. This is to avoid the instability when switching lanes
    if (lane_number != prev_lane_number) {
        closest_index++;
    }


    //Get the map coordinates of the next few points
    //Assuming staying on the second lane at the moment

    int num_next_index = 3;

    vector<double> next_map_x;
    vector<double> next_map_y;

    next_map_x.push_back(next_path_start_x);
    next_map_x.push_back(next_path_start2_x);

    next_map_y.push_back(next_path_start_y);
    next_map_y.push_back(next_path_start2_y);


    // Construct the array for spline fitting
    for (int i = 0; i < num_next_index; i++) {
        int waypoints_index = closest_index + i;

        next_map_x.push_back(map_waypoints_x[waypoints_index] +
                             ((lane_number + 0.5) * lane_width) * map_waypoints_dx[waypoints_index]);
        next_map_y.push_back(map_waypoints_y[waypoints_index] +
                             ((lane_number + 0.5) * lane_width) * map_waypoints_dy[waypoints_index]);

    }

    map_to_car_coords_array(cstate, next_map_x, next_map_y);


    // Find next points
    int points_to_generate = total_future_points - prev_points;

    fill_spline(next_map_x, next_map_y, next_x_vals, next_y_vals, points_to_generate,
                desired_speed, acceleration, cstate.car_speed);


    // Convert the coordinates back to the map coordinates

    car_to_map_cords_array(cstate, next_x_vals, next_y_vals);

    if (prev_points > 2) {
        next_x_vals.insert(next_x_vals.begin(), previous_path_x.begin(), previous_path_x.end());
        next_y_vals.insert(next_y_vals.begin(), previous_path_y.begin(), previous_path_y.end());
    }

}


void gen_traj_from_spline(car_state &cstate,
                          vector<double> &previous_path_x,
                          vector<double> &previous_path_y,
                          const vector<vector<double>> &sensor_fusion, const vector<double> &map_waypoints_x,
                          const vector<double> &map_waypoints_y, const vector<double> &map_waypoints_dx,
                          const vector<double> &map_waypoints_dy, vector<double> &next_x_vals,
                          vector<double> &next_y_vals) {

    const int lane_width = 4;
    int prev_lane_number = floor(cstate.car_d / lane_width);
    cout << "lane num:" << prev_lane_number << endl;

    assert(prev_lane_number < 3 && prev_lane_number >= 0);

    vector<double> lane_cost(3);

    get_lane_cost(cstate.car_s, sensor_fusion, lane_width, prev_lane_number, lane_cost);

    auto result_lane = min_element(lane_cost.begin(), lane_cost.end());

    // find the lane number with the lowest cost
    int lane_number = (result_lane - lane_cost.begin());

    cout << "selected lane number:" << lane_number << endl;
    assert(lane_number >= 0 && lane_number < 3);

    gen_traj_from_spline_x(cstate,
                           lane_number,
                           previous_path_x,
                           previous_path_y,
                           sensor_fusion, map_waypoints_x,
                           map_waypoints_y, map_waypoints_dx,
                           map_waypoints_dy, next_x_vals,
                           next_y_vals);


    cout << "total next points:" << next_x_vals.size() << endl;

    //print_map(next_x_vals, next_y_vals, 10);
}

void car_to_map_cords_array(const car_state &cstate, vector<double> &next_x_vals, vector<double> &next_y_vals) {
    for (int i = 0; i < next_x_vals.size(); i++) {
        vector<double> nc = car_to_map_coords(next_x_vals[i], next_y_vals[i],
                                              cstate.car_x, cstate.car_y, cstate.car_yaw);

        next_x_vals[i] = nc[0];
        next_y_vals[i] = nc[1];
    }
}

void map_to_car_coords_array(const car_state &cstate, vector<double> &next_map_x, vector<double> &next_map_y) {
    for (int i = 0; i < next_map_x.size(); i++) {
        //convert the map points to car coordinates
        vector<double> nc = map_to_car_coords(next_map_x[i], next_map_y[i],
                                              cstate.car_x, cstate.car_y, cstate.car_yaw);

        next_map_x[i] = nc[0];
        next_map_y[i] = nc[1];

    }
}

void
get_lane_cost(double car_s, const vector<vector<double>> &sensor_fusion, const int lane_width, int prev_lane_number,
              vector<double> &lane_cost) {

    for (int i = 0; i < 3; i++) {
        lane_cost[i] = 0;
    }


    // Sense other cars on the road
    cout << "reading sensor fusion data" << endl;
    for (int i = 0; i < sensor_fusion.size(); i++) {
        double neighbor_car_d = sensor_fusion[i][6];
        int on_lane = floor(neighbor_car_d / lane_width);

        double neighbor_car_s = sensor_fusion[i][5];
        double car_dist = neighbor_car_s - car_s;

        double cost = 1 / car_dist;
        if (lane_cost[on_lane] < cost) {
            lane_cost[on_lane] = cost;
        }

    }

    // "incumbent" advantage
    lane_cost[prev_lane_number] -= 0.01;

    cout << "lane cost:" << endl;
    for (int i = 0; i < 3; i++) {
        cout << lane_cost[i] << endl;
    }

}


