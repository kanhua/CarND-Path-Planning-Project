//
// Created by Kan-Hua Lee on 2017/11/09.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "plstream.h"
#include "vehicle_traj.h"


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

void switch_lane() {
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


int main() {

    test_s_only();
    test_s_d();
    switch_lane();
    return 0;
}