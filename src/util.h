//
// Created by Kan-Hua Lee on 2017/10/21.
//

#ifndef PATH_PLANNING_UTIL_H
#define PATH_PLANNING_UTIL_H

#endif //PATH_PLANNING_UTIL_H

#include <vector>
#include <cmath>
#include "spline.h"

using namespace std;

// For converting back and forth between radians and degrees.
constexpr double pi() { return 3.14159265359; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }


vector<double> transform_coords(double global_map_x, double global_map_y,
                                double global_car_x, double global_car_y, double car_yaw)
{
    double local_map_x;
    double local_map_y;

    car_yaw=deg2rad(car_yaw);

    double dx=global_map_x-global_car_x;
    double dy=global_map_y-global_car_y;

    local_map_x=cos(car_yaw)*dx+sin(car_yaw)*dy;
    local_map_y=-sin(car_yaw)*dx+cos(car_yaw)*dy;

    return {local_map_x,local_map_y};
}

vector<double> inv_transform_coords(double local_map_x, double local_map_y,
                                    double global_car_x, double global_car_y, double car_yaw)
{

    double global_map_x;
    double global_map_y;

    car_yaw=deg2rad(car_yaw);

    global_map_x=cos(car_yaw)*local_map_x-sin(car_yaw)*local_map_y+global_car_x;
    global_map_y=sin(car_yaw)*local_map_x+cos(car_yaw)*local_map_y+global_car_y;

    return {global_map_x,global_map_y};

}

void fill_spline(vector<double> &map_x, vector<double> &map_y,
                 vector<double> &traj_x,vector<double> &traj_y)
{

    vector<double> new_map_x=map_x;
    vector<double> new_map_y=map_y;


    // determine the required x

    int end_index=new_map_x.size()-1;

    double xd=new_map_x[end_index]-new_map_x[0];
    double delta_xd=0.25;

    // generate the spline trajectory
    tk::spline s;
    s.set_points(new_map_x,new_map_y);


    // select the points
    // select the first 50 points
    int num_points=50;


    vector<double> new_traj_x;
    vector<double> new_traj_y;


    for (int i=0;i<num_points;i++)
    {
        double test_x=new_map_x[0]+delta_xd*(i+1);
        new_traj_x.push_back(test_x);
        new_traj_y.push_back(s(test_x));
    }

    traj_x=new_traj_x;
    traj_y=new_traj_y;

}

void gen_traj(double start_x,double start_y,vector<double> &map_x, vector<double> &map_y,
              vector<double> &traj_x,vector<double> &traj_y)
{

    vector<double> new_map_x=map_x;
    vector<double> new_map_y=map_y;

    //generate the vector for feeding into spline
    new_map_x.insert(new_map_x.begin(),start_x);
    new_map_y.insert(new_map_y.begin(),start_y);

    cout << "map:value"<<endl;
    for (int i=0;i<new_map_x.size();i++)
    {
        cout << new_map_x[i]<<endl;
    }

    fill_spline(new_map_x,new_map_y,traj_x,traj_y);

}




