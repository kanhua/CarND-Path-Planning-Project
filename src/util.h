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


vector<double> map_to_car_coords(double global_map_x, double global_map_y,
                                 double global_car_x, double global_car_y, double car_yaw)
{
    double local_map_x;
    double local_map_y;

    double dx=global_map_x-global_car_x;
    double dy=global_map_y-global_car_y;

    local_map_x=cos(car_yaw)*dx+sin(car_yaw)*dy;
    local_map_y=-sin(car_yaw)*dx+cos(car_yaw)*dy;

    return {local_map_x,local_map_y};
}

vector<double> car_to_map_coords(double local_map_x, double local_map_y,
                                 double global_car_x, double global_car_y, double car_yaw)
{

    double global_map_x;
    double global_map_y;

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

void fill_spline_2(vector<double> &map_x, vector<double> &map_y,
                   vector<double> &traj_x,vector<double> &traj_y,
                   int points_to_generate,double desired_speed)
{

    // x-coordinate of the target point
    double max_x_shift=30;

    const double read_in_interval=0.02; // time interval between each reading of the simulator
    //double desired_speed=20; // desired speed in m/s


    tk::spline s;
    s.set_points(map_x,map_y);

    double start_x=map_x[1];
    double start_y=map_y[1];

    double max_end_x=start_x+max_x_shift;
    double max_end_y=s(max_end_x);

    double dx=max_end_x;
    double dy=max_end_y-start_y;

    double delta_d=sqrt(dx*dx+dy*dy);

    double grid_dx=desired_speed*read_in_interval;

    cout << "delta_d:"<<delta_d<<endl;

    double N=delta_d/grid_dx; // required number of segments
    cout <<"N:"<<N<<endl;

    vector<double> new_traj_x;
    vector<double> new_traj_y;

    for (int i=0;i<points_to_generate;i++)
    {
        double x=start_x+(dx/N)*(i+1);
        double y=s(x);

        new_traj_x.push_back(x);
        new_traj_y.push_back(y);
    }

    traj_x=new_traj_x;
    traj_y=new_traj_y;

}

void gen_traj(double start_x,double start_y,vector<double> &map_x, vector<double> &map_y,
              vector<double> &traj_x,vector<double> &traj_y)
{

    vector<double> new_map_x;
    vector<double> new_map_y;

    //new_map_x.push_back(start_x);
    //new_map_y.push_back(start_y);

    int start_index=0;
    //while(map_x[start_index]< start_x)
    //{
    //    start_index++;
    //}

    for (int i=start_index;i<map_x.size();i++)
    {
        new_map_x.push_back(map_x[i]);
        new_map_y.push_back(map_y[i]);
    }


    fill_spline_2(new_map_x,new_map_y,traj_x,traj_y,0,10);

}

void print_map(const vector<double> &map_x, const vector<double> &map_y, int number)
{
    int number_to_print;
    if (number==-1)
    {
        number_to_print=map_x.size();
    }
    else{
        number_to_print=number;
    }

    cout << "map:value"<<endl;
    for (int i=0;i<number_to_print;i++)
    {
        cout << map_x[i] <<","<< map_y[i]<<endl;
    }

}




