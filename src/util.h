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




