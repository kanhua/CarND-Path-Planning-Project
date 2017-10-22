//
// Created by Kan-Hua Lee on 2017/10/21.
//


#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "util.h"
#include "plstream.h"


void plot2d(const vector<double> &x, const vector<double> &y,const char *filename)
{
    plstream *pls;
    pls = new plstream();

    const int NSIZE=x.size();

    PLFLT plt_x[NSIZE], plt_y[NSIZE];

    for (int i=0;i<NSIZE;i++)
    {
        plt_x[i]=x[i];
        plt_y[i]=y[i];
    }


    // Parse and process command line arguments
    //pls->parseopts( &argc, argv, PL_PARSE_FULL );

    pls->sfnam( filename );       // file name
    pls->sdev( "png" );

    // Initialize plplot
    pls->init();

    // Create a labelled box to hold the plot.

    auto xmin_it=min_element(x.begin(),x.end());
    auto xmax_it=max_element(x.begin(),x.end());

    auto ymin_it=min_element(y.begin(),y.end());
    auto ymax_it=max_element(y.begin(),y.end());

    pls->col0(3);
    pls->env( *xmin_it, *xmax_it, *ymin_it, *ymax_it, 0, 0 );
    pls->lab( "x", "y=100 x#u2#d", "Simple PLplot demo of a 2D line plot" );


    pls->col0(2);
    // Plot the data that was prepared above.
    //pls->line( x.size(), plt_x, plt_y );
    pls->poin( x.size(), plt_x, plt_y ,9);


    // In C++ we don't call plend() to close PLplot library
    // this is handled by the destructor
    delete pls;
}



int main(){
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

    int plot_index_range=5;
    int start_index=5;

    vector<double> plot_map_x(plot_index_range);
    vector<double> plot_map_y(plot_index_range);

    for(int i=0;i<plot_index_range;i++)
    {
        plot_map_x[i]=map_waypoints_x[start_index+i];
        plot_map_y[i]=map_waypoints_y[start_index+i];
    }

    //plot2d(plot_map_x,plot_map_y,"global.png");

    double start_map_x=909.48;
    double start_map_y=1128.67;
    double start_yaw=0;

    for(int i=0;i<plot_index_range;i++)
    {
        vector<double> nc=transform_coords(plot_map_x[i],plot_map_y[i],
        start_map_x,start_map_y,start_yaw);
        plot_map_x[i]=nc[0];
        plot_map_y[i]=nc[1];
    }

    //plot2d(plot_map_x,plot_map_y,"local.png");

    //try spline

    vector<double> fitted_x;
    vector<double> fitted_y;

    //fill_spline(plot_map_x,plot_map_y,fitted_x,fitted_y);

    gen_traj(0,0,plot_map_x,plot_map_y,fitted_x,fitted_y);


    plot2d(fitted_x,fitted_y,"fitted.png");

    for(int i=0;i<fitted_x.size();i++)
    {
        vector<double> nc=inv_transform_coords(fitted_x[i],fitted_y[i],
                                           start_map_x,start_map_y,start_yaw);
        fitted_x[i]=nc[0];
        fitted_y[i]=nc[1];
    }

    plot2d(fitted_x,fitted_y,"inv_fitted.png");
    return 0;
}

