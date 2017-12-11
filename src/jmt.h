//
// Created by Kan-Hua Lee on 2017/12/11.
//

#ifndef PATH_PLANNING_JMT_H
#define PATH_PLANNING_JMT_H

#endif //PATH_PLANNING_JMT_H

#include <vector>
#include "Eigen-3.3/Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

std::vector<double> JMT(std::vector<double> start, std::vector<double> end, double T);

std::vector<double> fill_poly_traj(std::vector<double> a_vec, std::vector<double> t_vec);

