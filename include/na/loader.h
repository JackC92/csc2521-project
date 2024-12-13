#ifndef NA_LOADER_H
#define NA_LOADER_H

#include "Eigen/Core"
#include <string>

namespace na
{
    void load(const std::string &file, Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::MatrixXi &F);
} // namespace na

#endif // !NA_LOADER_H
