#include "na/loader.h"
#include "igl/readMESH.h"

void na::load(const std::string &file, Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::MatrixXi &F)
{
    igl::readMESH(file, V, T, F);
}
