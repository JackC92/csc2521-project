#include "na/loader.h"
#include "Eigen/Core"
#include "igl/readMESH.h"
#include <string>
#include <vector>

void load_mesh(const std::string &file, std::vector<Eigen::Vector3d> &verts, std::vector<Eigen::Vector4i> &tets)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
    Eigen::MatrixXi F;
    igl::readMESH(file, V, T, F);

    verts.resize(V.rows());
    for (Eigen::Index idx = 0; idx < V.rows(); ++idx)
    {
        verts[idx] = V.row(idx);
    }
    tets.resize(T.rows());
    for (Eigen::Index idx = 0; idx < T.rows(); ++idx)
    {
        tets[idx] = T.row(idx);
    }
}
