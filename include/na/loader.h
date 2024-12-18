#ifndef NA_LOADER_H
#define NA_LOADER_H

#include "Eigen/Core"
#include <string>
#include <vector>

void load_mesh(const std::string &file, std::vector<Eigen::Vector3d> &verts, std::vector<Eigen::Vector4i> &tets);

#endif // !NA_LOADER_H
