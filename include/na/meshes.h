#ifndef NA_MESHES_H
#define NA_MESHES_H

#include "Eigen/Core"
#include <vector>

void generate_bar(const double xmin,
                  const double xmax,
                  const double ymin,
                  const double ymax,
                  const double zmin,
                  const double zmax,
                  const int nx,
                  const int ny,
                  const int nz,
                  std::vector<Eigen::Vector3d> &verts,
                  std::vector<Eigen::Vector4i> &tets);

#endif // !NA_MESHES_H
