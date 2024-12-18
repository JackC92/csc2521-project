#include "na/meshes.h"
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
                  std::vector<Eigen::Vector4i> &tets)
{
    verts.clear();
    tets.clear();

    double dx = (xmax - xmin) / (nx - 1);
    double dy = (ymax - ymin) / (ny - 1);
    double dz = (zmax - zmin) / (nz - 1);

    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                verts.push_back(Eigen::Vector3d(xmin + i * dx, ymin + j * dy, zmin + k * dz));
            }
        }
    }

    for (int k = 0; k < nz - 1; ++k)
    {
        for (int j = 0; j < ny - 1; ++j)
        {
            for (int i = 0; i < nx - 1; ++i)
            {
                int idx000 = (i + 0) + (j + 0) * nx + (k + 0) * nx * ny;
                int idx001 = (i + 0) + (j + 0) * nx + (k + 1) * nx * ny;
                int idx010 = (i + 0) + (j + 1) * nx + (k + 0) * nx * ny;
                int idx011 = (i + 0) + (j + 1) * nx + (k + 1) * nx * ny;
                int idx100 = (i + 1) + (j + 0) * nx + (k + 0) * nx * ny;
                int idx101 = (i + 1) + (j + 0) * nx + (k + 1) * nx * ny;
                int idx110 = (i + 1) + (j + 1) * nx + (k + 0) * nx * ny;
                int idx111 = (i + 1) + (j + 1) * nx + (k + 1) * nx * ny;
                tets.emplace_back(idx000, idx010, idx011, idx111);
                tets.emplace_back(idx000, idx011, idx001, idx111);
                tets.emplace_back(idx000, idx001, idx101, idx111);
                tets.emplace_back(idx000, idx101, idx100, idx111);
                tets.emplace_back(idx000, idx100, idx110, idx111);
                tets.emplace_back(idx000, idx110, idx010, idx111);
            }
        }
    }
}