#include "na/loader.h"
#include "na/neohookean.h"
#include "na/tet_mesh.h"
#include "Eigen/Core"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/volume_mesh.h"
#include <vector>

TetMesh mesh;
NeoHookean material;
polyscope::PointCloud *vert_ptr;
polyscope::VolumeMesh *mesh_ptr;

void main_callback()
{
}

int main(int argc, char *argv[])
{
    std::vector<Eigen::Vector3d> verts;
    std::vector<Eigen::Vector4i> tets;
    na::load_mesh("../data/bar990.mesh", verts, tets);
    mesh = TetMesh(verts, verts, tets);
    material = NeoHookean(1.0e6, 1.0e7);

    polyscope::options::maxFPS = -1;
    polyscope::options::giveFocusOnShow = true;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::state::userCallback = main_callback;

    polyscope::init();
    vert_ptr = polyscope::registerPointCloud("Vertices", verts);
    mesh_ptr = polyscope::registerTetMesh("Mesh", verts, tets);
    polyscope::show();

    return 0;
}
