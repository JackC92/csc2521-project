#include "na/loader.h"
#include "na/meshes.h"
#include "na/neohookean.h"
#include "na/tet_mesh.h"
#include "Eigen/Dense"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/volume_mesh.h"
#include <chrono>
#include <utility>
#include <vector>

TetMesh mesh;
NeoHookean material;
std::vector<std::vector<std::pair<int, int>>> adj;
polyscope::PointCloud *vert_ptr;
polyscope::VolumeMesh *mesh_ptr;

std::vector<Eigen::Vector3d> x, xold;
std::vector<Eigen::Vector3d> v, vold;

struct
{
    int S = 5;
    double h = 1.0 / (60.0 * S);
    int nmax = 10;

    Eigen::Vector3d aext = Eigen::Vector3d(0.0, -9.81, 0.0);

    double density;
    double mu, lambda;
    double kd;

    bool is_start = false;
    bool with_damping = false;
} settings;

void simulate()
{
    static bool print = true;

    double h = settings.h;
    int nmax = settings.nmax;
    double density = settings.density;
    double kd = settings.kd;
    Eigen::Vector3d aext = settings.aext;
    Eigen::Vector3d ahat = aext.normalized();
    double amag = aext.norm();

    // Output of the VBD time stepping algorithm
    std::vector<Eigen::Vector3d> xnew(x.size()), vnew(x.size());

    std::vector<Eigen::Vector3d> inertia(x.size());
    for (std::size_t i = 0; i < x.size(); i++)
    {
        inertia[i] = x[i] + h * v[i] + h * h * aext;
    }

    if (vold.empty())
    {
        // Initialize guess with option (c) Inertia and acceleration
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            xnew[i] = inertia[i];
        }
    }
    else
    {
        // Initialize guess with adaptive initialization
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            Eigen::Vector3d at = (v[i] - vold[i]) / h;
            double acomp = at.dot(ahat);
            if (acomp > amag)
            {
                xnew[i] = x[i] + h * v[i] + h * h * aext;
            }
            else if (acomp < 0.0)
            {
                xnew[i] = x[i] + h * v[i];
            }
            else
            {
                xnew[i] = x[i] + h * v[i] + h * h * (acomp * ahat);
            }
        }
    }

    std::vector<Eigen::Vector3d> buf(x.size());
    for (int n = 1; n <= nmax; ++n)
    {
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            if (mesh.is_vertex_fixed(i))
            {
                xnew[i] = x[i];
                continue;
            }

            // Gauss-Seidel iteration
            double mass = density * mesh.get_mass(i);

            Eigen::Vector3d f = Eigen::Vector3d::Zero();
            Eigen::Matrix3d H = Eigen::Matrix3d::Zero();

            f -= mass * (xnew[i] - inertia[i]);
            H += mass * Eigen::Matrix3d::Identity();

            for (auto &pair : adj[i])
            {
                int tet_idx = pair.first;
                int local_idx = pair.second;
                const Eigen::Vector4i &tet = mesh.get_tet(tet_idx);

                Eigen::Matrix3d Dm_inv = mesh.get_Dm_inv(tet_idx);
                Eigen::Vector3d pFpxi;

                if (local_idx == 0)
                {
                    pFpxi = -Dm_inv.colwise().sum();
                }
                else if (local_idx == 1)
                {
                    pFpxi = Dm_inv.row(0);
                }
                else if (local_idx == 2)
                {
                    pFpxi = Dm_inv.row(1);
                }
                else if (local_idx == 3)
                {
                    pFpxi = Dm_inv.row(2);
                }

                Eigen::MatrixXd perturbation = Eigen::MatrixXd::Zero(3, 4);
                perturbation.col(0) = xnew[tet[0]] - x[tet[0]];
                perturbation.col(1) = xnew[tet[1]] - x[tet[1]];
                perturbation.col(2) = xnew[tet[2]] - x[tet[2]];
                perturbation.col(3) = xnew[tet[3]] - x[tet[3]];
                Eigen::Matrix3d F = mesh.compute_F(tet_idx, perturbation);
                Eigen::Vector9d gradF = material.gradF(F);
                Eigen::Matrix9d hessF = material.hessF(F);

                Eigen::Vector3d dEdxi = gradF.reshaped(3, 3) * pFpxi;
                Eigen::Matrix3d d2Edxi2 = Eigen::Matrix3d::Zero();
                for (Eigen::Index k = 0; k < 3; ++k)
                {
                    for (Eigen::Index j = 0; j < 3; ++j)
                    {
                        d2Edxi2 += hessF.block<3, 3>(3 * j, 3 * k) * pFpxi(j) * pFpxi(k);
                    }
                }

                double vol = mesh.get_rest_volume(tet_idx);
                f -= h * h * vol * dEdxi;
                H += h * h * vol * d2Edxi2;

                if (settings.with_damping)
                {
                    f -= kd * h * d2Edxi2 * (xnew[i] - x[i]);
                    H += kd * h * d2Edxi2;
                }
            }

            if (f.squaredNorm() < 1e-12)
            {
                buf[i] = xnew[i];
                continue;
            }
            else
            {
                Eigen::Vector3d dx = H.fullPivLu().solve(f);
                buf[i] = xnew[i] + dx;
            }
            xnew[i] = buf[i];
        }
    }

    for (std::size_t i = 0; i < x.size(); ++i)
    {
        mesh.set_vertex(i, xnew[i]);
        vnew[i] = (xnew[i] - x[i]) / h;
    }

    xold = x;
    vold = v;
    x = xnew;
    v = vnew;
}

void main_callback()
{
    ImGui::Checkbox("With Damping", &settings.with_damping);
    ImGui::Checkbox("Start Simulation", &settings.is_start);

    if (settings.is_start)
    {
        simulate();
    }

    vert_ptr->updatePointPositions(x);
    mesh_ptr->updateVertexPositions(x);
}

int main(int argc, char *argv[])
{
    settings.density = 1000.0;
    settings.mu = 1e6;
    settings.lambda = 1e7;
    settings.kd = 1e-6;

    std::vector<Eigen::Vector3d> verts;
    std::vector<Eigen::Vector4i> tets;
    load_mesh("../data/bar_test.mesh", verts, tets);

    mesh = TetMesh(verts, verts, tets);
    material = NeoHookean(settings.mu, settings.lambda);
    adj = TetMesh::compute_adjacency(tets);
    x = verts;
    v = std::vector<Eigen::Vector3d>(x.size(), Eigen::Vector3d::Zero());

    for (std::size_t i = 0; i < x.size(); ++i)
    {
        if (x[i][0] == -1.5)
        {
            mesh.set_vertex_fixed(i, true);
        }
    }

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
