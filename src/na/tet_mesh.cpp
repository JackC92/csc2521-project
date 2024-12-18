#include "na/tet_mesh.h"
#include "Eigen/Core"
#include "Eigen/Geometry"
#include <cassert>
#include <vector>

double compute_tet_volume(const Eigen::Vector3d &v0,
                          const Eigen::Vector3d &v1,
                          const Eigen::Vector3d &v2,
                          const Eigen::Vector3d &v3)
{
    return (v3 - v0).dot((v1 - v0).cross(v2 - v0)) / 6.0;
}

std::vector<double> compute_rest_volumes(const std::vector<Eigen::Vector4i> &tets,
                                         const std::vector<Eigen::Vector3d> &verts)
{
    std::vector<double> rest_volumes(tets.size());
    for (std::size_t idx = 0; idx < tets.size(); ++idx)
    {
        const Eigen::Vector4i &tet = tets[idx];
        const Eigen::Vector3d &v0 = verts[tet[0]];
        const Eigen::Vector3d &v1 = verts[tet[1]];
        const Eigen::Vector3d &v2 = verts[tet[2]];
        const Eigen::Vector3d &v3 = verts[tet[3]];
        rest_volumes[idx] = compute_tet_volume(v0, v1, v2, v3);
        assert(rest_volumes[idx] > 0.0);
    }
    return rest_volumes;
}

std::vector<Eigen::Matrix3d> compute_Dm_inv(const std::vector<Eigen::Vector4i> &tets,
                                            const std::vector<Eigen::Vector3d> &verts)
{
    std::vector<Eigen::Matrix3d> Dm_inv(tets.size());
    for (std::size_t idx = 0; idx < tets.size(); ++idx)
    {
        const Eigen::Vector4i &tet = tets[idx];
        Eigen::Matrix3d Dm;
        Dm.col(0) = verts[tet[1]] - verts[tet[0]];
        Dm.col(1) = verts[tet[2]] - verts[tet[0]];
        Dm.col(2) = verts[tet[3]] - verts[tet[0]];
        Dm_inv[idx] = Dm.fullPivLu().inverse();
        assert((Dm * Dm_inv[idx] - Eigen::Matrix3d::Identity()).lpNorm<Eigen::Infinity>() <= 1.0e-6);
    }
    return Dm_inv;
}

TetMesh::TetMesh(const std::vector<Eigen::Vector3d> &rest_verts,
                 const std::vector<Eigen::Vector3d> &verts,
                 const std::vector<Eigen::Vector4i> &tets)
{
    m_rest_vertices = rest_verts;
    m_vertices = verts;
    m_tets = tets;
    m_rest_volumes = compute_rest_volumes(tets, rest_verts);
    m_Dm_inv = compute_Dm_inv(tets, rest_verts);
    m_vert_is_fixed.resize(verts.size(), false);

    m_masses = std::vector<double>(verts.size(), 0.0);
    for (std::size_t idx = 0; idx < tets.size(); ++idx)
    {
        for (std::size_t i = 0; i < 4; ++i)
        {
            m_masses[tets[idx][i]] += 0.25 * m_rest_volumes[idx];
        }
    }
}

const std::vector<Eigen::Vector3d> &TetMesh::get_vertices() const
{
    return m_vertices;
}

const Eigen::Vector3d &TetMesh::get_vertex(const std::size_t idx) const
{
    return m_vertices[idx];
}

void TetMesh::set_vertex(const std::size_t idx, const Eigen::Vector3d &vert)
{
    m_vertices[idx] = vert;
}

const Eigen::Vector4i &TetMesh::get_tet(const std::size_t tet_idx) const
{
    return m_tets[tet_idx];
}

void TetMesh::set_vertex_fixed(const std::size_t idx, const bool is_fixed)
{
    m_vert_is_fixed[idx] = is_fixed;
}

bool TetMesh::is_vertex_fixed(const std::size_t idx) const
{
    return m_vert_is_fixed[idx];
}

Eigen::Matrix3d TetMesh::compute_F(const std::size_t tet_idx) const
{
    const Eigen::Vector4i &tet = m_tets[tet_idx];
    Eigen::Matrix3d Ds;
    Ds.col(0) = m_vertices[tet[1]] - m_vertices[tet[0]];
    Ds.col(1) = m_vertices[tet[2]] - m_vertices[tet[0]];
    Ds.col(2) = m_vertices[tet[3]] - m_vertices[tet[0]];
    return Ds * m_Dm_inv[tet_idx];
}

Eigen::Matrix3d TetMesh::compute_F(const std::size_t tet_idx, const Eigen::MatrixXd &perturbation) const
{
    const Eigen::Vector4i &tet = m_tets[tet_idx];
    Eigen::Matrix3d Ds;
    Ds.col(0) = (m_vertices[tet[1]] + perturbation.col(1)) - (m_vertices[tet[0]] + perturbation.col(0));
    Ds.col(1) = (m_vertices[tet[2]] + perturbation.col(2)) - (m_vertices[tet[0]] + perturbation.col(0));
    Ds.col(2) = (m_vertices[tet[3]] + perturbation.col(3)) - (m_vertices[tet[0]] + perturbation.col(0));
    return Ds * m_Dm_inv[tet_idx];
}

const double TetMesh::get_mass(const std::size_t idx) const
{
    return m_masses[idx];
}

const double TetMesh::get_rest_volume(const std::size_t tet_idx) const
{
    return m_rest_volumes[tet_idx];
}

const Eigen::Matrix3d &TetMesh::get_Dm_inv(const std::size_t tet_idx) const
{
    return m_Dm_inv[tet_idx];
}

std::vector<std::vector<std::pair<int, int>>> TetMesh::compute_adjacency(const std::vector<Eigen::Vector4i> &tets)
{
    int max_idx = 0;
    for (const Eigen::Vector4i &tet : tets)
    {
        max_idx = std::max(max_idx, tet.maxCoeff());
    }

    std::vector<std::vector<std::pair<int, int>>> adj(max_idx + 1);
    for (std::size_t idx = 0; idx < tets.size(); ++idx)
    {
        const Eigen::Vector4i &tet = tets[idx];
        for (int i = 0; i < 4; ++i)
        {
            adj[tet[i]].emplace_back(idx, i);
        }
    }
    return adj;
}
