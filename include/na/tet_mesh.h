#ifndef NA_TET_MESH_H
#define NA_TET_MESH_H

#include "Eigen/Core"
#include <vector>

class TetMesh
{
public:
    TetMesh() = default;
    TetMesh(const std::vector<Eigen::Vector3d> &rest_verts,
            const std::vector<Eigen::Vector3d> &verts,
            const std::vector<Eigen::Vector4i> &tets);

    const std::vector<Eigen::Vector3d> &get_vertices() const;
    const Eigen::Vector3d &get_vertex(const std::size_t idx) const;
    void set_vertex(const std::size_t idx, const Eigen::Vector3d &vert);

    const std::vector<Eigen::Vector4i> &get_tets() const;
    const Eigen::Vector4i &get_tet(const std::size_t tet_idx) const;

    void set_vertex_fixed(const std::size_t idx, const bool is_fixed);
    bool is_vertex_fixed(const std::size_t idx) const;

    Eigen::Matrix3d compute_F(const std::size_t tet_idx) const;
    Eigen::Matrix3d compute_F(const std::size_t tet_idx, const Eigen::MatrixXd &perturbation) const;

    const double get_mass(const std::size_t idx) const;
    const double get_rest_volume(const std::size_t tet_idx) const;
    const Eigen::Matrix3d &get_Dm_inv(const std::size_t tet_idx) const;

    static std::vector<std::vector<std::pair<int, int>>> compute_adjacency(const std::vector<Eigen::Vector4i> &tets);

private:
    std::vector<Eigen::Vector3d> m_vertices;
    std::vector<Eigen::Vector3d> m_rest_vertices;
    std::vector<Eigen::Vector4i> m_tets;
    std::vector<double> m_masses;
    std::vector<double> m_rest_volumes;
    std::vector<Eigen::Matrix3d> m_Dm_inv;

    std::vector<bool> m_vert_is_fixed;
};

#endif // !NA_TET_MESH_H
