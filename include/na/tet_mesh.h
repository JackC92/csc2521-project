#ifndef NA_TET_MESH_H
#define NA_TET_MESH_H

#include "Eigen/Core"
#include <vector>

double compute_tet_volume(const Eigen::Vector3d &v0,
                          const Eigen::Vector3d &v1,
                          const Eigen::Vector3d &v2,
                          const Eigen::Vector3d &v3);

std::vector<double> compute_rest_volumes(const std::vector<Eigen::Vector4i> &tets,
                                         const std::vector<Eigen::Vector3d> &verts);

std::vector<Eigen::Matrix3d> compute_Dm_inv(const std::vector<Eigen::Vector4i> &tets,
                                            const std::vector<Eigen::Vector3d> &verts);

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

    bool is_vertex_fixed(const std::size_t idx) const;

    Eigen::Matrix3d compute_F(const std::size_t tet_idx) const;
    Eigen::Matrix3d compute_F(const std::size_t tet_idx, const Eigen::MatrixXd& perturbation) const;

private:
    std::vector<Eigen::Vector3d> m_vertices;
    std::vector<Eigen::Vector3d> m_rest_vertices;
    std::vector<Eigen::Vector4i> m_tets;
    std::vector<double> m_rest_volumes;
    std::vector<Eigen::Matrix3d> m_Dm_inv;

    std::vector<bool> m_vert_is_fixed;
};

#endif // !NA_TET_MESH_H
