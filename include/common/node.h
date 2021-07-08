#ifndef SBS_COMMON_NODE_H
#define SBS_COMMON_NODE_H

#include <string>
#include <vector>

namespace sbs {
namespace common {

/**
 * @brief
 */
class renderable_node_t
{
  public:
    void set_id(std::string const& id);
    std::string const& id() const;

    void set_vao(unsigned int vao);
    void set_vbo(unsigned int vbo);
    void set_ebo(unsigned int ebo);

    unsigned int const& VAO() const;
    unsigned int const& VBO() const;
    unsigned int const& EBO() const;

    unsigned int& VAO();
    unsigned int& VBO();
    unsigned int& EBO();

    void mark_vertices_dirty();
    void mark_indices_dirty();
    void mark_should_render_wireframe();

    void mark_vertices_clean();
    void mark_indices_clean();
    void mark_should_render_triangles();

    bool should_transfer_vertices() const;
    bool should_transfer_indices() const;
    bool should_render_wireframe() const;

    bool is_environment_body() const;
    bool is_physically_simulated_body() const;

    void set_as_environment_body();
    void set_as_physically_simulated_body();

    std::vector<float> const& get_cpu_vertex_buffer() const;
    std::vector<std::uint32_t> const& get_cpu_index_buffer() const;

    virtual void prepare_vertices_for_rendering() = 0;
    virtual void prepare_indices_for_rendering()  = 0;

  protected:
    void transfer_vertices_for_rendering(std::vector<float>&& vertices);
    void transfer_indices_for_rendering(std::vector<std::uint32_t>&& indices);

  private:
    std::string id_;

    struct render_state_t
    {
        bool should_transfer_vertices = true;
        bool should_transfer_indices  = true;
        bool should_render_wireframe  = false;
    } render_state_;

    enum class body_type_t { environment, physical } body_type_;

    std::vector<float> cpu_vertex_buffer_;        ///< (x,y,z,nx,ny,nz,r,g,b)
    std::vector<std::uint32_t> cpu_index_buffer_; ///< (v1, v2, v3)

    unsigned int VBO_;
    unsigned int VAO_;
    unsigned int EBO_;
};

} // namespace common
} // namespace sbs

#endif // SBS_COMMON_NODE_H