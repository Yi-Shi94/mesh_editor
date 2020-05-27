//
//  Helpers.h
//  Assignment3
//
//  Created by syd on 2019/11/17.
//

#ifndef Helpers_h
#define Helpers_h

#include <string>
#include <vector>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Geometry>
using namespace std;
using namespace Eigen;

#ifdef _WIN32
#  include <windows.h>
#  undef max
#  undef min
#  undef DrawText
#endif

#ifndef __APPLE__
#  define GLEW_STATIC
#  include <GL/glew.h>
#endif

#ifdef __APPLE__
#   include <OpenGL/gl3.h>
#   define __gl_h_ /* Prevent inclusion of the old gl.h */
#else
#   ifdef _WIN32
#       include <windows.h>
#   endif
#   include <GL/gl.h>
#endif

class VertexArrayObject
{
public:
    unsigned int id;
    
    VertexArrayObject(unsigned int i) : id(i) {}
    
    // Create a new VAO
    void init();
    
    // Select this VAO for subsequent draw calls
    void bind();
    
    // Release the id
    void free();
};

class VertexBufferObject
{
public:
    typedef unsigned int GLuint;
    typedef int GLint;
    
    GLuint id;
    GLuint rows;
    GLuint cols;
    
    VertexBufferObject() : id(0), rows(0), cols(0) {}
    
    // Create a new empty VBO
    void init();
    
    // Updates the VBO with a matrix M
    void update(const Eigen::MatrixXf& M);
    
    // Select this VBO for subsequent draw calls
    void bind();
    
    // Release the id
    void free();
};
class ElementBufferObject
{
public:
    typedef unsigned int GLuint;
    typedef int GLint;
    
    GLuint id;
    GLuint rows;
    GLuint cols;
    
    ElementBufferObject() : id(0), rows(0), cols(0) {}
    
    // Create a new empty VBO
    void init();
    
    // Updates the VBO with a matrix M
    void update(const Eigen::MatrixXi& M);
    
    // Select this VBO for subsequent draw calls
    void bind();
    
    // Release the id
    void free();
};

// This class wraps an OpenGL program composed of two shaders
class Program
{
public:
    typedef unsigned int GLuint;
    typedef int GLint;
    
    GLuint vertex_shader;
    GLuint fragment_shader;
    GLuint program_shader;
    
    Program() : vertex_shader(0), fragment_shader(0), program_shader(0) { }
    
    // Create a new shader from the specified source strings
    bool init(const std::string &vertex_shader_string,
              const std::string &fragment_shader_string,
              const std::string &fragment_data_name);
    
    // Select this shader for subsequent draw calls
    void bind();
    // Release all OpenGL objects
    void free();
    // Return the OpenGL handle of a named shader attribute (-1 if it does not exist)
    GLint attrib(const std::string &name) const;
    // Return the OpenGL handle of a uniform attribute (-1 if it does not exist)
    GLint uniform(const std::string &name) const;
    // Bind a per-vertex array attribute
    GLint bindVertexAttribArray(const std::string &name, VertexBufferObject& VBO) const;
    GLuint create_shader_helper(GLint type, const std::string &shader_string);
    
};
boold 
class TriMesh{
public:
    vector<Vector3d> vertices;
    vector<Vector3d> faces;
    vector<Vector3d> normals;
    
    //MatrixXf indices;
    
    double scale = 1;
    double angle_x = 0;
    double angle_y = 0;
    double angle_z = 0;
    
    double x = 0;
    double y = 0;
    double z = 0;
    
    Vector3d barycenter;
    Vector3d color;
    vector<Vector3d> normals_phong;
    
    int start = 0;
    int tri_num = 0;
    
    int type = -1;
    int render_type;
    
    void set_trans_mat(float x,float y,float z,float angle_x,float angle_y,float angle_z,float scale);
    void set_color(Vector3d color);
    void set_render_type(int r_type);
    
    int get_render_type();
    MatrixXf get_trans_mat();
    MatrixXf get_matrix();
    void cal_normal_matrix();
    MatrixXf get_normal_matrix();
    
    void cal_phong_normal_matrix();
    MatrixXf get_phong_normal_matrix();
    Vector3d get_color();
    
    void get_bary();
    void compute_bary_n_max(int i);
    
    void set_normal();
    
    double is_hit(Vector3d ray_origin,Vector3d ray_direction);
    TriMesh(int target_type,int start);
    
};



void _check_gl_error(const char *file, int line);

/*
class Mesh_list{
public:
    std::vector<Mesh_object> obj_vec;
    std::vector<int> vertice_draw;
    int tri_num;
    int vertice_num;
    Mesh_list(): tri_num(0),vertice_num(0){};
    int hit(glm::vec3 ray_o,glm::vec4 ray_d);
    void add_object(Mesh_object o);
    void delete_object(int ind);
    Eigen::MatrixXf get_normal_matrix();
    Eigen::MatrixXf get_matrix();
    
};
 
///
 */
/// Usage
/// [... some opengl calls]
/// glCheckError();
///
#define check_gl_error() _check_gl_error(__FILE__,__LINE__)

#endif
