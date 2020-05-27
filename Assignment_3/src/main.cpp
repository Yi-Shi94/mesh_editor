// This example is heavily based on the tutorial at https://open.gl

// OpenGL Helpers to reduce the clutter
#include "Helpers.h"

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>
#else
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>
#endif

// Linear Algebra Library
#include <Eigen/Core>
// Timer
#include <chrono>

#include <iostream>
#include <fstream>
#include <math.h>
#include "Helpers.h"

using namespace std;
using namespace Eigen;

// VertexBufferObject wrapper
VertexBufferObject VBO;
VertexBufferObject VBO_C;
VertexBufferObject VBO_N;

// Contains the vertex positions
vector<TriMesh> Objects;

MatrixXf V(3,3);
MatrixXf V_color(3,3);
MatrixXf V_normal(3,3);

MatrixXf V_trans(4,4);
Matrix4f V_view(4,4);
Matrix4f V_proj(4,4);
Matrix4f V_scope(4,4);
int selected_index = -1;
int ball_flag = 0;
# define pi 3.1415926
/*CONFIGURE*/

//int camera_pers_type = 0; //0,1 //perspective,ortho
int camera_pers_type = 0 ;
int screenWidth = 640;
int screenHeight = 480;
int ball = 0;

Vector3d light_pos(0.,0.,5.);
Vector3d light_color(1.,1.,1.);
Vector3d cam_pos(0.,0.,5.);
/************/

MatrixXf scope_transform(float aspect_ratio){
    Matrix4f vmat;
    vmat<<aspect_ratio, 0, 0, 0,
          0,            1, 0, 0,
          0,            0, 1, 0,
          0,            0, 0, 1;
    return vmat;
}

MatrixXf view_transform(Vector3d eye, Vector3d center){
    Matrix4f camera_lookat;
    Vector3d up(0.0f,1.0f,0.0f);
    Vector3d X,Y,Z;
    Z = (eye-center).normalized();
    Y = up;
    X = Y.cross(Z);
    Y = Z.cross(X);
    
    X = X.normalized();
    Y = Y.normalized();
    
    camera_lookat<<X(0),X(1),X(2),-X.dot(eye),
                   Y(0),Y(1),Y(2),-Y.dot(eye),
                   Z(0),Z(1),Z(2),-Z.dot(eye),
                   0   ,0   ,0   ,1.0f;
    return camera_lookat;
}

MatrixXf projection_transform(int cam_type, float fov, int width, int height , float near, float far){
    Matrix4f pmat;
    
    if (cam_type==0){
        float yScale = 1./(tan(fov*3.1415926535/360));
        float xScale = yScale;
        //cout<<"xcale"<<xScale<<' '<<yScale<<endl;
        
        pmat << xScale, 0, 0, 0,
                0, yScale, 0, 0,
                0, 0,  -(near+far)/(far-near),-2*near*far/(far-near),
                0, 0, -1, 0;
    } else {
        float left  = -3;
        float right =  3;
        float up    =  3;
        float down  = -3;
        
        float xScale = 2./(right-left);
        float yScale = 2./(up-down);
        //cout<<"xcale"<<xScale<<' '<<yScale<<endl;
        pmat << xScale, 0, 0, -(right+left)/(right-left),
                 0, yScale, 0,-(up+down)/(up-down),
                 0, 0,  -2./(far-near),-(near+far)*1.0/(far-near),
                 0, 0, 0, 1;
    }
    return pmat;
}


void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);

    // Get the size of the window
    int width, height;
    glfwGetWindowSize(window, &width, &height);

    Vector4f p_screen(xpos,height-1-ypos,0,1);
    Vector4f p_canonical((p_screen[0]/width)*2-1,(p_screen[1]/height)*2-1,0,1);
    
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS){
        Vector4f p_cl = V_scope.inverse()*p_canonical;
        Vector4f ray_cl(p_cl(0),p_cl(1), -1.0, 1.0);
        Vector4f ray_view= V_proj.inverse() * ray_cl;

        Vector4f ray_view_d(ray_view(0),ray_view(1),-1.,0.);
        
        cout<<'a'<<ray_view_d<<'b'<<endl;
        
        Vector4f ray_world_3d = (V_view.inverse() * ray_view_d).normalized();
    
        cout<<'c'<<ray_world_3d<<'d'<<endl;
        
        Vector3d ray_wor_3d;
        if (camera_pers_type==0)
        {
            ray_wor_3d = Vector3d(ray_world_3d(0),ray_world_3d(1),ray_world_3d(2));
        } else {
            Vector4f ray_dir(0.,0., -1.,0.);
            Vector3d center(0,0,0);
            ray_dir = V_view.inverse() * ray_dir;
            cout<<ray_dir<<'l'<<endl;
            ray_wor_3d = Vector3d(ray_dir(0),ray_dir(1),ray_dir(2));
          
        }
        //cout<<ray_wor_3d<<endl;
        double min_dis = 100000;
        int min_index = -1;
        for (int i=0;i<Objects.size();i++){
            double t;
            if (camera_pers_type==0){
                t = Objects[i].is_hit(cam_pos, ray_wor_3d);
            } else {
                Vector3d pos(ray_view_d(0),ray_view_d(1),ray_view_d(2));
                t = Objects[i].is_hit(pos, ray_wor_3d);
            }
            if (t<0) continue;
            else {
                if (min_dis > t) {
                    min_dis = t;
                    min_index = i;
                }
            }
        }
        if (min_index>=0) selected_index = min_index;
    }
}

Vector3d fromCart2Sphere(float x, float y, float z){
    //r, theta, phi
    float r     = sqrt(x*x+y*y+z*z);
    float theta = acos(z/(r+0.000001));
    float phi   = atan2(y,x);
    Vector3d spherical(r,theta,phi);
    
    
    return spherical;
}

Vector3d fromSphere2Cart(float r,float theta,float phi){
    float x = r*sin(theta)*cos(phi);
    float y = r*sin(theta)*sin(phi);
    float z = r*cos(theta);
    Vector3d cart(x,y,z);
    return cart;
}

void ball_move(float dr,float dtheta,float dphi){
    
    Vector3d spher = fromCart2Sphere(cam_pos(0),cam_pos(1),cam_pos(2));
    float r = dr+spher(0);
    
    float addon = dtheta*pi/180;
    float theta = spher(1)+addon;
    float phi = dphi*pi/180+spher(2);
    if (r<0) r=0;
    cout<<theta/pi*180<<' '<<phi/pi*180<<endl;
    
    Vector3d cart = fromSphere2Cart(r,theta,phi);
    cam_pos(0) = cart(0);
    cam_pos(1) = cart(1);
    cam_pos(2) = cart(2);
}

void append_mesh(TriMesh obj){
    MatrixXf V_tmp = obj.get_matrix();
    MatrixXf V_ver_show(V.rows(),V.cols()+V_tmp.cols());
    V_ver_show<<V,V_tmp;
    V = V_ver_show;
    
    MatrixXf V_normal_tmp = obj.get_normal_matrix();
    MatrixXf V_normal_show(3,V_normal.cols()+V_normal_tmp.cols());
    //cout<<'s'<<V_normal_tmp.rows()<<' '<<V_normal_tmp.cols()<<' '<<V_normal.rows()<<' '<<V_normal.cols()<<endl;
    V_normal_show<<V_normal,V_normal_tmp;
    V_normal = V_normal_show;
    
    MatrixXf V_color_show(3,V.cols());
    for (int i=0;i<V.cols();i++){
        V_color_show.col(i)<<0.4,0.7,0.2;
    }
    V_color = V_color_show;
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    // Update the position of the first vertex if the keys 1,2, or 3 are pressed
    switch (key)
    {
        case GLFW_KEY_1:
        {
            if (action == GLFW_PRESS)
            {
                TriMesh obj(0,V.cols());
                Objects.push_back(obj);
                append_mesh(obj);
                cout << obj.get_trans_mat();
                cout <<"Total num obj: "<<Objects.size()<<" Verices: "<<V.cols()<<endl;
                VBO.update(V);
                VBO_C.update(V_color);
                VBO_N.update(V_normal);
            }
            break;
                
        }
        case GLFW_KEY_2:
        {
            if (action == GLFW_PRESS)
            {
                TriMesh obj1(1,V.cols());
                Objects.push_back(obj1);
                append_mesh(obj1);
                cout <<"Total num obj: "<<Objects.size()<<endl;
                VBO.update(V);
                VBO_C.update(V_color);
                VBO_N.update(V_normal);
            }
            break;
        }
            
        case GLFW_KEY_3:
        {
            if (action == GLFW_PRESS)
            {
                TriMesh obj2(2,V.cols());
                Objects.push_back(obj2);
                append_mesh(obj2);
                cout <<"Total num obj: "<<Objects.size()<<endl;
                VBO.update(V);
                VBO_C.update(V_color);
                VBO_N.update(V_normal);
            }
            break;
        }
            
        case GLFW_KEY_4:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                
                if (Objects[selected_index].render_type==2){
                    MatrixXf V_normal_tmp = Objects[selected_index].get_normal_matrix();
                    int start = Objects[selected_index].start;
                    int num_last = Objects[selected_index].tri_num*3;
                    V_normal.block(0,start,3,num_last) = V_normal_tmp;
                    VBO_N.update(V_normal);
                }
                Objects[selected_index].set_render_type(0);
               
            }
            break;
        }
        
        case GLFW_KEY_5:
        {
            if (action == GLFW_PRESS && selected_index>=0){

                MatrixXf V_normal_tmp = Objects[selected_index].get_normal_matrix();
                int start = Objects[selected_index].start;
                int length = Objects[selected_index].tri_num*3;
                V_normal.block(0,start,3,length) = V_normal_tmp;
                VBO_N.update(V_normal);
                 
                Objects[selected_index].set_render_type(1);
            }
            break;
        }
            
        case GLFW_KEY_6:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                if (Objects[selected_index].render_type!=2){
                    MatrixXf V_normal_tmp = Objects[selected_index].get_phong_normal_matrix();
                    int start = Objects[selected_index].start;
                    int length = Objects[selected_index].tri_num*3;
                    V_normal.block(0,start,3,length) = V_normal_tmp;
                    VBO_N.update(V_normal);
                }
                Objects[selected_index].set_render_type(2);
            }
            break;
        }
        case GLFW_KEY_7:
        {
            MatrixXf V_normal_tmp = Objects[selected_index].get_normal_matrix();
            int start = Objects[selected_index].start;
            int length = Objects[selected_index].tri_num*3;
            V_normal.block(0,start,3,length) = V_normal_tmp;
            VBO_N.update(V_normal);
            Objects[selected_index].set_render_type(-1);
            break;
        }
            
        case GLFW_KEY_9:
        {
            if (action == GLFW_PRESS){
               if (ball==0) ball=1;
               else ball= 0;
            }
            break;
        }
            
        
        case GLFW_KEY_0:
        {
            if (action == GLFW_PRESS){
                if (camera_pers_type==0) camera_pers_type=1;
                else camera_pers_type=0;
            }
            break;
        }
        
        case GLFW_KEY_L:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                Objects[selected_index].scale *= 1.2;
            }
            break;
        }
        
        case GLFW_KEY_K:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                Objects[selected_index].scale /= 1.2;
            }
            break;
        }
        
        case GLFW_KEY_T:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                Objects[selected_index].angle_z+=10;
            }
            break;
        }
            
        case GLFW_KEY_Y:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                Objects[selected_index].angle_z-=10;
            }
            break;
        }
            
            
        case GLFW_KEY_U:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                Objects[selected_index].angle_x+=10;
            }
            break;
        }
            
        case GLFW_KEY_I:
        {
            if (action == GLFW_PRESS && selected_index >= 0){
                Objects[selected_index].angle_x-=10;
            }
            break;
        }
            
            
        case GLFW_KEY_P:
        {
            if (action == GLFW_PRESS && selected_index >= 0){
                Objects[selected_index].angle_y+=10;
            }
            break;
        }
        
        case GLFW_KEY_O:
        {
            if (action == GLFW_PRESS && selected_index >= 0){
                Objects[selected_index].angle_y-=10;
            }
            break;
        }
        
        case GLFW_KEY_W:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                Objects[selected_index].y+=0.2;
            }
            break;
        }
        
        case GLFW_KEY_S:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                    Objects[selected_index].y-=0.2;
            }
            break;
        }
        
        case GLFW_KEY_A:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                
                Objects[selected_index].x-=0.2;
                
            }
            break;
        }
        case GLFW_KEY_D:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                 Objects[selected_index].x+=0.2;
            }
            break;
        }
        
        case GLFW_KEY_Q:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                Objects[selected_index].z+=0.2;
            }
            break;
        }
        
        case GLFW_KEY_E:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                Objects[selected_index].z-=0.2;
            }
            break;
        }
        
        case GLFW_KEY_F:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                TriMesh obj = Objects[selected_index];
                int start = obj.start;
                int length = obj.tri_num*3;
                
                MatrixXf V_normal_tmp(3,V_normal.cols()-length);
                MatrixXf V_tmp(3,V.cols()-length);
    
                V_normal_tmp<<V_normal.block(0,0,3,start),
                V_normal.block(0,start+length,3,V_normal.cols()-(start+length));
                V_normal = V_normal_tmp;
                
                V_tmp<<V.block(0,0,3,start),
                V.block(0,start+length,3, V.cols()-(start+length));
                V = V_tmp;
                
                for (int i=selected_index+1;i<Objects.size();i++){
                    Objects[i].start-=Objects[selected_index].tri_num*3;
                }
                                              
                Objects.erase (Objects.begin()+selected_index);
                VBO.update(V);
                VBO_N.update(V_normal);
                selected_index =-1;
            }
            break;
        }
            
        case GLFW_KEY_R:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                selected_index=-1;
            }
            break;
        }
            
        case GLFW_KEY_Z:
        {
            if (action == GLFW_PRESS){
                if (ball==0){
                    cam_pos[0] -= 0.2;
                } else {
                    ball_move(0.1, 0, 0);
                }
            }
            break;
        }
        
        case GLFW_KEY_X:
        {
            if (action == GLFW_PRESS){
                if (ball ==0){
                    cam_pos[0] += 0.2;
                } else {
                    ball_move(-0.1, 0, 0);
                }
            }
            break;
        }
            
        case GLFW_KEY_C:
        {
            if (action == GLFW_PRESS){
                if (ball ==0){
                    cam_pos[1] -= 0.2;
                } else {
                    ball_move(0, 10, 0);
                }
            }
            break;
        }
        
        case GLFW_KEY_V:
        {
            if (action == GLFW_PRESS){
                if (ball ==0){
                    cam_pos[1] += 0.2;
                } else {
                    ball_move(0, -10, 0);
                }
            }
            break;
        }
        
        case GLFW_KEY_B:
        {
            if (action == GLFW_PRESS){
                if (ball ==0){
                   cam_pos[2] += 0.2;
                } else {
                    ball_move(0, 0, 10);
                }
            }
            break;
        }
        
        case GLFW_KEY_N:
        {
            if (action == GLFW_PRESS){
                if (ball ==0){
                    cam_pos[2] -= 0.2;
                } else {
                     ball_move(0, 0, -10);
                }
            }
            break;
        }
            
        default:
        {
            break;
        }
        //VBO.update(V);
        //VBO_C.update(V_color);
        //VBO_N.update(V_normal);
    }
}

int main(void)
{
    GLFWwindow* window;
    
    // Initialize the library
    if (!glfwInit())
        return -1;
    
    // Activate supersampling
    glfwWindowHint(GLFW_SAMPLES, 8);

    // Ensure that we get at least a 3.2 context
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);

    // On apple we have to load a core profile with forward compatibility
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // Create a windowed mode window and its OpenGL context
    window = glfwCreateWindow(screenWidth,screenHeight, "Hello World", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);
    
    #ifndef __APPLE__
      glewExperimental = true;
      GLenum err = glewInit();
      if(GLEW_OK != err)
      {
        /* Problem: glewInit failed, something is seriously wrong. */
       fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
      }
      glGetError(); // pull and savely ignonre unhandled errors like GL_INVALID_ENUM
      fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
    #endif
    
    //glDepthFunc(GL_LESS);
    int major, minor, rev;
    major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
    minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
    rev = glfwGetWindowAttrib(window, GLFW_CONTEXT_REVISION);
    printf("OpenGL version recieved: %d.%d.%d\n", major, minor, rev);
    printf("Supported OpenGL is %s\n", (const char*)glGetString(GL_VERSION));
    printf("Supported GLSL is %s\n", (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));
    
    // Initialize the VAO
    // A Vertex Array Object (or VAO) is an object that describes how the vertex
    // attributes are stored in a Vertex Buffer Object (or VBO). This means that
    // the VAO is not the actual object storing the vertex data,
    // but the descriptor of the vertex data.
    //glEnable(GL_DEPTH_TEST);
    VertexArrayObject VAO(0);
    VAO.init();
    VAO.bind();

    // Initialize the VBO with the vertices data
    // A VBO is a data container that lives in the GPU memory
    VBO.init();
    VBO_C.init();
    VBO_N.init();
    
    load_meshes();
    
    V.resize(3,3);
    V <<0,0,0,
        0,0,0,
        0,0,0;
    VBO.update(V);
    
    V_color.resize(3,3);
    
    V_color <<0.0f,0.0f,0.0f,
              0.0f,0.0f,0.0f,
              1.0f,1.0f,1.0f;
    VBO_C.update(V_color);

    V_normal.resize(3,3);
    V_normal<<0,0,0,
              0,0,0,
              0,0,0;
    VBO_N.update(V_normal);
    
    V_trans<<1,0,0,0,
             0,1,0,0,
             0,0,1,0,
             0,0,0,1;
    
    // Initialize the OpenGL Program
    // A program controls the OpenGL pipeline and it must contains
    // at least a vertex shader and a fragment shader to be valid
    Program program;
    const GLchar* vertex_shader =
            "#version 150 core\n"
                    "in vec3 position;"
                    "in vec3 normal_local;"

                    "out vec3 normal;"
                    "out vec3 frag_pos;"
    
                    "uniform mat4 scope;"
                    "uniform mat4 trans;"
                    "uniform mat4 view;"
                    "uniform mat4 proj;"
    
                    "void main()"
                    "{"
                    "    frag_pos = vec3(trans * vec4(position, 1.0));"
                    "    gl_Position = scope*proj*view*trans*vec4(position, 1.0);"
                    "    normal = mat3(transpose(inverse(trans))) * normal_local;"
                    "}";
    
    const GLchar* fragment_shader =
        "#version 150 core\n"
        "in vec3 normal;"
        "in vec3 frag_pos;"
        "in vec3 f_color;"

        "out vec4 outColor;"
        "uniform vec3 mesh_color;"
        "uniform vec3 light_color;"
        "uniform vec3 light_pos;"
        "uniform vec3 cam_pos;"
    
        "void main()"
        "{"
    
        "    float ambient_weight =  0.2;"
        "    float specular_weight = 0.9;"
        "    float diff_weight =     1.1;"
        "    vec3 ambient_cont = ambient_weight * light_color;"
    
        "    vec3 normed_normal = normalize(normal);"
        "    vec3 light_dir = normalize(light_pos - frag_pos);"
        "    float diffuse = max(dot(normed_normal, light_dir), 0.0);"
        "    vec3 diffuse_cont = diff_weight * diffuse * light_color;"
    
        "    vec3 sight_dir = normalize(cam_pos - frag_pos);"
        "    vec3 ref_dir   = reflect(-light_dir, normed_normal);"
        "    float specular = pow(max(dot(sight_dir, ref_dir), 0.0), 123);"
        "    vec3 specular_cont = specular_weight * specular * light_color;"
    
        "    outColor = vec4(mesh_color * (ambient_cont + diffuse_cont + specular_cont), 1.0);"
        "}";
    
    // Compile the two shaders and upload the binary to the GPU
    // Note that we have to explicitly specify that the output "slot" called outColor
    // is the one that we want in the fragment buffer (and thus on screen)
    program.init(vertex_shader,fragment_shader,"outColor");
    program.bind();
    
    // The vertex shader wants the position of the vertices as an input.
    // The following line connects the VBO we defined above with the position "slot"
    // in the vertex shader
    program.bindVertexAttribArray("position",VBO);
    program.bindVertexAttribArray("normal_local",VBO_N);
    
    glUniform3f(program.uniform("light_color"),
                light_color(0),light_color(1),light_color(2));
    glUniform3f(program.uniform("light_pos"),
                light_pos(0),light_pos(1),light_pos(2));
    
    // Register the keyboard callback
    glfwSetKeyCallback(window, key_callback);

    // Register the mouse callback
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    // Update viewport
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    
    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window))
    {
        // Set the uniform value depending on the time difference
        // Clear the framebuffer
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
        //glClear(GL_COLOR_BUFFER_BIT);
        
        // Enable depth test
        glEnable(GL_DEPTH_TEST);
        // Bind your program
        VAO.bind();
        program.bind();
        
        //int width, height;
        glfwGetWindowSize(window, &screenWidth, &screenHeight);
        float aspect_ratio = float(screenHeight)/float(screenWidth); // corresponds to the necessary width scaling
        
        V_scope = scope_transform(aspect_ratio);
        glUniformMatrix4fv(program.uniform("scope"), 1, GL_FALSE, V_scope.data());
        
        Vector3d center(0.0f,0.0f,0.0f);
        V_view = view_transform(cam_pos,center);
        glUniformMatrix4fv(program.uniform("view"), 1, GL_FALSE,  V_view.data());
        //cout<<"view: "<<V_view<<endl;
        V_proj = projection_transform(camera_pers_type, 60, screenWidth,screenHeight, 0.1, 100);
        glUniformMatrix4fv(program.uniform("proj"), 1, GL_FALSE,  V_proj.data());
        
        glUniform3f(program.uniform("cam_pos"),
                    cam_pos(0),cam_pos(1),cam_pos(2));
        
        for (int i=0;i<Objects.size();i++){
            TriMesh cur_obj = Objects[i];
            MatrixXf trans = cur_obj.get_trans_mat();
            Vector3d mesh_color = cur_obj.get_color();
            glUniformMatrix4fv(program.uniform("trans"), 1, GL_FALSE,trans.data());
            
            if(i==selected_index){
                glUniform3f(program.uniform("mesh_color"),0.7,0.2,0.2);
            } else{
                glUniform3f(program.uniform("mesh_color"),
                    mesh_color(0),mesh_color(1),mesh_color(2));
            }
            
            if(cur_obj.render_type==0){
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                
                //glUniform3f(program.uniform("mesh_color"),0.1,0.1,0.1);
                glDrawArrays(GL_TRIANGLES,cur_obj.start,cur_obj.tri_num*3);
            }
            
            else if(cur_obj.render_type==1){
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                glDrawArrays(GL_TRIANGLES,cur_obj.start,cur_obj.tri_num*3);
                glUniform3f(program.uniform("mesh_color"),0.1,0.1,0.1);
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                glDrawArrays(GL_TRIANGLES,cur_obj.start,cur_obj.tri_num*3);
            }
            
            else if(cur_obj.render_type==2){
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                glDrawArrays(GL_TRIANGLES,cur_obj.start,cur_obj.tri_num*3);
            }
            
            else{
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                glDrawArrays(GL_TRIANGLES, cur_obj.start, cur_obj.tri_num*3);
            }
            //cout<<"tri:"<<Objects[i].start;
        }
        // Swap front and back buffers
        glfwSwapBuffers(window);
        // Poll for and process events
        glfwPollEvents();
    }
    
    // Deallocate opengl memory
    program.free();
    
    VAO.free();
    VBO.free();
    VBO_C.free();
    VBO_N.free();
    
    // Deallocate glfw internals
    glfwTerminate();
    return 0;
}
