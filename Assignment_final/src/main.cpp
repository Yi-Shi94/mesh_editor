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
# define pi 3.1415926

using namespace std;
using namespace Eigen;

// VertexBufferObject wrapper
VertexBufferObject VBO;
VertexBufferObject VBO_N;
VertexBufferObject VBO_T;

// Contains the vertex positions
vector<TriMesh> Objects;

MatrixXf V(3,3);
MatrixXf V_color(3,3);
MatrixXf V_normal(3,3);
MatrixXf V_T(2,3);

MatrixXf V_trans(4,4);
Matrix4f V_view(4,4);

Matrix4f V_proj(4,4);
Matrix4f V_scope(4,4);

int selected_index = -1;
int ball_flag = 0;
int flag_theta = 0;
int flag_phi = 0;
/*CONFIGURE*/

//int camera_pers_type = 0; //0,1 //perspective,ortho
int camera_pers_type = 0 ;
int screenWidth = 640;
int screenHeight = 480;
int ball = 0;
int mesh_edit = 0;
float cam_mesh_range = 1.5;
Vector3d light_pos(0.,0.,5.);
Vector3d light_color(1.,1.,1.);

Vector3d center(0.0f,0.0f,0.0f);
Vector3d cam_pos(0.,0.,5.);
Vector3d cam_pos_store(0.,0.,5.);
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


void update_inidicator_cube(){
    MatrixXf V_tmp = Objects[selected_index].populate_vertex();
    MatrixXf V_ver_show(V.rows(),V.cols()+V_tmp.cols());
    
    V_ver_show<<V,V_tmp;
    VBO.update(V_ver_show);
    
    MatrixXf V_normal_tmp = Objects[selected_index].populate_normal_vertex();
    MatrixXf V_normal_show(V_normal.rows(),
                           V_normal.cols()+V_normal_tmp.cols());
    
    V_normal_show<<V_normal,V_normal_tmp;
    VBO_N.update(V_normal_show);
    
    if(V_tmp.cols()>V_T.cols()){
        MatrixXf V_T_tmp(V_T.rows(),V_tmp.cols());
        for (int i=0;i<V_tmp.cols();i=i+3){
            V_T_tmp.col(i)<<V_T.col(0);
            V_T_tmp.col(i+1)<<V_T.col(1);
            V_T_tmp.col(i+2)<<V_T.col(2);
        }
        VBO_T.update(V_T_tmp);
    }
    
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
        Vector4f ray_cl(p_cl(0),p_cl(1), -1., 1.);
        Vector4f ray_view= V_proj.inverse() * ray_cl;
        Vector4f ray_view_d(ray_view(0),ray_view(1),-1.,0.);
        Vector4f ray_world_3d = V_view.inverse() * ray_view_d;
        //[u v pp w] =  V_scope * V_proj * V_view * V_model* [x y z 1] T
        Vector3d ray_wor_3d;
        Vector3d ray_origin_wor_3d;
        if (camera_pers_type==0)
        {
            ray_origin_wor_3d = cam_pos;
            ray_wor_3d = Vector3d(ray_world_3d(0),ray_world_3d(1),ray_world_3d(2));
        } else {
            Vector4f ray_dir(0., 0., -1.,0.);
            ray_dir = (V_view.inverse() * ray_dir).normalized();
            Vector4f ray_origin(ray_view(0),ray_view(1),2,0);
            Vector4f ray_origin_wor = V_view.inverse()*ray_origin;
            ray_origin_wor_3d = Vector3d(ray_origin_wor(0),ray_origin_wor(1),ray_origin_wor(2));
            ray_wor_3d = Vector3d(ray_dir(0),ray_dir(1),ray_dir(2));

        }
        //cout<<ray_wor_3d<<endl;
        double min_dis = 100000;
        int min_index = -1;
        for (int i=0;i<Objects.size();i++){
            double t = Objects[i].is_hit(ray_origin_wor_3d, ray_wor_3d);
            if (t<0) continue;
            else {
                if (min_dis > t) {
                    min_dis = t;
                    min_index = i;
                }
            }
        }
        if (min_index>=0) {
            selected_index = min_index;
            if (mesh_edit>0){
                TriMesh cur_obj =Objects[selected_index];
                //Objects[selected_index].select_vertex(ray_origin_wor_3d, ray_wor_3d, min_dis);
                
                Objects[selected_index].select_vertex_2(p_cl(0),p_cl(1),V_proj * V_view);
                update_inidicator_cube();
            }
        }
    }
}

Vector3d fromCart2Sphere(float x_c, float y_c, float z_c){
    float x = z_c;
    float y = x_c;
    float z = y_c;
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
    Vector3d cart(y,z,x);
    return cart;
}

void ball_move(float dr,float dtheta,float dphi){
    
    float x = cam_pos(0)-center(0);
    float y = cam_pos(1)-center(1);
    float z = cam_pos(2)-center(2);
    
    Vector3d spher = fromCart2Sphere(x,y,z);
    
    float r = dr+spher(0);
    //float addon = dtheta*pi/180;
    float theta;
    if (flag_theta==0){
        theta = spher(1)/pi*180+dtheta;
    } else {
        theta = spher(1)/pi*180-dtheta;
    }
    float phi;
    
    phi= spher(2)/pi*180+dphi;
    
    
    if (r<0.8) r=0;
    else if (r>5) r=5;
    
    if (theta>180) {
        theta = 180-(theta-180);
        phi += 180;
        phi = int(phi) % 180;
        if(flag_theta==0)
           flag_theta = 1;
        else flag_theta = 0;
    } else if (theta<0){
        theta = -theta;
        phi += 180;
        phi = int(phi) % 180;
        if(flag_theta==0)
            flag_theta = 1;
        else flag_theta = 0;
    }
    
    if (phi>=180) {
        phi = -180+(phi-180);
        
    } else if (phi<-180){
        phi = 180-(phi+180);
    }
    
    cout<<"x: "<<x<<" y: "<<y<<" z: "<<z<<endl;
    cout<<"theta_before: "<<spher(1)/pi*180<<" phi: "<<spher(2)/pi*180<<endl;
    cout<<"add_theta: "<<dtheta<<" add_phi: "<<dphi<<endl;
    cout<<"theta: "<<theta<<" phi: "<<phi<<endl;
    
    theta = theta/180*pi;
    phi = phi/180*pi;
    
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
    V_normal_show<<V_normal,V_normal_tmp;
    V_normal = V_normal_show;
    
    MatrixXf V_T_tmp(V_T.rows(),V.cols());
    for (int i=0;i<V.cols();i=i+3){
        V_T_tmp.col(i)<<V_T.col(0);
        V_T_tmp.col(i+1)<<V_T.col(1);
        V_T_tmp.col(i+2)<<V_T.col(2);
    }
    V_T = V_T_tmp;
}

void delete_mesh(){
    if(selected_index<0) return;
    TriMesh obj = Objects[selected_index];
    int start = obj.start;
    int length = obj.tri_num*3;
    
    MatrixXf V_normal_tmp(3,V_normal.cols()-length);
    MatrixXf V_tmp(3,V.cols()-length);
    MatrixXf V_T_tmp(V_T.rows(),V_T.cols()-length);
    
    V_normal_tmp<<V_normal.block(0,0,3,start),
        V_normal.block(0,start+length,3,V_normal.cols()-(start+length));
    V_normal = V_normal_tmp;
    
    V_tmp<<V.block(0,0,3,start),
        V.block(0,start+length,3, V.cols()-(start+length));
    V = V_tmp;
    
    V_T_tmp<<V_T.block(0,0,2,start),
        V_T.block(0,start+length,2, V_T.cols()-(start+length));
    V_T = V_T_tmp;
    
    for (int i=selected_index+1;i<Objects.size();i++){
        Objects[i].start-=Objects[selected_index].tri_num*3;
    }
    
    Objects.erase (Objects.begin()+selected_index);
    VBO.update(V);
    VBO_N.update(V_normal);
    VBO_T.update(V_T);
    selected_index =-1;
    return;
}

void update_mesh(){
    int start = Objects[selected_index].start;
    MatrixXf V_normal_tmp = V_normal.block(0,0,3,start);
    MatrixXf V_tmp = V.block(0,0,3,start);
    
    for(int i=selected_index;i<Objects.size();i++){
        Objects[i].start = V_tmp.cols();
        Objects[i].tri_num = Objects[i].faces.size();
        MatrixXf V_tmp_store= V_tmp;
        MatrixXf V_normal_tmp_store= V_normal_tmp;
        
        MatrixXf V_cur = Objects[i].get_matrix();
        MatrixXf V_normal_cur = Objects[i].get_normal_matrix();
        
        V_tmp.resize(3,V_tmp.cols()+V_cur.cols());
        V_tmp<<V_tmp_store,V_cur;
        
        V_normal_tmp.resize(3, V_normal_tmp.cols()+V_normal_cur.cols());
        V_normal_tmp<<V_normal_tmp_store,V_normal_cur;
        
    }
    
    V = V_tmp;
    V_normal = V_normal_tmp;
    
    MatrixXf V_T_tmp(2,V.cols());
    for (int i=0;i<V.cols();i=i+3){
        V_T_tmp.col(i)<<V_T.col(0);
        V_T_tmp.col(i+1)<<V_T.col(1);
        V_T_tmp.col(i+2)<<V_T.col(2);
    }
    V_T = V_T_tmp;
    
    VBO.update(V);
    VBO_N.update(V_normal);
    VBO_T.update(V_T);
    //cout<<V.cols()<<"after deleting"<<Objects[selected_index].selected_vertex_index.size()<<endl;
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    // Update the position of the first vertex if the keys 1,2, or 3 are pressed
    switch (key)
    {
        case GLFW_KEY_1:
        {
            /*
            if (action == GLFW_PRESS)
            {
                if (mesh_edit>0) break;
                TriMesh obj(0,V.cols());
                Objects.push_back(obj);
                append_mesh(obj);
                //cout << obj.get_trans_mat();
                cout <<"Total num obj: "<<Objects.size()<<" Verices: "<<V.cols()<<endl;
                VBO.update(V);
                VBO_T.update(V_T);
                VBO_N.update(V_normal);
            }
             */
            break;
        }
        case GLFW_KEY_2:
        {
            if (action == GLFW_PRESS)
            {
                if (mesh_edit>0) break;
                TriMesh obj1(1,V.cols());
                Objects.push_back(obj1);
                append_mesh(obj1);
                cout <<"Total num obj: "<<Objects.size()<<endl;
                VBO.update(V);
                VBO_T.update(V_T);
                VBO_N.update(V_normal);
            }
            break;
        }
            
        case GLFW_KEY_3:
        {
            if (action == GLFW_PRESS)
            {
                if (mesh_edit>0) break;
                TriMesh obj2(2,V.cols());
                Objects.push_back(obj2);
                append_mesh(obj2);
                cout <<"Total num obj: "<<Objects.size()<<endl;
                VBO.update(V);
                VBO_T.update(V_T);
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
        
        case GLFW_KEY_8:
        {
            if (action == GLFW_PRESS){
                if (selected_index >=0 && mesh_edit==0){
                    mesh_edit = 1;
                    ball=1;
                    cam_pos_store = cam_pos;
                    TriMesh o_cur = Objects[selected_index];
                    
                    cam_pos(0) = o_cur.x;
                    cam_pos(1) = o_cur.y;
                    cam_pos(2) = o_cur.z+cam_mesh_range;
                    
                    center(0) = o_cur.x;
                    center(1) = o_cur.y;
                    center(2) = o_cur.z;
                    
                } else if(selected_index >=0 && mesh_edit==1){
                    Objects[selected_index].save_off();
                    mesh_edit = 0;
                    ball = 0;
                    cam_pos = cam_pos_store;
                    center(0) = 0.0;
                    center(1) = 0.0;
                    center(2) = 0.0;
                }
            }
            break;
        }
            
        case GLFW_KEY_9:
        {
            if (action == GLFW_PRESS){
                if (mesh_edit>0) break;
                if (ball==0) {
                    ball=1;
                    if(selected_index>0){
                        TriMesh o_cur = Objects[selected_index];
                        
                        cam_pos(0) = o_cur.x;
                        cam_pos(1) = o_cur.y;
                        cam_pos(2) = o_cur.z+cam_mesh_range;
                        
                        center(0) = o_cur.x;
                        center(1) = o_cur.y;
                        center(2) = o_cur.z;
                    } else {
                        cam_pos(0) = 0;
                        cam_pos(1) = 0;
                        cam_pos(2) = cam_mesh_range;
                        
                        center(0) = 0;
                        center(1) = 0;
                        center(2) = 0;
                    }
                    
                }
                else ball= 0;
            }
            break;
        }
            
        
        case GLFW_KEY_0:
        {
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
                if(mesh_edit>0){
                    break;
                }
                Objects[selected_index].angle_z+=10;
            }
            break;
        }
            
        case GLFW_KEY_Y:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                if(mesh_edit>0){
                    break;
                }
                Objects[selected_index].angle_z-=10;
            }
            break;
        }
            
            
        case GLFW_KEY_U:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                if(mesh_edit>0){
                    break;
                }
                Objects[selected_index].angle_x+=10;
            }
            break;
        }
            
        case GLFW_KEY_I:
        {
            if (action == GLFW_PRESS && selected_index >= 0){
                if(mesh_edit>0){
                    break;
                }
                Objects[selected_index].angle_x-=10;
            }
            break;
        }
            
        case GLFW_KEY_P:
        {
            if (action == GLFW_PRESS && selected_index >= 0){
                if(mesh_edit>0){
                    break;
                }
                Objects[selected_index].angle_y+=10;
            }
            break;
        }
        
        case GLFW_KEY_O:
        {
            if (action == GLFW_PRESS && selected_index >= 0){
                if(mesh_edit>0){
                    break;
                }
                Objects[selected_index].angle_y-=10;
            }
            break;
        }
        
        case GLFW_KEY_W:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                Objects[selected_index].y+=0.2;
                if (mesh_edit>0){
                    cam_pos(0) = Objects[selected_index].x;
                    cam_pos(1) = Objects[selected_index].y;
                    cam_pos(2) = Objects[selected_index].z+cam_mesh_range;
                    center(0) = Objects[selected_index].x;
                    center(1) = Objects[selected_index].y;
                    center(2) = Objects[selected_index].z;
                }
            }
            break;
        }
        
        case GLFW_KEY_S:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                Objects[selected_index].y-=0.2;
                if (mesh_edit>0){
                    cam_pos(0) = Objects[selected_index].x;
                    cam_pos(1) = Objects[selected_index].y;
                    cam_pos(2) = Objects[selected_index].z+cam_mesh_range;
                    center(0) = Objects[selected_index].x;
                    center(1) = Objects[selected_index].y;
                    center(2) = Objects[selected_index].z;
                }
            }
            break;
        }
        
        case GLFW_KEY_A:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                
                Objects[selected_index].x-=0.2;
                if (mesh_edit>0){
                    cam_pos(0) = Objects[selected_index].x;
                    cam_pos(1) = Objects[selected_index].y;
                    cam_pos(2) = Objects[selected_index].z+cam_mesh_range;
                    center(0) = Objects[selected_index].x;
                    center(1) = Objects[selected_index].y;
                    center(2) = Objects[selected_index].z;
                }
                
            }
            break;
        }
        case GLFW_KEY_D:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                Objects[selected_index].x+=0.2;
                if (mesh_edit>0){
                    cam_pos(0) = Objects[selected_index].x;
                    cam_pos(1) = Objects[selected_index].y;
                    cam_pos(2) = Objects[selected_index].z+cam_mesh_range;
                    center(0) = Objects[selected_index].x;
                    center(1) = Objects[selected_index].y;
                    center(2) = Objects[selected_index].z;
                }
            }
            break;
        }
        
        case GLFW_KEY_Q:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                Objects[selected_index].z+=0.2;
                if (mesh_edit>0){
                    cam_pos(0) = Objects[selected_index].x;
                    cam_pos(1) = Objects[selected_index].y;
                    cam_pos(2) = Objects[selected_index].z+cam_mesh_range;
                    center(0) = Objects[selected_index].x;
                    center(1) = Objects[selected_index].y;
                    center(2) = Objects[selected_index].z;
                }
            }
            break;
        }
        
        case GLFW_KEY_E:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                Objects[selected_index].z-=0.2;
                if (mesh_edit>0){
                    cam_pos(0) = Objects[selected_index].x;
                    cam_pos(1) = Objects[selected_index].y;
                    cam_pos(2) = Objects[selected_index].z+cam_mesh_range;
                    center(0) = Objects[selected_index].x;
                    center(1) = Objects[selected_index].y;
                    center(2) = Objects[selected_index].z;
                }
            }
            break;
        }
        
        case GLFW_KEY_F:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                if (mesh_edit>0){
                    
                    if(Objects[selected_index].delete_vertex()){
                        cout<<'t'<<Objects[selected_index].tri_num<<endl;
                        update_mesh();
                    }
                } else{
                    delete_mesh();
                }
            }
            break;
        }
            
        case GLFW_KEY_R:
        {
            if (action == GLFW_PRESS && selected_index>=0){
                if (mesh_edit>0){
                    if (Objects[selected_index].selected_vertex_index.size()>0){
                        Objects[selected_index].selected_vertex_index.pop_back();
                        update_inidicator_cube();
                    }
                    
                } else{
                    selected_index=-1;
                }
            }
            break;
        }
            
        case GLFW_KEY_Z:
        {
            if (action == GLFW_PRESS){
                if (ball==0){
                    cam_pos[0] -= 0.2;
                } else {
                    ball_move(-0.05,0,0);
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
                    ball_move(0.05,0,0);
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
                    ball_move(0,3,0);
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
                    ball_move(0,-3,0);
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
                    ball_move(0,0,3);
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
                    ball_move(0,0,-3);
                }
            }
            break;
        }
            
        case GLFW_KEY_J:
        {
            if (action == GLFW_PRESS){
                if(mesh_edit>0){
                    if(Objects[selected_index].selected_vertex_index.size()!=2){
                        break;
                    } else{
                        if(Objects[selected_index].merge_vertex()){
                            update_mesh();
                            update_inidicator_cube();
                        }
                    }
                }
            }
            break;
        }
            
        case GLFW_KEY_G:
        {
            if (action == GLFW_PRESS){
                if(mesh_edit>0){
                    if(Objects[selected_index].selected_vertex_index.size()!=1){
                        break;
                    } else{
                        if(Objects[selected_index].translate_vertex(2, 0.004)){
                            update_mesh();
                            update_inidicator_cube();
                        }
                    }
                }
            }
            break;
        }
            
        case GLFW_KEY_H:
        {
            if(mesh_edit>0){
                if(Objects[selected_index].selected_vertex_index.size()!=1){
                    break;
                } else{
                    if(Objects[selected_index].translate_vertex(2, -0.004)){
                        update_mesh();
                        update_inidicator_cube();
                    }
                }
            }
            break;
        }
            
        case GLFW_KEY_UP:
        {
            if (action == GLFW_PRESS){
                if(Objects[selected_index].selected_vertex_index.size()!=1){
                    break;
                } else{
                    if(Objects[selected_index].translate_vertex(1, 0.004)){
                        update_mesh();
                        update_inidicator_cube();
                    }
                }
            }
            break;
        }
            
        case GLFW_KEY_DOWN:
        {
            if (action == GLFW_PRESS){
                if(Objects[selected_index].selected_vertex_index.size()!=1){
                    break;
                } else{
                    if(Objects[selected_index].translate_vertex(1, -0.004)){
                        update_mesh();
                        update_inidicator_cube();
                    }
                }
            }
            break;
        }
        
        case GLFW_KEY_LEFT:
        {
            if (action == GLFW_PRESS){
                if(Objects[selected_index].selected_vertex_index.size()!=1){
                    break;
                } else{
                    if(Objects[selected_index].translate_vertex(0, 0.004)){
                        update_mesh();
                        update_inidicator_cube();
                    }
                }
            }
            break;
        }
        
        case GLFW_KEY_RIGHT:
        {
            if (action == GLFW_PRESS){
                if(Objects[selected_index].selected_vertex_index.size()!=1){
                    break;
                } else{
                    if(Objects[selected_index].translate_vertex(0, -0.004)){
                        update_mesh();
                        update_inidicator_cube();
                    }
                }
            }
            break;
        }
            
        default:
        {
            break;
        }
    }
}

GLuint loadBMP_custom(const char * imagepath)
{
    // Data read from the header of the BMP file
    // Each BMP file begins by a 54-bytes header
    unsigned char header[54];
    // Position in the file where the actual data begins
    unsigned int dataPos;
    unsigned int width, height;
    unsigned int imageSize; // = width*height*3
    // Actual RGB data
    unsigned char * data;
    
    // Open the file
    printf("filename: %s\n",imagepath);
    FILE * file = fopen(imagepath,"rb");
    if (!file)
    {
        printf("Image could not be opened\n");
        return 0;
    }
    
    if ( fread(header, 1, 54, file)!=54 )
    {
        // If not 54 bytes read : problem
        printf("Not a correct BMP file\n");
        return false;
    }
    
    if ( header[0]!='B' || header[1]!='M' ){
        printf("Not a correct BMP file\n");
        return 0;
    }
    
    // Read ints from the byte array
    dataPos    = *(int*)&(header[0x0A]);
    imageSize  = *(int*)&(header[0x22]);
    width      = *(int*)&(header[0x12]);
    height     = *(int*)&(header[0x16]);
    
    // Create a buffer
    data = new unsigned char [imageSize];
    
    // Read the actual data from the file into the buffer
    fread(data,1,imageSize,file);
    
    //Everything is in memory now, the file can be closed
    fclose(file);
    
    // Create one OpenGL texture
    GLuint textureID;
    glGenTextures(1, &textureID);
    
    // "Bind" the newly created texture:
    // all future texture functions will modify this texture
    glBindTexture(GL_TEXTURE_2D, textureID);
    
    // Give the image to OpenGL
    glTexImage2D(GL_TEXTURE_2D, 0,GL_RGB, width, height, 0, \
                 GL_BGR, GL_UNSIGNED_BYTE, data);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    
    cout<<"texture: "<<imageSize<<"width: "<<width<<"height: "<<height<<endl;
    return textureID;
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
    
   
    int major, minor, rev;
    major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
    minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
    rev = glfwGetWindowAttrib(window, GLFW_CONTEXT_REVISION);
    printf("OpenGL version recieved: %d.%d.%d\n", major, minor, rev);
    printf("Supported OpenGL is %s\n", (const char*)glGetString(GL_VERSION));
    printf("Supported GLSL is %s\n", (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));
    
    glEnable(GL_TEXTURE_2D);
    GLuint textureID = loadBMP_custom("../data/grass.bmp");
    if (!textureID){
        return 0;
    }
    GLfloat vertices_texture[] = {
        0.05f, 0.05f,
        0.95f, 0.05f,
        0.5f, 0.95f,
        
    };
    
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
    VBO_N.init();
    VBO_T.init();
    
    load_meshes();
    
    V.resize(3,3);
    V <<0,0,0,
        0,0,0,
        0,0,0;
    VBO.update(V);
    
    V_T.resize(2,3);

    V_T <<0.1f,0.1f,0.9f,
          0.1f,0.9f,0.9f,
    
    VBO_T.update(V_T);

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
                    "in vec2 texCoord;"
    
                    "out vec2 TexCoord;"
                    "out vec3 normal;"
                    "out vec3 frag_pos;"
    
                    "uniform mat4 scope;"
                    "uniform mat4 trans;"
                    "uniform mat4 view;"
                    "uniform mat4 proj;"
    
                    "void main()"
                    "{"
                    "    frag_pos = vec3(trans * vec4(position, 1.0));"
                    "    gl_Position = scope*proj*view*trans*vec4(position,1.0);"
                    "    TexCoord = texCoord;"
                    "    normal = mat3(transpose(inverse(trans))) * normal_local;"
                    "}";
    
    const GLchar* fragment_shader =
        "#version 150 core\n"
        "in vec3 normal;"
        "in vec3 frag_pos;"
        "in vec2 TexCoord;"
    
        "uniform sampler2D ourTexture;"
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
        "    vec3 ambient_cont = ambient_weight * light_color * light_color;"
    
        "    vec3 normed_normal = normalize(normal);"
        "    vec3 light_dir = normalize(light_pos - frag_pos);"
        "    float diffuse = max(dot(normed_normal, light_dir), 0.0);"
        "    vec3 diffuse_cont = diff_weight * diffuse * light_color * light_color;"
    
        "    vec3 sight_dir = normalize(cam_pos - frag_pos);"
        "    vec3 ref_dir   = reflect(-light_dir, normed_normal);"
        "    float specular = pow(max(dot(sight_dir, ref_dir), 0.0), 123);"
        "    vec3 specular_cont = specular_weight * specular * light_color * light_color;"
        //"    vec4 texture_color = texture(ourTexture, TexCoord); "
        "    vec4 texture_color = vec4(1.0,1.0,1.0,1.0); "
        "    outColor = texture_color * vec4(mesh_color * (ambient_cont + diffuse_cont + specular_cont), 1.0);"
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
    program.bindVertexAttribArray("texCoord",VBO_T);
    
    glUniformMatrix3fv(program.uniform("scope"), 1, GL_FALSE, V_scope.data());
    
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
        
        glUniform3f(program.uniform("light_color"),
                    light_color(0),light_color(1),light_color(2));
        
        //if (mesh_edit>0){
        //    glUniform3f(program.uniform("light_pos"),
        //            cam_pos(0),cam_pos(1),cam_pos(2));
        // } else{
        
        glUniform3f(program.uniform("light_pos"),
                    cam_pos(0),cam_pos(1),cam_pos(2));
        // }
        
        glUniform3f(program.uniform("cam_pos"),
                    cam_pos(0),cam_pos(1),cam_pos(2));
        
        V_view = view_transform(cam_pos,center);
        glUniformMatrix4fv(program.uniform("view"), 1, GL_FALSE,  V_view.data());

        V_proj = projection_transform(camera_pers_type, 60, screenWidth,screenHeight, 0.1, 100);
        glUniformMatrix4fv(program.uniform("proj"), 1, GL_FALSE,  V_proj.data());
        
        for (int i=0;i<Objects.size();i++){
            TriMesh cur_obj = Objects[i];
            //cout<<"godam:"<<cur_obj.tri_num<<endl;
            MatrixXf trans = cur_obj.get_trans_mat();
            Vector3d mesh_color = cur_obj.get_color();
            glUniformMatrix4fv(program.uniform("trans"), 1, GL_FALSE,trans.data());
            
            if(i==selected_index){
                glUniform3f(program.uniform("mesh_color"),0.1,0.6,0.1);
                
            } else{
                glUniform3f(program.uniform("mesh_color"),
                    mesh_color(0),mesh_color(1),mesh_color(2));
            }
            
            if(cur_obj.render_type==0){
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                //glUniform3f(program.uniform("mesh_color"),0.1,0.1,0.1);
                glDrawArrays(GL_TRIANGLES,cur_obj.start,cur_obj.tri_num*3);
            }
            
            else if((i==selected_index && mesh_edit>0)||cur_obj.render_type==1){
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
        }
        
        if (mesh_edit>0 && selected_index>=0){
            vector<int> edit_obj = Objects[selected_index].selected_vertex_index;
            //cout<<'q'<<edit_obj.size()<<' '<<endl;
            if (Objects[selected_index].selected_vertex_index.size()>0){
                int start = V.cols();
                int num = 36;
                MatrixXf trans = Objects[selected_index].get_trans_mat();
                glUniform3f(program.uniform("mesh_color"),0.2,0.2,0.2);
                glUniformMatrix4fv(program.uniform("trans"),
                                   1,GL_FALSE,trans.data());
                
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                glDrawArrays(GL_TRIANGLES, start, num*Objects[selected_index].selected_vertex_index.size());
                /*
                for (int i=0;i<edit_obj.size();i++){
                    //cout<<'d'<<edit_obj.size()<<'i'<<i<<endl;
                    //cout<<'s'<<start<<endl;
                    //cout<<V.cols()<<endl;
                    glUniformMatrix4fv(program.uniform("trans"),
                               1,GL_FALSE,trans.data());
                    glUniform3f(program.uniform("mesh_color"),0.2,0.2,0.2);
                    //cout<<start<<' '<<num<<endl;
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    glDrawArrays(GL_TRIANGLES, start, num);
                    start += num;
                }
                 */
            }
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
    VBO_N.free();
    VBO_T.free();
    
    // Deallocate glfw internals
    glfwTerminate();
    return 0;
}
