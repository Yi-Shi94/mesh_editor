// This example is heavily based on the tutorial at https://open.gl

// OpenGL Helpers to reduce the clutter
#include "Helpers.h"
#include <iostream>
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>
#else
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>
#endif
//#include <paint.h>
// Linear Algebra Library
#include <Eigen/Core>
#include <Eigen/Dense>

// Timer
#include <chrono>
using namespace std;
using namespace Eigen;

// VertexBufferObject wrapper
VertexBufferObject VBO;
VertexBufferObject VBO_C;

// Contains the vertex positions
MatrixXf V(2,3);
MatrixXf V_curve(2,4);
MatrixXf V_color(3,3);
Matrix4f V_view(4,4);
Vector3f color_rgb(0.2f, 0.4f, 0.8f);

int selected_tri;
int selected_vertex;

MatrixXf coord_rela_tri(2,3);
double a_lx,a_ly,b_lx,b_ly,c_lx,c_ly;

int MODE=-1;
int NUM_TRIANGLE_CREATED = 0;
int CLICK_COUNTER = 0;

int ani_mv_remain = 0;
double ani_mv_interval =0;
double val_per_mv = 0;
bool ani_dir = true;
auto t_start = chrono::high_resolution_clock::now();

bool is_cur_in_triangle(Vector2f A,Vector2f B,Vector2f C,Vector2f P){
    //if p is in tri abc
    Vector2f v0 = C - A ; Vector2f v1 = B - A ; Vector2f v2 = P - A ;
    float val_00 = v0.dot(v0) ; float val_11 = v1.dot(v1) ;
    float inverDeno = 1 / (val_00 * val_11 -  v0.dot(v1) *  v0.dot(v1)) ;
    float u = (val_11 * v0.dot(v2) -  v0.dot(v1) * v1.dot(v2)) * inverDeno ;
    if (u < 0 || u > 1) // if u out of range, return directly
    {
        return false ;
    }
    float v = (val_00 * v1.dot(v2) - v0.dot(v1) * v0.dot(v2)) * inverDeno ;
    if (v < 0 || v > 1) // if v out of range, return directly
    {
        return false ;
    }
    return u + v <= 1 ;
}


int index_of_closest_vertex(double x, double y){
    Vector2f cur_pos(x,y);
    int cls_vert_index;
    double dist = 1000000.0;
    for (int i=3;i<=3*NUM_TRIANGLE_CREATED+2;i++){
        Vector2f vert(V(0,i),V(1,i));
        double cur_dist = (cur_pos-vert).norm();
        if (dist > cur_dist){
            dist = cur_dist;
            cls_vert_index = i;
        }
    }
    return cls_vert_index;
}

int index_of_selected_triangle(double xworld,double yworld){
    Vector2f cur(xworld,yworld);
    for (int i=1;i<=NUM_TRIANGLE_CREATED;i++){
        Vector2f A(V(0,i*3+0),V(1,i*3+0));
        Vector2f B(V(0,i*3+1),V(1,i*3+1));
        Vector2f C(V(0,i*3+2),V(1,i*3+2));
        if(is_cur_in_triangle(A,B,C,cur)){
            return i;
        }
    }
    return -1;
}

double triangle_area(double Ax,double Bx,double Cx,double Ay,double By,double Cy){
    return (Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))/2;
}

MatrixXf get_bar_center(MatrixXf Vert){
    return (Vert(seq(0,1),0)+Vert(seq(0,1),1)+Vert(seq(0,1),2))/3.0;
}

vector<double> cur_loc_get(GLFWwindow* window)
{
    double xpos, ypos;
    vector<double> coord;
    glfwGetCursorPos(window, &xpos, &ypos);
    
    // Get the size of the window
    int width, height;
    glfwGetWindowSize(window, &width, &height);
    
    Vector4f p_screen(xpos,height-1-ypos,0,1);
    Vector4f p_canonical((p_screen[0]/width)*2-1,(p_screen[1]/height)*2-1,0,1);
    Vector4f p_world = V_view.inverse()*p_canonical;
    // Convert screen position to world coordinates
    coord.push_back(p_world[0]);
    coord.push_back(p_world[1]);

    return coord;
}


void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void set_rotate_tri_intp(bool is_clockwise,double degree, double time_durance, int mv_num){
    t_start = chrono::high_resolution_clock::now();
    ani_mv_interval = time_durance/mv_num;
    val_per_mv = degree/(mv_num*1.0);
    ani_mv_remain = mv_num;
    ani_dir = is_clockwise;
    //auto t_start = chrono::high_resolution_clock::now();
}
    
void rotate_tri(int tri_target,bool is_clockwise,double degree){
    double theta_rad;
    if (is_clockwise) theta_rad = degree/180*3.1415926;
    else theta_rad = -degree/180*3.1415926;
    
    int i = tri_target;
    MatrixXf Vert(3,3);

    Vert<< V(0,i*3),V(0,i*3+1),V(0,i*3+2),
           V(1,i*3),V(1,i*3+1),V(1,i*3+2),
           1,1,1;
    
    MatrixXf bar_center(2,1);
    bar_center << get_bar_center(Vert);
    cout<<"bar: "<<bar_center<<endl;
    
    MatrixXf R(3,3);
    MatrixXf T(3,3);
    MatrixXf Tm(3,3);
    
    R<<cos(theta_rad),sin(theta_rad),0,
      -sin(theta_rad),cos(theta_rad),0,
       0,0,1;
    T<<1,0,bar_center(0,0),
       0,1,bar_center(1,0),
       0,0,1;
    Tm<<1,0,-bar_center(0,0),
        0,1,-bar_center(1,0),
        0,0,1;
    
    MatrixXf Vert_out(3,3);
    Vert_out << T*(R*(Tm*Vert));
    cout<<"out "<<Vert_out;
    V(0,i*3) = Vert_out(0,0); V(0,i*3+1) = Vert_out(0,1);V(0,i*3+2) = Vert_out(0,2);
    V(1,i*3) = Vert_out(1,0); V(1,i*3+1) = Vert_out(1,1);V(1,i*3+2) = Vert_out(1,2);
    
}

void scale_tri(int tri_target,bool is_enlarge){
    double scale;
    if (is_enlarge){
        scale = 1.25;
        cout<<scale<<endl;
    }
    else {
        scale = 0.75;
        cout<<scale<<endl;
    }
    
    int i = tri_target;
    MatrixXf Vert(2,3);
    
    Vert<< V(0,i*3),V(0,i*3+1),V(0,i*3+2),
           V(1,i*3),V(1,i*3+1),V(1,i*3+2);
    
    MatrixXf bar_center(2,1);
    bar_center << get_bar_center(Vert);
    cout<<"bar: "<<bar_center<<endl;
    
    V(0,i*3) = bar_center(0,0) + scale*(Vert(0,0)-bar_center(0,0));
    V(1,i*3) = bar_center(1,0) + scale*(Vert(1,0)-bar_center(1,0));
    V(0,i*3+1) = bar_center(0,0) + scale*(Vert(0,1)-bar_center(0,0));
    V(1,i*3+1) = bar_center(1,0) + scale*(Vert(1,1)-bar_center(1,0));
    V(0,i*3+2) = bar_center(0,0) + scale*(Vert(0,2)-bar_center(0,0));
    V(1,i*3+2) = bar_center(1,0) + scale*(Vert(1,2)-bar_center(1,0));
}


void drawBezier(MatrixXf control_points, double step) {
    double t = 0;
    //MatrixXf points_mat(int(1.0/step),2)
    V_curve.resize(V_curve.rows()+int(1.0/step),2);
    while (t < 1){
        V_curve(int(t/step),0) = pow((1 - t), 3) * control_points(0,0) + 3 * t * pow((1 -t), 2) * control_points(1,0) + 3 * (1-t) * pow(t, 2)* control_points(2,0) + pow (t, 3)* control_points(3,0);
        V_curve(int(t/step),1) = pow((1 - t), 3) * control_points(0,0) + 3 * t * pow((1 -t), 2) * control_points(1,0) + 3 * (1-t) * pow(t, 2)* control_points(2,0) + pow (t, 3)* control_points(3,0);
        t += step;
    }
}

void ini_color_vertex(){
    V_color.resize(V_color.rows(),V.cols());
    cout<<"V: "<<V.rows()<<' '<<V.cols()<<' '<<NUM_TRIANGLE_CREATED;
    cout<<"V_c: "<<V.rows()<<' '<<V.cols()<<' '<<NUM_TRIANGLE_CREATED;
    for (int i=0;i<V_color.cols();i++){
        V_color.col(i)<<0.2f, 0.2f, 0.8f;
    }
}
void highlight_color_triangle(int triangle_ind){
    V_color.col(3*triangle_ind) <<1.0f,0.0f,0.0f;
    V_color.col(3*triangle_ind+1) <<1.0f,0.0f,0.0f;
    V_color.col(3*triangle_ind+2) <<1.0f,0.0f,0.0f;
}

void paint_color_vertex(int vertex_ind,Vector3f color){
    V_color.col(vertex_ind)<<color;
}

void remove_tri(int tri_target)
{
    MatrixXf V_display(2,3);
    
    if (tri_target!= NUM_TRIANGLE_CREATED){
        MatrixXf before(2,tri_target*3);
        before << V(seqN(0,2),seq(0,tri_target*3-1));
        MatrixXf after(2,NUM_TRIANGLE_CREATED*3-tri_target*3);
        after << V(seqN(0,2),seq(tri_target*3+3,last));
        V_display.resize(before.rows(),before.cols()+after.cols());
        V_display << before,after;
        
    } else{
        MatrixXf before(2,tri_target*3);
        before << V(all, seq(0,tri_target*3-1));
        V_display.resize(before.rows(),before.cols());
        V_display << before;
    }
    V.resize(V_display.rows(),V_display.cols());
    V<<V_display;
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    // Get the position of the mouse in the window
     // NOTE: y axis is flipped in glfw
    vector<double> coord = cur_loc_get(window);
    double xworld = coord[0];
    double yworld = coord[1];
    // Update the position of the first vertex if the left button is pressed
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
    {
        if(MODE==1){
            cout<<"releases"<<endl;
            CLICK_COUNTER = 0;
        }
    }
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        
        if (MODE==0){
            selected_tri = -1;
            //cout<<CLICK_COUNTER;
            MatrixXf V_display(2,3);
            MatrixXf V_display_c(3,3);
            
            MatrixXf V_temp(2,3);
            MatrixXf V_temp_c(3,3);
            
            if (CLICK_COUNTER>=3) {
                a_lx = 0; a_ly = 0;
                b_lx = 0; b_ly = 0;
                c_lx = 0; c_ly = 0;
                CLICK_COUNTER=0;
            }
            
            if(CLICK_COUNTER==0) {
                a_lx = xworld;
                a_ly = yworld;
            }else if (CLICK_COUNTER==1){
                b_lx = xworld;
                b_ly = yworld;
                
            }else if (CLICK_COUNTER==2) {
                c_lx = xworld;
                c_ly = yworld;
                V_temp << a_lx,b_lx,c_lx,a_ly,b_ly,c_ly;
                V_display.resize(V.rows(),V.cols()+V_temp.cols());
                V_display<<V,V_temp;
                V.resize(V_display.rows(),V_display.cols());
                V<<V_display;
                
                V_temp_c<<  0.0,0.0,0.0,
                             0.0,0.0,0.0,
                             1.0,1.0,1.0;
                
                V_display_c.resize(V_color.rows(),V_color.cols()+V_temp_c.cols());
                V_display_c<<V_color,V_temp_c;
                V_color.resize(V_display_c.rows(),V_display_c.cols());
                V_color<<V_display_c;
                
                NUM_TRIANGLE_CREATED ++;
                selected_tri = NUM_TRIANGLE_CREATED;
            }
            CLICK_COUNTER ++;
            
        }  else if (MODE ==1){
            int index = index_of_selected_triangle(xworld,yworld);
            if (index>0 && index<=NUM_TRIANGLE_CREATED){
                selected_tri = index;
                CLICK_COUNTER ++;
                coord_rela_tri(0,0) = xworld - V(0,index*3);
                coord_rela_tri(1,0) = yworld - V(1,index*3);
                coord_rela_tri(0,1) = xworld - V(0,index*3+1);
                coord_rela_tri(1,1) = yworld - V(1,index*3+1);
                coord_rela_tri(0,2) = xworld - V(0,index*3+2);
                coord_rela_tri(1,2) = yworld - V(1,index*3+2);
                ini_color_vertex();
                
            }
        } else if (MODE==2){
            selected_tri = -1;
            int index = index_of_selected_triangle(xworld,yworld);
            cout<<"i: "<<index<<" ; num: "<<NUM_TRIANGLE_CREATED<<endl;
            if (index>0 && index<=NUM_TRIANGLE_CREATED){
                remove_tri(index);
                NUM_TRIANGLE_CREATED--;
            }
            
        } else if (MODE == 3){
            CLICK_COUNTER  = 1;
            selected_vertex = index_of_closest_vertex(xworld,yworld);
            
        } else if (MODE == 5){
            CLICK_COUNTER  = 1;
            int index = index_of_selected_triangle(xworld,yworld);
            selected_tri = index;
        } else if (MODE == 7){
            if (CLICK_COUNTER == 3){
                //V_curve.resize()
                
            }
                
            CLICK_COUNTER++;
        }
    }
    // Upload the change to the GPU
    VBO_C.update(V_color);
    VBO.update(V);
    
}

void view_transform(double scale,double x_trans,double y_trans){
    MatrixXf V_view_new(4,4);
    V_view_new <<
        scale, 0,      0,  x_trans,
        0,     scale,  0,  y_trans,
        0,     0,      1,        0,
        0,     0,      0,        1;
    V_view = V_view*V_view_new;
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    // Update the position of the first vertex if the keys 1,2, or 3 are pressed
    
    switch (key)
    {
        case GLFW_KEY_EQUAL:
            view_transform(1.2,0,0);
            break;
        case GLFW_KEY_MINUS:
            view_transform(0.8,0,0);
            break;
        case GLFW_KEY_W:
            view_transform(1,0,-0.2);
            break;
        case GLFW_KEY_S:
            view_transform(1,0,0.2);
            break;
        case GLFW_KEY_A:
            view_transform(1,0.2,0);
            break;
        case GLFW_KEY_D:
            view_transform(1,-0.2,0);
            break;
            
        case GLFW_KEY_I:
            MODE=0;
            CLICK_COUNTER = 0;
            break;
        case GLFW_KEY_O:
            MODE=1;
            CLICK_COUNTER = 0;
            break;
        case GLFW_KEY_P:
            MODE=2;
            CLICK_COUNTER = 0;
            break;
        case GLFW_KEY_Y:
            MODE=5;
            CLICK_COUNTER = 0;
            break;

        case GLFW_KEY_H:
            if(action == GLFW_PRESS && MODE==5 && selected_tri>0){
                set_rotate_tri_intp(true,50.0,0.3,200);
            }
            if (action == GLFW_PRESS && MODE==1 && selected_tri>0){
                rotate_tri(selected_tri,true,10.0);
            }
            break;
        case GLFW_KEY_J:
            if(action == GLFW_PRESS && MODE==5 && selected_tri>0){
                set_rotate_tri_intp(false,50.0,0.3,200);
            }
            
            if (action == GLFW_PRESS && MODE==1 && selected_tri>0){
                rotate_tri(selected_tri,false,10.0);
            }
            break;
        case GLFW_KEY_K:
            if (action == GLFW_PRESS && MODE==1 and selected_tri>0){
                scale_tri(selected_tri,true);
            }
            break;
        case GLFW_KEY_L:
            if (action == GLFW_PRESS && MODE==1  and selected_tri>0){
                scale_tri(selected_tri,false);
            }
            break;
        
        case GLFW_KEY_C:
            MODE=3;
            ini_color_vertex();
            CLICK_COUNTER = 0;
            break;
            
        case GLFW_KEY_1:
            if (MODE ==3 && CLICK_COUNTER)
                paint_color_vertex(selected_vertex,Vector3f(1.0f, 0.0f, 0.0f));
            break;
        case GLFW_KEY_2:
            if (MODE ==3 && CLICK_COUNTER)
                paint_color_vertex(selected_vertex,Vector3f(0.8f, 0.2f, 0.0f));
            break;
        case GLFW_KEY_3:
            if (MODE ==3 && CLICK_COUNTER)
                paint_color_vertex(selected_vertex,Vector3f(0.6f, 0.4f, 0.2f));
            break;
        case GLFW_KEY_4:
            if (MODE ==3 && CLICK_COUNTER)
                paint_color_vertex(selected_vertex,Vector3f(0.4f, 0.6f, 0.2f));
            break;
        case GLFW_KEY_5:
            if (MODE ==3 && CLICK_COUNTER)
                paint_color_vertex(selected_vertex,Vector3f(0.2f, 0.8f, 0.2f));
            break;
        case GLFW_KEY_6:
            if (MODE ==3 && CLICK_COUNTER)
                paint_color_vertex(selected_vertex,Vector3f(0.0f, 1.0f, 0.0f));
            break;
        case GLFW_KEY_7:
            if (MODE ==3 && CLICK_COUNTER)
                paint_color_vertex(selected_vertex,Vector3f(0.0f, 0.4f, 0.4f));
            break;
        case GLFW_KEY_8:
            if (MODE ==3 && CLICK_COUNTER)
                paint_color_vertex(selected_vertex,Vector3f(0.0f, 0.2f, 0.6f));
            break;
        case GLFW_KEY_9:
            if (MODE ==3 && CLICK_COUNTER)
                paint_color_vertex(selected_vertex,Vector3f(0.0f, 0.2f, 0.6f));
            break;
            
        default:
            break;
    }
    // Upload the change to the GPU
    VBO.update(V);
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
    window = glfwCreateWindow(640, 480, "Hello World", NULL, NULL);
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

    // Initialize the VAO
    // A Vertex Array Object (or VAO) is an object that describes how the vertex
    // attributes are stored in a Vertex Buffer Object (or VBO). This means that
    // the VAO is not the actual object storing the vertex data,
    // but the descriptor of the vertex data.
    VertexArrayObject VAO;
    VAO.init();
    VAO.bind();

    // Initialize the VBO with the vertices data
    // A VBO is a data container that lives in the GPU memory
    VBO.init();
    VBO_C.init();
    
    V.resize(2,3);
    V <<0,0,0,0,0,0;
    VBO.update(V);

    V_color.resize(3,3);
    V_color <<0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,1.0f,1.0f,1.0f;
    VBO_C.update(V_color);
    
    V_view <<
    1,     0,      0,        0,
    0,     1,      0,        0,
    0,     0,      1,        0,
    0,     0,      0,        1;
    
    // Initialize the OpenGL Program
    // A program controls the OpenGL pipeline and it must contains
    // at least a vertex shader and a fragment shader to be valid
    Program program;
    
    const GLchar* vertex_shader =
    "#version 150 core\n"
    "in vec2 position;"
    "uniform mat4 view;"
    "in vec3 color;"
    "out vec3 f_color;"
    "void main()"
    "{"
    "    gl_Position = view * vec4(position, 0.0, 1.0);"
    "    f_color = color;"
    "}";
    
    const GLchar* fragment_shader =
    "#version 150 core\n"
    "in vec3 f_color;"
    "out vec4 outColor;"
    "uniform vec3 triangleColor;"
    "void main()"
    "{"
    "    outColor = vec4(f_color, 1.0);"
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
    program.bindVertexAttribArray("color",VBO_C);
    // Save the current time --- it will be used to dynamically change the triangle color

    // Register the keyboard callback
    glfwSetKeyCallback(window, key_callback);

    // Register the mouse callback
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    // Update viewport
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    

    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window))
    {
        // Bind your VAO (not necessary if you have only one)
        VAO.bind();

        // Bind your program
        program.bind();
      
        glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        // Enable blending test
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        //glUniform3f(program.uniform("triangleColor"),0.2f, 0.4f, 0.8f);

        vector<double> coord = cur_loc_get(window);
        int cur_tri = index_of_selected_triangle(coord[0], coord[1]);
        
        cout<<"tri: "<<cur_tri<<endl;
        
        if (MODE==0){
            MatrixXf V_temp(2,3);
            MatrixXf V_temp_c(3,3);
            if (CLICK_COUNTER==0){
                VBO.update(V);
                VBO_C.update(V_color);
                glDrawArrays(GL_TRIANGLES, 0, V.cols());
                
            }else if (CLICK_COUNTER==1){
                MatrixXf V_L(2,2);
                MatrixXf V_L_c(3,2);
                
                V_L<<coord[0],a_lx,coord[1],a_ly;
                V_L_c<<0.0f,0.0f,0.0f,0.0f,0.0f,0.0f;
                
                MatrixXf V_d(V.rows(),V.cols()+V_L.cols());
                MatrixXf V_c(V_color.rows(),V_color.cols()+V_L_c.cols());
                
                V_d<<V,V_L;
                V_c<<V_color,V_L_c;
                
                VBO.update(V_d);
                VBO_C.update(V_c);
                
                glDrawArrays(GL_TRIANGLES, 0, V.cols());
                //glUniform3f(program.uniform("triangleColor"),0.0f, 0.0f, 0.0f);
                glDrawArrays(GL_LINES,V.cols(),2);
                
            } else if (CLICK_COUNTER==2){
                V_temp<<coord[0],a_lx,b_lx,coord[1],a_ly,b_ly;
                V_temp_c<<0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f;
                
                MatrixXf V_d(V.rows(),V.cols()+V_temp.cols());
                MatrixXf V_c(V_color.rows(),V_color.cols()+V_temp_c.cols());
                
                V_d<<V,V_temp;
                V_c<<V_color,V_temp_c;
                
                VBO.update(V_d);
                VBO_C.update(V_c);
                
                glDrawArrays(GL_TRIANGLES, 0, V.cols());
                glDrawArrays(GL_LINE_LOOP, V.cols(),3);
                
            } else if (CLICK_COUNTER==3){
                cout<<"triangle created"<<endl;
                cout<<selected_tri<<endl;
                ini_color_vertex();
                highlight_color_triangle(selected_tri);
                VBO.update(V);
                VBO_C.update(V_color);
                
                glDrawArrays(GL_TRIANGLES, 0, V.cols());
            }
            
        } else if (MODE==1){
            
            if (selected_tri>0 && CLICK_COUNTER==1){
                //ini_color_vertex();
                V(0,selected_tri*3) = coord[0] - coord_rela_tri(0,0);
                V(1,selected_tri*3) = coord[1] - coord_rela_tri(1,0);
                V(0,selected_tri*3+1) = coord[0] - coord_rela_tri(0,1);
                V(1,selected_tri*3+1) = coord[1] - coord_rela_tri(1,1);
                V(0,selected_tri*3+2) = coord[0] - coord_rela_tri(0,2);
                V(1,selected_tri*3+2) = coord[1] - coord_rela_tri(1,2);
                ini_color_vertex();
                highlight_color_triangle(selected_tri);
            }
            
            VBO.update(V);
            VBO_C.update(V_color);
            glDrawArrays(GL_TRIANGLES, 0, V.cols());
            
        } else if (MODE==2){
            cout<<V.rows()<<' '<<V.cols()<<endl;
            VBO.update(V);
            ini_color_vertex();
            VBO_C.update(V_color);
            glDrawArrays(GL_TRIANGLES, 0, V.cols());
            
        } else if (MODE==3){
            cout<<"closest vertex ind: "<<selected_vertex<<endl;
            //ini_color_vertex();
            VBO.update(V);
            VBO_C.update(V_color);
            glDrawArrays(GL_TRIANGLES, 0, V.cols());
        } else if (MODE==5){
            if (selected_tri>0){
                //ini_color_vertex();
                ini_color_vertex();
                highlight_color_triangle(selected_tri);
                auto t_now = chrono::high_resolution_clock::now();
                float time = chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();
                if (time>=ani_mv_interval && ani_mv_remain>0){
                    rotate_tri(selected_tri,ani_dir,val_per_mv);
                    ani_mv_remain--;
                    t_start = chrono::high_resolution_clock::now();
                }
            
            }
            
            VBO.update(V);
            VBO_C.update(V_color);
            glDrawArrays(GL_TRIANGLES, 0, V.cols());
        }
        glUniformMatrix4fv(program.uniform("view"), 1, GL_FALSE, V_view.data());
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

    // Deallocate glfw internals
    glfwTerminate();
    return 0;
}
