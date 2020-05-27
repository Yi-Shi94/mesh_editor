// C++ include
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <string>
#include <stdlib.h>
//#include "tbb/task_scheduler_init.h"
//#include "tbb/parallel_for.h"
//#include "tbb/blocked_range.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"
#include "utils.h"
#include "obj.h"

#define size_img_1D 800*800

// Shortcut to avoid Eigen:: and std:: everywhere, DO NOT USE IN .h
using namespace std;
using namespace Eigen;
//using namespace tbb;

//#include "tbb/parallel_for.h"
//#include "tbb/blocked_range.h"

static const size_t N = 9;


void part1()
{
    std::cout << "Part 1" << std::endl;
    const std::string filename("part1.png");
    MatrixXd CRGB[3] = {MatrixXd::Zero(800,800),MatrixXd::Zero(800,800),MatrixXd::Zero(800,800)};
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask
    MatrixXd C = MatrixXd::Zero(800,800);
    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(-1,1,1);
    Vector3d x_displacement(2.0/A.cols(),0,0);
    Vector3d y_displacement(0,-2.0/A.rows(),0);

    //mt, rad, x,y,z, r,g,b
    
    Sphere sphere_0(0, 0.5, 0.1,-0.3,-1,  1,1,1);
    Sphere *sphere_ptr_0 = &sphere_0;
    
    Sphere sphere_1(0, 0.3, 0.4,0.4,-2, 1,1,1);
    Sphere *sphere_ptr_1 = &sphere_1;
    //x,y,z
    //CGLightSc light_0(5,3,-3);
    CGLightSc light_0(1,1,3);

    vector<CGLightSc> light_lst;
    light_lst.push_back(light_0);
    
    vector<CGObject*> obj_lst;
    obj_lst.push_back(sphere_ptr_0);
    obj_lst.push_back(sphere_ptr_1);
    
    Vector3d ray_direction = RowVector3d(0,0,-1);
    
    int cols = C.cols();
    int rows = C.rows();
    
    
    for (unsigned ind=0;ind<size_img_1D;ind++)
    {
        int j = ind % rows;
        int i = ind / rows;
            // Prepare the ray
        Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            
        double min_dis = A_LARGE_NUMBER;
        int min_obj_ind = -1;
        double t = -1;

        for (int ind=0; ind<obj_lst.size(); ++ind)
        {
                
            CGObject *obj = obj_lst[ind];
            double t = obj->t_hit(ray_origin,ray_direction);
            //cout<<ind<<' ' <<t<<endl;
            if (t < min_dis && t>0){
                min_dis = t;
                min_obj_ind = ind;
            }
        }
            
        if(min_dis<A_LARGE_NUMBER){
                
            CGObject *closest_obj = obj_lst[min_obj_ind];
            Vector3d ray_intersection = ray_origin+min_dis*ray_direction;
            Vector3d c_vec = closest_obj->get_shading(ray_origin,ray_direction,ray_intersection,obj_lst,light_lst,0,min_obj_ind);
                
            for (int c=0;c<3;c++) CRGB[c](i,j)= c_vec[c];
            A(i,j) = 1;
        }
    }
    
    write_matrix_to_png(CRGB[0],CRGB[1],CRGB[2],A,filename);
}

void part2()
{
    std::cout << "Part 2" << std::endl;
    const std::string filename("part2.png");
    MatrixXd CRGB[3] = {MatrixXd::Zero(800,800),MatrixXd::Zero(800,800),MatrixXd::Zero(800,800)};
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask
    //MatrixXd C = MatrixXd::Zero(800,800);
    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(-1,1,1);
    Vector3d x_displacement(2.0/A.cols(),0,0);
    Vector3d y_displacement(0,-2.0/A.rows(),0);
    
    //mt, rad, x,y,z, r,g,b
    
    Sphere sphere_0(2, 0.5, 0.1,-0.3,-1,  1,0.2,0.6);
    Sphere *sphere_ptr_0 = &sphere_0;
    
    Sphere sphere_1(0, 0.3, 0.4,0.4,-2, 0.1,1,0.3);
    Sphere *sphere_ptr_1 = &sphere_1;
    //x,y,z
    CGLightSc light_0(5,3,-3);
    CGLightSc light_1(1,1,3);
    
    vector<CGLightSc> light_lst;
    light_lst.push_back(light_0);
    light_lst.push_back(light_1);
    
    vector<CGObject*> obj_lst;
    obj_lst.push_back(sphere_ptr_0);
    obj_lst.push_back(sphere_ptr_1);
    
    Vector3d ray_direction = RowVector3d(0,0,-1);
    
    int cols = A.cols();
    int rows = A.rows();
    
    
    for (unsigned ind=0;ind<size_img_1D;ind++)
    {
        int j = ind % rows;
        int i = ind / rows;
        // Prepare the ray
        Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
        
        double min_dis = A_LARGE_NUMBER;
        int min_obj_ind = -1;
        double t = -1;
        
        for (int ind=0; ind<obj_lst.size(); ++ind)
        {
            
            CGObject *obj = obj_lst[ind];
            double t = obj->t_hit(ray_origin,ray_direction);
            //cout<<ind<<' ' <<t<<endl;
            if (t < min_dis && t>0){
                min_dis = t;
                min_obj_ind = ind;
            }
        }
        
        if(min_dis<A_LARGE_NUMBER){
            
            CGObject *closest_obj = obj_lst[min_obj_ind];
            Vector3d ray_intersection = ray_origin+min_dis*ray_direction;
            Vector3d c_vec = closest_obj->get_shading(ray_origin,ray_direction,ray_intersection,obj_lst,light_lst,0,min_obj_ind);
            
            for (int c=0;c<3;c++) CRGB[c](i,j)= c_vec[c];
            A(i,j) = 1;
        }
    }

    write_matrix_to_png(CRGB[0],CRGB[1],CRGB[2],A,filename);
}


void part3(){
    std::cout << "Part 3" << std::endl;
    
    const std::string filename("part3.png");
    //MatrixXd C = MatrixXd::Zero(800,800); // Store the color
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask
    MatrixXd CRGB[3] = {MatrixXd::Zero(800,800),MatrixXd::Zero(800,800),MatrixXd::Zero(800,800)};
    
    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d screen_origin(-1,1,1);
    Vector3d camera_origin(0,0,4);
    Vector3d x_displacement(2.0/A.cols(),0,0);
    Vector3d y_displacement(0,-2.0/A.rows(),0);
    
    //mt, rad, x,y,z, r,g,b
    
    Sphere sphere_0(2, 0.5, 0.1,-0.3,-1,  1,0.2,0.6);
    Sphere *sphere_ptr_0 = &sphere_0;
    
    Sphere sphere_1(0, 0.3, 0.4,0.4,-2, 0.1,1,0.3);
    Sphere *sphere_ptr_1 = &sphere_1;
    //x,y,z
    CGLightSc light_0(5,3,-3);
    CGLightSc light_1(1,1,3);
    
    vector<CGLightSc> light_lst;
    light_lst.push_back(light_0);
    light_lst.push_back(light_1);
    
    vector<CGObject*> obj_lst;
    obj_lst.push_back(sphere_ptr_0);
    obj_lst.push_back(sphere_ptr_1);
    
    
    int cols = A.cols();
    int rows = A.rows();
    
    
    for (unsigned ind=0;ind<size_img_1D;ind++)
    {
        int j = ind % rows;
        int i = ind / rows;
        // Prepare the ray
        
        Vector3d screen_coord = screen_origin + double(i)*x_displacement + double(j)*y_displacement;
        Vector3d ray_direction = (screen_coord-camera_origin).normalized();
        Vector3d ray_origin = camera_origin;
        
        double min_dis = A_LARGE_NUMBER;
        int min_obj_ind = -1;
        double t = -1;
        
        for (int ind=0; ind<obj_lst.size(); ++ind)
        {
            CGObject *obj = obj_lst[ind];
            double t = obj->t_hit(ray_origin,ray_direction);
            //cout<<ind<<' ' <<t<<endl;
            if (t < min_dis && t>0){
                min_dis = t;
                min_obj_ind = ind;
            }
        }
        
        if(min_dis<A_LARGE_NUMBER){
            
            CGObject *closest_obj = obj_lst[min_obj_ind];
            Vector3d ray_intersection = ray_origin+min_dis*ray_direction;
            Vector3d c_vec = closest_obj->get_shading(ray_origin,ray_direction,ray_intersection,obj_lst,light_lst,0,min_obj_ind);
            
            for (int c=0;c<3;c++) CRGB[c](i,j)= c_vec[c];
            A(i,j) = 1;
        }
    }
    
    write_matrix_to_png(CRGB[0],CRGB[1],CRGB[2],A,filename);

}


void part4_1(){
    cout << "Part 4_1" << std::endl;
    const std::string filename("part4_1.png");
    MatrixXd Amap = MatrixXd::Zero(200,200); // Store the alpha mask
    vector<MatrixXd> Cmap = {MatrixXd::Zero(200,200),MatrixXd::Zero(200,200),MatrixXd::Zero(200,200)};
    //MatrixXd CRGB[3] = {MatrixXd::Zero(800,800),MatrixXd::Zero(800,800),MatrixXd::Zero(800,800)};
    Vector3d screen_origin(-1.5,2,1);
    Vector3d camera_origin(-0.5,1,4);
    double scale = 10;
    Vector3d x_displacement(2.0/Amap.cols(),0,0);
    Vector3d y_displacement(0,-2.0/Amap.rows(),0);
    
    const int num_obj = 1;
    
    TriMesh trims_0(2,0.3,0.4,0.7);
    trims_0.readMesh("../data/bunny.off");
    Vector3d rgb = trims_0.rgb;
    
    TriMesh obj_lst[1] = {trims_0};
    const int num_light = 2;
    const Vector3d light_position_0(-1,1,2);
    const Vector3d light_position_1(2,-2,12);
    
    Vector3d light_position_lst[num_light] = {light_position_0,light_position_1};
    
    // Single light sourcesphere_material_0
    // for (unsigned c=0;c<3;c++)
    //{
    for (unsigned i=0;i<Amap.cols();i++)
    {
        for (unsigned j=0;j<Amap.rows();j++)
        {
            //cout<<' '<<i<<' '<<j<<endl;
            
            Vector3d screen_coord = screen_origin + double(i)*x_displacement + double(j)*y_displacement;
            Vector3d ray_direction = (screen_coord-camera_origin).normalized();
            Vector3d ray_origin = camera_origin;
            
            double min_dis = A_LARGE_NUMBER;
            int min_obj_ind = -1;
            double t = -1;
            
            
            for (int p_ind=0;p_ind<trims_0.get_plain_num();p_ind++)
            {
                Vector3d vertices_indeces = trims_0.plains[p_ind];
                //3*1
                Vector3d vert_a = scale*trims_0.vertices[vertices_indeces[0]];
                Vector3d vert_b = scale*trims_0.vertices[vertices_indeces[1]];
                Vector3d vert_c = scale*trims_0.vertices[vertices_indeces[2]];
                
                double a_x = vert_a[0];
                double a_y = vert_a[1];
                double a_z = vert_a[2];
                //cout<<a_x<<a_y<<a_z<<endl;
                double b_x = vert_b[0];
                double b_y = vert_b[1];
                double b_z = vert_b[2];
                
                double c_x = vert_c[0];
                double c_y = vert_c[1];
                double c_z = vert_c[2];
                
                double e_x = ray_origin[0];
                double e_y = ray_origin[1];
                double e_z = ray_origin[2];
                
                double d_x = ray_direction[0];
                double d_y = ray_direction[1];
                double d_z = ray_direction[2];
                
                double A11 = a_x-b_x;
                double A12 = a_x-c_x;
                double A21 = a_y-b_y;
                double A22 = a_y-c_y;
                double A31 = a_z-b_z;
                double A32 = a_z-c_z;
                
                double Xae = a_x-e_x;
                double Yae = a_y-e_y;
                double Zae = a_z-e_z;
                
                double Axd_sub = A21*A32-A22*A31;
                double Ayd_sub = A11*A32-A12*A31;
                double Azd_sub = A11*A22-A21*A12;
                
                double A = d_x*Axd_sub - d_y*Ayd_sub + d_z*Azd_sub;
                double t = (Xae*Axd_sub - Yae*Ayd_sub +  Zae*Azd_sub)/A;
                
                
                if (t<0){
                    continue;
                }
                
                double gamma = (d_x * (A21*Zae-Yae*A31) -
                                d_y * (A11*Zae-Xae*A31) +
                                d_z * (A11*Yae-Xae*A21))/A;
                
                if (gamma<0 || gamma>1){
                    continue;
                }
                
                
                double beta = (d_x * (Yae*A32-Zae*A22) -
                               d_y * (Xae*A32-Zae*A12) +
                               d_z * (Xae*A22-Yae*A12))/A;
                
                if(beta<0 || beta>1-gamma){
                    continue;
                }
                
                if (t<min_dis){
                    min_dis = t;
                    min_obj_ind = p_ind;
                }
                
            }
            
            Vector3d vertices_indeces = trims_0.plains[min_obj_ind];
            Vector3d a = scale*trims_0.vertices[vertices_indeces[0]];
            Vector3d b = scale*trims_0.vertices[vertices_indeces[1]];
            Vector3d c = scale*trims_0.vertices[vertices_indeces[2]];
            
            Vector3d intersection = ray_origin + min_dis*ray_direction;
            Vector3d inter_norm = ((a-c).cross(b-c)).normalized();
            for (int c=0;c<3;c++) Cmap[c](i,j) +=  0.2;
            for (int ind=0; ind<num_light; ++ind){
                Vector3d l = light_position_lst[ind] - intersection;
                Vector3d v = ray_origin - intersection;
                double lambertian  = l.normalized().transpose() * inter_norm;
                double specular = inter_norm.transpose()* ((v+l).normalized());
                if (trims_0.material == 0)
                {
                    for (int c=0;c<3;c++) Cmap[c](i,j) += 0.4*rgb[c]*max(lambertian,0.);
                }
                else if (trims_0.material == 1){
                    for (int c=0;c<3;c++) Cmap[c](i,j) += 0.4*rgb[c]*pow(max(specular,0.),1000);
                } else{
                    for (int c=0;c<3;c++) Cmap[c](i,j) += 0.4*rgb[c]*(max(lambertian,0.) + pow(max(specular,0.),1000));
                }
                    
            }
                
                // Disable the alpha mask for this pixel
            Amap(i,j) = 1;

        }
    }
    //}
    write_matrix_to_png(Cmap[0],Cmap[1],Cmap[2],Amap,filename);
}


void part4_2(){
    cout << "Part 4_2" << std::endl;
    const std::string filename("part4_2.png");
    MatrixXd Amap = MatrixXd::Zero(200,200); // Store the alpha mask
    //MatrixXd Cmap = MatrixXd::Zero(400,400);
    MatrixXd Cmap[3] = {MatrixXd::Zero(200,200),MatrixXd::Zero(200,200),MatrixXd::Zero(200,200)};
    Vector3d screen_origin(-6,6,100);
    Vector3d camera_origin(-6,6,120);
    Vector3d x_displacement(2./Amap.cols(),0,0);
    Vector3d y_displacement(0,-2./Amap.rows(),0);
    
    const int num_obj = 1;
    
    TriMesh trims_0(2,0.3,0.4,0.7);
    trims_0.readMesh("../data/bumpy_cube.off");
    Vector3d rgb = trims_0.rgb;
    
    TriMesh obj_lst[1] = {trims_0};
    const int num_light = 2;
    const Vector3d light_position_0(0,10,40);
    const Vector3d light_position_1(2,-2,6);
    
    Vector3d light_position_lst[num_light] = {light_position_0,light_position_1};
    
    // Single light sourcesphere_material_0
    // for (unsigned c=0;c<3;c++)
    //{
    for (unsigned i=0;i<Amap.cols();i++)
    {
        //cout<<' '<<i<<endl;
        for (unsigned j=0;j<Amap.rows();j++)
        {
            
            Vector3d screen_coord = screen_origin + double(i)*x_displacement + double(j)*y_displacement;
            Vector3d ray_direction = (screen_coord-camera_origin).normalized();
            Vector3d ray_origin = camera_origin;
            
            double min_dis = A_LARGE_NUMBER;
            int min_obj_ind = -1;
            double t = -1;
            
            for (int p_ind=0;p_ind<trims_0.get_plain_num();p_ind++)
            {
                Vector3d vertices_indeces = trims_0.plains[p_ind];
                
                //3*1
                Vector3d vert_a = trims_0.vertices[vertices_indeces[0]];
                Vector3d vert_b = trims_0.vertices[vertices_indeces[1]];
                Vector3d vert_c = trims_0.vertices[vertices_indeces[2]];
                
                double a_x = vert_a[0];
                double a_y = vert_a[1];
                double a_z = vert_a[2];
                //cout<<a_x<<a_y<<a_z<<endl;
                double b_x = vert_b[0];
                double b_y = vert_b[1];
                double b_z = vert_b[2];
                
                double c_x = vert_c[0];
                double c_y = vert_c[1];
                double c_z = vert_c[2];
                
                double e_x = ray_origin[0];
                double e_y = ray_origin[1];
                double e_z = ray_origin[2];
                
                double d_x = ray_direction[0];
                double d_y = ray_direction[1];
                double d_z = ray_direction[2];
                
                double A11 = a_x-b_x;
                double A12 = a_x-c_x;
                double A21 = a_y-b_y;
                double A22 = a_y-c_y;
                double A31 = a_z-b_z;
                double A32 = a_z-c_z;
                
                double Xae = a_x-e_x;
                double Yae = a_y-e_y;
                double Zae = a_z-e_z;
                
                double Axd_sub = A21*A32-A22*A31;
                double Ayd_sub = A11*A32-A12*A31;
                double Azd_sub = A11*A22-A21*A12;
                
                double A = d_x*Axd_sub - d_y*Ayd_sub + d_z*Azd_sub;
                double t = (Xae*Axd_sub - Yae*Ayd_sub +  Zae*Azd_sub)/A;
                
                if (t<0){
                    continue;
                }
                
                double gamma = (d_x * (A21*Zae-Yae*A31) -
                                d_y * (A11*Zae-Xae*A31) +
                                d_z * (A11*Yae-Xae*A21))/A;
                
                if (gamma<0 || gamma>1){
                    continue;
                }
                
                
                double beta = (d_x * (Yae*A32-Zae*A22) -
                               d_y * (Xae*A32-Zae*A12) +
                               d_z * (Xae*A22-Yae*A12))/A;
                
                if(beta<0 || beta>1-gamma){
                    continue;
                }
                
                if (t<min_dis){
                    min_dis = t;
                    min_obj_ind = p_ind;
                }
            }
            
            Vector3d vertices_indeces = trims_0.plains[min_obj_ind];
            Vector3d a = trims_0.vertices[vertices_indeces[0]];
            Vector3d b = trims_0.vertices[vertices_indeces[1]];
            Vector3d c = trims_0.vertices[vertices_indeces[2]];
            
            Vector3d intersection = ray_origin + min_dis*ray_direction;
            Vector3d inter_norm = ((a-c).cross(b-c)).normalized();
                
            for (int c=0;c<3;c++) Cmap[c](i,j) +=  0.2;
            for (int ind=0; ind<num_light; ++ind){
                Vector3d l = light_position_lst[ind] - intersection;
                Vector3d v = ray_origin - intersection;
                double lambertian  = l.normalized().transpose() * inter_norm;
                double specular = inter_norm.transpose()* ((v+l).normalized());
                if (trims_0.material == 0)
                {
                    for (int c=0;c<3;c++) Cmap[c](i,j) += 0.5*rgb[c]*max(lambertian,0.);
                }
                else if (trims_0.material == 1){
                    for (int c=0;c<3;c++) Cmap[c](i,j) += 0.5*rgb[c]*pow(max(specular,0.),100);
                } else{
                    for (int c=0;c<3;c++) Cmap[c](i,j) += 0.5*rgb[c]*( max(lambertian,0.) +pow(max(specular,0.),100));
                }
                
            }
            Amap(i,j) = 1;
        }
            
    }
   
    write_matrix_to_png(Cmap[0],Cmap[1],Cmap[2],Amap,filename);
}


void part5n6(){
    
    std::cout << "Part 5" << std::endl;
    
    const std::string filename("part5n6.png");
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask
    MatrixXd CRGB[3] = {MatrixXd::Zero(800,800),MatrixXd::Zero(800,800),MatrixXd::Zero(800,800)};
    
    Vector3d screen_origin(-3,3,2);
    Vector3d camera_origin(0,0,12);
    Vector3d x_displacement(6.0/A.cols(),0,0);
    Vector3d y_displacement(0,-6.0/A.rows(),0);

    //mt,y, r,g,b
    Surface sf(3, 2.21, 1, 1, 1);
    Surface *sf_ptr = &sf;
    //mt, rad, x,y,z, r,g,b
    Sphere sphere_0(0, 1.9, -2.0,-0.4, -3,  0.8,0.1,0.1);
    Sphere *sphere_ptr_0 = &sphere_0;
    
    Sphere sphere_1(2, 1.4, 1.3,-0.9, -3,   0.2,0.2,0.9);
    Sphere *sphere_ptr_1 = &sphere_1;
    
    //x,y,z
    CGLightSc light_0(-4.5,4,-3);
    CGLightSc light_1(-5,3,-3);
    
    vector<CGLightSc> light_lst;
    light_lst.push_back(light_0);
    light_lst.push_back(light_1);
    
    vector<CGObject*> obj_lst;
    obj_lst.push_back(sphere_ptr_0);
    obj_lst.push_back(sphere_ptr_1);
    obj_lst.push_back(sf_ptr);
    
    // Single light sourcesphere_material_0
    for (unsigned i=0;i<A.cols();i++)
    {
        for (unsigned j=0;j<A.rows();j++)
        {
            Vector3d screen_coord = screen_origin + double(i)*x_displacement + double(j)*y_displacement;
            Vector3d ray_direction = (screen_coord-camera_origin).normalized();
            Vector3d ray_origin = camera_origin;
            
            double min_dis = A_LARGE_NUMBER;
            int min_obj_ind= -1;
            
            for (int ind=0; ind<obj_lst.size(); ++ind)
            {
                CGObject *obj = obj_lst[ind];
                double t = obj->t_hit(ray_origin,ray_direction);
                if (t > 0 && t < min_dis){
                    min_dis = t;
                    min_obj_ind = ind;
                }
            }
            
            if(min_obj_ind != -1){
                CGObject *closest_obj = obj_lst[min_obj_ind];
                Vector3d ray_intersection = ray_origin+min_dis*ray_direction;
                Vector3d c_vec = closest_obj->get_shading(ray_origin,ray_direction,ray_intersection,obj_lst,light_lst,1,min_obj_ind);
                for (int c=0;c<3;c++) {
                    CRGB[c](i,j)=c_vec[c];
                }
            }
            A(i,j) = 1;
        }
    }
    
    write_matrix_to_png(CRGB[0],CRGB[1],CRGB[2],A,filename);
}



int main()
{   part1();
    part2();
    part3();
    part4_1();
    part4_2();
    part5n6();ioff=
