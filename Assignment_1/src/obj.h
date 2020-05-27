//
//  obj.h
//  
//
//  Created by syd on 2019/10/2.
//
#ifndef obj_h
#define obj_h

#include<Eigen/Geometry>
#include<math.h>

#define A_LARGE_NUMBER 1000
#define A_SMALL_NUMBER 0.0001
#define ambientStrength 0.4
#define lambertianStrength 0.4
#define specularStrength 1

#define SP_PW 45
using namespace std;
using namespace Eigen;

class CGLightSc{
public:
    
    double intens;
    Vector3d rgb;
    Vector3d source_coord;
    CGLightSc(double x,double y,double z);
    void set_color_intens(double red,double green,double blue,double intens);
    Vector3d get_coord(){return source_coord;};
    Vector3d get_color(){return rgb;};
};
CGLightSc::CGLightSc(double x,double y,double z):
source_coord(Vector3d(x,y,z)){
    rgb = Vector3d(1.,1.,1.);
    intens=1;
    cout<<"Light Source: "<<x<<' '<<y<<' '<<z<<endl;
}

void  CGLightSc::set_color_intens(double red,double green,double blue,double intens){
    intens = intens;
    rgb = Vector3d(red,green,blue);
}

class CGObject{
public:
    CGObject(int material,double red,double green,double blue);
    int material;
    Vector3d rgb;
    Vector3d get_color(){return rgb;};
    virtual bool is_hit(Vector3d ray_origin,Vector3d ray_direction)=0;
    virtual double t_hit(Vector3d ray_origin, Vector3d ray_direction)=0;
    virtual Vector3d get_normal(Vector3d intersection)=0;
    virtual Vector3d get_shading(Vector3d ray_origin,Vector3d ray_direction, Vector3d ray_intersection,
                                 vector<CGObject *>obj_lst, vector<CGLightSc> light_lst, bool display_shadow,int index_self)=0;
    
};

CGObject::CGObject(int material_type, double red, double green, double blue):
material(material_type){
    rgb = Vector3d(red,green,blue);
}


class Surface:public CGObject
{
public:
    double y_level;
    Surface(int material_type,double y_level,double red,double green,double blue);
    //double t_hit2(Vector3d screen_coord,Vector3d ray_origin,Vector3d ray_direction);
    double t_hit(Vector3d ray_origin,Vector3d ray_direction);
    bool is_hit(Vector3d ray_origin,Vector3d ray_direction);
    Vector3d get_normal(Vector3d ray_direction);
    Vector3d get_shading(Vector3d ray_origin,Vector3d ray_direction, Vector3d ray_intersection,
                         vector<CGObject *>obj_lst, vector<CGLightSc> light_lst, bool display_shadow,int index_self);
};

Surface::Surface(int material_type,double y_level,double red,double green,double blue):
CGObject(material_type,red,green,blue),y_level(y_level){}
/*
 double Surface::t_hit2(Vector3d screen_coord,Vector3d ray_origin,Vector3d ray_direction){
 if (screen_coord(1) < ray_origin(1) && screen_coord(1)>y_level) return (ray_origin(1)-y_level)/ray_direction(1);
 else if (screen_coord(1) > ray_origin(1) && screen_coord(1)<y_level) return (-ray_origin(1)+y_level)/ray_direction(1);
 else return A_LARGE_NUMBER;
 }
 */
Vector3d Surface::get_normal(Vector3d ray_direction)
{
    //if (ray_direction[1] <= 0) {
    //    ray_normal = - ray_normal;
    //}
    return Vector3d(0,1.,0);
}

double Surface::t_hit(Vector3d ray_origin,Vector3d ray_direction){
    if (is_hit(ray_origin,ray_direction)){
        return (ray_origin(1)-y_level)/ray_direction(1);
    }
    return -1;
}

bool Surface::is_hit(Vector3d ray_origin,Vector3d ray_direction){
    if((ray_origin(1)-y_level)/ray_direction(1)>0) return true;
    else return false;
}

Vector3d Surface::get_shading(Vector3d ray_origin,Vector3d ray_direction, Vector3d ray_intersection,
                              vector<CGObject *>obj_lst, vector<CGLightSc> light_lst, bool display_shadow, int index_self){
    
    Vector3d vec_RGB(0.,0.,0.);
    Vector3d ray_normal = get_normal(ray_intersection);
    
    for (int ind=0;ind<light_lst.size();ind++){
        Vector3d l = light_lst[ind].get_coord() - ray_intersection;
        Vector3d v = ray_origin - ray_intersection;
        bool shadow_hit_flag = false;
        if (display_shadow){
            Vector3d new_p = ray_intersection+A_SMALL_NUMBER*l.normalized();
            for(int o=0;o < obj_lst.size();o++){
                if(o!=index_self && obj_lst[o]->is_hit(new_p,(light_lst[ind].get_coord()-new_p).normalized())){
                    shadow_hit_flag = true;
                    break;
                }
            }
        }
        
        double lambertian  = l.normalized().transpose() * ray_normal;
        double specular = ray_normal.transpose()* ((v+l).normalized());
        
        if (material==-1) continue;
        else if (material == 0) for (int am=0;am<3;am++) vec_RGB[am] +=  rgb[am]*ambientStrength + lambertianStrength * rgb[am] * max(lambertian,0.);
        else if (material == 1) for (int am=0;am<3;am++) vec_RGB[am] +=  rgb[am]*ambientStrength + specularStrength * rgb[am] * pow(max(specular,0.),SP_PW);
        else if (material == 2) for (int am=0;am<3;am++) vec_RGB[am] +=  rgb[am]*ambientStrength + rgb[am]*(specularStrength * pow(max(specular,0.),SP_PW)+ 0.8*lambertianStrength * max(lambertian,0.));
        if (material == 3) {
            if (display_shadow && shadow_hit_flag){
                
            } else {
                for (int am=0;am<3;am++) vec_RGB[am] += 0.3* rgb[am]*( ambientStrength+ specularStrength * pow(max(specular,0.),10000)+lambertianStrength * max(lambertian,0.));
            }
            double length = 2 * ray_direction.transpose() * ray_normal;
            Vector3d r= (ray_direction - length * ray_normal).normalized();
            Vector3d p = ray_intersection + r*A_SMALL_NUMBER;
            double t = -1;
            double min_dis = A_LARGE_NUMBER;
            int min_obj_ind = -1;
            for(int o=0;o<obj_lst.size();o++){
                if(o!=index_self && obj_lst[o]->is_hit(p,r)){
                    double t = obj_lst[o]-> t_hit(p, r);
                    if (t<min_dis) {
                        min_dis = t;
                        min_obj_ind = o;
                    }
                }
            }
            if (min_obj_ind>=0){
                Vector3d new_ray_intersection = p+r*min_dis;
                //cout<<"from   "<<material<<" <-> "<<obj_lst[min_obj_ind]->material<<endl;
                vec_RGB += obj_lst[min_obj_ind] -> get_shading(p, r, new_ray_intersection, obj_lst, light_lst, 1, min_obj_ind);
            }
        }
        else continue;
    }
    return vec_RGB;
}

class Sphere:public CGObject
{
public:
    double Radius;
    Vector3d center;
    Sphere(int material_type,double Radius, double x, double y,
           double z,double red,double green,double blue);
    
    double get_R(){return Radius;};
    Vector3d get_center(){return center;};
    double t_hit(Vector3d ray_origin, Vector3d ray_direction);
    bool is_hit(Vector3d ray_origin,Vector3d ray_direction);
    Vector3d get_normal(Vector3d intersection);
    Vector3d get_shading(Vector3d ray_origin,Vector3d ray_direction,Vector3d ray_intersection,
                         vector<CGObject *>obj_lst, vector<CGLightSc> light_lst, bool is_shadow,int index_self);
};

Sphere::Sphere(int material_type,double Radius,double x,double y,
               double z,double red,double green,double blue):
CGObject(material_type,red,green,blue){
    this->Radius = Radius;
    this->center = Vector3d(x,y,z);
    cout<<"Sphere:x->"<<x<<" y->"<<y<<" z->"<<z<<" R->"<<Radius<<endl;
    cout<<"color: "<<red<<" : "<<green<<" : "<<blue<<endl;
}

bool Sphere::is_hit(Vector3d ray_origin, Vector3d ray_direction){
    Vector3d e_minus_c = ray_origin-center;
    double bb = 2*e_minus_c.transpose()*ray_direction;
    double aa = ray_direction.transpose()*ray_direction;
    double cc = e_minus_c.transpose()*e_minus_c - Radius*Radius;
    double dd = bb*bb-4*aa*cc;
    if (dd<0) return false;
    else return true;
}

double Sphere::t_hit(Vector3d ray_origin, Vector3d ray_direction){
    Vector3d e_minus_c = ray_origin-center;
    
    double bb = 2*e_minus_c.transpose()*ray_direction;
    double aa = ray_direction.transpose()*ray_direction;
    double cc = e_minus_c.transpose()*e_minus_c - Radius*Radius;
    double dd = bb*bb-4*aa*cc;
    
    if (dd<0){
        return A_LARGE_NUMBER;
    }else{
        return (-bb-sqrt(dd))/(2*aa);
    }
}
Vector3d Sphere::get_normal(Vector3d ray_intersection){
    return (ray_intersection-center).normalized();
}

Vector3d Sphere::get_shading(Vector3d ray_origin,Vector3d ray_direction, Vector3d ray_intersection,vector<CGObject *>obj_lst, vector<CGLightSc> light_lst, bool display_shadow, int index_self)
{
    Vector3d vec_RGB(0,0,0);
    
    for (int am=0;am<3;am++) vec_RGB[am] = rgb[am] * ambientStrength;
    for (int ind=0;ind<light_lst.size();ind++){
        Vector3d l = light_lst[ind].get_coord() - ray_intersection;
        
        if (display_shadow){
            Vector3d new_p = ray_intersection + A_SMALL_NUMBER*l;
            bool hit_flag = false;
            for(int o=0;o<obj_lst.size();o++){
                if(o!=index_self && obj_lst[o]->is_hit(new_p,(light_lst[ind].get_coord()-new_p).normalized())){
                    hit_flag = true;
                    break;
                }
            }
            if (hit_flag){
                continue;
            }
        }
        
        Vector3d ray_normal = get_normal(ray_intersection);
        Vector3d v = ray_origin - ray_intersection;
        double lambertian  = l.normalized().transpose() * ray_normal;
        double specular = ray_normal.transpose()* ((v+l).normalized());
        
        if(material == -1)continue;
        
        else if (material == 0) for (int am=0;am<3;am++) vec_RGB[am] += 1.8 * lambertianStrength * rgb[am] * max(lambertian,0.);
        else if (material == 1) for (int am=0;am<3;am++) vec_RGB[am] += specularStrength * rgb[am]  * pow(max(specular,0.),7000);
        else if (material == 2) for (int am=0;am<3;am++) vec_RGB[am] += 1.8 * rgb[am] * (specularStrength * pow(max(specular,0.),7000) + lambertianStrength * max(lambertian,0.));
        else if (material == 3) {
            double length = 2 * ray_direction.transpose() * ray_normal;
            Vector3d r= (ray_direction - length * ray_normal).normalized();
            Vector3d p = ray_intersection + r*A_SMALL_NUMBER;
            double t = -1;
            double min_dis = A_LARGE_NUMBER;
            int min_obj_ind = -1;
            bool hit_flag = false;
            for(int o=0;o<obj_lst.size();o++){
                if(o!=index_self && obj_lst[o]->is_hit(p,r)){
                    double t = obj_lst[o]-> t_hit(ray_origin, ray_direction);
                    if (t<min_dis) {
                        min_dis = t;
                        min_obj_ind = o;
                        hit_flag = true;
                    }
                }
            }
            if (hit_flag){
                
                Vector3d new_ray_intersection = p + r * min_dis;
                Vector3d reflect_RGB = obj_lst[min_obj_ind] -> get_shading(p, r, new_ray_intersection, obj_lst, light_lst, 1, min_obj_ind);
                //cout<<reflect_RGB<<endl;
                //cout<<vec_RGB<<endl;
                vec_RGB = reflect_RGB;
            }
        }
        else continue;
    }
    return vec_RGB;
}

/*
 Vector3d mirror_helper(CGObject * obj, Vector3d ray_origin,Vector3d ray_direction, Vector3d ray_intersection,
 vector<CGObject*> obj_lst, vector<CGLightSc> light_lst, bool display_shadow, int index_self, int jump_remain)
 {
 Vector3d vec_blank(0,0,0);
 if (recurse_remain<=0) {
 return vec_blank;
 }
 Vector3d obj->get_shading(ray_origin,ray_direction,ray_intersection,
 obj_lst, light_lst, 0, index_self);
 
 }
 */


class Triang:public CGObject
{
public:
    Vector3d a,b,c;
    Triang(int material_type, Vector3d a, Vector3d b,
           Vector3d c,double red,double green,double blue);
    
    double t_hit(Vector3d ray_origin, Vector3d ray_direction);
    bool is_hit(Vector3d ray_origin,Vector3d ray_direction);
    Vector3d get_normal(Vector3d intersection);
    Vector3d get_shading(Vector3d ray_origin ,Vector3d ray_direction,Vector3d ray_intersection,
                         vector<CGObject *>obj_lst, vector<CGLightSc> light_lst, bool is_shadow,int index_self);
};

Triang::Triang(int material_type, Vector3d a, Vector3d b, Vector3d c,double red,double green,double blue):
CGObject(material_type,red,green,blue){
    this->a = a;
    this->b = b;
    this->c = c;
    //cout<<"Triangle:a->"<<a<<" b->"<<b<<" c->"<<c<<endl;
    //cout<<"color: "<<red<<" : "<<green<<" : "<<blue<<endl;
}

Vector3d Triang::get_normal(Vector3d intersection){
    return ((a-c).cross(b-c)).normalized();
}

Vector3d Triang::get_shading(Vector3d ray_origin, Vector3d ray_direction
                             , Vector3d ray_intersection,vector<CGObject *>obj_lst, vector<CGLightSc> light_lst, bool is_shadow, int index_self){
    
    Vector3d inter_norm = get_normal(ray_intersection);
    Vector3d vec_RGB(0,0,0);
    
    for (int am=0;am<3;am++) vec_RGB[am] = rgb[am]*ambientStrength;
    for (int ind=0; ind<light_lst.size(); ++ind){
        Vector3d l = light_lst[ind].get_coord() - ray_intersection;
        Vector3d v = ray_origin - ray_intersection;
        double lambertian  = l.normalized().transpose() * inter_norm;
        double specular = inter_norm.transpose()* ((v+l).normalized());
        
        if(material==-1) continue;
        else if (material == 0) for (int am=0;am<3;am++) vec_RGB[am] += 0.4 * lambertianStrength * rgb[am] *max(lambertian,0.);
        else if (material == 1) for (int am=0;am<3;am++) vec_RGB[am] += specularStrength * rgb[am] * pow(max(specular,0.),SP_PW);
        else if (material == 2) for (int am=0;am<3;am++) vec_RGB[am] += rgb[am]*(specularStrength * pow(max(specular,0.),SP_PW) + 0.4 * lambertianStrength * max(lambertian,0.));
        else continue;
    }
    return vec_RGB;
}

bool Triang::is_hit(Vector3d ray_origin, Vector3d ray_direction)
{
    double a_x = a[0]; double a_y = a[1]; double a_z = a[2];
    double b_x = b[0]; double b_y = b[1]; double b_z = b[2];
    double c_x = c[0]; double c_y = c[1]; double c_z = c[2];
    double e_x = ray_origin[0];    double e_y = ray_origin[1];    double e_z = ray_origin[2];
    double d_x = ray_direction[0]; double d_y = ray_direction[1]; double d_z = ray_direction[2];
    double A11 = a_x-b_x; double A12 = a_x-c_x;
    double A21 = a_y-b_y; double A22 = a_y-c_y;
    double A31 = a_z-b_z; double A32 = a_z-c_z;
    double Xae = a_x-e_x; double Yae = a_y-e_y; double Zae = a_z-e_z;
    double Axd_sub = A21*A32-A22*A31;
    double Ayd_sub = A11*A32-A12*A31;
    double Azd_sub = A11*A22-A21*A12;
    double A = d_x*Axd_sub - d_y*Ayd_sub + d_z*Azd_sub;
    double t = (Xae*Axd_sub - Yae*Ayd_sub +  Zae*Azd_sub)/A;
    if (t<0){
        return false;
    }
    double gamma = (d_x * (A21*Zae-Yae*A31) - d_y * (A11*Zae-Xae*A31) +  d_z * (A11*Yae-Xae*A21))/A;
    if (gamma<0 || gamma>1){
        return false;
    }
    double beta = (d_x * (Yae*A32-Zae*A22)-d_y * (Xae*A32-Zae*A12)+d_z * (Xae*A22-Yae*A12))/A;
    if(beta<0 || beta>1-gamma){
        return false;
    }
    return true;
}

double Triang::t_hit(Vector3d ray_origin, Vector3d ray_direction)
{
    double a_x = a[0]; double a_y = a[1]; double a_z = a[2];
    double b_x = b[0]; double b_y = b[1]; double b_z = b[2];
    double c_x = c[0]; double c_y = c[1]; double c_z = c[2];
    double e_x = ray_origin[0];    double e_y = ray_origin[1];    double e_z = ray_origin[2];
    double d_x = ray_direction[0]; double d_y = ray_direction[1]; double d_z = ray_direction[2];
    
    double A11 = a_x-b_x; double A12 = a_x-c_x;
    double A21 = a_y-b_y; double A22 = a_y-c_y;
    double A31 = a_z-b_z; double A32 = a_z-c_z;
    double Xae = a_x-e_x; double Yae = a_y-e_y; double Zae = a_z-e_z;
    double Axd_sub = A21*A32-A22*A31;
    double Ayd_sub = A11*A32-A12*A31;
    double Azd_sub = A11*A22-A21*A12;
    
    double A = d_x*Axd_sub - d_y*Ayd_sub + d_z*Azd_sub;
    double t = (Xae*Axd_sub - Yae*Ayd_sub +  Zae*Azd_sub)/A;
    if (t<0){
        return -1;
    }
    
    double gamma = (d_x * (A21*Zae-Yae*A31) - d_y * (A11*Zae-Xae*A31) +  d_z * (A11*Yae-Xae*A21))/A;
    if (gamma<0 || gamma>1){
        return -1;
    }
    
    double beta = (d_x * (Yae*A32-Zae*A22)-d_y * (Xae*A32-Zae*A12)+d_z * (Xae*A22-Yae*A12))/A;
    if(beta<0 || beta>1-gamma){
        return -1;
    }
    return t;
}


class TriMesh{
public:
    vector<Vector3d> vertices;
    vector<Vector3d> plains;
    vector<Triang> triangles;
    Vector3d rgb;
    int material;
    TriMesh(int material_type,double red,double green,double blue);
    void readMesh(string fname);
    int get_plain_num(){return plains.size();};
    int get_triangle_num(){return triangles.size();};
    int get_vertices_num(){return vertices.size();};
};


TriMesh::TriMesh(int material_type,double red,double green,double blue)
{
    rgb = Vector3d(red,green,blue);
    material = material_type;
}

void TriMesh::readMesh(string fname){
    ifstream infile;
    infile.open(fname);
    cout<<"Reading: "<<fname<<endl;
    if (!infile) cerr << "Could not open the file!" << endl;
    
    string line_str;
    int index_line = 0;
    
    while (getline(infile, line_str))
    {
        if (index_line==0) {
            if(line_str!="OFF"){
                cerr<<"Only off file allowed!"<<endl;
                return;
            }
        } else {
            istringstream iss(line_str);
            
            if (index_line==1) {
                int a1,a2,a3;
                if (!(iss >> a1 >> a2 >> a3)) {
                    cerr<<"off corruted!"<<endl;
                    break;
                } else {
                    //cout<<"num:"<<a1<<' '<<a2<<' '<<endl;
                }
            } else {
                vector<string> data_vec;
                string tmp;
                while (iss >> tmp) {
                    data_vec.push_back(tmp);
                }
                if (data_vec.size()==3) {
                    Vector3d vertex(stod(data_vec[0]),stod(data_vec[1]),stod(data_vec[2]));
                    vertices.push_back(vertex);
                } else if (data_vec.size()==4){
                    Vector3d plain(stoi(data_vec[1]),stoi(data_vec[2]),stoi(data_vec[3]));
                    plains.push_back(plain);
                    //cout<<data_vec[1]<<' '<<data_vec[2]<<' '<<data_vec[3];
                    Vector3d a = vertices[stoi(data_vec[1])];
                    Vector3d b = vertices[stoi(data_vec[2])];
                    Vector3d c = vertices[stoi(data_vec[3])];
                    Triang trg(material, a, b, c, rgb[0], rgb[1], rgb[2]);
                    triangles.push_back(trg);
                    //cout<<a<<'|'<<b<<'|'<<c<<' '<<endl;
                } else{
                    
                    cerr<<"double/int size:"<<data_vec.size()<<endl;
                }
            }
        }
        index_line ++;
    }
    cout<<"success: "<<triangles.size()<<' '<<vertices.size()<<endl;
}





#endif /* obj_h */
