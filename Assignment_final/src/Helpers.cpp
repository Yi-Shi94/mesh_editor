#include "Helpers.h"

#include <iostream>
#include <fstream>
#include <algorithm>
//#include <math.h>
using namespace std;
using namespace Eigen;

vector<Vector3d> bumpy_vertices;
//vector<Vector3d> bumpy_faces;
vector<Vector3d> bumpy_tris;
vector<Vector3d> bumpy_faces;
Vector3d bumpy_bary(0.0292562,0.000497742,0.00340984);
//Vector3d bumpy_bary(0,0,0);

vector<Vector3d> bunny_vertices;
//vector<Vector3d> bunny_faces;
vector<Vector3d> bunny_tris;
vector<Vector3d> bunny_faces;
Vector3d bunny_bary(-0.0281731,0.0941791,0.00829992);
//Vector3d bunny_bary(0,0,0);

vector<Vector3d> cube_vertices;
//vector<Vector3d> cube_faces;
vector<Vector3d> cube_tris;
vector<Vector3d> cube_faces;
Vector3d cube_bary(0,0,0);

vector<Vector3d> bary;
Vector3d ntris(12,1000,1000);

void VertexArrayObject::init()
{
  glGenVertexArrays(1, &id);
  check_gl_error();
}

void VertexArrayObject::bind()
{
  glBindVertexArray(id);
  check_gl_error();
}

void VertexArrayObject::free()
{
  glDeleteVertexArrays(1, &id);
  check_gl_error();
}

void VertexBufferObject::init()
{
  glGenBuffers(1,&id);
  check_gl_error();
}

void VertexBufferObject::bind()
{
  glBindBuffer(GL_ARRAY_BUFFER,id);
  check_gl_error();
}

void VertexBufferObject::free()
{
  glDeleteBuffers(1,&id);
  check_gl_error();
}

void VertexBufferObject::update(const Eigen::MatrixXf& M)
{
  assert(id != 0);
  glBindBuffer(GL_ARRAY_BUFFER, id);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*M.size(), M.data(), GL_DYNAMIC_DRAW);
  rows = M.rows();
  cols = M.cols();
  check_gl_error();
}

bool Program::init(
  const std::string &vertex_shader_string,
  const std::string &fragment_shader_string,
  const std::string &fragment_data_name)
{
  using namespace std;
  vertex_shader = create_shader_helper(GL_VERTEX_SHADER, vertex_shader_string);
  fragment_shader = create_shader_helper(GL_FRAGMENT_SHADER, fragment_shader_string);

  if (!vertex_shader || !fragment_shader)
    return false;

  program_shader = glCreateProgram();

  glAttachShader(program_shader, vertex_shader);
  glAttachShader(program_shader, fragment_shader);

  glBindFragDataLocation(program_shader, 0, fragment_data_name.c_str());
  glLinkProgram(program_shader);

  GLint status;
  glGetProgramiv(program_shader, GL_LINK_STATUS, &status);

  if (status != GL_TRUE)
  {
    char buffer[512];
    glGetProgramInfoLog(program_shader, 512, NULL, buffer);
    cerr << "Linker error: " << endl << buffer << endl;
    program_shader = 0;
    return false;
  }

  check_gl_error();
  return true;
}

void Program::bind()
{
  glUseProgram(program_shader);
  check_gl_error();
}

GLint Program::attrib(const std::string &name) const
{
  return glGetAttribLocation(program_shader, name.c_str());
}

GLint Program::uniform(const std::string &name) const
{
  return glGetUniformLocation(program_shader, name.c_str());
}

GLint Program::bindVertexAttribArray(
        const std::string &name, VertexBufferObject& VBO) const
{
  GLint id = attrib(name);
  if (id < 0)
    return id;
  if (VBO.id == 0)
  {
    glDisableVertexAttribArray(id);
    return id;
  }
  VBO.bind();
  glEnableVertexAttribArray(id);
  glVertexAttribPointer(id, VBO.rows, GL_FLOAT, GL_FALSE, 0, 0);
  check_gl_error();

  return id;
}

void Program::free()
{
  if (program_shader)
  {
    glDeleteProgram(program_shader);
    program_shader = 0;
  }
  if (vertex_shader)
  {
    glDeleteShader(vertex_shader);
    vertex_shader = 0;
  }
  if (fragment_shader)
  {
    glDeleteShader(fragment_shader);
    fragment_shader = 0;
  }
  check_gl_error();
}

GLuint Program::create_shader_helper(GLint type, const std::string &shader_string)
{
  using namespace std;
  if (shader_string.empty())
    return (GLuint) 0;

  GLuint id = glCreateShader(type);
  const char *shader_string_const = shader_string.c_str();
  glShaderSource(id, 1, &shader_string_const, NULL);
  glCompileShader(id);

  GLint status;
  glGetShaderiv(id, GL_COMPILE_STATUS, &status);

  if (status != GL_TRUE)
  {
    char buffer[512];
    if (type == GL_VERTEX_SHADER)
      cerr << "Vertex shader:" << endl;
    else if (type == GL_FRAGMENT_SHADER)
      cerr << "Fragment shader:" << endl;
    else if (type == GL_GEOMETRY_SHADER)
      cerr << "Geometry shader:" << endl;
    cerr << shader_string << endl << endl;
    glGetShaderInfoLog(id, 512, NULL, buffer);
    cerr << "Error: " << endl << buffer << endl;
    return (GLuint) 0;
  }
  check_gl_error();

  return id;
}

void _check_gl_error(const char *file, int line)
{
  GLenum err (glGetError());

  while(err!=GL_NO_ERROR)
  {
    std::string error;

    switch(err)
    {
      case GL_INVALID_OPERATION:      error="INVALID_OPERATION";      break;
      case GL_INVALID_ENUM:           error="INVALID_ENUM";           break;
      case GL_INVALID_VALUE:          error="INVALID_VALUE";          break;
      case GL_OUT_OF_MEMORY:          error="OUT_OF_MEMORY";          break;
      case GL_INVALID_FRAMEBUFFER_OPERATION:  error="INVALID_FRAMEBUFFER_OPERATION";  break;
    }
    std::cerr << "GL_" << error.c_str() << " - " << file << ":" << line << std::endl;
    err = glGetError();
  }
}

MatrixXf ini_color_vertex(MatrixXf V_Shape){
    MatrixXf V_Color(3,V_Shape.size());
    for (int i=0;i < V_Shape.cols();i++){
        V_Color.col(i)<< 0.1f, 0.1f, 0.9f;
    }
    return V_Color;
}


void readMesh(string fname,int type){
    ifstream infile;
    infile.open(fname);
    if (!infile) cerr << "Could not open the file!" << endl;
    string line_str;
    int index_line = 0;
    
    vector<Vector3d> vertices;
    vector<Vector3d> face_vertex;
    vector<Vector3d> face_index;
    
    while (getline(infile, line_str))
    {
        if (index_line==0) {
            if(line_str!="OFF"){
                cerr<<"Only off file allowed!"<<endl;
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
                    Vector3d vertex(stod(data_vec[0]),
                                    stod(data_vec[1]),
                                    stod(data_vec[2]));
                    vertices.push_back(vertex-bary[type]);
                } else if (data_vec.size()==4){
                    Vector3d plain(stoi(data_vec[1]),
                                   stoi(data_vec[2]),
                                   stoi(data_vec[3]));
                    //cout<<data_vec[1]<<' '<<data_vec[2]<<' '<<data_vec[3];
                    
                    face_index.push_back(plain);
                    
                    Vector3d a = vertices[stoi(data_vec[1])];
                    Vector3d b = vertices[stoi(data_vec[2])];
                    Vector3d c = vertices[stoi(data_vec[3])];
                    face_vertex.push_back(a);
                    face_vertex.push_back(b);
                    face_vertex.push_back(c);
                } else{
                    cerr<<"double/int size:"<<data_vec.size()<<endl;
                }
            }
        }
        if (type==0) {
            cube_vertices = vertices;
            cube_faces = face_index;
            cube_tris = face_vertex;
        }
        else if (type==1) {
            bunny_vertices = vertices;
            bunny_faces = face_index;
            bunny_tris = face_vertex;
        }
        else {
            bumpy_vertices = vertices;
            bumpy_faces = face_index;
            bumpy_tris = face_vertex;
        }
        index_line ++;
    }
    cout<<"loading success! "<<fname<<' '<<vertices.size()<<' '<<face_vertex.size()<<endl;
    //return ;
}

void load_meshes(){
    srand(static_cast<unsigned>(time(0)));
    //load_cube();
    bary.push_back(cube_bary);
    bary.push_back(bunny_bary);
    bary.push_back(bumpy_bary);
    readMesh("../data/cube.off",0);
    readMesh("../data/bunny.off",1);
    readMesh("../data/bumpy_cube.off",2);
}

vector<Vector3d> TriMesh::load_vertices(int i){
    if(i==0) return cube_vertices;
    else if (i==1) return bunny_vertices;
    else return bumpy_vertices;
}

vector<Vector3d> TriMesh::load_tris(int i){
    if(i==0) return cube_tris;
    else if (i==1) return bunny_tris;
    else return bumpy_tris;
}

vector<Vector3d> TriMesh::load_faces(int i){
    if(i==0) return cube_faces;
    else if (i==1) return bunny_faces;
    else return bumpy_faces;
}

vector<Vector3d> TriMesh::get_vertices(){
    return this->vertices;
}

vector<Vector3d> TriMesh::get_tris(){
    return this->tris;
}

vector<Vector3d> TriMesh::get_faces(){
    return this->faces;
}


Vector3d replace_ind_in_tri(int replace, int target, Vector3d tri){
    Vector3d res = tri;
    for(int i=0;i<3;i++){
        if(target==int(res(i))){
            res(i)=replace;
        }
    }
    return res;
}


bool if_has_vert_ind(int target, Vector3d tri){
    for(int i=0;i<3;i++){
        if(int(tri(i))==target){
            return true;
        }
    }
    return false;
}

bool if_has_vec_ind(int target, vector<int> vec){
    if (vec.size()==0){
        //cout<<"empty vec"<<endl;
        return false;
    }
    for(int i=0;i<vec.size();i++){
        if(vec[i]==target){
            return true;
        }
    }
    return false;

}

int get_pivot(int vertex_index,Vector3d tri){
    for(int i=0;i<3;i++){
        if(int(tri(i))!=vertex_index){
            return int(tri(i));
        }
    }
    return -1;
}

vector<int> insert_neighbor_vert(vector<int> before, int a, int b){
    vector<int> after = before;
    if(!if_has_vec_ind(a,before)) after.push_back(a);
    if(!if_has_vec_ind(b,after)) after.push_back(b);
    return after;
}


vector<vector<int>> generate_neighbor_vert(vector<Vector3d> faces, int vert_num){
    vector<vector<int>> neighbors;
    for (int i=0;i<vert_num;i++){
        vector<int> tmp;
        neighbors.push_back(tmp);
    }
    
    for (int i=0;i<faces.size();i++){
        int a = int(faces[i](0));
        int b = int(faces[i](1));
        int c = int(faces[i](2));
        vector<int> neib_a = neighbors[a];
        vector<int> neib_b = neighbors[b];
        vector<int> neib_c = neighbors[c];
        neighbors[a] = insert_neighbor_vert(neib_a,b,c);
        neighbors[b] = insert_neighbor_vert(neib_b,a,c);
        neighbors[c] = insert_neighbor_vert(neib_c,a,b);
        
    }
    return neighbors;
}

vector<vector<int>> generate_neighbor_face(vector<Vector3d> faces, int vert_num){
    vector<vector<int>> neighbors;
    for (int i=0;i<vert_num;i++){
        vector<int> tmp;
        neighbors.push_back(tmp);
    }
    
    for (int i=0;i<faces.size();i++){
        int a = int(faces[i](0));
        int b = int(faces[i](1));
        int c = int(faces[i](2));
        neighbors[a].push_back(i);
        neighbors[b].push_back(i);
        neighbors[c].push_back(i);
    }
    return neighbors;
}


void print_vec_vec(vector<vector<int>> tar){
    for (int i=0;i<tar.size();i++){
        cout<<i<<"-> ";
        for(int j=0;j<tar[i].size();j++){
            cout<<tar[i][j]<<' ';
        }
        cout<<endl;
    }
}
TriMesh::TriMesh(int target_type,int start){
    
    this->type = target_type;
    this->vertices = load_vertices(this->type);
    this->faces = load_faces(this->type);
    this->tris = load_tris(this->type);
    this->vert_neigbors = generate_neighbor_vert(this->faces,this->vertices.size());
    this->face_neigbors = generate_neighbor_face(this->faces,this->vertices.size());
    print_vec_vec(this->vert_neigbors);
    cout<<"---------------------------"<<endl;
    print_vec_vec(this->face_neigbors);
    compute_bary_n_max(target_type);
    if(target_type == 0){
        set_trans_mat(0,0,0,0,0,0,1);
    }else if(target_type == 1){
        set_trans_mat(0,0,0,0,0,0,1/0.156141);
    }else if (target_type == 2){
        set_trans_mat(0,0,0,0,0,0,1/8.75694);
    }
    Vector3d ini_color(0.1,0.5,0.5);
    
    
    this->render_type = -1;
    this->set_color(ini_color);
    this->start = start;
    this->get_bary();
    this->tri_num = int(ntris(target_type));
    this->cal_normal_matrix();
    this->cal_phong_normal_matrix();
}


MatrixXf TriMesh::populate_vertex(){
    MatrixXf out(3,cube_tris.size()*this->selected_vertex_index.size());
    double scale_adj = 1;
    if (this->type==0) scale_adj = 40.0;
    else if (this->type==1) scale_adj = 200.0;
    else scale_adj = 5.0;
    for (int i=0;i<this->selected_vertex_index.size();i++){
        Vector3d bary_indi = this->vertices[this->selected_vertex_index[i]];
        for (int j=0;j<cube_tris.size();j++){
            Vector3d ver_cube =
                cube_tris[j]/scale_adj+bary_indi;
            out.col(i*cube_tris.size()+j)<<ver_cube(0),ver_cube(1),ver_cube(2);
        }
    }
    //cout<<'m'<<out.cols()<<endl;
    return out;
}

MatrixXf TriMesh::populate_normal_vertex(){
    MatrixXf out(3,cube_tris.size()*this->selected_vertex_index.size());
    for (int i=0;i<this->selected_vertex_index.size();i++){
        Vector3d bary_indi = this->vertices[this->selected_vertex_index[i]];
        for (int j=0;j<cube_tris.size();j++){
            out.col(i*cube_tris.size()+j)<<0,0,1;
        }
    }
    return out;
}

vector<float> world_2_screen(Matrix4f trans, Vector3d xyz){
    Vector4f xyz1(xyz(0),xyz(1),xyz(2),1);
    Vector4f uv11= trans*xyz1;
    vector<float> res;
    res.push_back(uv11(0));
    res.push_back(uv11(1));
    return res;
    
}

void TriMesh::select_vertex_2(float x,float y,MatrixXf trans){
    int index_vertex = -1;
    int face_hit_index = this->selected_face;
    Matrix4f model = this->get_trans_mat();
    Vector3d vert_indeces = this->faces[face_hit_index];
    
    Vector3d vert_a = this->vertices[int(vert_indeces(0))];
    Vector3d vert_b = this->vertices[int(vert_indeces(1))];
    Vector3d vert_c = this->vertices[int(vert_indeces(2))];
    
    vector<float> uv_a = world_2_screen(trans*model,vert_a);
    vector<float> uv_b = world_2_screen(trans*model,vert_b);
    vector<float> uv_c = world_2_screen(trans*model,vert_c);
    
    double dis_a = abs(uv_a[0]-x)+abs(uv_a[1]-y);
    double dis_b = abs(uv_b[0]-x)+abs(uv_b[1]-y);
    double dis_c = abs(uv_c[0]-x)+abs(uv_c[1]-y);
    
    if (dis_a>dis_b){
        if(dis_b>dis_c){
            index_vertex = int(vert_indeces(2));
        } else{
            index_vertex = int(vert_indeces(1));
        }
    } else {
        if(dis_a>dis_c){
            index_vertex = int(vert_indeces(2));
        } else{
            index_vertex = int(vert_indeces(0));
        }
    }
    
    cout<<"selected vert ind: "<<index_vertex<<endl;
    if (index_vertex>=0){
        if(find(this->selected_vertex_index.begin(),        this->selected_vertex_index.end(),index_vertex) ==
           this->selected_vertex_index.end())
            this->selected_vertex_index.push_back(index_vertex);
    }
    
}

void TriMesh::select_vertex(Vector3d ray_origin,Vector3d ray_direction,double t){
    
    Vector3d target = ray_origin+t*ray_direction;
    int index_vertex = -1;
    int face_hit_index = this->selected_face;
    Vector3d vert_indeces = this->faces[face_hit_index];
    
    int vert_a_ind = int(vert_indeces(0));
    int vert_b_ind = int(vert_indeces(1));
    int vert_c_ind = int(vert_indeces(2));
    
    //cout<<vert_a_ind<<' '<<vert_b_ind<<' '<<vert_c_ind<<endl;
    
    //cout<<"vert: "<<this->vertices.size()<<endl;
    //cout<<"face: "<<this->faces.size()<<endl;
    //cout<<"tris: "<<this->tris.size()<<endl;
    
    Vector3d vert_a = this->vertices[vert_a_ind];
    Vector3d vert_b = this->vertices[vert_b_ind];
    Vector3d vert_c = this->vertices[vert_c_ind];
    
    double dis_a = (vert_a-target).norm();
    double dis_b = (vert_b-target).norm();
    double dis_c = (vert_c-target).norm();
    
    if (dis_a>dis_b){
        if(dis_b>dis_c){
            index_vertex = vert_c_ind;
        } else{
            index_vertex = vert_b_ind;
        }
    } else {
        if(dis_a>dis_c){
            index_vertex = vert_c_ind;
        } else{
            index_vertex = vert_a_ind;
        }
    }
    cout<<"selected vert ind: "<<index_vertex<<endl;
    if (index_vertex>=0){
        if(find(this->selected_vertex_index.begin(),        this->selected_vertex_index.end(),index_vertex) ==
            this->selected_vertex_index.end())
            this->selected_vertex_index.push_back(index_vertex);
    }
}


void TriMesh::update_tris(){
    vector<Vector3d> face_vertex;
    for (int i=0;i<this->faces.size();i++){
        Vector3d plain = this->faces[i];
        face_vertex.push_back(this->vertices[plain(0)]);
        face_vertex.push_back(this->vertices[plain(1)]);
        face_vertex.push_back(this->vertices[plain(2)]);
    }
    
    this->tris = face_vertex;
    this->tri_num = this->faces.size();
    this->face_neigbors = generate_neighbor_face(this->faces,this->vertices.size());
    this->vert_neigbors = generate_neighbor_vert(this->faces,this->vertices.size());
    this->cal_normal_matrix();
    //this->cal_phong_normal_matrix();
}

void print_ev(vector<Vector3d> s,vector<Vector3d> s2){
    cout<<"vertices "<<s.size()<<endl;
    for (int i=0;i<s.size();i++){
        cout<<s[i](0)<<' '<<s[i](1)<<' '<<s[i](2)<<' '<<endl;
    }
    cout<<"faces "<<s2.size()<<endl;
    for (int i=0;i<s2.size();i++){
        cout<<s2[i](0)<<' '<<s2[i](1)<<' '<<s2[i](2)<<' '<<endl;
    }
}
vector<vector<int>> insert_vec_in_vec(vector<int> edge,vector<vector<int>>edges){
    int flag = 0;
    for(int i=0;i<edges.size();i++){
        if((edge[0]==edges[i][0]&&edge[1]==edges[i][1])
           ||(edge[1]==edges[i][0]&&edge[0]==edges[i][1]))
            flag =1;
            break;
    }
    if(flag==0){
        vector<vector<int>> new_edges = edges;
        new_edges.push_back(edge);
        return new_edges;
    }
    return edges;
}

bool if_vec_in_vec(int a,int b,vector<vector<int>>edges){
    if(edges.size()==0){
        return false;
    }
    int flag = 0;
    for(int i=0;i<edges.size();i++){
        if((a==edges[i][0]&&b==edges[i][1])
           ||(b==edges[i][0]&&a==edges[i][1]))
            flag =1;
        break;
    }
    if(flag==0){
        return false;
    }
    return true;
}

vector<vector<int>> get_edges(vector<Vector3d> tris){
    vector<vector<int>>edges;
    for (int i=0;i<tris.size();i++){
        int a = int(tris[i](0));
        int b = int(tris[i](1));
        int c = int(tris[i](2));
        if(!if_vec_in_vec(a,b,edges)) {
            vector<int> edge_ab;
            edge_ab.push_back(a);
            edge_ab.push_back(b);
            edges.push_back(edge_ab);
        }
        if(!if_vec_in_vec(b,c,edges)) {
            vector<int> edge_bc;
            edge_bc.push_back(b);
            edge_bc.push_back(c);
            edges.push_back(edge_bc);
        }
        if(!if_vec_in_vec(c,a,edges)) {
            vector<int> edge_ca;
            edge_ca.push_back(c);
            edge_ca.push_back(a);
            edges.push_back(edge_ca);
        }
    }
    return edges;
}


int num_edge_vec3d(vector<int> edge,vector<Vector3d> tri_to_be_inspected){
    int counter = 0;
    for(int i=0;i<tri_to_be_inspected.size();i++){
        int a = int(tri_to_be_inspected[i](0));
        int b = int(tri_to_be_inspected[i](1));
        int c = int(tri_to_be_inspected[i](2));
        if(edge[0]==a){
            if(edge[1]==b||edge[1]==c){
                counter+=1;
            }
        } else if (edge[0]==b){
            if(edge[1]==a||edge[1]==c){
                counter+=1;
            }
        } else if (edge[0]==c){
            if(edge[1]==a||edge[1]==b){
                counter+=1;
            }
        }
    }
    return counter;
}


bool check_local(vector<Vector3d> tris_candidate,vector<Vector3d> tris_all){
    return true;
    vector<vector<int>> edges = get_edges(tris_candidate);
    for (int i=0;i<edges.size();i++){
        int num = num_edge_vec3d(edges[i],tris_all);
        if(num!=1) {
            cout<<num<<endl;
            return false;
        }
    }
    return true;
}



bool TriMesh::delete_vertex(){
    if(this->selected_vertex_index.size()!=1){
        cout<<"only support deleting one at a time."<<this->selected_vertex_index.size()<<endl;
        return false;
    }
    int vertex_index = this->selected_vertex_index[0];
    
    vector<Vector3d> tri_to_be_inspected;
    vector<int> ind_tri_to_delete;
    
    vector<int> neighbors_of_target = this->vert_neigbors[vertex_index];
    vector<int> faces_candidates = this->face_neigbors[vertex_index];
    
    int pivot_index = 0;
    int pivot = neighbors_of_target[pivot_index];
    bool checked=false;
    
    
        for (int i=0;i<this->faces.size();i++){
            Vector3d cur_face = this->faces[i];
            if (if_has_vert_ind(vertex_index,cur_face)) {
                if(if_has_vert_ind(pivot,cur_face)){
                    ind_tri_to_delete.push_back(i);
                }
                else{
                    this->faces[i] = replace_ind_in_tri(pivot,vertex_index,cur_face);
                    tri_to_be_inspected.push_back(this->faces[i]);
                }
            }
            
        }
    
    
   // print_ev(this->vertices,this->faces);
    vector<Vector3d> face_tmp;
    for(int j=0;j<this->faces.size();j++){
        bool flag = true;
        for (int i=0;i<ind_tri_to_delete.size();i++){
            if(j==ind_tri_to_delete[i]){
                flag=false;
               //break;
            }
        }
        if (flag==true){
            face_tmp.push_back(this->faces[j]);
        }
    }
    
    this->faces = face_tmp;
    this->vertices.erase(this->vertices.begin()+vertex_index);
    
    for (int i=0;i<this->faces.size();i++){
        for(int j=0;j<3;j++){
            if(this->faces[i](j)>vertex_index){
                this->faces[i](j)-=1;
            }
        }
    }
    //print_ev(this->vertices,this->faces);
    cout<<"deleted vertexi index: "<<this->selected_vertex_index[0]<<endl;
    this->selected_vertex_index.clear();
    this->update_tris();
    return true;
}

bool TriMesh::translate_vertex(int axis,double amount){
    if (this->selected_vertex_index.size()==1){
        int vertex_index = this->selected_vertex_index[0];
        this->vertices[vertex_index](axis) += amount;
        if(!this->check()){
            this->vertices[vertex_index](axis) -= amount;
            return false;
        } else {
            this->update_tris();
            return true;
        }
    }
    return false;
}

vector<int> TriMesh::is_edge(int vert_ind_a, int vert_ind_b){
    vector<int> tri_with_edge;
    for (int i=0;i<this->faces.size();i++){
        if(if_has_vert_ind(vert_ind_a,this->faces[i]) &&
           if_has_vert_ind(vert_ind_b,this->faces[i]))
            tri_with_edge.push_back(i);
    }
    return tri_with_edge;
}

float random_interpo_float(float a,float b){
    float LO,HI;
    if (a<b){
        LO = a;
        HI = b;
    } else {
        LO = b;
        HI = a;
    }
    float r3 = LO+static_cast<float>(rand())/(static_cast<float>(RAND_MAX/(HI-LO)));
    return r3;
}

Vector3d random_interpo_vec3d(Vector3d tri_ver_edge_a,Vector3d tri_ver_edge_b){
    Vector3d res(0,0,0);
    for(int i=0;i<3;i++)
        res(i) = random_interpo_float(tri_ver_edge_a(i),tri_ver_edge_b(i));
    return res;
}

bool TriMesh::merge_vertex(){
    if(this->selected_vertex_index.size()!=2){
        return false;
    }
    sort(this->selected_vertex_index.begin(),
              this->selected_vertex_index.begin()+2);
    
    vector<int> tri_ind_edge = this->is_edge(this->selected_vertex_index[0],
                                              this->selected_vertex_index[1]);
    if (tri_ind_edge.size()==0){
        return false;
    }
    
    Vector3d tri_ver_edge_a = this->vertices[this->selected_vertex_index[0]];
    Vector3d tri_ver_edge_b = this->vertices[this->selected_vertex_index[1]];
    Vector3d tri_ver_edge_interpo = random_interpo_vec3d(tri_ver_edge_a,tri_ver_edge_b);
    
    vector<Vector3d> face_tmp;
    for(int j=0;j<this->faces.size();j++){
        bool flag = true;
        for (int i=0;i<tri_ind_edge.size();i++){
            if(j==tri_ind_edge[i]){
                flag=false;
                break;
            }
        }
        if (flag==true){
            face_tmp.push_back(this->faces[j]);
        }
    }
    
    this->vertices[this->selected_vertex_index[0]] = tri_ver_edge_interpo;
    this->vertices.erase(this->vertices.begin()+this->selected_vertex_index[1]);
    
    for (int i=0;i<face_tmp.size();i++){
        Vector3d cur_face = face_tmp[i];
        int vertex_index_0 = this->selected_vertex_index[0];
        int vertex_index_1 = this->selected_vertex_index[1];
        if(if_has_vert_ind(vertex_index_1,cur_face)){
            face_tmp[i] = replace_ind_in_tri(vertex_index_0,vertex_index_1,cur_face);
        }
        
        for(int j=0;j<3;j++){
            if(face_tmp[i](j)>vertex_index_1){
                face_tmp[i](j)-=1;
            }
        }
    }
        
    this->faces = face_tmp;
    this->selected_vertex_index.clear();
    this->update_tris();
    return true;
}

bool TriMesh::check(){
    return true;
}

/*
bool TriMesh::check_local(int vertex_index_1,int vertex_index_2){
    int count = 0;
    for (int i=0;i<face_tmp.size();i++){
        if(if_has_vert_ind(vertex_index_1,cur_face) && if_has_vert_ind(vertex_index_2,cur_face)){
            count+=1;
        }
    }
    if (count!=2) return false;
    else return true;
    
}
*/
int TriMesh::save_off(){
    ofstream outfile ("/Users/syd/Documents/CG/Assignment_3/data/mesh_out.off",ofstream::binary);
    // allocate memory for file content
    outfile << "OFF\n";
    outfile << to_string(this->vertices.size())<<
            ' '<<to_string(this->faces.size())<<
            ' '<<'0'<<"\n";
    for (const auto &e : this->vertices){
        outfile << to_string(e(0))<<
            ' '<<to_string(e(1))<<
            ' '<<to_string(e(2))<<
            "\n";
    }
    for (const auto &f : this->faces){
        outfile <<"3 "<<to_string(int(f(0)))<<
                  ' '<<to_string(int(f(1)))<<
                  ' '<<to_string(int(f(2)))<<
                    "\n";
    }
    cout<<"save complete!"<<endl;
    outfile.close();
    return 0;
}

void TriMesh::get_bary(){
    barycenter = bary[this->type]+Vector3d(x,y,z);
    //cout<<barycenter<<endl;
}

void TriMesh::compute_bary_n_max(int type){
     vector<Vector3d> vertices = this->get_tris();

     double x=0,y=0,z=0;
     for(int i=0;i<vertices.size();i++){
         x+=vertices[i](0);
         y+=vertices[i](1);
         z+=vertices[i](2);
     }
     //cout<<x_max<<" x "<<x_min<<' '<<x_max-x_min <<endl;
     //cout<<y_max<<" y "<<y_min<<' '<<y_max-y_min <<endl;
     //cout<<z_max<<" z "<<z_min<<' '<<z_max-z_min <<endl;
     x=x/vertices.size();
     y=y/vertices.size();
     z=z/vertices.size();
     Vector3d res(x,y,z);
     barycenter=res;
     cout<<type<<" bary "<<res<<endl;
}

MatrixXf TriMesh::get_matrix(){
    vector<Vector3d> vertices = this->tris;
    MatrixXf v_matrix(3,vertices.size());
    for(int i=0;i<vertices.size();i++){
        v_matrix.col(i)<<vertices[i](0),vertices[i](1),vertices[i](2);
    }
    return v_matrix;
}

MatrixXf TriMesh::get_phong_normal_matrix(){
    {
        MatrixXf v_matrix(3,tri_num*3);
        for(int i=0;i<normals_phong.size();i=i+1){
            Vector3d normal = normals_phong[i];
            v_matrix.col(i)<<normal(0),normal(1),normal(2);
        }
        return v_matrix;
    }
}

void TriMesh:: cal_phong_normal_matrix(){
    vector<Vector3d> face_index = this->get_faces();
    vector<Vector3d> sum_vec;
    vector<int> num_vec;
    
    for (int i=0;i<normals.size();i++){
        Vector3d ini(0,0,0);
        sum_vec.push_back(ini);
        num_vec.push_back(0);
    }
    
    for (int i=0;i<face_index.size();i++){
        sum_vec[face_index[i](0)]+=normals[3*i];
        num_vec[face_index[i](0)]++;
        
        sum_vec[face_index[i](1)]+=normals[3*i+1];
        num_vec[face_index[i](1)]++;
        
        sum_vec[face_index[i](2)]+=normals[3*i+2];
        num_vec[face_index[i](2)]++;
    }
    
    for(int i=0;i<face_index.size();i++){
        normals_phong.push_back(sum_vec[face_index[i](0)]
                                /float(num_vec[face_index[i](0)]));
        normals_phong.push_back(sum_vec[face_index[i](1)]
                                /float(num_vec[face_index[i](1)]));
        normals_phong.push_back(sum_vec[face_index[i](2)]
                                /float(num_vec[face_index[i](2)]));
    }
}

MatrixXf TriMesh::get_normal_matrix()
{
    MatrixXf v_matrix(3,tri_num*3);
    for(int i=0;i<this->normals.size();i=i+3){
        Vector3d normal = normals[i];
        v_matrix.col(i)  <<normal(0),normal(1),normal(2);
        v_matrix.col(i+1)<<normal(0),normal(1),normal(2);
        v_matrix.col(i+2)<<normal(0),normal(1),normal(2);
    }
     return v_matrix;
}

void TriMesh::cal_normal_matrix(){
    vector<Vector3d> normal_per_ver;
    Matrix4f trans = this->get_trans_mat();
    //cout<<tri_num<<' '<<vertices.size()<<endl;
    for(int i=0;i<this->tris.size();i=i+3){
        MatrixXf points(4,3);
        points.col(0) << this->tris[i](0),this->tris[i](1),this->tris[i](2),1;
        points.col(1) << this->tris[i+1](0),this->tris[i+1](1),this->tris[i+1](2),1;
        points.col(2) << this->tris[i+2](0),this->tris[i+2](1),this->tris[i+2](2),1;
        points = trans*points;
        
        Vector3d x1(points(0,0),points(1,0),points(2,0));
        Vector3d x2(points(0,1),points(1,1),points(2,1));
        Vector3d x3(points(0,2),points(1,2),points(2,2));
        
        Vector3d a = x2 - x1;
        Vector3d b = x3 - x1;
        Vector3d normal = a.cross(b).normalized();
        normal_per_ver.push_back(normal);
        normal_per_ver.push_back(normal);
        normal_per_ver.push_back(normal);
    }
    this->normals = normal_per_ver;
}

MatrixXf TriMesh::get_trans_mat(){
    
    MatrixXf rot_x(4,4);
    MatrixXf rot_y(4,4);
    MatrixXf rot_z(4,4);
    MatrixXf scal(4,4);
    MatrixXf trans(4,4);
    //MatrixXf trans_final(4,4);
    
    float rad_x = angle_x*3.1415926535/180.;
    float rad_y = angle_y*3.1415926535/180.;
    float rad_z = angle_z*3.1415926535/180.;
    
    float cos_x = cos(rad_x);
    float cos_y = cos(rad_y);
    float cos_z = cos(rad_z);
    
    float sin_x = sin(rad_x);
    float sin_y = sin(rad_y);
    float sin_z = sin(rad_z);
    float s = scale;
    scal<<s,0,0,0,
          0,s,0,0,
          0,0,s,0,
          0,0,0,1;
    
    trans<< 1,0,0,x,
            0,1,0,y,
            0,0,1,z,
            0,0,0,1;
    
    rot_x<<  1,0,0,0,
            0,cos_x,-sin_x,0,
            0,sin_x,cos_x,0,
            0,0,0,1;
    
    rot_y<<cos_y,    0,     sin_y,    0,
           0,        1,     0,        0,
          -sin_y,    0,     cos_y,    0,
           0,        0,     0,        1.;
    
    rot_z<<cos_z,-sin_z,0,0,
          sin_z,cos_z,  0,0,
          0,    0,      1,0,
          0,    0,      0,1;
    
    //trans_final<< 1,0,0,barycenter(0),
    //            0,1,0,barycenter(1),
     //           0,0,1,barycenter(2),
     //           0,0,0,1;
    
    //compute_bary_n_max(type);
    
    return trans*rot_z*rot_y*rot_x*scal;
}

void TriMesh::set_trans_mat(float x,float y, float z, float angle_x, float angle_y, float angle_z,float s){
    this->x = x;
    this->y = y;
    this->z = z;
    this->angle_x = angle_x;
    this->angle_y = angle_y;
    this->angle_z = angle_z;
    this->scale = s;
}


double TriMesh::is_hit(Vector3d ray_origin, Vector3d ray_direction)
{
    vector<Vector3d> tris = this->get_tris();
    Matrix4f trans = this->get_trans_mat();
    double t_min_valid = 1000000000.0;
    int face_min_valid = -1;
    
    for (int i=0; i<tris.size(); i=i+3){
        MatrixXf points(4,3);
        points.col(0) << tris[i](0),tris[i](1),tris[i](2),1;
        points.col(1) << tris[i+1](0),tris[i+1](1),tris[i+1](2),1;
        points.col(2) << tris[i+2](0),tris[i+2](1),tris[i+2](2),1;
        MatrixXf res = trans*points;
        
        double a_x = res(0,0); double a_y = res(1,0); double a_z = res(2,0);
        double b_x = res(0,1); double b_y = res(1,1); double b_z = res(2,1);
        double c_x = res(0,2); double c_y = res(1,2); double c_z = res(2,2);
        double e_x = ray_origin[0];
        double e_y = ray_origin[1];
        double e_z = ray_origin[2];
        
        double d_x = ray_direction[0];
        double d_y = ray_direction[1];
        double d_z = ray_direction[2];
    
        double A11 = a_x-b_x; double A12 = a_x-c_x;
        double A21 = a_y-b_y; double A22 = a_y-c_y;
        double A31 = a_z-b_z; double A32 = a_z-c_z;
        double Xae = a_x-e_x; double Yae = a_y-e_y; double Zae = a_z-e_z;
        double Axd_sub = A21*A32-A22*A31;
        double Ayd_sub = A11*A32-A12*A31;
        double Azd_sub = A11*A22-A21*A12;
    
        double A = d_x*Axd_sub - d_y*Ayd_sub + d_z*Azd_sub;
        double t = (Xae*Axd_sub - Yae*Ayd_sub +  Zae*Azd_sub)/A;
        if (t<=0){
            continue;
        }
    
        double gamma = (d_x * (A21*Zae-Yae*A31) - d_y * (A11*Zae-Xae*A31) + d_z * (A11*Yae-Xae*A21))/A;
        if (gamma<0 || gamma>1){
            continue;
        }
    
        double beta = (d_x * (Yae*A32-Zae*A22)-d_y * (Xae*A32-Zae*A12) + d_z * (Xae*A22-Yae*A12))/A;
        if(beta<0 || beta>1-gamma){
            continue;
        }
        if (t_min_valid>t){
            t_min_valid = t;
            face_min_valid = i;
        }
    }
    if(face_min_valid<0){
        return -1;
    } else {
        this->selected_face = face_min_valid/3;
        return t_min_valid;
    }
}


Vector3d TriMesh::get_color(){
    return color;
}

void TriMesh::set_color(Vector3d c){
    color = c;
}


void TriMesh::set_render_type(int t){
    render_type = t;
}


int TriMesh::get_render_type(){
    return render_type;
}

