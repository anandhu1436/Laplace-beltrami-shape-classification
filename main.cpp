
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/eigs.h>
#include <igl/vertex_triangle_adjacency.h>
#include <Eigen/Core>
#include <random>



#include <tgmath.h>
#include <math.h>
#include <fstream>
#include <filesystem>
#include <vector>




namespace fs = std::__fs::filesystem;
using namespace std;

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd N_vertices;
std::vector<std::vector<int>> neighbour_map;
string file;



float distance(Eigen::Vector3d a,Eigen::Vector3d b){
    return sqrt(pow(a.x()-b.x(),2)+pow(a.y()-b.y(),2)+pow(a.z()-b.z(),2));
}

void build_neighbour(){
  neighbour_map.clear();
  std::vector<std::vector<int>> VF;
  std::vector<std::vector<int>> VFi;
  igl::vertex_triangle_adjacency(V.rows(), F, VF, VFi);
  for (int i = 0; i < V.rows(); ++i) {
      std::vector<int> neigh;
      std::vector<double> dist;
      neigh.clear();
      dist.clear();
      for (int j = 0; j < VF[i].size(); ++j) {
          int f = VF[i][j];
          for (int k = 0; k < 3; ++k) {
              int v_idx = F(f, k);
              if (v_idx != i && std::find(neigh.begin(), neigh.end(), v_idx) == neigh.end()) {
                  neigh.push_back(v_idx);
              }
          }
      }

  neighbour_map.push_back(neigh);
  }
}

void add_noise(Eigen::MatrixXd& V, double noise_scale) {
    // Create a random number generator
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, noise_scale);

    // Add noise to each vertex
    for (int i = 0; i < V.rows(); i++) {
        Eigen::Vector3d noise(distribution(generator), distribution(generator), distribution(generator));
        V.row(i) += noise;
    }
}

void add_noise_along_normals(Eigen::MatrixXd& V, double noise_mag)
{
    // Compute the vertex normals if they are not given
    Eigen::MatrixXd N_vertices;
    igl::per_vertex_normals(V, F, N_vertices);

    // Perturb each vertex along its normal direction
    // for(int i = 0; i < V.rows(); i++) {
    //     Eigen::Vector3d noise = noise_mag * N_vertices.row(i);
    //     V.row(i) += noise;
    // }

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, noise_mag);
    for (int i = 0; i < V.rows(); i++) {
        // Eigen::Vector3d n = N_vertices[i];
        double magnitude = distribution(generator);
        V.row(i) += magnitude * N_vertices.row(i);
    }
}


void bilateral(Eigen::MatrixXd& V, Eigen::MatrixXi F) {
    Eigen::MatrixXd N_vertices;
    igl::per_vertex_normals(V, F, N_vertices);

    
    std::vector<std::vector<double>> dist_map;

    for(int num=0;num<1;++num){
      dist_map.clear();
      for (int i = 0; i < V.rows(); ++i) {
        std::vector<double> dist;
        dist.clear();
        for (int j = 0; j < neighbour_map[i].size(); ++j) {
              dist.push_back(distance(V.row(i),V.row(neighbour_map[i][j])));
                }
            
        dist_map.push_back(dist);
        }

    
    for (int i = 0; i < V.rows(); ++i) {
        Eigen::Vector3d n = N_vertices.row(i);
        double min_dist = *std::min_element(dist_map[i].begin(), dist_map[i].end());
        float sigma_c = min_dist;

        double average_offset = 0;
        int neighbour_count = 0;
        std::vector<double> offset_dist;
        offset_dist.clear();
        for (int j = 0; j < neighbour_map[i].size(); ++j) {
            float offset_d=0;
            if (dist_map[i][j] < 2 * sigma_c) {
              Eigen::Vector3d vec = V.row(i) - V.row(neighbour_map[i][j]);
                // float t = dist_map[i][j];
                float d = n.x() * vec.x() + n.y() * vec.y() + n.z() * vec.z();
                offset_d=sqrt(d*d);
                average_offset += offset_d;
                neighbour_count++;
            }
            offset_dist.push_back(offset_d);
        }
        if (neighbour_count > 0) {
            average_offset /= neighbour_count;
        }

        float off_set_dis = 0;
        for (int j = 0; j < neighbour_map[i].size(); ++j) {
            if (dist_map[i][j] < 2 * sigma_c) {
                off_set_dis += (offset_dist[j] - average_offset) * (offset_dist[j] - average_offset);
            }
        }
        float sigma_s ;
        off_set_dis= off_set_dis / neighbour_count;
        if (sqrt(off_set_dis) < 1.0e-12)sigma_s=sqrt(off_set_dis) + 1.0e-12;
        else sigma_s=sqrt(off_set_dis);
        sigma_s=2*sigma_s;
        float sum = 0;
        float normalizer = 0;
        for (int j = 0; j < neighbour_map[i].size(); ++j) {
            if (dist_map[i][j] < 2 * sigma_c) {
                Eigen::Vector3d vec =  V.row(neighbour_map[i][j])-V.row(i);
                float t = dist_map[i][j];
                float h = n.x() * vec.x() + n.y() * vec.y() + n.z() * vec.z();
                double wc = exp(-1 * t * t / (2 * sigma_c * sigma_c));
                double ws = exp(-1 * h * h / (2 * sigma_s * sigma_s));
                sum += wc * ws * h;
                normalizer += wc * ws;
            }
        }
        if (normalizer != 0) {
          // std::cout<<(sum / normalizer)<<std::endl;
            V.row(i) = V.row(i) + (sum / normalizer) * N_vertices.row(i);
        }
    }
    }
}






bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
  std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;
  if (key == 'B')
  {
    bilateral(V,F);
  }
  else if (key == 'R')
  {
    add_noise(V,0.001);
  }
  else if (key == 'N')
  {
    add_noise_along_normals(V,0.001);
  }
  else if (key == 'S')
  {
    string output_file="../denoised/"+file.substr(file.find_last_of("/") + 1);

    igl::writeOBJ(output_file, V, F);
  }
  else if (key == 'W')
  {
    string output_file="../noise/"+file.substr(file.find_last_of("/") + 1);

    igl::writeOBJ(output_file, V, F);
  }

  else
    return false;

  viewer.data(0).clear();
  viewer.data(0).set_mesh(V, F);

  return false;
}




// ------------ main program ----------------
int main(int argc, char *argv[])
{
  if(argc<2) {
    std::cout << "Error: input file required (.ply)" << std::endl;
    return 0;
  }
  std::cout << "reading input file: " << argv[1] << std::endl;
  file = argv[1];
  if(file.substr(file.find_last_of(".") + 1) == "off") {
        igl::readOFF(file,V,F);
    }else if(file.substr(file.find_last_of(".") + 1) == "obj"){
        igl::readOBJ(file,V,F);
    }else {
        std::cout << "ELSE..." << std::endl;
    }
  std::cout << "Points: " << V.rows() << std::endl;

  build_neighbour();

  igl::opengl::glfw::Viewer viewer; // create the 3d viewer
  viewer.callback_key_down = &key_down;
  viewer.data(0).set_mesh(V, F);
  viewer.launch(); // run the viewer
}




