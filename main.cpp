#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>

#include <igl/gaussian_curvature.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>

#include "HalfedgeBuilder.cpp"



using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V;
MatrixXi F;


MatrixXd N_faces;   //computed calling pre-defined functions of LibiGL
MatrixXd N_vertices; //computed calling pre-defined functions of LibiGL


MatrixXd lib_N_vertices;  //computed using face-vertex structure of LibiGL
MatrixXi lib_Deg_vertices;//computed using face-vertex structure of LibiGL

MatrixXd he_N_vertices; //computed using the HalfEdge data structure


// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier){
    switch(key){
        case '1':
            viewer.data().set_normals(N_faces);
            return true;
        case '2':
            viewer.data().set_normals(N_vertices);
            return true;
        case '3':
            viewer.data().set_normals(lib_N_vertices);
            return true;
        case '4':
            viewer.data().set_normals(he_N_vertices);
            return true;
        default: break;
    }
    return false;
}

    

MatrixXd surfaceLaplace(HalfedgeDS he){
    int nV = he.sizeOfVertices();
    
    MatrixXd laplace;
    laplace.setZero(he.sizeOfVertices(), he.sizeOfVertices());
    
    for (int vIndex = 0; vIndex < nV; vIndex++) {
        int first_he = he.getEdge(vIndex);
        int current_he = first_he;
        
        float sum_angles = 0;
        float area_dual_cell = 0;
        
        Vector3d  v1 = V.row(vIndex);
        
        //traverse around 1-ring
        do {
            //get vertices of the face
            int v2_index = he.getTarget(he.getNext(current_he));
            int v3_index = he.getTarget(he.getPrev(current_he));
            Vector3d  v2 = V.row(v2_index);
            Vector3d  v3 = V.row(v3_index);
            
            //tip angle
            float v1angle = acos((v2 - v1).dot(v3 - v1)/((v2 - v1).norm() * (v3 - v1).norm()));
            sum_angles += v1angle;
            
            //dual cell area
            //           float v2angle = F(i,j);
            float alpha = acos((v1 - v3).dot(v2 - v3) / ((v1 - v3).norm() * (v2 - v3).norm()));
            
            // find A_voronoi
            float v1v2 = (v1 - v2).norm(); //v1v2 length
            
            int v2_next_index = he.getTarget(he.getNext(he.getOpposite(he.getNext(current_he)))); //corresponding v2 of next triangle
            Vector3d v2_next = V.row(v2_next_index);
            float beta = acos((v1 - v2_next).dot(v2 - v2_next) / ((v1 - v2_next).norm() * (v2 - v2_next).norm()));
            
            float cot_alpha = cos(alpha) / sin(alpha);
            float cot_beta = cos(beta) / sin(beta);
            laplace(vIndex,v2_index) = -0.5 * (cot_alpha + cot_beta);
            laplace(vIndex,vIndex) += 0.5 * (cot_alpha + cot_beta);
            
//            std::cout << "i = " << vIndex << " j = " << v2_index << std::endl;//
//            std::cout << "alpha = " << alpha << " beta = " << beta << std::endl;
//            std::cout << "cot_alpha = " << cot_alpha << " cot_beta = " << cot_beta << std::endl;
//            std::cout << "(1/2) * (cot_alpha + cot_beta) = " << 0.5 * (cot_alpha + cot_beta) << std::endl;
//            std::cout << laplace << std::endl;
            
            current_he = he.getOpposite(he.getNext(current_he));
        }
        while (current_he != first_he); // for each face
        
        
    }
    return laplace;
}

    
    int vertex_neighbours(HalfedgeDS he, int v) {
        int result=0;
        int e=he.getEdge(v);
        MatrixXd u_,v_;
        RowVector3d v1=V.row(he.getTarget(e));
        RowVector3d v2=V.row(he.getTarget(he.getOpposite(e)));
        MatrixXd neigh(V.rows(),3);
        neigh.row(0)=v2-v1;//v1-v2;
        result++;
        int pe=he.getOpposite(he.getNext(e));
        
        MatrixXd adjucencyM;
        adjucencyM.setZero(he.sizeOfVertices(), he.sizeOfVertices());
        
        
        while(pe!=e){
            v2=V.row(he.getTarget(he.getOpposite(pe)));
            neigh.row(result)=v2-v1;//v1-v2;
            result++;
            pe=he.getOpposite(he.getNext(pe));

        }
        MatrixXd neighbours(result,3);

        for (int i=0;i<result;i++){
            neighbours.row(i)=neigh.row(i);
        }


        // std::cout<<neighbours<<std::endl;

        MatrixXd cov = neighbours.transpose()*neighbours;
        SelfAdjointEigenSolver<MatrixXd> eigen_solver((cov.transpose())*cov);
        u_ = eigen_solver.eigenvalues();
        v_ = eigen_solver.eigenvectors();

        std::cout<<u_<<std::endl;
        std::cout<<v_<<std::endl;
        return result;

    
    }

    
    /**
     * Return the number of boundaries of the mesh
   * Exercice
     **/
int countBoundaries(HalfedgeDS he){
    // initialize boundary count to zero
    int boundaries = 0;
    // set bool list to track visited halfedges
    bool visited[he.sizeOfHalfedges()];
    //initialize with false
    memset(visited,false,sizeof(visited));


    for(int i=0;i<he.sizeOfHalfedges();i++){
        if(!visited[i]){
            int last = he.getNext(he.getNext(he.getNext(i)));

            if(i != last){
                visited[i] = true;
                last = he.getNext(i);

                while(last != i){
                    visited[last] = true;
                    last = he.getNext(last);
                }

                boundaries++;
            }
        }
    }


    return boundaries;
}





// ------------ main program ----------------
int main(int argc, char *argv[]) {

//    if(argc<2) {
//        std::cout << "Error: input file required (.OFF)" << std::endl;
//        return 0;
//    }
//    std::cout << "reading input file: " << argv[1] << std::endl;

//    igl::readOFF(argv[1], V, F);
        igl::readOFF("../data/cube_tri.off",V,F);

    //print the number of mesh elements
        //std::cout << "Points: " << V.rows() << std::endl;


       HalfedgeBuilder* builder=new HalfedgeBuilder();  //

       HalfedgeDS he=builder->createMesh(V.rows(), F);  //

    // compute vertex degrees
        // vertexDegreeStatistics(he);  //
        // lib_vertexDegrees();

    // compute number of boundaries
        int B=countBoundaries(he);  //
        std::cout << "The mesh has " << B << " boundaries" << std::endl << std::endl;//
    // find_neighbours(he);
//    for (int i=0;i<V.rows();++i){
//
//        std::cout << "number of degrees for vertex-"<<i<<" = "<< vertex_neighbours(he,i)<<std::endl;
//    }
//
    
    MatrixXd L = surfaceLaplace(he);
    std::cout << "========== Laplace matrix =========== " << std::endl;
    std::cout<< L <<std::endl;
    
    SelfAdjointEigenSolver<MatrixXd> eigen_solver(L);
    MatrixXd u_,v_;
    u_ = eigen_solver.eigenvalues();
    v_ = eigen_solver.eigenvectors();
    
    std::cout << "========== eigenvalues =========== " << std::endl;
    std::cout<< u_ <<std::endl;
    std::cout << "========== eigenvectors =========== " << std::endl;
    std::cout<< v_ <<std::endl;
    


 ///////////////////////////////////////////

  // Plot the mesh with pseudocolors
  igl::opengl::glfw::Viewer viewer; // create the 3d viewer

  viewer.callback_key_down = &key_down;
  viewer.data().show_lines = false;
  viewer.data().set_mesh(V, F);  //
  viewer.data().set_normals(N_faces);  //
  std::cout<<
    "Press '1' for per-face normals calling pre-defined functions of LibiGL."<<std::endl<<
    "Press '2' for per-vertex normals calling pre-defined functions of LibiGL."<<std::endl<<
    "Press '3' for lib_per-vertex normals using face-vertex structure of LibiGL ."<<std::endl<<
    "Press '4' for HE_per-vertex normals using HalfEdge structure."<<std::endl;


    VectorXd Z;
    Z.setZero(V.rows(),1);

//  Z colors
   // Use the z coordinate as a scalar field over the surface
    for(int i ; i<V.rows();++i){
    Z[i]= V.row(i).z();
    }

    MatrixXd C;
  // Assign per-vertex colors
    igl::jet(Z,true,C);
    viewer.data().set_colors(C);  // Add per-vertex colors

  //viewer.core(0).align_camera_center(V, F);  //not needed
    viewer.launch(); // run the viewer
}
