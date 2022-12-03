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
#include <tgmath.h>
#include <math.h>
#include<fstream>

#include "HalfedgeBuilder.cpp"



using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V;
MatrixXi F;


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


MatrixXd* surfaceLaplace(HalfedgeDS he){
    int nV = he.sizeOfVertices();
    MatrixXd weights, diagonal;
    weights.setZero(he.sizeOfVertices(), he.sizeOfVertices());
    diagonal.setZero(he.sizeOfVertices(), he.sizeOfVertices());
    
    for (int vIndex = 0; vIndex < nV; vIndex++) {
        int first_he = he.getEdge(vIndex);
        int current_he = first_he;
        
        float sum_angles = 0;
        float area_dual_cell = 0;
        
        Vector3d  v1 = V.row(vIndex); //i
        
        //traverse around 1-ring
        do {
            //get vertices of the face
            int v2_index = he.getTarget(he.getNext(current_he));
            int v3_index = he.getTarget(he.getPrev(current_he));
            Vector3d  v2 = V.row(v2_index); //j
            Vector3d  v3 = V.row(v3_index); //vertex with alpha
            
            
            
            int v2_next_index = he.getTarget(he.getNext(he.getOpposite(he.getNext(current_he))));
            Vector3d v2_next = V.row(v2_next_index); //vertex with beta
            
            float alpha = acos((v1 - v3).dot(v2 - v3) / ((v1 - v3).norm() * (v2 - v3).norm()));
            float beta = acos((v1 - v2_next).dot(v2 - v2_next) / ((v1 - v2_next).norm() * (v2 - v2_next).norm()));
            
            float cot_alpha = cos(alpha) / sin(alpha);
            float cot_beta = cos(beta) / sin(beta);
            
            weights(vIndex,v2_index) = -0.5 * (cot_alpha + cot_beta);
            weights(vIndex,vIndex) += 0.5 * (cot_alpha + cot_beta);
            
            // find A_voronoi
            float v1v2 = (v1 - v2).norm(); //v1v2 length
            
            float v3v1v2 = acos((v3 - v1).dot(v2 - v1) / ((v3 - v1).norm() * (v2 - v1).norm()));
            float v1v2v3 = acos((v1 - v2).dot(v3 - v2) / ((v1 - v2).norm() * (v3 - v2).norm()));
            float _90 = M_PI/2;
            if (alpha < _90 && v3v1v2 < _90 && v1v2v3 < _90){
                area_dual_cell = (1.0/8) * (cot_alpha + cot_beta) * v1v2;
            } else {
                if (v3v1v2 >= _90 ){
                    area_dual_cell = (1.0/2) * (v3 - v1).norm() * (v2 - v1).norm();
                } else {
                    area_dual_cell = (1.0/4) * (v3 - v1).norm() * (v2 - v1).norm();
                }
            }
            
            diagonal(vIndex,vIndex) += area_dual_cell; //Ai
            
            //            std::cout << "i = " << vIndex << " j = " << v2_index << std::endl;//
            //            std::cout << "alpha = " << alpha << " beta = " << beta << std::endl;
            //            std::cout << "cot_alpha = " << cot_alpha << " cot_beta = " << cot_beta << std::endl;
            //            std::cout << "(1/2) * (cot_alpha + cot_beta) = " << 0.5 * (cot_alpha + cot_beta) << std::endl;
            //            std::cout << laplace << std::endl;
            
            current_he = he.getOpposite(he.getNext(current_he));
        }
        while (current_he != first_he); // for each face
        
    }
    
    
    MatrixXd* matrices = new MatrixXd[2];
    matrices[0] = weights;
    matrices[1] = diagonal;
    
    return matrices;
}

void saveData(string fileName, MatrixXd  matrix)
{
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
 
    ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}

MatrixXd findGPS(MatrixXd V, MatrixXi F){
    
    HalfedgeBuilder* builder=new HalfedgeBuilder();  //
    
    HalfedgeDS he=builder->createMesh(V.rows(), F);  //
    
    // compute vertex degrees
    // vertexDegreeStatistics(he);  //
    // lib_vertexDegrees();
    
    // compute number of boundaries
    int B=countBoundaries(he);  //
    std::cout << "The mesh has " << B << " boundaries" << std::endl << std::endl;//
    
    MatrixXd* L = surfaceLaplace(he);
    std::cout << "========== weights matrix =========== " << std::endl;
    std::cout<< L[0].rows() <<std::endl;
//    std::cout<< L[0] <<std::endl;
    std::cout << "========== area matrix =========== " << std::endl;
    std::cout<< L[1].rows() <<std::endl;
//    std::cout<< L[1] <<std::endl;
//    std::cout << "========== det of weights =========== " << std::endl;
//    std::cout<< L[0].determinant() <<std::endl;
//    std::cout << "========== det of area =========== " << std::endl;
//    std::cout<< L[1].determinant() <<std::endl;
    
    SelfAdjointEigenSolver<MatrixXd> eigen_solver(L[1].inverse()*L[0]);
    //    std::cout << "========== eigenvalues =========== " << std::endl;
    //    std::cout<< u_ <<std::endl;
    //    std::cout << "========== eigenvectors =========== " << std::endl;
    //    std::cout<< v_ <<std::endl;
    //    std::cout << "========== eigenvector size =========== " << std::endl;
    //    std::cout<< v_.row(0) << std::endl;
    //
    
    MatrixXd u_,v_;
    u_ = eigen_solver.eigenvalues();
    v_ = eigen_solver.eigenvectors();
    
//    GeneralizedEigenSolver<MatrixXd> eigen_solver;
//    eigen_solver.compute(L[0], L[1]);
//    cout << "The (complex) numerators of the generalzied eigenvalues are: " << eigen_solver.alphas().transpose() << endl;
//    cout << "The (real) denominatore of the generalzied eigenvalues are: " << eigen_solver.betas().transpose() << endl;
//    cout << "The (complex) generalzied eigenvalues are (alphas./beta): " << eigen_solver.eigenvalues() << endl;
//    cout << "The (complex) generalzied eigenvectors are (alphas./beta): " << eigen_solver.eigenvectors().col(0) << endl;
    
    saveData("/Users/kulakshifernando/Edu/INF574/Project/eigenvalues.csv",u_);
    saveData("/Users/kulakshifernando/Edu/INF574/Project/eigenvectors.csv",v_);
    
        int cols = v_.cols();
        int d = std::min(cols - 1, 15);
    
    MatrixXd gps;
    gps.setZero(V.rows(), d);
    for (int i = 0; i < V.rows(); i++) {
        for (int j = 1; j <= d; j++) {
            double eigenVal = u_.row(j)[0];
            gps(i,j-1) = v_(i,j)/sqrt(eigenVal);
        }
    }
    saveData("/Users/kulakshifernando/Edu/INF574/Project/gps.csv",gps);
    std::cout << "========== gps =========== " << std::endl;
    std::cout<< gps <<std::endl;
    
    return gps;
}


// ------------ main program ----------------
int main(int argc, char *argv[]) {
    
    MatrixXd N;
//    igl::readOFF("../data/cube_tri.off",V,F);
        igl::readOBJ("/Users/kulakshifernando/Edu/INF574/Project/TD4/data/horse-poses/horse-01.obj",V,F);
    
    
    //    MatrixXd V1;
    //    MatrixXi F1;
    //    igl::readOBJ("/Users/kulakshifernando/Edu/INF574/Project/TD4/data/horse-poses/horse-01.obj",V1,F1);
    
    
    MatrixXd gps = findGPS(V, F);
    //    MatrixXd gps1 = findGPS(V1, F1);
    std::cout<< gps.rowwise().sum() <<std::endl;
    //    std::cout<< gps1.rowwise().sum() <<std::endl;
    
    
    // Plot the mesh with pseudocolors
    igl::opengl::glfw::Viewer viewer; // create the 3d viewer
    
    //  viewer.callback_key_down = &key_down;
    viewer.data().show_lines = false;
    viewer.data().set_mesh(V, F);  //
    //  viewer.data().set_normals(N_faces);  //
    
    viewer.launch(); // run the viewer
}
