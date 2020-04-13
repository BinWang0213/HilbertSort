// SFCSort.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>


#include "hilbertSort.h"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>

#include <Eigen/Dense>


//---------------------------Main function----------------------------
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatD;

void writeConnectivityLine2VTK(const std::string& fileName, Eigen::Ref<MatD> vertices) {

    std::ofstream out(fileName);
    int NumSLs = 1;
    int numVertices = vertices.rows();

    out << "# vtk DataFile Version 4.1\n";
    out << "vtk output\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "\n";

    out << "POINTS " << numVertices << " float\n";
    for (unsigned int i = 0; i < numVertices; ++i)
            out << vertices(i,0) << " " << vertices(i, 1) << " " << vertices(i, 2) << std::endl;
    
    out << "\n";

    out << "VERTICES " << numVertices << " " << numVertices*2 << std::endl;
    for (int i = 0; i < numVertices; i++)
        out << 1 << " " << i << std::endl;

    out << "\n";

    out << "LINES " << NumSLs << " " << numVertices + NumSLs << std::endl;
    int vertexID = 0;
    out << numVertices;
    for (int i = 0; i < numVertices; i++) {
        out << " " << vertexID;
        vertexID += 1;
    }
    out << "\n";
    
    out << "\n\n";

    out << "POINT_DATA " << numVertices << "\n";
    out << "FIELD FieldData 1\n"; //Now we only set one field here

    out << "VertID 1 " << numVertices << " int\n";
    for (int i = 0; i < numVertices; ++i)
        out << i << " " << "\n";;

    out.close();
}



void appendPtsInBox(MatD& vertices, int NumPts,
    double BoundingBoxXmin, double BoundingBoxXmax,
    double BoundingBoxYmin, double BoundingBoxYmax,
    double BoundingBoxZmin, double BoundingBoxZmax) {

    MatD vert_old(vertices);
    MatD vert_new(NumPts, 3);

    for (unsigned int i = 0; i < NumPts; ++i) {

        double xr = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 1.0));//Random number between 0,1
        double yr = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 1.0));//Random number between 0,1
        double zr = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 1.0));//Random number between 0,1

        xr = BoundingBoxXmin + xr * (BoundingBoxXmax - BoundingBoxXmin);
        yr = BoundingBoxYmin + yr * (BoundingBoxYmax - BoundingBoxYmin);
        zr = BoundingBoxZmin + zr * (BoundingBoxZmax - BoundingBoxZmin);

        vert_new.row(i) << xr, yr, zr;
        //printf("(%d)%lf,%lf,%lf\n", i, xr, yr, zr);
    }

    vertices.resize(vertices.rows() + NumPts, 3);
    vertices << vert_old, vert_new;
}


std::tuple< MatD,std::vector<int> > Zoltan_LocalHSFC_Order(int dim, Eigen::Ref<MatD> pts) {
    //Modified from zoltan/src/order/hsfcOrder.c
    int n = pts.rows();
    int numGeomDims = dim;
    std::cout << "[Info]\tSorting " << n << " points......."<<std::endl;

    /******************************************************************/
    /* Scale Pts Into Unit Box [0,1] x [0,1] x [0,1]
    /******************************************************************/
    MatD pts_scale=pts;

    Eigen::ArrayXd bbox_min(3), bbox_max(3), bbox_width(3);
    for (int i = 0; i < 3; ++i) {//Find bounding box
        bbox_min(i) = pts.col(i).minCoeff();  
        bbox_max(i) = pts.col(i).maxCoeff();
    }
    bbox_width = bbox_max - bbox_min;
    std::cout <<"\tBbox="<< bbox_min.transpose() << " " << bbox_max.transpose() << std::endl;
    //printf("Bbox width=%lf,%lf,%lf\n", bbox_width[0], bbox_width[1], bbox_width[2]);

    //enlarge bbox slightly to include all points
    bbox_min -= 0.01 * bbox_width;
    bbox_max += 0.01 * bbox_width;
    bbox_width = bbox_max - bbox_min;
    //printf("Bbox="); std::cout << bbox_min.transpose() << " " << bbox_max.transpose() << std::endl;
    //printf("Bbox width=%lf,%lf,%lf\n", bbox_width[0], bbox_width[1], bbox_width[2]);


    //scale pts into unit box
    for (int i = 0; i < n; ++i) 
        for (int j = 0; j < 3; ++j) {
            pts_scale(i,j) -= bbox_min(j);
            pts_scale(i, j) /= bbox_width(j);
        }

    for (int i = 0; i < 3; ++i) {//Find bounding box
        bbox_min(i) = pts_scale.col(i).minCoeff();
        bbox_max(i) = pts_scale.col(i).maxCoeff();
        bbox_width(i) = bbox_max(i) - bbox_min(i);
    }
    //printf("Bbox_scale="); std::cout << bbox_min.transpose() << " " << bbox_max.transpose() << std::endl;
    //printf("Bbox_scale width=%lf,%lf,%lf\n", bbox_width[0], bbox_width[1], bbox_width[2]);
    //printf("(%lf,%lf,%lf)->(%lf,%lf,%lf)\n", pts[i]->at(0), pts[i]->at(1), pts[i]->at(2),
    //    pts_scale[i][0], pts_scale[i][1], pts_scale[i][2]);
    //std::cout << pts_scale << std::endl;

    /******************************************************************/
    /* Generate hsfc keys and indices to be sorted                    */
    /******************************************************************/

    /* space filling curve functions */
    double (*fhsfc)(double*) = NULL;
    if (numGeomDims == 1) fhsfc = Zoltan_HSFC_InvHilbert1d;
    else if (numGeomDims == 2) fhsfc = Zoltan_HSFC_InvHilbert2d;
    else if (numGeomDims == 3) fhsfc = Zoltan_HSFC_InvHilbert3d;


    std::vector<double> hsfcKey(n);
    std::vector<int> sort_idx(n);
    for (int i = 0; i < n; ++i) {
        hsfcKey[i] = fhsfc(&pts_scale.data()[i * 3]);
        sort_idx[i] = i;
        //printf("%d %.15lf\n", coordIndx[i], hsfcKey[i]);
    }

    /******************************************************************/
    /* Sort indices based on keys                                     */
    /******************************************************************/

    auto start = std::chrono::system_clock::now();
    Zoltan_quicksort_pointer_inc_double(sort_idx.data(), hsfcKey.data(), 0, n - 1);
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "\tSortTime=" << elapsed.count() << " ms\n";


    //Update pts based on sorted order 
    pts_scale = pts;
    for (int i = 0; i < n; ++i) {
        pts_scale.row(i) = pts.row(sort_idx[i]);
        //printf("%d->%d\n", i, coordIndx[i]);
    }

    //printf("Before=\n");
    //std::cout << pts_scale << std::endl;

    //printf("After=\n");
    //std::cout << pts << std::endl;
    std::cout << "\tDone!\n";

    return std::make_tuple(pts_scale, sort_idx);
}


int main()
{
    int numPts = 500;

    //-----------------Zoltan-------------------
    std::cout<<"[Zoltan] Sorting.....\n";

    MatD vertices;
    std::vector<int> sort_idx;

    //Append random pts on a large box
    appendPtsInBox(vertices, numPts,         //Container, NumPts
        0.0, 10.0, 0.0, 10.0, 0.0, 10.0);//Box: Xmin,Xmax,Ymin,Ymax,Zmin,Zmax

    //Append random pts on a small box
    appendPtsInBox(vertices, numPts,        //Container, NumPts
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0);//Box: Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
    //std::cout << vertices<<std::endl;

    std::cout<<"Random rows\n";
    //Random the input pts again
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(vertices.rows());
    perm.setIdentity();
    std::random_shuffle(perm.indices().data(), perm.indices().data() + perm.indices().size());
    vertices = perm * vertices; // permute rows

    //std::cout << vertices << std::endl;

    writeConnectivityLine2VTK("BeforeSort_Zontal.vtk", vertices);

    auto start = std::chrono::system_clock::now();
    std::tie(vertices, sort_idx) = Zoltan_LocalHSFC_Order(3,vertices);
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "RunTime=" << elapsed.count() << " ms\n";
    //std::cout << vertices << std::endl;

    std::cout<<"\t Done !\n";
    writeConnectivityLine2VTK("AfterSort_Zontal.vtk", vertices);
}

// ----------------
// Python interface
// ----------------
namespace py = pybind11;

std::tuple<MatD, py::array> pyHilbertSort(int Dim, Eigen::Ref<MatD> Pts_Nx3, bool writeConnectivityMap = false) {
    //non-copy vector https://github.com/pybind/pybind11/issues/1042
    py::scoped_ostream_redirect stream(
        std::cout,                               // std::ostream&
        py::module::import("sys").attr("stdout") // Python output
    );

    std::vector<int> idx_sort;
    MatD pts_sort(Pts_Nx3.rows(),Pts_Nx3.rows());
    std::tie(pts_sort, idx_sort) = Zoltan_LocalHSFC_Order(3, Pts_Nx3);

    if (writeConnectivityMap){
        std::string fname="OrigPts_"+std::to_string(Pts_Nx3.rows())+".vtk";
        writeConnectivityLine2VTK(fname, Pts_Nx3);
        std::cout<<"[Info]\tWrite original point ordering to " << fname <<std::endl;

        fname="SorttedPts_"+std::to_string(Pts_Nx3.rows())+".vtk";
        writeConnectivityLine2VTK(fname, pts_sort);
        std::cout<<"[Info]\tWrite sortted point ordering to " << fname <<std::endl;
    } 

    return std::make_tuple(pts_sort, py::array(idx_sort.size(), idx_sort.data()));
}

PYBIND11_MODULE(pyHilbertSort, m) {
    m.doc() = R"pbdoc(
        pyHilbertSort library
    )pbdoc";

    m.def("hilbertSort", &pyHilbertSort, R"pbdoc(
        Sort 1/2/3D points based on Hilbert curve
        
        Input Args:
        ----------- 
        1. [Int] Dimension 
        2. [Numpy.Array Nx3] InputPoints 
        3. [Bool] writeConnectivityMap, if output original/sorted points connectivity map into VTK

        Author:Bin Wang (binwang.0213@gmail.com)
        Date: April, 2020
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

