// SFCSort.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>

#include <Eigen/Dense>

#include "Hilbert_Zoltan.h"


void writeConnectivityLine2VTK(const std::string& fileName, MatD vertices) {

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

int main()
{
    size_t numPts = 1000;

    //-----------------Zoltan-------------------
    {
        printf("[Zoltan] Sorting.....\n");

        MatD vertices;

        //Append random pts on a large box
        appendPtsInBox(vertices, numPts,         //Container, NumPts
            0.0, 10.0, 0.0, 10.0, 0.0, 10.0);//Box: Xmin,Xmax,Ymin,Ymax,Zmin,Zmax

        //Append random pts on a small box
        appendPtsInBox(vertices, numPts,        //Container, NumPts
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0);//Box: Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
        //std::cout << vertices<<std::endl;

        printf("Random rows\n");
        //Random the input pts again
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(vertices.rows());
        perm.setIdentity();
        std::random_shuffle(perm.indices().data(), perm.indices().data() + perm.indices().size());
        vertices = perm * vertices; // permute rows

        //std::cout << vertices << std::endl;

        writeConnectivityLine2VTK("BeforeSort_Zontal.vtk", vertices);

        auto start = std::chrono::system_clock::now();
        Zoltan_LocalHSFC_Order(3,vertices);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "RunTime=" << elapsed.count() << "ms ";
        //std::cout << vertices << std::endl;

        printf("\t Done !\n");
        writeConnectivityLine2VTK("AfterSort_Zontal.vtk", vertices);
    }
    
  

}
