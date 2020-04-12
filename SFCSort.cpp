// SFCSort.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>

#include "Hilbert_gmsh.h"
#include "Hilbert_Zoltan.h"


void writeConnectivityLine2VTK(const std::string& fileName, std::vector<SPoint3*> verts) {

    std::ofstream out(fileName);
    int NumSLs = 1;
    int numVertices = verts.size();

    out << "# vtk DataFile Version 4.1\n";
    out << "vtk output\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "\n";

    out << "POINTS " << numVertices << " float\n";
    for (auto& v : verts)
            out << v->x() << " " << v->y() << " " << v->z() << std::endl;
    
    out << "\n";

    out << "VERTICES " << numVertices << " " << numVertices*2 << std::endl;
    for (int i = 0; i < verts.size(); i++)
        out << 1 << " " << i << std::endl;

    out << "\n";

    out << "LINES " << NumSLs << " " << numVertices + NumSLs << std::endl;
    int vertexID = 0;
    out << numVertices;
    for (int i = 0; i < verts.size(); i++) {
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



void appendPtsInBox(std::vector<SPoint3*>& vertices, int NumPts,
    double BoundingBoxXmin, double BoundingBoxXmax,
    double BoundingBoxYmin, double BoundingBoxYmax,
    double BoundingBoxZmin, double BoundingBoxZmax) {

    for (unsigned int i = 0; i < NumPts; ++i) {

        double xr = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 1.0));//Random number between 0,1
        double yr = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 1.0));//Random number between 0,1
        double zr = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / 1.0));//Random number between 0,1

        xr = BoundingBoxXmin + xr * (BoundingBoxXmax - BoundingBoxXmin);
        yr = BoundingBoxYmin + yr * (BoundingBoxYmax - BoundingBoxYmin);
        zr = BoundingBoxZmin + zr * (BoundingBoxZmax - BoundingBoxZmin);

        SPoint3* v = new SPoint3(xr, yr, zr);
        vertices.push_back(v);
        //printf("(%d)%lf,%lf,%lf\n", i, xr, yr, zr);
    }

}

int main()
{
    //-----------------Gmsh-------------------
    size_t numPts = 10000;
    {
        printf("[Gmsh] Sorting.....\n");

        //Generate random points
        std::vector<SPoint3*> vertices;

        //Append random pts on a large box
        appendPtsInBox(vertices, numPts,         //Container, NumPts
            0.0, 10.0, 0.0, 10.0, 0.0, 10.0);//Box: Xmin,Xmax,Ymin,Ymax,Zmin,Zmax

        //Append random pts on a small box
        appendPtsInBox(vertices, numPts,        //Container, NumPts
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0);//Box: Xmin,Xmax,Ymin,Ymax,Zmin,Zmax

        //Random the input pts again
        auto rng = std::default_random_engine{};
        std::shuffle(std::begin(vertices), std::end(vertices), rng);

        //writeConnectivityLine2VTK("BeforeSort.vtk", vertices);

        // HilbertSort h;
        auto start = std::chrono::system_clock::now();
        HilbertSort h(1000);
        h.Apply(vertices);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "RunTime=" << elapsed.count() << "ms ";

        printf("\t Done !\n");
        writeConnectivityLine2VTK("AfterSort.vtk", vertices);
    }
    

    //-----------------Zoltan-------------------
    {
        printf("[Zoltan] Sorting.....\n");

        std::vector<SPoint3*> vertices;

        //Append random pts on a large box
        appendPtsInBox(vertices, numPts,         //Container, NumPts
            0.0, 10.0, 0.0, 10.0, 0.0, 10.0);//Box: Xmin,Xmax,Ymin,Ymax,Zmin,Zmax

        //Append random pts on a small box
        appendPtsInBox(vertices, numPts,        //Container, NumPts
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0);//Box: Xmin,Xmax,Ymin,Ymax,Zmin,Zmax

        //Random the input pts again
        auto rng = std::default_random_engine{};
        std::shuffle(std::begin(vertices), std::end(vertices), rng);

        //writeConnectivityLine2VTK("BeforeSort_Zontal.vtk", vertices);

        auto start = std::chrono::system_clock::now();
        Zoltan_LocalHSFC_Order(3,vertices);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "RunTime=" << elapsed.count() << "ms ";

        printf("\t Done !\n");
        writeConnectivityLine2VTK("AfterSort_Zontal.vtk", vertices);
    }
    
  

}
