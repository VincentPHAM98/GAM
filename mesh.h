#ifndef MESH_H
#define MESH_H

#include <QGLWidget>
#include <vector>
#include <string>

using namespace std;

// TO MODIFY
class Point
{
public:
    double _x;
    double _y;
    double _z;

    Point():_x(),_y(),_z() {}
    Point(float x_, float y_, float z_):_x(x_),_y(y_),_z(z_) {}
};

class Vertex{
public:
    int pointIndex;
    int faceIndex;

    Vertex(int _pointIndex):pointIndex(_pointIndex){}
    Vertex(int _pointIndex, int _faceIndex):pointIndex(_pointIndex), faceIndex(_faceIndex) {}
};

class Face{
public:
    int vertices[3];
    int adjFaces[3];

    Face(int _v1, int _v2, int _v3){
        vertices[0] = _v1;
        vertices[1] = _v2;
        vertices[2] = _v3;
    }

    Face(int _v1, int _v2, int _v3, int _f1, int _f2, int _f3){
        vertices[0] = _v1;
        vertices[1] = _v2;
        vertices[2] = _v3;

        adjFaces[0] = _f1;
        adjFaces[1] = _f2;
        adjFaces[2] = _f3;
    }
};

//** TP : TO MODIFY

class Mesh
{
  // (Q ou STL)Vector of vertices
  // (Q ou STL)Vector of faces
  // Those who do not know about STL Vectors should have a look at cplusplus.com examples
public:
    vector<Point> points;
    vector<Vertex> vertices;
    vector<Face> faces;
    // map <arrete>, <face, indice sommet opposé>
    map<pair<int, int>, pair<int, int>> topology;
    Mesh(); // Constructors automatically called to initialize a Mesh (default strategy)
    //~Mesh(); // Destructor automatically called before a Mesh is destroyed (default strategy)
    void drawMesh();
    void drawMeshWireFrame();
    void initTetrahedron();
    void initPyramid();
    void init2dBBox();
    void initFile(string filepath);
    void handleFace(int indVertex1, int indVertex2, int indVertex3, int faceIndex);
    void handleEdge(int indVertex1, int indVertex2, int faceIndex, int opposedVertex);
    void clearData();

    void splitTriangle(int indFace, int indVertex);
    void edgeFlip(int indFace1, int indFace2);
    float orientaion2D(Point p1, Point p2, Point p3);
    bool isInside(Point p, Face f);
    void insertPoint(Point p);


};

class GeometricWorld //Generally used to create a singleton instance
{
  QVector<Point> _bBox;  // Bounding box // ou std::vector
public :
  GeometricWorld();
  void draw();
  void drawWireFrame();
  // ** TP Can be extended with further elements;
  Mesh _mesh;
};


#endif // MESH_H
