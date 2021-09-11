#ifndef MESH_H
#define MESH_H

#include <QGLWidget>
#include <fstream>

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

class Triangle;

//** TP : TO MODIFY
class Vertex
{
public:
    Point p;
    int triangleIdx = -1;
};

class Triangle
{
public:
    std::array<uint, 3> vertices; // Trigonometric order
    std::array<uint, 3> adjacent;
    // Constraint: vertex i facing adjacent triangle i
    Triangle();
//    ~Triangle();
};

class Mesh
{
    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;
public:
    Mesh();
    //~Mesh();
    void readOffFile (std::string path);
    void findTopology(const std::vector<Point> &points, const std::vector<std::array<uint, 3> > &faces);
    void drawMesh(bool wireframe = false);
    void drawMeshWireFrame();
//    void loadOff();
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
