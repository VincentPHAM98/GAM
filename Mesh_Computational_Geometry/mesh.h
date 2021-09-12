#ifndef MESH_H
#define MESH_H

#include <QGLWidget>
#include <fstream>
#include <iterator>

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
    struct Iterator_on_vertices
    {

        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = Vertex;
        using pointer           = Vertex*;
        using reference         = Vertex&;
        Iterator_on_vertices(pointer ptr): m_ptr(ptr) {}

        reference operator*() const { return *m_ptr; }
        pointer operator->() { return m_ptr; }

        // Prefix increment
        Iterator_on_vertices& operator++() { currentIdx++; m_ptr++; return *this; }

        // Postfix increment
        Iterator_on_vertices operator++(int) { currentIdx++; Iterator_on_vertices tmp = *this; ++(*this); return tmp; }

        std::size_t getIdx() { return currentIdx; }

        friend bool operator== (const Iterator_on_vertices& a, const Iterator_on_vertices& b) { return a.m_ptr == b.m_ptr; }
        friend bool operator!= (const Iterator_on_vertices& a, const Iterator_on_vertices& b) { return a.m_ptr != b.m_ptr; }

    private:
        pointer m_ptr;
        std::size_t currentIdx = 0;
    };
    Iterator_on_vertices vertices_begin() {
        return Iterator_on_vertices(&(*vertices.begin()));
    }
    Iterator_on_vertices vertices_past_the_end() {
        return Iterator_on_vertices(&(*vertices.end()));
    }

    struct Iterator_on_faces
    {
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = Triangle;
        using pointer           = Triangle*;
        using reference         = Triangle&;
        Iterator_on_faces(pointer ptr): m_ptr(ptr) {}

        reference operator*() const { return *m_ptr; }
        pointer operator->() { return m_ptr; }

        // Prefix increment
        Iterator_on_faces& operator++() { currentIdx++; m_ptr++; return *this; }

        // Postfix increment
        Iterator_on_faces operator++(int) { currentIdx++; Iterator_on_faces tmp = *this; ++(*this); return tmp; }

        std::size_t getIdx() { return currentIdx; }

        friend bool operator== (const Iterator_on_faces& a, const Iterator_on_faces& b) { return a.m_ptr == b.m_ptr; }
        friend bool operator!= (const Iterator_on_faces& a, const Iterator_on_faces& b) { return a.m_ptr != b.m_ptr; }

    private:
        pointer m_ptr;
        std::size_t currentIdx = 0;
    };
    Iterator_on_faces faces_begin() {
        return Iterator_on_faces(&(*triangles.begin()));
    };
    Iterator_on_faces faces_past_the_end() {
        return Iterator_on_faces(&(*triangles.end()));
    };

    void loadOFF (std::string path);
    void findTopology();
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
