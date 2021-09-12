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
        Iterator_on_vertices() = delete;
        Iterator_on_vertices(const Iterator_on_vertices &copy): m_ptr(copy.m_ptr), currentIdx(copy.currentIdx) {}
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

    struct Iterator_on_faces
    {
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = Triangle;
        using pointer           = Triangle*;
        using reference         = Triangle&;
        Iterator_on_faces() = delete;
        Iterator_on_faces(const Iterator_on_faces &copy): m_ptr(copy.m_ptr), currentIdx(copy.currentIdx) {}
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

    struct Circulator_on_faces
    {
//        using iterator_category = std::forward_iterator_tag;
//        using difference_type   = std::ptrdiff_t;
        using value_type        = Triangle;
        using pointer           = Triangle*;
        using reference         = Triangle&;
        Circulator_on_faces() = delete;
        Circulator_on_faces(const Circulator_on_faces &c): m_mesh(c.m_mesh), m_ptr(c.m_ptr), m_center_idx(c.m_center_idx) {}
//        Circulator_on_faces(Mesh &mesh, pointer ptr, const Vertex *center): m_mesh(mesh), m_ptr(ptr), m_center(center){}
        Circulator_on_faces(Mesh &mesh, pointer ptr, std::size_t center_idx): m_mesh(mesh), m_ptr(ptr), m_center_idx(center_idx){}
        Circulator_on_faces(Mesh &mesh, pointer ptr, Iterator_on_vertices vIt): m_mesh(mesh), m_ptr(ptr){
            m_center_idx = vIt.getIdx();
        }

        reference operator*() const { return *m_ptr; }
        pointer operator->() { return m_ptr; }

        // Prefix increment
        Circulator_on_faces& operator++() {
            // In the vertex array of the triangle, find the index corresponding to the
            // index of the center vertex of the circulator. The index after this one will
            // be the vertex opposite to the next triangle in counter-clockwise order.
            int idx;
            for (idx = 0 ; idx < 3 ; ++idx)
                if (m_ptr->vertices[idx] == m_center_idx) break;
            int next = idx < 2 ? idx+1 : 0;
            m_ptr = &(m_mesh.triangles[m_ptr->adjacent[next]]);
            return *this;

        }
        friend bool operator== (const Circulator_on_faces& a, const Circulator_on_faces& b) { return a.m_ptr == b.m_ptr; }
        friend bool operator!= (const Circulator_on_faces& a, const Circulator_on_faces& b) { return a.m_ptr != b.m_ptr; }
    private:
        Mesh &m_mesh;
        pointer m_ptr;
//        const Vertex *m_center;
        std::size_t m_center_idx;
//        std::size_t tIdx;
    };

    struct Circulator_on_vertices
    {
//        using iterator_category = std::forward_iterator_tag;
//        using difference_type   = std::ptrdiff_t;
        using value_type        = Vertex;
        using pointer           = Vertex*;
        using reference         = Vertex&;
        Circulator_on_vertices(Mesh &mesh, pointer ptr, const Vertex *center): m_mesh(mesh), m_ptr(ptr), m_center(center){}
        Circulator_on_vertices(Mesh &mesh, pointer ptr, std::size_t center_idx): m_mesh(mesh), m_ptr(ptr), m_center_idx(center_idx){}
        Circulator_on_vertices(Mesh &mesh, pointer ptr, Iterator_on_vertices vIt): m_mesh(mesh), m_ptr(ptr){
            m_center_idx = vIt.getIdx();
        }

        reference operator*() const { return *m_ptr; }
        pointer operator->() { return m_ptr; }

        // Prefix increment
        Circulator_on_vertices& operator++() {

        }
    private:
        Mesh &m_mesh;
        pointer m_ptr;
        const Vertex *m_center;
        std::size_t m_center_idx;
        std::size_t tIdx;
    };

    Iterator_on_vertices vertices_begin() {
        return Iterator_on_vertices(&(*vertices.begin()));
    }
    Iterator_on_vertices vertices_past_the_end() {
        return Iterator_on_vertices(&(*vertices.end()));
    }
    Iterator_on_faces faces_begin() {
        return Iterator_on_faces(&(*triangles.begin()));
    }
    Iterator_on_faces faces_past_the_end() {
        return Iterator_on_faces(&(*triangles.end()));
    }
//    Circulator_on_faces incident_faces(const Vertex &v) {
//        return Circulator_on_faces(*this, &(triangles[v.triangleIdx]), &v);
//    }

    Circulator_on_faces incident_faces(std::size_t vIdx) {
        return Circulator_on_faces(*this, &(triangles[vertices[vIdx].triangleIdx]), vIdx);
    }

//    Circulator_on_vertices adjacent_vertices(std::size_t vIdx) {
//        const auto &inc_triangle = triangles[vertices[vIdx].triangleIdx];

//        return Circulator_on_vertices(*this, vertex, &vertices[vIdx]);
//    }

    void loadOFF (std::string path);
    void findTopology();
    void findTopology(const std::vector<Point> &points, const std::vector<std::array<uint, 3> > &faces);
    void drawMesh(bool wireframe = false);
    void drawMeshWireFrame();
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
