#ifndef MESH_H
#define MESH_H

#include <QGLWidget>
#include <QDebug>
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
    QVector3D operator -(const Point & p);
    friend std::ostream& operator <<(std::ostream &os, const Point &p);
};


class Triangle;

//** TP : TO MODIFY
class Vertex
{
public:
    Vertex(): p() {}
    Vertex(const Point &_p): p(_p) {}
    Vertex(const Point &_p, int id): p(_p), triangleIdx(id) {}
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
    Triangle(const Triangle &t): vertices(t.vertices), adjacent(t.adjacent) {}
    Triangle(std::array<uint, 3> _vertices): vertices(_vertices) {}
    Triangle(std::array<uint, 3> _vertices, std::array<uint, 3> _adjacent): vertices(_vertices), adjacent(_adjacent) {}
//    ~Triangle();
    // Gives internal index of given vertex index
    size_t getInternalIdx(size_t vertexIdx, int shift = 0) const;
    size_t getExternalIdx(size_t vertexIdx, int shift = 0) const;
};

class Mesh
{
    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;
    std::vector<QVector3D> laplacians;
    std::vector<double> curvature;
public:
    // Arête de 2 sommets
    using Edge = std::pair<int, int>;
    Mesh();
    //~Mesh();
    struct Iterator_on_vertices
    {

        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = Vertex;
        using pointer           = Vertex*;
        using reference         = Vertex&;
        Iterator_on_vertices() : m_ptr(nullptr), currentIdx(-1) {}
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
        Iterator_on_faces(): m_ptr(nullptr), currentIdx(-1) {}
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
        Circulator_on_faces(const Circulator_on_faces &c)
            : m_mesh(c.m_mesh), m_ptr(c.m_ptr), m_center_idx(c.m_center_idx), m_index(c.m_index) {}
//        Circulator_on_faces(Mesh &mesh, pointer ptr, const Vertex *center): m_mesh(mesh), m_ptr(ptr), m_center(center){}
        Circulator_on_faces(Mesh &mesh, pointer ptr, std::size_t center_idx, int index)
            : m_mesh(mesh), m_ptr(ptr), m_center_idx(center_idx), m_index(index){}
//        Circulator_on_faces(Mesh &mesh, pointer ptr, Iterator_on_vertices vIt): m_mesh(mesh), m_ptr(ptr){
//            m_center_idx = vIt.getIdx();
//        }

        reference operator*() const { return *m_ptr; }
        pointer operator->() { return m_ptr; }

        // Prefix increment
        Circulator_on_faces& operator++() {
            // In the vertex array of the triangle, find the index corresponding to the
            // index of the center vertex of the circulator. The index after this one will
            // be the vertex opposite to the next triangle in counter-clockwise order.
            int next = m_ptr->getInternalIdx(m_center_idx, 1);
            m_index = m_ptr->adjacent[next];
            m_ptr = &(m_mesh.triangles[m_index]);
            return *this;
        }
        friend bool operator== (const Circulator_on_faces& a, const Circulator_on_faces& b) { return a.m_ptr == b.m_ptr; }
        friend bool operator!= (const Circulator_on_faces& a, const Circulator_on_faces& b) { return a.m_ptr != b.m_ptr; }

        int globalIdx() { return m_index; }
    private:
        Mesh &m_mesh;
        pointer m_ptr;
        std::size_t m_center_idx;
        int m_index;
    };

    struct Circulator_on_vertices
    {
//        using iterator_category = std::forward_iterator_tag;
//        using difference_type   = std::ptrdiff_t;
        using value_type        = Vertex;
        using pointer           = Vertex*;
        using reference         = Vertex&;
        Circulator_on_vertices(Mesh &mesh, pointer ptr, std::size_t center_idx, Circulator_on_faces cof, uint vertexIdx): m_mesh(mesh),
            m_ptr(ptr), m_center_idx(center_idx), m_cof(cof), m_vertexIdx(vertexIdx){
        }

        reference operator*() const { return *m_ptr; }
        pointer operator->() { return m_ptr; }

        // Prefix increment
        Circulator_on_vertices& operator++() {
            ++m_cof;
            int next = m_cof->getInternalIdx(m_center_idx, 1);
            m_vertexIdx = m_cof->vertices[next];
            m_ptr = &(m_mesh.vertices[m_vertexIdx]);
            return *this;
        }

        friend bool operator== (const Circulator_on_vertices& a, const Circulator_on_vertices& b) {return a.m_ptr == b.m_ptr; }
        friend bool operator!= (const Circulator_on_vertices& a, const Circulator_on_vertices& b) {return a.m_ptr != b.m_ptr; }

        const Triangle &getCurrentFace() const { return *m_cof; }

        int globalVertexIdx() { return m_vertexIdx; }
        int globalFaceIdx() { return m_cof.globalIdx(); }
    private:
        Mesh &m_mesh;
        pointer m_ptr;
        std::size_t m_center_idx;
        Circulator_on_faces m_cof;
        int m_vertexIdx;
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

    Circulator_on_faces incident_faces(std::size_t vIdx) {
        return Circulator_on_faces(*this, &(triangles[vertices[vIdx].triangleIdx]), vIdx, vertices[vIdx].triangleIdx);
    }

    Circulator_on_vertices adjacent_vertices(std::size_t vIdx) {
        auto cof = incident_faces(vIdx);
        int next = cof->getInternalIdx(vIdx, 1);
        uint globalVertexIdx = cof->vertices[next];
        Vertex *adjVert = &vertices[globalVertexIdx];
        return Circulator_on_vertices(*this, adjVert, vIdx, cof, globalVertexIdx);
    }

    float area(const Triangle &t);
    void loadOFF (std::string path);
    void findTopology();
    void findTopology(const std::vector<Point> &points, const std::vector<std::array<uint, 3> > &faces);
    void calculateLaplacian();
    uint splitTriangle(uint indFace, const Point &p);
    uint splitTriangleMiddle(int indFace);
    void edgeFlip(int indFace1, int indFace2);
    float orientation2D(Point p1, Point p2, Point p3) const;
    float orientation2D(int i1, int i2, int i3) const;
    float orientation2D(const Triangle &t) const;
    bool isInside(const Point &p, const Triangle &t) const;
    void insertPoint2D(Point p);
    bool is2D(int indF);

    void drawMesh(bool wireframe = false);
    void drawMeshLaplacian(bool wireframe = false);
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
