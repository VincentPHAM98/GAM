#ifndef MESH_H
#define MESH_H

#include <math.h>

#include <QDebug>
#include <QGLWidget>
#include <QVector3D>
#include <fstream>
#include <iterator>
#include <queue>
#include <set>
#include <string>
#include <vector>

using namespace std;

class Point {
   public:
    double _x = 0.;
    double _y = 0.;
    double _z = 0.;

    Point() : _x(0.), _y(0.), _z(0.) {}
    Point(const Point& p) : _x(p._x), _y(p._y), _z(p._z) {}
    Point& operator=(const Point& p) {
        _x = p._x;
        _y = p._y;
        _z = p._z;
        return *this;
    }
    Point(double x_, double y_, double z_) : _x(x_), _y(y_), _z(z_) {}
    Point(const QVector3D& v) : _x(v.x()), _y(v.y()), _z(v.z()) {}

    /**
     * @brief cross / produit vectoriel
     * @param v
     * @return résultat du produit vectoriel avec v
     */
    Point cross(Point v);

    /**
     * @brief dot / produit scalaire
     * @param v
     * @return résultat du produit sclaire avec v
     */
    double dot(Point v);

    /**
     * @brief normalize / normalisation du vecteur
     * @return vecteur normalisé
     */
    Point normalize();

    /**
     * @brief length
     * @return norme du vecteur
     */
    double length();

    /**
     * @brief tangente
     * @param a
     * @param b
     * @param c
     * @return tangente de ABC
     */
    double tangente(Point a, Point b, Point c);

    Point operator+(const Point& p) {
        double x = _x + p._x;
        double y = _y + p._y;
        double z = _z + p._z;
        return Point(x, y, z);
    }

    Point operator/(int i) {
        double x = _x / i;
        double y = _y / i;
        double z = _z / i;
        return Point(x, y, z);
    }

    QVector3D operator-(const Point& p);

    Point operator*(const float f) {
        Point res;
        res._x = _x * f;
        res._y = _y * f;
        res._z = _z * f;
        return res;
    }

    friend std::ostream& operator<<(std::ostream& os, const Point& p);
};

class Vertex {
   public:
    Point p;
    int triangleIdx = -1;
    bool isDeleted = false;

    Vertex() : p() {}
    Vertex(const Point& _p) : p(_p) {}
    Vertex(const Point& _p, int id) : p(_p), triangleIdx(id) {}

    void remove() { isDeleted = true; }
};

class Triangle {
   public:
    bool isDeleted = false;
    std::array<uint, 3>
        vertices;  // Trigonometric order
    std::array<uint, 3> adjacent;
    // Constraint: vertex i facing adjacent triangle i
    Triangle();
    Triangle(const Triangle& t) : vertices(t.vertices), adjacent(t.adjacent) {}
    Triangle(std::array<uint, 3> _vertices) : vertices(_vertices) {}
    Triangle(std::array<uint, 3> _vertices, std::array<uint, 3> _adjacent) : vertices(_vertices), adjacent(_adjacent) {}
    //    ~Triangle();
    /**
     * @brief getInternalIdx
     * @param vertexIdx
     * @param shift
     * @return Gives internal index of given vertex index
     */
    int getInternalIdx(size_t vertexIdx, int shift = 0) const;

    /**
     * @brief getAdjacentFaceFromGlobalVertex
     * @param vertexIdx
     * @return Gives adjacent face from global vertex
     */
    int getAdjacentFaceFromGlobalVertex(int vertexIdx) const;

    /**
     * @brief getVertexFromAdjacentFace
     * @param faceIdx
     * @return Gives vertex from adjacent face
     */
    uint& getVertexFromAdjacentFace(uint faceIdx);
    void remove() { isDeleted = true; }
};

class Mesh {
    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;
    std::vector<QVector3D> laplacians;
    std::vector<double> curvature;

    std::map<uint, Point> voronoiCenter;

   public:
    Mesh();
    //~Mesh();
    struct Iterator_on_vertices {
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = Vertex;
        using pointer = Vertex*;
        using reference = Vertex&;
        Iterator_on_vertices(Mesh& m) : m_mesh(m), m_ptr(nullptr), currentIdx(-1) {}
        Iterator_on_vertices(const Iterator_on_vertices& copy) : m_mesh(copy.m_mesh), m_ptr(copy.m_ptr), currentIdx(copy.currentIdx) {}
        Iterator_on_vertices(Mesh& m, pointer ptr) : m_mesh(m), m_ptr(ptr) {}

        reference operator*() const { return *m_ptr; }
        pointer operator->() { return m_ptr; }

        // Prefix increment
        Iterator_on_vertices& operator++() {
            do {
                m_ptr++;
                currentIdx++;
            } while (m_ptr->isDeleted && currentIdx < m_mesh.vertices.size());
            return *this;
        }

        // Postfix increment
        Iterator_on_vertices operator++(int) {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        std::size_t getIdx() { return currentIdx; }

        friend bool operator==(const Iterator_on_vertices& a, const Iterator_on_vertices& b) { return a.m_ptr == b.m_ptr; }
        friend bool operator!=(const Iterator_on_vertices& a, const Iterator_on_vertices& b) { return a.m_ptr != b.m_ptr; }

       private:
        Mesh& m_mesh;
        pointer m_ptr;
        std::size_t currentIdx = 0;
    };

    struct Iterator_on_faces {
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = Triangle;
        using pointer = Triangle*;
        using reference = Triangle&;
        Iterator_on_faces(Mesh& m) : m_mesh(m), m_ptr(nullptr), currentIdx(-1) {}
        Iterator_on_faces(const Iterator_on_faces& copy) : m_mesh(copy.m_mesh), m_ptr(copy.m_ptr), currentIdx(copy.currentIdx) {}
        Iterator_on_faces(Mesh& m, pointer ptr) : m_mesh(m), m_ptr(ptr) {}

        reference operator*() const { return *m_ptr; }
        pointer operator->() { return m_ptr; }

        // Prefix increment
        Iterator_on_faces& operator++() {
            do {
                m_ptr++;
                currentIdx++;
            } while (m_ptr->isDeleted && currentIdx < m_mesh.triangles.size());

            return *this;
        }

        // Postfix increment
        Iterator_on_faces operator++(int) {
            Iterator_on_faces tmp = *this;
            ++(*this);
            return tmp;
        }

        std::size_t getIdx() { return currentIdx; }

        friend bool operator==(const Iterator_on_faces& a, const Iterator_on_faces& b) { return a.m_ptr == b.m_ptr; }
        friend bool operator!=(const Iterator_on_faces& a, const Iterator_on_faces& b) { return a.m_ptr != b.m_ptr; }

       private:
        Mesh& m_mesh;
        pointer m_ptr;
        std::size_t currentIdx = 0;
    };

    struct Circulator_on_faces {
        using value_type = Triangle;
        using pointer = Triangle*;
        using reference = Triangle&;
        Circulator_on_faces() = delete;
        Circulator_on_faces(const Circulator_on_faces& c)
            : m_mesh(c.m_mesh), m_ptr(c.m_ptr), m_center_idx(c.m_center_idx), m_index(c.m_index) {}
        Circulator_on_faces(Mesh& mesh, pointer ptr, std::size_t center_idx, int index)
            : m_mesh(mesh), m_ptr(ptr), m_center_idx(center_idx), m_index(index) {}

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
        friend bool operator==(const Circulator_on_faces& a, const Circulator_on_faces& b) { return a.m_ptr == b.m_ptr; }
        friend bool operator!=(const Circulator_on_faces& a, const Circulator_on_faces& b) { return a.m_ptr != b.m_ptr; }

        int globalIdx() { return m_index; }

       private:
        Mesh& m_mesh;
        pointer m_ptr;
        std::size_t m_center_idx;
        int m_index;
    };

    struct Circulator_on_vertices {
        using value_type = Vertex;
        using pointer = Vertex*;
        using reference = Vertex&;
        Circulator_on_vertices(Mesh& mesh, pointer ptr, std::size_t center_idx, Circulator_on_faces cof, uint vertexIdx)
            : m_mesh(mesh),
              m_ptr(ptr),
              m_center_idx(center_idx),
              m_cof(cof),
              m_vertexIdx(vertexIdx) {
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

        friend bool operator==(const Circulator_on_vertices& a, const Circulator_on_vertices& b) { return a.m_ptr == b.m_ptr; }
        friend bool operator!=(const Circulator_on_vertices& a, const Circulator_on_vertices& b) { return a.m_ptr != b.m_ptr; }

        const Triangle& getCurrentFace() const { return *m_cof; }

        int globalVertexIdx() { return m_vertexIdx; }
        int globalFaceIdx() { return m_cof.globalIdx(); }

       private:
        Mesh& m_mesh;
        pointer m_ptr;
        std::size_t m_center_idx;
        Circulator_on_faces m_cof;
        int m_vertexIdx;
    };

    Iterator_on_vertices vertices_begin() {
        auto it = Iterator_on_vertices(*this, vertices.data());
        if (it->isDeleted)
            ++it;
        return it;
    }
    Iterator_on_vertices vertices_past_the_end() {
        return Iterator_on_vertices(*this, &(*vertices.end()));
    }
    Iterator_on_faces faces_begin() {
        auto it = Iterator_on_faces(*this, triangles.data());
        if (it->isDeleted)
            ++it;
        return it;
    }
    Iterator_on_faces faces_past_the_end() {
        return Iterator_on_faces(*this, &(*triangles.end()));
    }

    Circulator_on_faces incident_faces(std::size_t vIdx) {
        return Circulator_on_faces(*this, &(triangles[vertices[vIdx].triangleIdx]), vIdx, vertices[vIdx].triangleIdx);
    }

    Circulator_on_vertices adjacent_vertices(std::size_t vIdx) {
        auto cof = incident_faces(vIdx);
        int next = cof->getInternalIdx(vIdx, 1);
        uint globalVertexIdx = cof->vertices.at(next);
        Vertex* adjVert = &vertices.at(globalVertexIdx);
        return Circulator_on_vertices(*this, adjVert, vIdx, cof, globalVertexIdx);
    }

    float area(const Triangle& t);

    /**
     * @brief loads off file into data structure
     * @param path
     */
    void loadOFF(std::string path);

    /**
     * @brief findTopology
     */
    void findTopology();

    /**
     * @brief compute laplacian for the mesh
     */
    void calculateLaplacian();

    /**
     * @brief splitTriangle into 3 sub triangles at location of indVertex
     * @param indFace
     * @param indVertex
     * @param delaunay
     */
    void splitTriangle(int indFace, int indVertex, bool delaunay);
    uint splitTriangle(uint indFace, const Point& p, bool delaunay);
    uint splitTriangleMiddle(int indFace, bool delaunay);

    /**
     * @brief edgeFlip edge between indFace1 and inFace2 if possible, computes local delaunay if need
     * @param indFace1
     * @param indFace2
     * @param delaunay
     */
    void edgeFlip(int indFace1, int indFace2, bool delaunay);

    /**
     * @brief orientation2D
     * @param p1
     * @param p2
     * @param p3
     * @return orientation of triangle
     */
    float orientation2D(Point p1, Point p2, Point p3) const;
    float orientation2D(int i1, int i2, int i3) const;
    float orientation2D(const Triangle& t) const;

    /**
     * @brief getNormal
     * @param a
     * @param b
     * @param c
     * @return normal of triangle
     */
    Point getNormal(Point a, Point b, Point c);
    Point getNormal(Triangle f);
    Point getNormal(int idFace);

    /**
     * @brief isInside
     * @param p
     * @param t
     * @return bool p inside t
     */
    bool isInside(const Point& p, const Triangle& t) const;

    /**
     * @brief inserts p in mesh if possible, computes local delaunay if needed
     * @param p
     * @param delaunay
     */
    void insertPoint2D(const Point& p, bool delaunay);
    void insertRandPoint2D(int max, bool delaunay);

    // checks if elements are 2D (not defined on z axis)
    bool is2D(int indF);
    bool isVert2D(int indV);
    bool isFace2D(int indF);

    /**
     * @brief findFacesWithCommonEdge of idVert1 and idVert2
     * @param idVert1
     * @param idVert2
     * @return common edge
     */
    std::pair<int, int> findFacesWithCommonEdge(uint idVert1, uint idVert2);
    std::set<int> adjacentVerticesOfVertex(uint indV);
    void changeIncidentFacesOfFaceVertices(uint idFace, std::pair<int, int> deletedFaces);
    int collapseEdge(uint idVert1, uint idVert2);
    void collapseShortestEdge();

    // DRAW FUNCTIONS
    void drawMesh();
    void drawMeshWireFrame();
    void drawMeshLaplacian(bool wireframe = false);
    void drawVoronoi();

    void test();

    int currentFace = 0;
    bool highlightNeighbors;
    uint selectedVertex1 = 1, selectedVertex2 = 2;

    // complete l'enveloppe convexe avec des edge flip
    void completeConvexHull(int idFace, int idVert, bool delaunay);
    // utility functions to find indexes
    int findAdjFace(int idFace, int id2find);
    int vertIndexInFace(int idFace, int idVert);
    int infiniteInFace(int idFace);

    // retourne l'id du sommet dans idFace2 opposé à idFace1
    int opposedVert(int idFace1, int idFace2);

    /**
     * @brief isInCirconscrit
     * @param idFace
     * @param p
     * @return bool p inside face idFace
     */
    bool isInCirconscrit(int idFace, Point p);

    /**
     * @brief checkFaceDelaunayGlobal, chekcs one face and its neighbors if they are delaunay
     * @param nonDelaunay
     * @param idFace
     */
    void checkFaceDelaunayGlobal(queue<pair<int, int>>& nonDelaunay, int idFace);

    /**
     * @brief checkFaceDelaunayLocal, chekcs one face and its neighbors if they are delaunay and checks them if they are not
     * @param nonDelaunay
     * @param idFace
     */
    void checkFaceDelaunayLocal(queue<pair<int, int>>& nonDelaunay, int idFace);

    /**
     * @brief makeDelaunay, update to make mesh Delaunay
     */
    void makeDelaunay();

    /**
     * @brief localDelaunay, update to make face delaunay locally and its neighbors
     * @param idF
     */
    void localDelaunay(int idF);

    /**
     * @brief tangente
     * @param a
     * @param b
     * @param c
     * @return tangente of ABC
     */
    float tangente(Point a, Point b, Point c);

    /**
     * @brief centreCercleCirconscrit
     * @param idF
     * @return coordinates of CCC
     */
    Point centreCercleCirconscrit(int idF);

    /**
     * @brief computeVoronoi, computes voronoi cells and stores them
     */
    void computeVoronoi();

    friend class MainWindow;
};

class GeometricWorld  //Generally used to create a singleton instance
{
    QVector<Point> _bBox;  // Bounding box // ou std::vector
   public:
    GeometricWorld();
    void draw();
    void drawWireFrame();
    void drawVoronoi();
    // ** TP Can be extended with further elements;
    Mesh _mesh;
};

#endif  // MESH_H
