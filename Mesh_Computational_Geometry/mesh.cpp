#include "mesh.h"
#include <iostream>
#include <set>
#include <QFileDialog>
#include <QDebug>
#include <QDir>
#include <QVector3D>
#include <QtMath>
#include <sstream>

// The following functions could be displaced into a module OpenGLDisplayGeometricWorld that would include mesh.h

// Draw a Point
void glPointDraw(const Point & p) {
    glVertex3f(p._x, p._y, p._z);
}

Triangle::Triangle() {}

size_t Triangle::getInternalIdx(size_t vertexIdx, int shift = 0) const
{
    size_t i;
    for (i = 0; i < 3; ++i)
        if (vertices[i] == vertexIdx) break;
    return (i + shift) % 3;
}

float Mesh::area(const Triangle &t) {

    QVector3D ab = vertices[t.vertices[1]].p - vertices[t.vertices[0]].p;
    QVector3D ac = vertices[t.vertices[2]].p - vertices[t.vertices[0]].p;
    return 0.5 * QVector3D::crossProduct(ab, ac).length();
}

Mesh::Mesh() {
    loadOFF("../data/queen.off");
//    loadOFF("../data/sphere.off");
//    loadOFF("../data/cube.off");

    calculateLaplacian();
    // Circulator on vertex test
//    Circulator_on_vertices cov = adjacent_vertices(0);
//    Iterator_on_vertices its = vertices_begin();
//    for (; its != vertices_past_the_end(); ++its) {
//        int cmpt=0;
//        auto cv=adjacent_vertices(its.getIdx());
//        ++cv;
//        auto end = adjacent_vertices(its.getIdx());
//        for (; cv!=end; ++cv)
//            cmpt++;
//        std::cout<< "valence of the vertex : " << cmpt << std::endl;
//    }
    // Circulator on faces test
//    Iterator_on_vertices its = vertices_begin();
//    for (; its != vertices_past_the_end(); ++its) {
//        int cmpt=0;
//        auto cf=incident_faces(its.getIdx());
//        ++cf;
//        auto end = incident_faces(its.getIdx());
//        for (; cf!=end; ++cf)
//            cmpt++;
//        std::cout<< "valence of the vertex : " << cmpt << std::endl;
//    }

}


std::pair<int,int> make_ordered_pair(int a, int b) {
    return (a < b ? std::make_pair(a,b) : std::make_pair(b,a));
}

int findOppositeIdx(int a, int b) {
    std::set<int> s = {0,1,2};
    s.erase(a);
    s.erase(b);
    return *s.begin();
}

void Mesh::loadOFF(std::string path)
{
    vertices.clear();
    triangles.clear();
    std::ifstream off(path);
    if (off.is_open()) {
        std::string line;
        uint nbVertices, nbFaces;
        std::getline(off,line);
        // Getting vertex and face info from file
        std::istringstream ss(line);
        ss >> nbVertices;
        ss >> nbFaces;

        std::vector<std::string> strArr(5);

        // Iterating on each line of the file
        while (std::getline(off,line)) {
            strArr.clear();
            ss.clear();
            ss.str(line);
            std::string value;
            // populate array of strings
            while (ss >> value) {
                strArr.push_back(value);
            }

            if (strArr.size() == 3) { // it's a vertex description line
                Vertex v;
                v.p._x = std::stof(strArr[0]);
                v.p._y = std::stof(strArr[1]);
                v.p._z = std::stof(strArr[2]);
                v.triangleIdx = -1;
                vertices.push_back(v);
            }
            else if (strArr.size() == 4) { // it's a face index line
                Triangle t;
                t.vertices[0] = std::stoi(strArr[1]);
                t.vertices[1] = std::stoi(strArr[2]);
                t.vertices[2] = std::stoi(strArr[3]);
                triangles.push_back(t);
            }
        }
        off.close();

        findTopology();
    } else std::cout << "Ouverture du fichier impossible." << std::endl;
}

void Mesh::findTopology() {
    std::map<std::pair<int,int>, std::pair<int,int>> topo; // EDGE -> FACE

    for (std::size_t i = 0 ; i < triangles.size() ; ++i) {
        for (int j = 0 ; j < 3 ; ++j) {
            uint vertexIdx = triangles[i].vertices[j];
            if (vertices[vertexIdx].triangleIdx == -1)
                vertices[vertexIdx].triangleIdx = i;	// Affectation de l'id du triangle courant aux sommets n'étant pas attachés à un triangle
            // Créer une paire de sommets avec les id dans l'ordre croissant
            int next = j < 2 ? j+1 : 0;
            uint nextVertexIdx = triangles[i].vertices[next];
            int oppositeIdx = findOppositeIdx(j,next);
            const auto pair = topo.insert({make_ordered_pair(vertexIdx, nextVertexIdx),
                                         std::make_pair(i, oppositeIdx)});
            if (pair.second == false) { // Si l'insertion n'a pas pu être faite dans la map
                int otherFaceIdx = pair.first->second.first;
                triangles[i].adjacent[oppositeIdx] = otherFaceIdx;
                triangles[otherFaceIdx].adjacent[pair.first->second.second] = i;
            }
        }
    }
}

void Mesh::findTopology(const std::vector<Point> &points, const std::vector<std::array<uint,3>> &faces) {
    std::map<std::pair<int,int>, std::pair<int,int>> topo; // EDGE -> FACE
    for (const auto &point : points) {
        Vertex v;
        v.p = point;
        v.triangleIdx = -1;
        vertices.push_back(v);
    }

    // Parcours de chaque triangle
    for (std::size_t i = 0 ; i < faces.size(); ++i) {
        Triangle t;
        for (int j = 0 ; j < 3 ; ++j) {
            uint vertexIdx = faces[i][j];
            if (vertices[vertexIdx].triangleIdx == -1)
                vertices[vertexIdx].triangleIdx = i;	// Affectation de l'id du triangle courant aux sommets n'étant pas attachés à un triangle
            t.vertices[j] = vertexIdx;
            // Créer une paire de sommets avec les id dans l'ordre croissant
            int next = j < 2 ? j+1 : 0;
            uint nextVertexIdx = faces[i][next];
            int oppositeIdx = findOppositeIdx(j,next);
            const auto pair = topo.insert({make_ordered_pair(vertexIdx, nextVertexIdx),
                                         std::make_pair(i, oppositeIdx)});
            if (pair.second == false) { // Si l'insertion n'a pas pu être faite dans la map
                int otherFaceIdx = pair.first->second.first;
                t.adjacent[oppositeIdx] = otherFaceIdx;
                triangles[otherFaceIdx].adjacent[pair.first->second.second] = i;
            }
        }
        triangles.push_back(t);
    }
}

void Mesh::calculateLaplacian()
{
    for (auto vIt = vertices_begin(); vIt != vertices_past_the_end(); ++vIt) {
        const Point &ui = vIt->p;
        double areas = 0.;
        // Calculating area
        for (auto tCir = ++incident_faces(vIt.getIdx()); tCir != incident_faces(vIt.getIdx()); ++tCir ) {
            areas += area(*tCir);
        }
        areas /= 3;

        double laplacianX = 0., laplacianY = 0., laplacianZ = 0.;
        for (auto vCir = ++adjacent_vertices(vIt.getIdx()); vCir != adjacent_vertices(vIt.getIdx()); ++vCir) {
            double cotAij, cotBij;
            const Triangle &cF = vCir.getCurrentFace();
            int u0 = cF.getInternalIdx(vIt.getIdx(),1);
            int u1 = (u0 + 1) % 3;
            int v0 = u0;
            int v1 = (v0 + 2) % 3;
            QVector3D u = vertices[cF.vertices[u1]].p - vertices[cF.vertices[u0]].p;
            QVector3D v = vertices[cF.vertices[v1]].p - vertices[cF.vertices[v0]].p;
            cotAij = QVector3D::dotProduct(u,v) / QVector3D::crossProduct(u,v).length();

            const Triangle &nF = vCir.getNextFace();
            u0 = nF.getInternalIdx(vIt.getIdx(),2);
            u1 = (u0 + 2) % 3;
            v0 = u0;
            v1 = (v0 + 1) % 3;
            u = vertices[nF.vertices[u1]].p - vertices[nF.vertices[u0]].p;
            v = vertices[nF.vertices[v1]].p - vertices[nF.vertices[v0]].p;
            cotBij = QVector3D::dotProduct(u,v) / QVector3D::crossProduct(u,v).length();

            laplacianX += (cotAij + cotBij)*(vCir->p._x - ui._x);
            laplacianY += (cotAij + cotBij)*(vCir->p._y - ui._y);
            laplacianZ += (cotAij + cotBij)*(vCir->p._z - ui._z);
        }
        laplacianX = (1 / (2 * areas)) * laplacianX;
        laplacianY = (1 / (2 * areas)) * laplacianY;
        laplacianZ = (1 / (2 * areas)) * laplacianZ;
//        std::cout << "laplacianX: " << laplacianX << std::endl;
        QVector3D l(laplacianX,laplacianY,laplacianZ);
        laplacians.push_back(l);
        curvature.push_back(l.length() / 2);
    }
    for (const auto &vec:laplacians)
        qDebug() << vec.length() / 2;
    // normalize curvature
    const auto pair = std::minmax_element(curvature.begin(), curvature.end());
    double min = *pair.first;
    double max = *pair.second;
    for (auto &d:curvature) {
        std::cout << "dBefore: " << d << std::endl;
        d = (d-min) / max * 10;
        std::cout << "dAfter: " << d << std::endl;
    }
}

void Mesh::drawMesh(bool wireframe) {
    for (const auto & face : triangles) {
        if (wireframe) {
           glBegin(GL_LINE_STRIP);
        } else {
            glBegin(GL_TRIANGLES);
        }
        double v = curvature[face.vertices[0]];
        glColor3d(v,v,v);
        glPointDraw(vertices[face.vertices[0]].p);
        v = curvature[face.vertices[1]];
        glColor3d(v,v,v);
        glPointDraw(vertices[face.vertices[1]].p);
        v = curvature[face.vertices[2]];
        glColor3d(v,v,v);
        glPointDraw(vertices[face.vertices[2]].p);
        glEnd();
    }
}

GeometricWorld::GeometricWorld()
{
    double width=0.5, depth=0.6, height=0.8;
    _bBox.push_back(Point(-0.5*width,-0.5*depth,-0.5*height)); //0
    _bBox.push_back(Point(-0.5*width,0.5*depth,-0.5*height)); // 1
    _bBox.push_back(Point(0.5*width,-0.5*depth,-0.5*height)); // 2
    _bBox.push_back(Point(-0.5*width,-0.5*depth,0.5*height)); // 3
}


//Example with a bBox
void GeometricWorld::draw() {
//    _mesh.drawMesh(true);
    _mesh.drawMesh(false);
}

//Example with a wireframe bBox
void GeometricWorld::drawWireFrame() {
    glColor3d(0,1,0);
    glBegin(GL_LINE_STRIP);
    glPointDraw(_bBox[0]);
    glPointDraw(_bBox[1]);
    glEnd();
    glColor3d(0,0,1);
    glBegin(GL_LINE_STRIP);
    glPointDraw(_bBox[0]);
    glPointDraw(_bBox[2]);
    glEnd();
    glColor3d(1,0,0);
    glBegin(GL_LINE_STRIP);
    glPointDraw(_bBox[0]);
    glPointDraw(_bBox[3]);
    glEnd();
}


QVector3D Point::operator -(const Point &p)
{
    return QVector3D(p._x - _x, p._y - _y, p._z - _z);
}
