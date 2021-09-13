#include "mesh.h"
#include <iostream>
#include <set>
#include <QFileDialog>
#include <QDebug>
#include <QDir>
#include <sstream>

// The following functions could be displaced into a module OpenGLDisplayGeometricWorld that would include mesh.h

// Draw a Point
void glPointDraw(const Point & p) {
    glVertex3f(p._x, p._y, p._z);
}

Triangle::Triangle() {}

Mesh::Mesh() {
    loadOFF("../data/queen.off");
//    loadOFF("/home/mikail/Downloads/cube.off");
    Iterator_on_vertices its = vertices_begin();
    for (; its != vertices_past_the_end(); ++its) {
        int cmpt=0;
        auto cf=incident_faces(its.getIdx());
        ++cf;
        auto end = incident_faces(its.getIdx());
        for (; cf!=end; ++cf)
            cmpt++;
        std::cout<< "valence of the vertex : " << cmpt << std::endl;
    }

//    loadOFF("/home/mikail/Downloads/sphere.off");
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

void Mesh::drawMesh(bool wireframe) {
    int i = 0;
    for (const auto & face : triangles) {
        if (wireframe) {
           glBegin(GL_LINE_STRIP);
        } else {
            glColor3d(1,0,1);


            glBegin(GL_TRIANGLES);
        }

        glPointDraw(vertices[face.vertices[0]].p);
        glPointDraw(vertices[face.vertices[1]].p);
        glPointDraw(vertices[face.vertices[2]].p);
        glEnd();
        ++i;
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
    _mesh.drawMesh(true);
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

