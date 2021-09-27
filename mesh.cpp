#include "mesh.h"
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

// Draw a Point
void glPointDraw(const Point & p) {
    glVertex3f(p._x, p._y, p._z);
}

Mesh::Mesh(){
    initFile("../2D_mesh_test.off");
//    init2dBBox();
//    points.push_back(Point(.25, 0, .25));
//    vertices.push_back(4);
//    splitTriangle(0,4);
}

void Mesh::initTetrahedron(){
    clearData();
    points.push_back(Point());
    points.push_back(Point(0.0, 0.0, 1.0));
    points.push_back(Point(1.0, 0.0, 0.0));
    points.push_back(Point(0.5, 1.0, 0.5));

    for(int i=0; i<points.size(); i++){
        vertices.push_back(Vertex(i));
    }

    faces.push_back(Face(0, 2, 1));
    faces.push_back(Face(0, 3, 2));
    faces.push_back(Face(0, 1, 3));
    faces.push_back(Face(1, 2, 3));
}

void Mesh::initPyramid(){
    clearData();
    points.push_back(Point());
    points.push_back(Point(0.0, 0.0, 1.0));
    points.push_back(Point(1.0, 0.0, 0.0));
    points.push_back(Point(1.0, 0.0, 1.0));
    points.push_back(Point(0.5, 1.0, 0.5));

    for(int i=0; i<points.size(); i++){
        vertices.push_back(Vertex(i));
    }

    faces.push_back(Face(0, 2, 1));
    faces.push_back(Face(1, 2, 3));
    faces.push_back(Face(0, 1, 4));
    faces.push_back(Face(1, 3, 4));
    faces.push_back(Face(3, 2, 4));
    faces.push_back(Face(2, 0, 4));
}

void Mesh::init2dBBox(){
    clearData();
    points.push_back(Point());
    points.push_back(Point(0.0, 0.0, 1.0));
    points.push_back(Point(1.0, 0.0, 0.0));
    points.push_back(Point(1.0, 0.0, 1.0));

    for(int i=0; i<points.size(); i++){
        vertices.push_back(Vertex(i));
    }

    faces.push_back(Face(0, 2, 1, -1, -1, -1));
    faces.push_back(Face(1, 2, 3, -1, -1, -1));

    handleFace(0, 2, 1, 0);
    handleFace(1, 2, 3, 1);
}


void Mesh::initFile(string filepath){
    clearData();
    int nbVertices,nbFaces;
    ifstream file;
    string line, word1, word2, word3, word4;
    stringstream sline;

    cout << "initiating : " << filepath << endl;
    file.open(filepath);
    if(file.is_open()){
        getline(file, line);
        sline = stringstream(line);
        // reading number of faces
        getline(sline, word1, ' ');
        nbVertices = stoi(word1);
        // reading number of vertices
        getline(sline, word1, ' ');
        nbFaces = stoi(word1);

        // reading vertices coordinates
        for(int i=0; i<nbVertices; i++){
            getline(file, line);
            sline = stringstream(line);
            getline(sline, word1, ' ');
            getline(sline, word2, ' ');
            getline(sline, word3, ' ');

            points.push_back(Point(stof(word1), stof(word2), stof(word3)));
        }

        // creating vertices vector
        for(int i=0; i<points.size(); i++){
            vertices.push_back(Vertex(i, -1));
        }

        // reading description of faces (only triangles for the moment)
        for(int i=0; i<nbFaces; i++){
            getline(file, line);
            sline = stringstream(line);
            getline(sline, word1, ' ');
            getline(sline, word2, ' ');
            getline(sline, word3, ' ');
            getline(sline, word4, ' ');
            faces.push_back(Face(stoi(word2), stoi(word3), stoi(word4), -1, -1, -1));

            handleFace(stoi(word2), stoi(word3), stoi(word4), i);
        }
    }
}

void Mesh::handleFace(int indVertex1, int indVertex2, int indVertex3, int faceIndex){
    handleEdge(indVertex1, indVertex2, faceIndex, 2);
    handleEdge(indVertex2, indVertex3, faceIndex, 0);
    handleEdge(indVertex3, indVertex1, faceIndex, 1);
}

void Mesh::handleEdge(int indVertex1, int indVertex2, int faceIndex, int opposedVertex){
    // if edge in map
    if(topology.find(make_pair(indVertex2, indVertex1)) != topology.end()){
        faces[faceIndex].adjFaces[opposedVertex] = topology[make_pair(indVertex2, indVertex1)].first;
        faces[topology[make_pair(indVertex2, indVertex1)].first].adjFaces[topology[make_pair(indVertex2, indVertex1)].second] = faceIndex;
    }
    else{
        topology.insert(make_pair(make_pair(indVertex1, indVertex2), make_pair(faceIndex, opposedVertex)));
    }
}

void Mesh::clearData(){
    points.clear();
    vertices.clear();
    faces.clear();
    topology.clear();
}

void Mesh::splitTriangle(int indFace, int indVertex){
    int nbFace = faces.size();

    // Creating sub triangles
    faces.push_back(Face(faces[indFace].vertices[1], faces[indFace].vertices[2], indVertex,
            nbFace + 1, indFace, faces[indFace].adjFaces[0]));

    faces.push_back(Face(faces[indFace].vertices[2], faces[indFace].vertices[0], indVertex,
            indFace, nbFace, faces[indFace].adjFaces[1]));

    faces[indFace].vertices[2] = indVertex;
    faces[indFace].adjFaces[0] = nbFace;
    faces[indFace].adjFaces[1] = nbFace + 1;

    // Updating topology
    handleFace(faces[indFace].vertices[0], faces[indFace].vertices[1], faces[indFace].vertices[2], indFace);
    handleFace(faces[nbFace].vertices[0], faces[nbFace].vertices[1], faces[nbFace].vertices[2], nbFace);
    handleFace(faces[nbFace + 1].vertices[0], faces[nbFace + 1].vertices[1], faces[nbFace + 1].vertices[2], nbFace + 1);
}

void Mesh::edgeFlip(int indFace1, int indFace2){
    int sommetF1, sommetF2, indSommetF1, indSommetF2;
    for(int i=0; i<3; i++){
        if(faces[indFace1].adjFaces[i] == indFace2){
            sommetF1 = faces[indFace1].vertices[i];
            indSommetF1 = i;
        }
        if(faces[indFace2].adjFaces[i] == indFace1){
            sommetF2 = faces[indFace2].vertices[i];
            indSommetF2 = i;
        }
    }
    topology.erase(make_pair(faces[indFace1].vertices[(indSommetF1+1)%3], faces[indFace1].vertices[(indSommetF1+2)%3]));
    topology.erase(make_pair(faces[indFace1].vertices[(indSommetF1+2)%3], faces[indFace1].vertices[(indSommetF1+1)%3]));

    faces[indFace1].vertices[(indSommetF1+2) % 3] = sommetF2;
    faces[indFace2].vertices[(indSommetF2+2) % 3] = sommetF1;
    // update adjacence
    faces[indFace1].adjFaces[(indSommetF1+2) % 3] = indFace2;
    faces[indFace2].adjFaces[(indSommetF2+2) % 3] = indFace1;
    int temp = faces[indFace1].adjFaces[(indSommetF1+2) % 3];
    faces[indFace1].adjFaces[indSommetF1] = faces[indFace2].adjFaces[(indSommetF2+2) % 3];
    faces[indFace2].adjFaces[indSommetF2] = temp;
}

// cross / vectoriel
float Mesh::orientaion2D(Point p1, Point p2, Point p3){
    int x1 = p2._x - p1._x;
    int y1 = p2._y - p1._y;
    int x2 = p3._x - p1._x;
    int y2 = p3._y - p1._y;

    return x1*y2 - x2*y1;
}

bool Mesh::isInside(Point p, Face f){
    return (orientaion2D(points[f.vertices[0]], points[f.vertices[1]], p) > 0 &&
            orientaion2D(points[f.vertices[1]], points[f.vertices[2]], p) > 0 &&
            orientaion2D(points[f.vertices[2]], points[f.vertices[0]], p) > 0);
}

void insertPoint(Point p){

}

void Mesh::drawMesh(){
    glColor3d(1,1,1);

    for ( const auto & face : faces){
        glBegin(GL_TRIANGLES);
        glPointDraw(points[face.vertices[0]]);
        glPointDraw(points[face.vertices[1]]);
        glPointDraw(points[face.vertices[2]]);
        glEnd();
    }
}

void Mesh::drawMeshWireFrame(){
    glColor3d(1,1,1);

    for ( const auto & face : faces){
        glBegin(GL_LINE_STRIP);
        glPointDraw(points[face.vertices[0]]);
        glPointDraw(points[face.vertices[1]]);
        glEnd();

        glBegin(GL_LINE_STRIP);
        glPointDraw(points[face.vertices[0]]);
        glPointDraw(points[face.vertices[2]]);
        glEnd();

        glBegin(GL_LINE_STRIP);
        glPointDraw(points[face.vertices[1]]);
        glPointDraw(points[face.vertices[2]]);
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

// The following functions could be displaced into a module OpenGLDisplayGeometricWorld that would include mesh.h

//Example with a bBox
void GeometricWorld::draw() {
    _mesh.drawMesh();

}

//Example with a wireframe bBox
void GeometricWorld::drawWireFrame() {
    _mesh.drawMeshWireFrame();
}

