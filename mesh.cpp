#include "mesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <time.h>

using namespace std;

// Draw a Point
void glPointDraw(const Point & p) {
//    cout << p._x << " " << p._y << " " << p._z << endl;
    glVertex3f(p._x, p._y, p._z);
}

Mesh::Mesh(){
    srand (time(NULL));
    currentFace = 0;
    highlightNeighbors = 0;

    cout << orientation2D(Point(2,2,0), Point(1,1,0), Point(0,2,0)) << endl;
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
    points.push_back(Point(0.0, 1.0, 0.0));
    points.push_back(Point(1.0, 0.0, 0.0));
    points.push_back(Point(1.0, 1.0, 0.0));

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
        cout << "file opened succesfully" << endl;
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

        // reading description of faces (only triangles)
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
    cout << "init end" << endl;

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
    // update des voisins
    faces[faces[indFace].adjFaces[0]].adjFaces[findAdjFace(faces[indFace].adjFaces[0], indFace)] = nbFace;
    faces[faces[indFace].adjFaces[1]].adjFaces[findAdjFace(faces[indFace].adjFaces[1], indFace)] = nbFace + 1;

    // Creating sub triangles
    faces.push_back(Face(faces[indFace].vertices[1], faces[indFace].vertices[2], indVertex,
            nbFace + 1, indFace, faces[indFace].adjFaces[0]));

    faces.push_back(Face(faces[indFace].vertices[2], faces[indFace].vertices[0], indVertex,
            indFace, nbFace, faces[indFace].adjFaces[1]));

    faces[indFace].vertices[2] = indVertex;
    faces[indFace].adjFaces[0] = nbFace;
    faces[indFace].adjFaces[1] = nbFace + 1;
}

void Mesh::splitTriangle(int idFace, Point p){
    points.push_back(p);
    vertices.push_back(Vertex(vertices.size()));
    splitTriangle(idFace, vertices.back().pointIndex);
}

void Mesh::splitTriangleAtCenter(int idFace){
    Point p;
    Point a = points[faces[idFace].vertices[0]];
    Point b = points[faces[idFace].vertices[1]];
    Point c = points[faces[idFace].vertices[2]];

    p._x = (a._x + b._x + c._x)/3;
    p._y = (a._y + b._y + c._y)/3;
    p._z = (a._z + b._z + c._z)/3;
    splitTriangle(idFace, p);

}

void Mesh::edgeFlip(int indFace1, int indFace2){
    int sommetF1 = -1, sommetF2, indSommetF1, indSommetF2;
    // cherche l'arrête commune et stock des infos pour la suite
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
    if(sommetF1 == -1){
        cout << "Cannot flip this edge" << endl;
        return;
    }
    // update des sommets pour flipper l'arrête
    faces[indFace1].vertices[(indSommetF1+2) % 3] = sommetF2;
    faces[indFace2].vertices[(indSommetF2+2) % 3] = sommetF1;

    // update de l'adjacence des faces flippés
    int temp = faces[indFace1].adjFaces[(indSommetF1+1) % 3];
    int temp2 = faces[indFace2].adjFaces[(indSommetF2+1) % 3];
    faces[indFace1].adjFaces[indSommetF1] = temp2;
    faces[indFace2].adjFaces[indSommetF2] = temp;

    faces[indFace1].adjFaces[(indSommetF1+1) % 3] = indFace2;
    faces[indFace2].adjFaces[(indSommetF2+1) % 3] = indFace1;

    // update de l'adjacence des faces voisines flippés
    for(int i = 0; i < 3; i++){
        if (faces[temp].adjFaces[i] == indFace1){
            faces[temp].adjFaces[i] = indFace2;
        }
        if (faces[temp2].adjFaces[i] == indFace2)
            faces[temp2].adjFaces[i] = indFace1;
    }
}

// cross / vectoriel
double Mesh::orientation2D(Point p1, Point p2, Point p3){
    double uX = p2._x - p1._x;
    double uY = p2._y - p1._y;
    double vX = p3._x - p1._x;
    double vY = p3._y - p1._y;
    return uX*vY - vX*uY;
}

double Mesh::orientation2D(int v1, int v2, int v3){
    double uX = points[v2]._x - points[v1]._x;
    double uY = points[v2]._y - points[v1]._y;
    double vX = points[v3]._x - points[v1]._x;
    double vY = points[v3]._y - points[v1]._y;
    return uX*vY - vX*uY;
}

bool Mesh::isInside(Point p, Face f){
    return (orientation2D(points[f.vertices[0]], points[f.vertices[1]], p) > 0. &&
            orientation2D(points[f.vertices[1]], points[f.vertices[2]], p) > 0. &&
            orientation2D(points[f.vertices[2]], points[f.vertices[0]], p) > 0.);
}

bool Mesh::isVert2D(int indV){
    return points[indV]._z == 0.;
}

bool Mesh::isFace2D(int indF){
    return isVert2D(faces[indF].vertices[0]) &&
           isVert2D(faces[indF].vertices[1]) &&
           isVert2D(faces[indF].vertices[2]);
}

int Mesh::findAdjFace(int idFace, int id2find){
    for(int i = 0; i < 3; i++){
        if(faces[idFace].adjFaces[i] == id2find)
            return i;
    }
    return -1;
}

int Mesh::vertIndexInFace(int idFace, int idVert){
    for(int i = 0; i < 3; i++){
        if(faces[idFace].vertices[i] == idVert)
            return i;
    }
    return -1;
}

int Mesh::infiniteInFace(int idFace){
    for(int i = 0; i < 3; i++){
        if(!isVert2D(faces[idFace].vertices[i]))
            return i;
    }
    return -1;
}

void Mesh::insertPoint2D(Point p){
    cout << "début insertion" << endl;
    // on commence sur une face aléatoire pas reliée au point infini
    int indF = rand() % faces.size();
    while(!isFace2D(indF))
        indF = rand() % faces.size();

    // marche pour trouver la bonne face
    int n = 0;
    while(!isInside(p, faces[indF])){
        cout << indF << endl;
        if(orientation2D(points[faces[indF].vertices[0]], points[faces[indF].vertices[1]], p) < 0.)
            indF = faces[indF].adjFaces[2];
        else if(orientation2D(points[faces[indF].vertices[1]], points[faces[indF].vertices[2]], p) < 0.)
            indF = faces[indF].adjFaces[0];
        else if(orientation2D(points[faces[indF].vertices[2]], points[faces[indF].vertices[0]], p) < 0.)
            indF = faces[indF].adjFaces[1];
        // si on est hors de l'enveloppe
        if(!isFace2D(indF))
            break;
        if(n++ > faces.size()){
            cout << "Failed to insert point" << endl;
            return;
        }
    }
    cout << "fin marche" << endl;
    int indV = vertices.size();
    points.push_back(p);
    vertices.push_back(Vertex(indV));
    splitTriangle(indF, indV);
    cout << "split" << endl << endl;

    // TO DO finir l'enveloppe convexe si ajout en dehors
    if(!isFace2D(indF)){
        indF = faces[indF].adjFaces[infiniteInFace(indF)];
        completeConvexHull(indF, indV);
    }
    currentFace = indF;

}

void  Mesh::completeConvexHull(int idFace, int idVert){
    // Complete à droite (sens trigo)
    int tempFace = idFace;
    int face2Flip1 = faces[tempFace].adjFaces[(vertIndexInFace(tempFace, idVert)+2)%3];
    int face2Flip2 = faces[face2Flip1].adjFaces[(vertIndexInFace(face2Flip1, idVert))];
    int idV2 = faces[face2Flip1].vertices[(infiniteInFace(face2Flip1)+1)%3];
    int idV3 = faces[face2Flip2].vertices[(infiniteInFace(face2Flip2)+1)%3];

    while(orientation2D(idVert, idV2, idV3) < 0.){
        edgeFlip(face2Flip1, face2Flip2);
        tempFace = faces[tempFace].adjFaces[(vertIndexInFace(tempFace, idVert)+2)%3];
        face2Flip1 = faces[tempFace].adjFaces[(vertIndexInFace(tempFace, idVert)+2)%3];
        face2Flip2 = faces[face2Flip1].adjFaces[(vertIndexInFace(face2Flip1, idVert))];
        idV2 = faces[face2Flip1].vertices[(infiniteInFace(face2Flip1)+1)%3];
        idV3 = faces[face2Flip2].vertices[(infiniteInFace(face2Flip2)+1)%3];
    }

    // Complete à gauche (sens non trigo)
    tempFace = idFace;
    face2Flip1 = faces[tempFace].adjFaces[(vertIndexInFace(tempFace, idVert)+1)%3];
    face2Flip2 = faces[face2Flip1].adjFaces[(vertIndexInFace(face2Flip1, idVert))];
    idV2 = faces[face2Flip1].vertices[(infiniteInFace(face2Flip1)+2)%3];
    idV3 = faces[face2Flip2].vertices[(infiniteInFace(face2Flip2)+2)%3];

    while(orientation2D(idVert, idV2, idV3) > 0.){
        edgeFlip(face2Flip1, face2Flip2);
        tempFace = faces[tempFace].adjFaces[(vertIndexInFace(tempFace, idVert)+1)%3];
        face2Flip1 = faces[tempFace].adjFaces[(vertIndexInFace(tempFace, idVert)+1)%3];
        face2Flip2 = faces[face2Flip1].adjFaces[(vertIndexInFace(face2Flip1, idVert))];
        idV2 = faces[face2Flip1].vertices[(infiniteInFace(face2Flip1)+2)%3];
        idV3 = faces[face2Flip2].vertices[(infiniteInFace(face2Flip2)+2)%3];
    }
}

void Mesh::insertRandPoint2D(int max){
    insertPoint2D(Point(rand()%(2*max) - max, rand()%(2*max) - max, 0));
}

int Mesh::opposedVert(int idFace1, int idFace2){
    for(int i = 0; i < 3; i++){
        if(faces[idFace2].adjFaces[i] == idFace1)
            return i;
    }
    return -1;
}


// cherche parmi les voisins de idFace si leurs arrêtes sont Delaunay et les ajoute à la queue
void Mesh::checkFaceDelaunay(queue<pair<int, int>> & nonDelaunay, int idFace){
    cout << idFace << endl;
    for(int j = 0; j < 3; j++){
        int idAdjacent = faces[idFace].adjFaces[j];
        // ne regarde que des nouvelles arrête pas dans une face infinie
        if(idAdjacent > idFace && isFace2D(idAdjacent)){
            float determinant;
            Point a,b,c,d;
            a = points[faces[idFace].vertices[0]];
            b = points[faces[idFace].vertices[1]];
            c = points[faces[idFace].vertices[2]];
            d = points[faces[idAdjacent].vertices[opposedVert(idFace, idAdjacent)]];

            determinant = a._x - d._x * b._y - d._y * (c._x*c._x - d._x*d._x) + (c._y*c._y - d._y*d._y) +
                          b._x - d._x * c._y - d._y * (a._x*a._x - d._x*d._x) + (a._y*a._y - d._y*d._y) +
                          c._x - d._x * a._y - d._y * (b._x*b._x - d._x*d._x) + (b._y*b._y - d._y*d._y) -
                          c._x - d._x * b._y - d._y * (a._x*a._x - d._x*d._x) + (a._y*a._y - d._y*d._y) -
                          b._x - d._x * a._y - d._y * (c._x*c._x - d._x*d._x) + (c._y*c._y - d._y*d._y) -
                          a._x - d._x * c._y - d._y * (b._x*b._x - d._x*d._x) + (b._y*b._y - d._y*d._y);

            cout << "determinant : " << determinant << endl;
            // si non delaunay localement, on ajoute à la queue
            if(determinant > 0.){
                nonDelaunay.push({idFace, idAdjacent});
            }
        }
    }
}

void Mesh::makeDelaunay(){
    cout << "début delaunay" << endl;
    // triangles adjacents non delaunay / arrêtes non delaunay
    queue<pair<int, int>> nonDelaunay;

    // on cherche les arrêtes non delaunay
    for(int i = 0; i < faces.size(); i++){
        checkFaceDelaunay(nonDelaunay, i);
    }

    // flip des arrêtes qu'on a trouvé
    while(!nonDelaunay.empty()){
        edgeFlip(nonDelaunay.front().first, nonDelaunay.front().second);
        nonDelaunay.pop();
    }

    cout << "fin delaunay" << endl;
}

void Mesh::drawMesh(){
    if(!faces.empty()){
        glColor3d(1, 0, 0);
        glBegin(GL_TRIANGLES);
        glPointDraw(points[faces[currentFace].vertices[0]]);
        glPointDraw(points[faces[currentFace].vertices[1]]);
        glPointDraw(points[faces[currentFace].vertices[2]]);
        glEnd();

        if(highlightNeighbors){
            for(int i = 0; i < 3; i++){
                int idNeighbor = faces[currentFace].adjFaces[i];
                glColor3d(0, 0, 1);
                glBegin(GL_TRIANGLES);

                glPointDraw(points[faces[idNeighbor].vertices[0]]);
                glPointDraw(points[faces[idNeighbor].vertices[1]]);
                glPointDraw(points[faces[idNeighbor].vertices[2]]);
                glEnd();
            }
        }

        glColor3d(1, 1, 1);
        for ( const auto & face : faces){
            glBegin(GL_TRIANGLES);
            glPointDraw(points[face.vertices[0]]);
            glPointDraw(points[face.vertices[1]]);
            glPointDraw(points[face.vertices[2]]);
            glEnd();
        }
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

