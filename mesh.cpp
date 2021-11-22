#include "mesh.h"

#include <time.h>

#include <QDebug>
#include <QDir>
#include <QFileDialog>
#include <QVector3D>
#include <QtMath>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>

using namespace std;

// The following functions could be displaced into a module OpenGLDisplayGeometricWorld that would include mesh.h

std::ostream &operator<<(std::ostream &os, const Point &p) {
    os << "Point(" << p._x << "," << p._y << "," << p._z << ") ";
    return os;
}

// Draw a Point
void glPointDraw(const Point &p) {
    glVertex3f(p._x, p._y, p._z);
}

Triangle::Triangle() {}

// returns -1 if not present, internal vertex Id otherwise
int Triangle::getInternalIdx(size_t vertexIdx, int shift) const {
    bool found = false;
    int i;
    for (i = 0; i < 3; ++i)
        if (vertices[i] == vertexIdx) {
            found = true;
            break;
        }
    if (found)
        return (i + shift) % 3;
    else
        return -1;
}

int Triangle::getAdjacentFaceFromGlobalVertex(int vertexIdx) const {
    return adjacent[getInternalIdx(vertexIdx)];
}

uint &Triangle::getVertexFromAdjacentFace(uint faceIdx) {
    int i;
    for (i = 0; i < 3; ++i)
        if (adjacent[i] == faceIdx)
            return adjacent[i];
}

float Mesh::area(const Triangle &t) {
    QVector3D ab = vertices[t.vertices[1]].p - vertices[t.vertices[0]].p;
    QVector3D ac = vertices[t.vertices[2]].p - vertices[t.vertices[0]].p;
    return 0.5 * QVector3D::crossProduct(ab, ac).length();
}

Point Point::cross(Point v) {
    Point p;
    p._x = _y * v._z - _z * v._y;
    p._y = _z * v._x - _x * v._z;
    p._z = _x * v._y - _y * v._x;
    return p;
}

double Point::dot(Point v) {
    return _x * v._x + _y * v._y + _z * v._z;
}

Point Point::normalize() {
    float d = sqrt(_x * _x + _y * _y + _z * _z);
    return Point(_x / d, _y / d, _z / d);
}

double Point::length() {
    return sqrt(_x * _x + _y * _y + _z * _z);
}

Mesh::Mesh() {
    currentFace = 0;
    highlightNeighbors = 0;
    srand(time(NULL));
}

std::pair<int, int> make_ordered_pair(int a, int b) {
    return (a < b ? std::make_pair(a, b) : std::make_pair(b, a));
}

int findOppositeIdx(int a, int b) {
    std::set<int> s = {0, 1, 2};
    s.erase(a);
    s.erase(b);
    return *s.begin();
}

void Mesh::loadOFF(std::string path) {
    vertices.clear();
    triangles.clear();
    std::ifstream off(path);
    if (off.is_open()) {
        std::string line;
        uint nbVertices, nbFaces;
        std::getline(off, line);
        // Getting vertex and face info from file
        std::istringstream ss(line);
        ss >> nbVertices;
        ss >> nbFaces;

        std::vector<std::string> strArr(5);

        // Iterating on each line of the file
        while (std::getline(off, line)) {
            strArr.clear();
            ss.clear();
            ss.str(line);
            std::string value;
            // populate array of strings
            while (ss >> value) {
                strArr.push_back(value);
            }

            if (strArr.size() == 3) {  // it's a vertex description line
                Vertex v;
                v.p._x = std::stof(strArr[0]);
                v.p._y = std::stof(strArr[1]);
                v.p._z = std::stof(strArr[2]);
                v.triangleIdx = -1;
                vertices.push_back(v);
            } else if (strArr.size() == 4) {  // it's a face index line
                Triangle t;
                t.vertices[0] = std::stoi(strArr[1]);
                t.vertices[1] = std::stoi(strArr[2]);
                t.vertices[2] = std::stoi(strArr[3]);
                triangles.push_back(t);
            }
        }
        off.close();

        findTopology();
    } else
        std::cout << "Ouverture du fichier impossible." << std::endl;
}

void Mesh::findTopology() {
    std::map<std::pair<int, int>, std::pair<int, int>> topo;  // EDGE -> FACE

    for (std::size_t i = 0; i < triangles.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            uint vertexIdx = triangles[i].vertices[j];
            if (vertices[vertexIdx].triangleIdx == -1)
                vertices[vertexIdx].triangleIdx = i;  // Affectation de l'id du triangle courant aux sommets n'étant pas attachés à un triangle
            // Créer une paire de sommets avec les id dans l'ordre croissant
            int next = j < 2 ? j + 1 : 0;
            uint nextVertexIdx = triangles[i].vertices[next];
            int oppositeIdx = findOppositeIdx(j, next);
            const auto pair = topo.insert({make_ordered_pair(vertexIdx, nextVertexIdx),
                                           std::make_pair(i, oppositeIdx)});
            if (pair.second == false) {  // Si l'insertion n'a pas pu être faite dans la map
                int otherFaceIdx = pair.first->second.first;
                triangles[i].adjacent[oppositeIdx] = otherFaceIdx;
                triangles[otherFaceIdx].adjacent[pair.first->second.second] = i;
            }
        }
    }
}

void Mesh::calculateLaplacian() {
    for (auto vIt = vertices_begin(); vIt != vertices_past_the_end(); ++vIt) {
        const Point &ui = vIt->p;
        double areas = 0.;
        // Calculating area
        for (auto tCir = ++incident_faces(vIt.getIdx()); tCir != incident_faces(vIt.getIdx()); ++tCir) {
            areas += area(*tCir);
        }
        areas /= 3;

        double laplacianX = 0., laplacianY = 0., laplacianZ = 0.;
        for (auto vCir = ++adjacent_vertices(vIt.getIdx()); vCir != adjacent_vertices(vIt.getIdx()); ++vCir) {
            double cotAij, cotBij;
            const Triangle &cF = vCir.getCurrentFace();
            int u0 = cF.getInternalIdx(vIt.getIdx(), 1);
            int u1 = (u0 + 1) % 3;
            int v0 = u0;
            int v1 = (v0 + 2) % 3;
            QVector3D u = vertices[cF.vertices[u1]].p - vertices[cF.vertices[u0]].p;
            QVector3D v = vertices[cF.vertices[v1]].p - vertices[cF.vertices[v0]].p;
            cotAij = QVector3D::dotProduct(u, v) / QVector3D::crossProduct(u, v).length();

            //            const Triangle &nF = vCir.getNextFace();
            const Triangle &nF = vCir.getCurrentFace();
            u0 = nF.getInternalIdx(vIt.getIdx(), 2);
            u1 = (u0 + 2) % 3;
            v0 = u0;
            v1 = (v0 + 1) % 3;
            u = vertices[nF.vertices[u1]].p - vertices[nF.vertices[u0]].p;
            v = vertices[nF.vertices[v1]].p - vertices[nF.vertices[v0]].p;
            cotBij = QVector3D::dotProduct(u, v) / QVector3D::crossProduct(u, v).length();

            laplacianX += (cotAij + cotBij) * (vCir->p._x - ui._x);
            laplacianY += (cotAij + cotBij) * (vCir->p._y - ui._y);
            laplacianZ += (cotAij + cotBij) * (vCir->p._z - ui._z);
        }
        laplacianX = (1 / (2 * areas)) * laplacianX;
        laplacianY = (1 / (2 * areas)) * laplacianY;
        laplacianZ = (1 / (2 * areas)) * laplacianZ;
        //        std::cout << "laplacianX: " << laplacianX << std::endl;
        QVector3D l(laplacianX, laplacianY, laplacianZ);
        laplacians.push_back(l);
        curvature.push_back(l.length() / 2);
    }
    // normalize curvature
    const auto pair = std::minmax_element(curvature.begin(), curvature.end());
    double min = *pair.first;
    double max = *pair.second;
    for (double &d : curvature) {
        d = (d - min) / (max - min);
    }
}

// retourne l'id du nouveau sommet
uint Mesh::splitTriangle(uint indFace, const Point &newVertexPos, bool delaunay) {
    Triangle *t = &triangles[indFace];  // current working face

    // Mettre a jour les indices de triangles des sommets de la face splitée
    for (int i = 0; i < 3; i++) {
        auto cof = incident_faces(t->vertices[i]);
        auto begin = cof;
        do {
            if (cof.globalIdx() != indFace)
                // if (!isVert2D(t->vertices[i]) || isFace2D(cof.globalIdx()))
                break;
            ++cof;
        } while (true);
        vertices.at(t->vertices[i]).triangleIdx = cof.globalIdx();
    }

    // Adding new vertex to Mesh
    Vertex v(newVertexPos, indFace);
    vertices.push_back(v);
    // id du nouveau vertex
    int idVertex = vertices.size() - 1;

    // Adding new triangles into structure
    uint t1Idx = (int)triangles.size();
    uint t2Idx = t1Idx + 1;
    triangles.push_back(Triangle({t->vertices[1], t->vertices[2], (uint)idVertex},
                                 {t2Idx, indFace, t->adjacent[0]}));
    t = &triangles[indFace];  // current working face
    triangles.push_back(Triangle({t->vertices[2], t->vertices[0], (uint)idVertex},
                                 {indFace, t1Idx, t->adjacent[1]}));
    t = &triangles[indFace];  // current working face

    // Modify adjacent triangles' opposite vertex indices
    Triangle *currentT = &triangles[t->adjacent[0]];
    int currentTIdxToModify = currentT->getInternalIdx(t->vertices[1], 1);
    currentT->adjacent[currentTIdxToModify] = t1Idx;

    currentT = &triangles[t->adjacent[1]];
    currentTIdxToModify = currentT->getInternalIdx(t->vertices[2], 1);
    currentT->adjacent[currentTIdxToModify] = t2Idx;

    // finally modifying split triangle
    t->vertices[2] = idVertex;
    t->adjacent[0] = t1Idx;
    t->adjacent[1] = t2Idx;

    if (delaunay) {
        localDelaunay(indFace);
        localDelaunay(t1Idx);
        localDelaunay(t2Idx);
    }

    return idVertex;
}

uint Mesh::splitTriangleMiddle(int indFace, bool delaunay) {
    const Triangle &t = triangles[indFace];
    // find middle point of triangle
    double pX = 0., pY = 0., pZ = 0.;
    for (int i = 0; i < 3; ++i) {
        pX += vertices[t.vertices[i]].p._x;
        pY += vertices[t.vertices[i]].p._y;
        pZ += vertices[t.vertices[i]].p._z;
    }
    return splitTriangle(indFace, Point(pX / 3., pY / 3., pZ / 3.), delaunay);
}

void Mesh::edgeFlip(int indFace1, int indFace2, bool delaunay) {
    auto faces = std::make_pair(indFace1, indFace2);
    changeIncidentFacesOfFaceVertices(indFace1, faces);
    changeIncidentFacesOfFaceVertices(indFace2, faces);

    int sommetF1 = -1, sommetF2, indSommetF1, indSommetF2;
    // cherche l'arrête commune et stock des infos pour la suite
    for (int i = 0; i < 3; i++) {
        if (triangles[indFace1].adjacent[i] == indFace2) {
            sommetF1 = triangles[indFace1].vertices[i];
            indSommetF1 = i;
        }
        if (triangles[indFace2].adjacent[i] == indFace1) {
            sommetF2 = triangles[indFace2].vertices[i];
            indSommetF2 = i;
        }
    }
    if (sommetF1 == -1) {
        cout << "Impossible de flip" << endl;
        return;
    }
    // update des sommets pour flipper l'arrête
    triangles[indFace1].vertices[(indSommetF1 + 2) % 3] = sommetF2;
    triangles[indFace2].vertices[(indSommetF2 + 2) % 3] = sommetF1;

    // update de l'adjacence des faces flippés
    int temp = triangles[indFace1].adjacent[(indSommetF1 + 1) % 3];
    int temp2 = triangles[indFace2].adjacent[(indSommetF2 + 1) % 3];
    triangles[indFace1].adjacent[indSommetF1] = temp2;
    triangles[indFace2].adjacent[indSommetF2] = temp;

    triangles[indFace1].adjacent[(indSommetF1 + 1) % 3] = indFace2;
    triangles[indFace2].adjacent[(indSommetF2 + 1) % 3] = indFace1;

    // update de l'adjacence des faces voisines flippés
    for (int i = 0; i < 3; i++) {
        if (triangles[temp].adjacent[i] == indFace1) {
            triangles[temp].adjacent[i] = indFace2;
        }
        if (triangles[temp2].adjacent[i] == indFace2)
            triangles[temp2].adjacent[i] = indFace1;
    }

    if (delaunay) {
        localDelaunay(indFace1);
        localDelaunay(indFace2);
    }
}

float Mesh::orientation2D(Point p1, Point p2, Point p3) const {
    p1._z = 0.;
    p2._z = 0.;
    p3._z = 0.;
    QVector3D normale = QVector3D::crossProduct(p2 - p1, p3 - p1);
    return normale.z();
}

float Mesh::orientation2D(int i1, int i2, int i3) const {
    return orientation2D(vertices[i1].p, vertices[i2].p, vertices[i3].p);
}

float Mesh::orientation2D(const Triangle &t) const {
    return orientation2D(
        vertices[t.vertices[0]].p,
        vertices[t.vertices[1]].p,
        vertices[t.vertices[2]].p);
}

bool Mesh::isInside(const Point &p, const Triangle &t) const {
    return (orientation2D(vertices[t.vertices[0]].p, vertices[t.vertices[1]].p, p) > 0. &&
            orientation2D(vertices[t.vertices[1]].p, vertices[t.vertices[2]].p, p) > 0. &&
            orientation2D(vertices[t.vertices[2]].p, vertices[t.vertices[0]].p, p) > 0.);
}

bool Mesh::is2D(int indF) {
    return vertices[triangles[indF].vertices[0]].p._z == 0. &&
           vertices[triangles[indF].vertices[1]].p._z == 0. &&
           vertices[triangles[indF].vertices[2]].p._z == 0.;
}

void Mesh::insertPoint2D(const Point &p, bool delaunay) {
    // on commence sur une face aléatoire pas reliée au point infini
    int indF;
    do {
        indF = rand() % triangles.size();
    } while (triangles[indF].isDeleted || !is2D(indF));

    // marche pour trouver la bonne face
    int n = 0;
    while (!isInside(p, triangles[indF])) {
        if (orientation2D(vertices[triangles[indF].vertices[0]].p, vertices[triangles[indF].vertices[1]].p, p) < 0.)
            indF = triangles[indF].adjacent[2];
        else if (orientation2D(vertices[triangles[indF].vertices[1]].p, vertices[triangles[indF].vertices[2]].p, p) < 0.)
            indF = triangles[indF].adjacent[0];
        else if (orientation2D(vertices[triangles[indF].vertices[2]].p, vertices[triangles[indF].vertices[0]].p, p) < 0.)
            indF = triangles[indF].adjacent[1];
        // si on est hors de l'enveloppe
        if (!isFace2D(indF))
            break;
        if (n++ > triangles.size()) {
            cout << "Impossible d'insérer" << endl;
            return;
        }
    }
    int indV = splitTriangle(indF, p, delaunay);

    // complete l'enveloppe convexe si ajout en dehors
    if (!isFace2D(indF)) {
        indF = triangles[indF].adjacent[infiniteInFace(indF)];
        completeConvexHull(indF, indV, delaunay);
    }
    currentFace = indF;
}

bool Mesh::isVert2D(int indV) {
    return vertices[indV].p._z == 0.;
}

bool Mesh::isFace2D(int indF) {
    return isVert2D(triangles[indF].vertices[0]) &&
           isVert2D(triangles[indF].vertices[1]) &&
           isVert2D(triangles[indF].vertices[2]);
}

int Mesh::findAdjFace(int idFace, int id2find) {
    for (int i = 0; i < 3; i++) {
        if (triangles[idFace].adjacent[i] == id2find)
            return i;
    }
    return -1;
}

int Mesh::vertIndexInFace(int idFace, int idVert) {
    for (int i = 0; i < 3; i++) {
        if (triangles[idFace].vertices[i] == idVert)
            return i;
    }
    return -1;
}

int Mesh::infiniteInFace(int idFace) {
    for (int i = 0; i < 3; i++) {
        if (!isVert2D(triangles[idFace].vertices[i]))
            return i;
    }
    return -1;
}

void Mesh::completeConvexHull(int idFace, int idVert, bool delaunay) {
    // Complete à droite (sens trigo)
    int tempFace = idFace;
    int face2Flip1 = triangles[tempFace].adjacent[(vertIndexInFace(tempFace, idVert) + 2) % 3];
    int face2Flip2 = triangles[face2Flip1].adjacent[(vertIndexInFace(face2Flip1, idVert))];
    int idV2 = triangles[face2Flip1].vertices[(infiniteInFace(face2Flip1) + 1) % 3];
    int idV3 = triangles[face2Flip2].vertices[(infiniteInFace(face2Flip2) + 1) % 3];

    while (orientation2D(idVert, idV2, idV3) < 0.) {
        edgeFlip(face2Flip1, face2Flip2, delaunay);
        tempFace = triangles[tempFace].adjacent[(vertIndexInFace(tempFace, idVert) + 2) % 3];
        face2Flip1 = triangles[tempFace].adjacent[(vertIndexInFace(tempFace, idVert) + 2) % 3];
        face2Flip2 = triangles[face2Flip1].adjacent[(vertIndexInFace(face2Flip1, idVert))];
        idV2 = triangles[face2Flip1].vertices[(infiniteInFace(face2Flip1) + 1) % 3];
        idV3 = triangles[face2Flip2].vertices[(infiniteInFace(face2Flip2) + 1) % 3];
    }

    // Complete à gauche (sens non trigo)
    tempFace = idFace;
    face2Flip1 = triangles[tempFace].adjacent[(vertIndexInFace(tempFace, idVert) + 1) % 3];
    face2Flip2 = triangles[face2Flip1].adjacent[(vertIndexInFace(face2Flip1, idVert))];
    idV2 = triangles[face2Flip1].vertices[(infiniteInFace(face2Flip1) + 2) % 3];
    idV3 = triangles[face2Flip2].vertices[(infiniteInFace(face2Flip2) + 2) % 3];

    while (orientation2D(idVert, idV2, idV3) > 0.) {
        edgeFlip(face2Flip1, face2Flip2, delaunay);
        tempFace = triangles[tempFace].adjacent[(vertIndexInFace(tempFace, idVert) + 1) % 3];
        face2Flip1 = triangles[tempFace].adjacent[(vertIndexInFace(tempFace, idVert) + 1) % 3];
        face2Flip2 = triangles[face2Flip1].adjacent[(vertIndexInFace(face2Flip1, idVert))];
        idV2 = triangles[face2Flip1].vertices[(infiniteInFace(face2Flip1) + 2) % 3];
        idV3 = triangles[face2Flip2].vertices[(infiniteInFace(face2Flip2) + 2) % 3];
    }
}

void Mesh::insertRandPoint2D(int max, bool delaunay) {
    insertPoint2D(Point(rand() % (2 * max) - max, rand() % (2 * max) - max, 0), delaunay);
}

int Mesh::opposedVert(int idFace1, int idFace2) {
    for (int i = 0; i < 3; i++) {
        if (triangles[idFace2].adjacent[i] == idFace1)
            return i;
    }
    return -1;
}

Point Mesh::getNormal(Point a, Point b, Point c) {
    auto u = b - a;
    auto v = c - a;
    u.normalize();
    v.normalize();

    auto normale = QVector3D::crossProduct(u, v);
    return Point(normale.x(), normale.y(), normale.z());
}

Point Mesh::getNormal(Triangle f) {
    Point a, b, c;
    a = vertices[f.vertices[0]].p;
    b = vertices[f.vertices[1]].p;
    c = vertices[f.vertices[2]].p;

    return getNormal(a, b, c);
}

Point Mesh::getNormal(int idFace) {
    return getNormal(triangles[idFace]);
}

bool Mesh::isInCirconscrit(int idF, Point p) {
    Point _p[3];  // points de la face idF projetés dans le paraboloïde
    for (int i = 0; i < 3; i++) {
        _p[i]._x = vertices[triangles[idF].vertices[i]].p._x;
        _p[i]._y = vertices[triangles[idF].vertices[i]].p._y;
        _p[i]._z = _p[i]._x * _p[i]._x + _p[i]._y * _p[i]._y;
    }

    Point pTemp = Point(p._x, p._y, p._x * p._x + p._y * p._y);  // point p projeté dans le paraboloïde

    Point normale = getNormal(_p[0], _p[1], _p[2]);
    Point pa = Point(_p[0] - pTemp);

    return normale.dot(pa) < 0;
}

// cherche parmi les voisins de idFace si leurs arrêtes sont Delaunay et les ajoute à la queue
void Mesh::checkFaceDelaunayGlobal(queue<pair<int, int>> &nonDelaunay, int idFace) {
    for (int j = 0; j < 3; j++) {
        int idAdjacent = triangles[idFace].adjacent[j];
        // ne regarde que des nouvelles arrête pas dans une face infinie
        if (idAdjacent > idFace && isFace2D(idFace) && isFace2D(idAdjacent)) {
            Vertex d = vertices[triangles[idAdjacent].vertices[opposedVert(idFace, idAdjacent)]];
            // si non delaunay localement, on ajoute à la queue
            if (isInCirconscrit(idFace, d.p)) {
                nonDelaunay.push({idFace, idAdjacent});
            }
        }
    }
}

void Mesh::makeDelaunay() {
    // triangles adjacents non delaunay / arrêtes non delaunay
    queue<pair<int, int>> nonDelaunay;

    // on cherche les arrêtes non delaunay
    // for (int i = 0; i < triangles.size(); i++) {
    for (auto it = faces_begin(); it != faces_past_the_end(); ++it) {
        checkFaceDelaunayGlobal(nonDelaunay, it.getIdx());
    }
    // tant que la triangulation n'est pas de Delaunay
    while (!nonDelaunay.empty()) {
        // flip des arrêtes qu'on a trouvé
        while (!nonDelaunay.empty()) {
            edgeFlip(nonDelaunay.front().first, nonDelaunay.front().second, 0);
            nonDelaunay.pop();
        }
        // on re verifie s'il en reste encore
        for (auto it = faces_begin(); it != faces_past_the_end(); ++it) {
            checkFaceDelaunayGlobal(nonDelaunay, it.getIdx());
        }
    }
    cout << "global delaunay" << endl;
}

void Mesh::checkFaceDelaunayLocal(queue<pair<int, int>> &nonDelaunay, int idFace) {
    for (int j = 0; j < 3; j++) {
        int idAdjacent = triangles[idFace].adjacent[j];
        // ne regarde que des nouvelles arrête pas dans une face infinie
        if (isFace2D(idFace) && isFace2D(idAdjacent)) {
            Vertex d = vertices[triangles[idAdjacent].vertices[opposedVert(idFace, idAdjacent)]];
            // si non delaunay localement, on ajoute à la queue
            if (isInCirconscrit(idFace, d.p)) {
                nonDelaunay.push({idFace, idAdjacent});
            }
        }
    }
}

void Mesh::localDelaunay(int idF) {
    // triangles adjacents non delaunay / arrêtes non delaunay
    queue<pair<int, int>> nonDelaunay;

    checkFaceDelaunayLocal(nonDelaunay, idF);

    while (!nonDelaunay.empty()) {
        edgeFlip(nonDelaunay.front().first, nonDelaunay.front().second, 0);
        localDelaunay(nonDelaunay.front().first);
        localDelaunay(nonDelaunay.front().second);
        nonDelaunay.pop();
    }
}

std::pair<int, int> Mesh::findFacesWithCommonEdge(uint idVert1, uint idVert2) {
    std::vector<int> faces;
    auto cof = incident_faces(idVert1);
    auto begin = cof;
    do {
        if (cof->getInternalIdx(idVert2) != -1)
            faces.push_back(cof.globalIdx());
        ++cof;
    } while (cof != begin);

    return std::make_pair<int, int>((int)faces[0], (int)faces[1]);
}

std::set<int> Mesh::adjacentVerticesOfVertex(uint indV) {
    std::set<int> vertices;
    auto cov = adjacent_vertices(indV);
    auto begin = cov;
    do {
        vertices.insert(cov.globalVertexIdx());
        ++cov;
    } while (cov != begin);
    return vertices;
}

// Si les sommets d'un triangle ont pour indice de triangle lui même,
// change les sommets pour avoir pour indice les triangles voisins.
void Mesh::changeIncidentFacesOfFaceVertices(uint idFace, std::pair<int, int> deletedFaces) {
    Triangle &t = triangles.at(idFace);
    for (int i = 0; i < 3; i++) {
        Vertex &v = vertices.at(t.vertices[i]);
        // si le sommet pointe vers les triangles allant disparaitre
        auto cof = incident_faces(t.vertices[i]);
        while (v.triangleIdx == deletedFaces.first || v.triangleIdx == deletedFaces.second) {
            ++cof;
            v.triangleIdx = cof.globalIdx();
        }
    }
}

// idVert1 a garder et déplacer et idVert2 à supprimer
int Mesh::collapseEdge(uint idVert1, uint idVert2) {
    // Verifier precond: pas de sommet en commun sur les 2 arêtes a collapse
    auto adjacentVerticesOfVert1 = adjacentVerticesOfVertex(idVert1);
    auto adjacentVerticesOfVert2 = adjacentVerticesOfVertex(idVert2);
    std::vector<int> intersection;
    std::set_intersection(
        adjacentVerticesOfVert1.begin(), adjacentVerticesOfVert1.end(),
        adjacentVerticesOfVert2.begin(), adjacentVerticesOfVert2.end(),
        std::back_inserter(intersection));
    if (intersection.size() != 2) {
        // modal qt pour prevenir
        std::cout << "Err: les arêtes à collapse ont une arête en commun." << std::endl;
        return 1;
    }

    // Deplacer le premier sommet au milieu des 2
    Point middle = (vertices[idVert1].p + vertices[idVert2].p) / 2;
    // std::cout << "v1: " << vertices[idVert1].p << "v2: " << vertices[idVert2].p << "middle: " << middle << std::endl;
    vertices[idVert1].p = middle;
    // std::cout << "after: " << vertices[idVert1].p << std::endl;

    // Trouver les 2 faces qui vont etre supprimées
    auto faces = findFacesWithCommonEdge(idVert1, idVert2);

    // changer la face incidente des sommets des faces allant etre supprimees
    changeIncidentFacesOfFaceVertices(faces.first, faces);
    changeIncidentFacesOfFaceVertices(faces.second, faces);
    // Lister les faces qui vont perdre leur sommet
    // pour arranger leur sommet plus tard
    auto cof = incident_faces(idVert2);
    auto begin = cof;
    std::vector<int> patchUp;
    do {
        if (cof.globalIdx() != faces.first && cof.globalIdx() != faces.second) {
            patchUp.push_back(cof.globalIdx());
        }
        ++cof;
    } while (cof != begin);
    // Mettre en ordre l'adjacence des faces adjacentes aux faces a supprimer
    int adjTriangle1Idx, adjTriangle2Idx;
    //      Adjacence face 1:
    adjTriangle1Idx = triangles[faces.first].getAdjacentFaceFromGlobalVertex(idVert1);
    adjTriangle2Idx = triangles[faces.first].getAdjacentFaceFromGlobalVertex(idVert2);
    triangles[adjTriangle1Idx].getVertexFromAdjacentFace(faces.first) = adjTriangle2Idx;
    triangles[adjTriangle2Idx].getVertexFromAdjacentFace(faces.first) = adjTriangle1Idx;
    //      Adjacence face 2:
    adjTriangle1Idx = triangles[faces.second].getAdjacentFaceFromGlobalVertex(idVert1);
    adjTriangle2Idx = triangles[faces.second].getAdjacentFaceFromGlobalVertex(idVert2);
    triangles[adjTriangle1Idx].getVertexFromAdjacentFace(faces.second) = adjTriangle2Idx;
    triangles[adjTriangle2Idx].getVertexFromAdjacentFace(faces.second) = adjTriangle1Idx;

    // arranger les sommets des faces de tout a l'heure
    for (auto &&idFace : patchUp) {
        triangles[idFace].vertices[triangles[idFace].getInternalIdx(idVert2)] = idVert1;
    }

    // Supprimer l'arrete et les 2 faces
    triangles[faces.first].remove();
    triangles[faces.second].remove();
    vertices[idVert2].remove();
    return 0;
}

void Mesh::collapseShortestEdge() {
    // on cherche les arêtes les plus courtes
    using Edge = std::pair<int, int>;
    std::priority_queue<std::pair<double, Edge>,
                        std::vector<std::pair<double, Edge>>,
                        std::greater<std::pair<double, Edge>>>
        queue;

    for (auto it = vertices_begin(); it != vertices_past_the_end(); ++it) {
        auto cov = adjacent_vertices(it.getIdx());
        auto begin = cov;
        do {
            Edge edge = std::make_pair((int)it.getIdx(), cov.globalVertexIdx());
            auto edgeVector = cov->p - it->p;
            double length = edgeVector.lengthSquared();
            queue.push(std::make_pair(length, edge));
            ++cov;
        } while (cov != begin);
    }
    int success;
    do {
        success = collapseEdge(queue.top().second.first, queue.top().second.second);
        queue.pop();
    } while (success && queue.size());
}

double Point::tangente(Point a, Point b, Point c) {
    Point ba = a - b;
    Point bc = c - b;
    Point crossABC = bc.cross(ba);
    double dotABC = bc.dot(ba);

    // pour éviter d'avoir des nan avec un triangle rectangle
    if (dotABC != 0) {
        if (crossABC._z < 0)
            return -crossABC.length() / bc.dot(ba);
        return crossABC.length() / bc.dot(ba);
    }
    return 42;
}

Point Mesh::centreCercleCirconscrit(int idF) {
    // Methode des tangentes vu en cours
    Point _p[3];
    double tan[3];
    for (int i = 0; i < 3; i++) {
        _p[i] = vertices[triangles[idF].vertices[i]].p;
    }
    for (int i = 0; i < 3; i++) {
        tan[i] = _p[i].tangente(vertices[triangles[idF].vertices[(i + 2) % 3]].p,
                                vertices[triangles[idF].vertices[(i) % 3]].p,
                                vertices[triangles[idF].vertices[(i + 1) % 3]].p);
        // si triangle rectangle, centre du cercle au milieu de l'hypotenuse
        if (tan[i] == 42) {
            Point p = vertices[triangles[idF].vertices[(i + 2) % 3]].p * 0.5 + vertices[triangles[idF].vertices[(i + 1) % 3]].p * 0.5;
            return p;
        }
    }
    // ponderation avec les coordonnees barycentriques
    double alpha = (tan[1] + tan[2]);
    double beta = (tan[0] + tan[2]);
    double gamma = (tan[1] + tan[0]);
    double somme = alpha + beta + gamma;
    alpha /= somme;
    beta /= somme;
    gamma /= somme;
    Point p = vertices[triangles[idF].vertices[0]].p * alpha + vertices[triangles[idF].vertices[1]].p * beta + vertices[triangles[idF].vertices[2]].p * gamma;
    return p;
}

void Mesh::computeVoronoi() {
    cout << "computing Voronoi" << endl;
    voronoiCenter.clear();

    for (auto it = faces_begin(); it != faces_past_the_end(); ++it) {
        if (isFace2D(it.getIdx())) {
            voronoiCenter.insert({it.getIdx(), centreCercleCirconscrit(it.getIdx())});
        }
        // else
        //     voronoiCenter.push_back(Point(0, 0, 0));
    }
}

void Mesh::drawMesh() {
    if (!triangles.empty()) {
        glColor3d(1, 0, 0);
        glBegin(GL_TRIANGLES);
        glPointDraw(vertices[triangles[currentFace].vertices[0]].p);
        glPointDraw(vertices[triangles[currentFace].vertices[1]].p);
        glPointDraw(vertices[triangles[currentFace].vertices[2]].p);
        glEnd();

        if (highlightNeighbors) {
            for (int i = 0; i < 3; i++) {
                int idNeighbor = triangles[currentFace].adjacent[i];
                glColor3d(0, 0, 1);
                glBegin(GL_TRIANGLES);

                glPointDraw(vertices[triangles[idNeighbor].vertices[0]].p);
                glPointDraw(vertices[triangles[idNeighbor].vertices[1]].p);
                glPointDraw(vertices[triangles[idNeighbor].vertices[2]].p);
                glEnd();
            }
        }

        glColor3d(1, 1, 1);
        // for (int i = 0; i < triangles.size(); i++)  {
        for (auto it = faces_begin(); it != faces_past_the_end(); ++it) {
            glBegin(GL_TRIANGLES);
            glPointDraw(vertices[it->vertices[0]].p);
            glPointDraw(vertices[it->vertices[1]].p);
            glPointDraw(vertices[it->vertices[2]].p);
            glEnd();
        }
    }
}

void Mesh::drawMeshWireFrame() {
    // color selected vertices (for edge collapse)
    if (!triangles.empty()) {
        glColor3d(0, 1, 1);
        glPointSize(7.F);
        glBegin(GL_POINTS);
        glPointDraw(vertices[selectedVertex1].p);
        glColor3d(1, 0, 1);
        glPointDraw(vertices[selectedVertex2].p);
        glEnd();

        glColor3d(1, 1, 1);
        for (auto it = faces_begin(); it != faces_past_the_end(); ++it) {
            glBegin(GL_LINE_STRIP);
            glPointDraw(vertices[it->vertices[0]].p);
            glPointDraw(vertices[it->vertices[1]].p);
            glEnd();

            glBegin(GL_LINE_STRIP);
            glPointDraw(vertices[it->vertices[0]].p);
            glPointDraw(vertices[it->vertices[2]].p);
            glEnd();

            glBegin(GL_LINE_STRIP);
            glPointDraw(vertices[it->vertices[1]].p);
            glPointDraw(vertices[it->vertices[2]].p);
            glEnd();
        }
    }
}

void Mesh::drawMeshLaplacian() {
    // for (const auto &face : triangles) {
    for (auto it = faces_begin(); it != faces_past_the_end(); ++it) {
        // glBegin(GL_LINE_STRIP);
        glBegin(GL_TRIANGLES);
        double v = curvature[it->vertices[0]];
        QColor c(0, 0, 0);
        c.setHsv((std::log(1 + 100 * v)), 255, 255);
        glColor3d(c.red(), c.green(), c.blue());
        glPointDraw(vertices[it->vertices[0]].p);
        v = curvature[it->vertices[1]];
        glColor3d(c.red(), c.green(), c.blue());
        glPointDraw(vertices[it->vertices[1]].p);
        v = curvature[it->vertices[2]];
        glColor3d(c.red(), c.green(), c.blue());
        glPointDraw(vertices[it->vertices[2]].p);
        glEnd();
    }
}

void Mesh::drawVoronoi() {
    glColor3d(0, 1, 0);
    if (!voronoiCenter.empty()) {
        glBegin(GL_LINES);
        for (auto it = faces_begin(); it != faces_past_the_end(); ++it) {
            int i = it.getIdx();
            for (int j = 0; j < 3; j++) {
                if (isFace2D(i) && isFace2D(triangles[i].adjacent[j])) {
                    auto pointIt = voronoiCenter.find(i);
                    auto adjacentIt = voronoiCenter.find(triangles[i].adjacent[j]);
                    if (pointIt != voronoiCenter.end() && adjacentIt != voronoiCenter.end()) {
                        glPointDraw(pointIt->second);
                        glPointDraw(adjacentIt->second);
                    }
                }
            }
        }
        glEnd();
    }
}

GeometricWorld::GeometricWorld() {
    double width = 0.5, depth = 0.6, height = 0.8;
    _bBox.push_back(Point(-0.5 * width, -0.5 * depth, -0.5 * height));  //0
    _bBox.push_back(Point(-0.5 * width, 0.5 * depth, -0.5 * height));   // 1
    _bBox.push_back(Point(0.5 * width, -0.5 * depth, -0.5 * height));   // 2
    _bBox.push_back(Point(-0.5 * width, -0.5 * depth, 0.5 * height));   // 3
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

void GeometricWorld::drawVoronoi() {
    _mesh.drawVoronoi();
}

void GeometricWorld::drawLaplacien() {
    _mesh.drawMeshLaplacian();
}

QVector3D Point::operator-(const Point &p) {
    return QVector3D(p._x - _x, p._y - _y, p._z - _z);
}
