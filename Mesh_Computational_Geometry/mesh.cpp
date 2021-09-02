#include "mesh.h"

// The following functions could be displaced into a module OpenGLDisplayGeometricWorld that would include mesh.h

// Draw a Point
void glPointDraw(const Point & p) {
    glVertex3f(p._x, p._y, p._z);
}

Mesh::Mesh() {
    // Tetrahèdre
    // 	Vertices
    //		Base
    vertices.push_back(Point(0.,0.,0.));
    vertices.push_back(Point(0.,1.,0.));
    vertices.push_back(Point(1.,1.,0.));
    //		Sommet
    vertices.push_back(Point(0.,0.,1.));
    // 	Faces
    faces.push_back({0,1,2}); // Base
    faces.push_back({0,1,3});
    faces.push_back({1,2,3});
    faces.push_back({2,0,3});
}

void Mesh::drawMesh() {
    int i = 0;
    for (const auto & face : faces) {
        if (i % 3 == 0)
            glColor3d(1,0,0);
        else if (i % 3 == 1)
            glColor3d(0,1,0);
        else
            glColor3d(0,0,1);

        glBegin(GL_TRIANGLES);
        glPointDraw(vertices[face[0]]);
        glPointDraw(vertices[face[1]]);
        glPointDraw(vertices[face[2]]);
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
//    glColor3d(1,0,0);
//    glBegin(GL_TRIANGLES);
//    glPointDraw(_bBox[0]);
//    glPointDraw(_bBox[1]);
//    glPointDraw(_bBox[2]);
//    glEnd();

//    glColor3d(0,1,0);
//    glBegin(GL_TRIANGLES);
//    glPointDraw(_bBox[0]);
//    glPointDraw(_bBox[2]);
//    glPointDraw(_bBox[3]);
//    glEnd();

//    glColor3d(0,0,1);
//    glBegin(GL_TRIANGLES);
//    glPointDraw(_bBox[0]);
//    glPointDraw(_bBox[3]);
//    glPointDraw(_bBox[1]);
//    glEnd();

    _mesh.drawMesh();
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

