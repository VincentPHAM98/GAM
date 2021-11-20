#include "mainwindow.h"

#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent),
                                          ui(new Ui::MainWindow) {
    ui->setupUi(this);
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::on_pushButton_3_released() {
    ui->widget->toggleWireFrame = ui->widget->toggleWireFrame ^= true;
}

void MainWindow::on_pushButton_2_released() {
    ui->widget->_geomWorld._mesh.insertRandPoint2D(10, ui->incrementalDelaunayBox->isChecked());
}

void MainWindow::on_initButton_released() {
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open off file"), "", tr("Off files (*.off)"));
    ui->widget->_geomWorld._mesh.loadOFF(fileName.toStdString());
}

void MainWindow::on_splitCenterButton_released() {
    ui->widget->_geomWorld._mesh.splitTriangleMiddle(ui->widget->_geomWorld._mesh.currentFace, ui->incrementalDelaunayBox->isChecked());
}

void MainWindow::on_selectFaceButton_released() {
    int buttonValue = ui->faceSelectBox->value();
    if (buttonValue >= 0 && buttonValue < ui->widget->_geomWorld._mesh.triangles.size())
        ui->widget->_geomWorld._mesh.currentFace = buttonValue;
}

void MainWindow::on_highlightNeighborsBox_stateChanged(int arg1) {
    ui->widget->_geomWorld._mesh.highlightNeighbors = !ui->widget->_geomWorld._mesh.highlightNeighbors;
}

void MainWindow::on_edgeFlipButton_released() {
    int face1 = ui->flipFace1->value();
    int face2 = ui->flipFace2->value();
    ui->widget->_geomWorld._mesh.edgeFlip(face1, face2, 0);
}

void MainWindow::on_addPointButton_released() {
    ui->widget->_geomWorld._mesh.insertPoint2D(Point(ui->xBox->text().toFloat(), ui->yBox->text().toFloat(), 0)
                                               , ui->incrementalDelaunayBox->isChecked());
}

void MainWindow::on_DelaunayButton_released() {
    ui->widget->_geomWorld._mesh.makeDelaunay();
}

void MainWindow::on_collapseEdgeButton_released() {
    uint vertex1 = ui->collapseVertice1->value();
    uint vertex2 = ui->collapseVertice2->value();
    ui->widget->_geomWorld._mesh.collapseEdge((uint)vertex1, (uint)vertex2);
}

void MainWindow::on_collapseVertice1_valueChanged(int arg1) {
    if (arg1 == ui->collapseVertice2->value()) {
        ui->collapseVertice2->setValue(arg1 + 1);
    }
    ui->widget->_geomWorld._mesh.selectedVertex1 = arg1;
}

void MainWindow::on_collapseVertice2_valueChanged(int arg1) {
    if (arg1 == ui->collapseVertice1->value()) {
        ui->collapseVertice2->setValue(arg1 - 1);
    }
    ui->widget->_geomWorld._mesh.selectedVertex2 = arg1;
}

void MainWindow::on_edgeCollapseButton_released()
{
    ui->widget->_geomWorld._mesh.collapseShortestEdge();
}
