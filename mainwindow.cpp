#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_3_released()
{
    ui->widget->toggleWireFrame = ui->widget->toggleWireFrame ^= true;
}

void MainWindow::on_pushButton_2_released()
{
    ui->widget->_geomWorld._mesh.insertRandPoint2D(5);
}

void MainWindow::on_initButton_released()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open off file"), "", tr("Off files (*.off)"));
    ui->widget->_geomWorld._mesh.initFile(fileName.toStdString());
    ui->nbTriangleLabel->setText(QString::number(ui->widget->_geomWorld._mesh.faces.size()));

//    ui->widget->_geomWorld._mesh.splitTriangleAtCenter(0);
}

void MainWindow::on_splitCenterButton_released()
{
    ui->widget->_geomWorld._mesh.splitTriangleAtCenter(ui->widget->_geomWorld._mesh.currentFace);
}

void MainWindow::on_selectFaceButton_released()
{
    int buttonValue = ui->faceSelectBox->value();
    if(buttonValue >= 0 && buttonValue < ui->widget->_geomWorld._mesh.faces.size())
        ui->widget->_geomWorld._mesh.currentFace = buttonValue;
}

void MainWindow::on_highlightNeighborsBox_stateChanged(int arg1)
{
    ui->widget->_geomWorld._mesh.highlightNeighbors = !ui->widget->_geomWorld._mesh.highlightNeighbors;
}

void MainWindow::on_edgeFlipButton_released()
{
    int face1 = ui->flipFace1->value();
    int face2 = ui->flipFace2->value();
    cout << "flipping : " << face1 << " " << face2 << endl;
    ui->widget->_geomWorld._mesh.edgeFlip(face1, face2);
}
