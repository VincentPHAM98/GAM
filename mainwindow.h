#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <iostream>
using namespace std;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT

   public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

   private slots:
    void on_pushButton_3_released();

    void on_pushButton_2_released();

    void on_initButton_released();

    void on_splitCenterButton_released();

    void on_selectFaceButton_released();

    void on_highlightNeighborsBox_stateChanged(int arg1);

    void on_edgeFlipButton_released();

    void on_addPointButton_released();

    void on_DelaunayButton_released();

    void on_collapseEdgeButton_released();

    void on_collapseVertice1_valueChanged(int arg1);

    void on_collapseVertice2_valueChanged(int arg1);

    void on_edgeCollapseButton_released();

    void on_voronoiButton_released();

    void on_voronoiBox_stateChanged(int arg1);

    void on_faceSelectBox_valueChanged(int arg1);

    void on_laplacienBox_stateChanged(int arg1);

    void on_computeLaplacianButton_released();

private:
    Ui::MainWindow *ui;
};

#endif  // MAINWINDOW_H
