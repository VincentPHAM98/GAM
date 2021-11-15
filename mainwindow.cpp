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

void MainWindow::on_fileSelectBtn_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this,
        tr("Open .off"), QDir::homePath() + "/Downloads/", tr(".off Files (*.off)"), nullptr, QFileDialog::DontUseNativeDialog);
}
