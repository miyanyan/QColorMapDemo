#include "HotPlot.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    HotPlot w;
    w.resize(500, 400);
    w.show();
    return a.exec();
}
