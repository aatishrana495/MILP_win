#include "milp.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Milp w;
    w.show();
    return a.exec();
}
