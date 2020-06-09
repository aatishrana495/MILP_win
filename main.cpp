#include "milp.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Milp *w;
    w = new Milp();
    w->show();
    return a.exec();
}
