#ifndef MILP_H
#define MILP_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <complex>
#include <QMainWindow>
#include <QApplication>

#include <ui_milp.h>

#define MAX_N 1001
#define MAX_M 1001

typedef long long lld;
typedef unsigned long long llu;
using namespace std;


QT_BEGIN_NAMESPACE
namespace Ui {
class Milp;
}
QT_END_NAMESPACE

class Milp : public QMainWindow
{
  Q_OBJECT

public:
  explicit Milp(QWidget *parent = 0);
  ~Milp();

private slots:
  void update_all();
  void reset_all();
  void solve();


private:
  Ui::Milp *ui;
  int n,m;
  double A[MAX_M][MAX_N], b[MAX_M], c[MAX_N], v;
  int N[MAX_N], B[MAX_M];
  pair<vector<double>, double> ret;

  void pivot(int x,  int y);
  int iterate_simplex();
  int initialise_simplex();
  pair<vector<double>, double> simplex();


  // ui related
  void read_equations();
  void reset_inputs();
  void reset_line_edits();
  void update_table_inputs();
  void print_constraint_equations();

  void update_optimize_table();
  void read_optimize_equation();
  void print_optimize_equation();

  void display_result_coefficients();
  void print_result_coefficients();
  void update_result_table();

  void update_inputs();

};


#endif // MILP_H
