#include "milp.h"
#include "ui_milp.h"
#include <QApplication>
/*
 The Simplex algorithm aims to solve a linear program - optimising a linear function subject
 to linear constraints. As such it is useful for a very wide range of applications.

 N.B. The linear program has to be given in *slack form*, which is as follows:
 maximise
     c_1 * x_1 + c_2 * x_2 + ... + c_n * x_n + v
 subj. to
     a_11 * x_1 + a_12 * x_2 + ... + a_1n * x_n + b_1 = s_1
     a_21 * x_1 + a_22 * x_2 + ... + a_2n * x_n + b_2 = s_2
     ...
     a_m1 * x_1 + a_m2 * x_2 + ... + a_mn * x_n + b_m = s_m
 and
     x_1, x_2, ..., x_n, s_1, s_2, ..., s_m >= 0

 Every linear program can be translated into slack form; the parameters to specify are:
     - the number of variables, n, and the number of constraints, m;
     - the matrix A = [[A_11, A_12, ..., A_1n], ..., [A_m1, A_m2, ..., A_mn]];
     - the vector b = [b_1, b_2, ..., b_m];
     - the vector c = [c_1, c_2, ..., c_n] and the constant v.

 Complexity:    O(m^(n/2)) worst case
                O(n + m) average case (common)
*/

Milp::Milp(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::Milp){
      ui->setupUi(this);
      // initialise variables
      n = m = 0;

      // update
      reset_all();
      update_all();

      // connect slots
      connect(ui->update_inputs_btn, SIGNAL(pressed()),this,SLOT(update_all()));
      connect(ui->reset_inputs_btn, SIGNAL(pressed()),this,SLOT(reset_all()));
      //connect(ui->read_equations_btn, SIGNAL(pressed()),this,SLOT(read_equations()));
      connect(ui->solve_btn, SIGNAL(pressed()),this,SLOT(solve()));
  }

Milp::~Milp(){}

/*
 pivot yth variable around xth constraint
*/

void Milp::pivot(int x, int y){
  // printf("Pivoting variable %d around constraint %d.\n", y, x);

  // first rearrange the x-th row
  for (int j=0;j<n;j++)
  {
      if (j != y)
      {
          A[x][j] /= -A[x][y];
      }
  }
  b[x] /= -A[x][y];
  A[x][y] = 1.0 / A[x][y];

  // now rearrange the other rows
  for (int i=0;i<m;i++)
  {
      if (i != x)
      {
          for (int j=0;j<n;j++)
          {
              if (j != y)
              {
                  A[i][j] += A[i][y] * A[x][j];
              }
          }
          b[i] += A[i][y] * b[x];
          A[i][y] *= A[x][y];
      }
  }

  // now rearrange the objective function
  for (int j=0;j<n;j++)
  {
      if (j != y)
      {
          c[j] += c[y] * A[x][j];
      }
  }
  v += c[y] * b[x];
  c[y] *= A[x][y];

  // finally, swap the basic & nonbasic variable
  swap(B[x], N[y]);
}

/*
 Run a single iteration of the simplex algorithm.
 Returns: 0 if OK, 1 if STOP, -1 if UNBOUNDED
*/

int Milp::iterate_simplex(){
  // printf("--------------------\n");
  // printf("State:\n");
  // printf("Maximise: ");
  // for (int j=0;j<n;j++) printf("%lfx_%d + ", c[j], N[j]);
  // printf("%lf\n", v);
  // printf("Subject to:\n");
  // for (int i=0;i<m;i++)
  // {
  //     for (int j=0;j<n;j++) printf("%lfx_%d + ", A[i][j], N[j]);
  //     printf("%lf = x_%d\n", b[i], B[i]);
  // }

  // getchar(); // uncomment this for debugging purposes!

  int ind = -1, best_var = -1;
  for (int j=0;j<n;j++)
  {
      if (c[j] > 0)
      {
          if (best_var == -1 || N[j] < ind)
          {
              ind = N[j];
              best_var = j;
          }
      }
  }
  if (ind == -1) return 1;

  double max_constr = INFINITY;
  int best_constr = -1;
  for (int i=0;i<m;i++)
  {
      if (A[i][best_var] < 0)
      {
          double curr_constr = -b[i] / A[i][best_var];
          if (curr_constr < max_constr)
          {
              max_constr = curr_constr;
              best_constr = i;
          }
      }
  }
  if (isinf(max_constr)) return -1;
  else pivot(best_constr, best_var);

  return 0;
}


/*
(Possibly) converts the LP into a slack form with a feasible basic solution.
Returns 0 if OK, -1 if INFEASIBLE
*/


int Milp::initialise_simplex(){
  int k = -1;
  double min_b = -1;
  for (int i=0;i<m;i++)
  {
      if (k == -1 || b[i] < min_b)
      {
          k = i;
          min_b = b[i];
      }
  }

  if (b[k] >= 0) // basic solution feasible!
  {
      for (int j=0;j<n;j++) N[j] = j;
      for (int i=0;i<m;i++) B[i] = n + i;
      return 0;
  }

  // generate auxiliary LP
  n++;
  for (int j=0;j<n;j++) N[j] = j;
  for (int i=0;i<m;i++) B[i] = n + i;

  // store the objective function
  double c_old[MAX_N];
  for (int j=0;j<n-1;j++) c_old[j] = c[j];
  double v_old = v;

  // aux. objective function
  c[n-1] = -1;
  for (int j=0;j<n-1;j++) c[j] = 0;
  v = 0;
  // aux. coefficients
  for (int i=0;i<m;i++) A[i][n-1] = 1;

  // perform initial pivot
  pivot(k, n - 1);

  // now solve aux. LP
  int code;
  while (!(code = iterate_simplex()));

  assert(code == 1); // aux. LP cannot be unbounded!!!

  if (v != 0) return -1; // infeasible!

  int z_basic = -1;
  for (int i=0;i<m;i++)
  {
      if (B[i] == n - 1)
      {
          z_basic = i;
          break;
      }
  }

  // if x_n basic, perform one degenerate pivot to make it nonbasic
  if (z_basic != -1) pivot(z_basic, n - 1);

  int z_nonbasic = -1;
  for (int j=0;j<n;j++)
  {
      if (N[j] == n - 1)
      {
          z_nonbasic = j;
          break;
      }
  }
  assert(z_nonbasic != -1);

  for (int i=0;i<m;i++)
  {
      A[i][z_nonbasic] = A[i][n-1];
  }
  swap(N[z_nonbasic], N[n - 1]);

  n--;
  for (int j=0;j<n;j++) if (N[j] > n) N[j]--;
  for (int i=0;i<m;i++) if (B[i] > n) B[i]--;

  for (int j=0;j<n;j++) c[j] = 0;
  v = v_old;

  for (int j=0;j<n;j++)
  {
      bool ok = false;
      for (int jj=0;jj<n;jj++)
      {
          if (j == N[jj])
          {
              c[jj] += c_old[j];
              ok = true;
              break;
          }
      }
      if (ok) continue;
      for (int i=0;i<m;i++)
      {
          if (j == B[i])
          {
              for (int jj=0;jj<n;jj++)
              {
                  c[jj] += c_old[j] * A[i][jj];
              }
              v += c_old[j] * b[i];
              break;
          }
      }
  }

  return 0;
}


/*
 Runs the simplex algorithm to optimise the LP.
 Returns a vector of -1s if unbounded, -2s if infeasible.
*/

pair<vector<double>, double> Milp::simplex(){
  if (initialise_simplex() == -1)
  {
      return make_pair(vector<double>(n + m, -2), INFINITY);
  }

  int code;
  while (!(code = iterate_simplex()));

  if (code == -1) return make_pair(vector<double>(n + m, -1), INFINITY);

  vector<double> ret;
  ret.resize(n + m);
  for (int j=0;j<n;j++)
  {
      ret[N[j]] = 0;
  }
  for (int i=0;i<m;i++)
  {
      ret[B[i]] = b[i];
  }

  return make_pair(ret, v);
}


void Milp::update_inputs(){
  n = atoi(ui->variable_input->text().toStdString().c_str());
  m = atoi(ui->constraints_input->text().toStdString().c_str());
}

void Milp::reset_inputs(){
  n = m = 0;
}

void Milp::reset_line_edits(){
  ui->variable_input->setText("");
  ui->constraints_input->setText("");
  ui->variable_input->setPlaceholderText("Enter the no of variables");
  ui->constraints_input->setPlaceholderText("Enter the no of contraints");
}

void Milp::update_table_inputs(){
  if(n != 0 && m != 0){
    ui->table_input->setShowGrid(true);
    ui->table_input->verticalHeader()->setVisible(false);
    ui->table_input->setRowCount(m+1);
    ui->table_input->setColumnCount(n+2);
    for(int i=0;i<n;i++){
      string str = "var" + to_string(i+1);
      QString s = QString::fromUtf8(str.c_str());
      ui->table_input->setItem(0, i, new QTableWidgetItem(s));
    }
    ui->table_input->setItem(0, n, new QTableWidgetItem("constant"));
    ui->table_input->setItem(0, n+1, new QTableWidgetItem("equate"));
    for(int i=1;i<m+1;i++){
      for(int j=0;j<n+2;j++){
          ui->table_input->setItem(i, j, new QTableWidgetItem(""));
      }
    }
  }else{
    ui->table_input->setShowGrid(true);
    ui->table_input->verticalHeader()->setVisible(false);
    ui->table_input->setRowCount(0);
    ui->table_input->setColumnCount(0);
  }
}

void Milp::read_equations(){
  for(int i=1;i<m+1;i++){
    for(int j=0;j<n+2;j++){
        if(j==n+1){
          // assuming that they will be greater than 0
          A[i-1][j] = (double)atoi(ui->table_input->item(i,j)->text().toStdString().c_str());
          if(A[i-1][j] >= 0){
            b[i-1] -= A[i-1][j];
          }else{
            b[i-1] += A[i-1][j];
          }
        }else if(j==n){
          b[i-1] = (double)atoi(ui->table_input->item(i,j)->text().toStdString().c_str());
        }else{
           if(ui->table_input->item(i,j)!=NULL){
             A[i-1][j] = (double)atoi(ui->table_input->item(i,j)->text().toStdString().c_str());
             //cout << A[i-1][j] << endl;
           }else{
             A[i-1][j] = 0;
           }
        }
    }
  }
  // print_constraint_equations();
  read_optimize_equation();
}

void Milp::print_constraint_equations(){
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
        cout << A[i][j] << " ";
    }
    cout << b[i] << endl;
  }
}

void Milp::update_optimize_table(){
  if(n != 0 && m != 0){
    ui->table_optimize->setShowGrid(true);
    ui->table_optimize->verticalHeader()->setVisible(false);
    ui->table_optimize->setRowCount(2);
    ui->table_optimize->setColumnCount(n+1);
    for(int i=0;i<n;i++){
      string str = "var" + to_string(i+1);
      QString s = QString::fromUtf8(str.c_str());
      ui->table_optimize->setItem(0, i, new QTableWidgetItem(s));
    }
    ui->table_optimize->setItem(0, n, new QTableWidgetItem("constant"));
    for(int i=0;i<n+1;i++){
      ui->table_optimize->setItem(1, i, new QTableWidgetItem(""));
    }
  }
  else{
    ui->table_optimize->setShowGrid(true);
    ui->table_optimize->verticalHeader()->setVisible(false);
    ui->table_optimize->setRowCount(0);
    ui->table_optimize->setColumnCount(0);
  }
}

void Milp::read_optimize_equation(){
  for(int i=0;i<n+1;i++){
    if(i==n){
        v = (double)atoi(ui->table_optimize->item(1,i)->text().toStdString().c_str());
    }else{
      if(ui->table_optimize->item(1,i)!=NULL){
        c[i] = (double)atoi(ui->table_optimize->item(1,i)->text().toStdString().c_str());
        //cout << A[i-1][j] << endl;
      }else{
        c[i] = 0;
      }
    }
  }
  // print_optimize_equation();
}

void Milp::print_optimize_equation(){
  for(int i=0;i<n;i++){
    cout << c[i] << " ";
  }
  cout << v << endl;
}

void Milp::solve(){
  read_equations();
  ret = simplex();
  //print_result_coefficients();
  display_result_coefficients();
}

void Milp::print_result_coefficients(){
  if (isinf(ret.second))
  {
      if (ret.first[0] == -1){
        printf("Objective function unbounded!\n");
      }
      else if (ret.first[0] == -2){
        printf("Linear program infeasible!\n");
      }
  }
  else
  {
    printf("Solution: (");
    for (int i=0;i<n;i++){

      printf("%lf%s", ret.first[i], (i < n + m - 1) ? ", " : ")\n");
    }
    printf("Optimal objective value: %lf\n", ret.second);
  }
}

void Milp::display_result_coefficients(){
  if (isinf(ret.second))
  {
      if (ret.first[0] == -1){
        ui->result_label->setText("The result is: Objective function is unbounded");
      }
      else if (ret.first[0] == -2){
        ui->result_label->setText("The result is: Linear program infeasible!");
      }
  }
  else
  {
    ui->result_label->setText("The result is: Linear program is feasible and solutions is in the table below");
    ui->table_result->setShowGrid(true);
    ui->table_result->verticalHeader()->setVisible(false);
    ui->table_result->setRowCount(2);
    ui->table_result->setColumnCount(n+1);
    for(int i=0;i<n;i++){
      string str = "var" + to_string(i+1);
      QString s = QString::fromUtf8(str.c_str());
      ui->table_result->setItem(0, i, new QTableWidgetItem(s));
    }
    ui->table_result->setItem(0, n, new QTableWidgetItem("OOV"));
    for (int i=0;i<n;i++){
      string str = to_string(ret.first[i]);
      QString s = QString::fromUtf8(str.c_str());
      ui->table_result->setItem(1, i, new QTableWidgetItem(s));
    }
    string str = to_string(ret.second);
    QString s = QString::fromUtf8(str.c_str());
    ui->table_result->setItem(1, n, new QTableWidgetItem(s));
    ui->oov_label->setText("The Optimised Objective Value(OOV) is " + s);
  }

}

void Milp::update_result_table(){
  ui->result_label->setText("The result is: NA");
  ui->table_result->setShowGrid(true);
  ui->table_result->verticalHeader()->setVisible(false);
  ui->table_result->setRowCount(0);
  ui->table_result->setColumnCount(0);
  ui->oov_label->setText("The Optimised Objective Value(OOV) is NA");
}


void Milp::update_all(){
  update_inputs();
  update_result_table();
  update_table_inputs();
  update_optimize_table();
}

void Milp::reset_all(){
  reset_inputs();
  reset_line_edits();
  update_all();
}


