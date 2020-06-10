#include "milp.h"
#include "ui_milp.h"
#include <QApplication>


Milp::Milp(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::Milp){
      printf("gi");
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



void Milp::canonicalize ( vector <vector <double> > & A,
            vector <double>& B,
            vector <double>& C,
            vector <int>& var_r,
            double & obj
            )
{
  int m = A.size(), n = C.size();
  for (int r=0; r<m; r++) {
    int bc = var_r[r];
    if ( fabs(A[r][bc] - 1.0) > EPS) {
      double p = A[r][bc];
      for (int c=0; c<n; c++) A[r][c] /= p;
      B[r] /= p;
    }
    if (fabs(C[bc]) > EPS) {
      double p = C[bc];
      for (int c=0; c<n; c++) C[c] -= A[r][c] * p;
      obj -= B[r] * p;
    }
  }
}

bool Milp::pivoting ( vector <vector <double> > & A,
        vector <double>& B,
        vector <double>& C,
        vector <int>& var_r,
        double & obj
        )
{
  int m = A.size(), n = C.size();
 while (1) {
    int ev = 0;                    // entering variable id
    for (ev=0; ev<n; ev++)
      if (C[ev] < -EPS) break;

    if (ev == n) break;           // optimum reached.
          int lvr = -1;                  // leaving variable
    double minRatio;
    for (int r=0; r<m; r++) {
      if (A[r][ev] > EPS) {
    if ( lvr < 0 || B[r]/A[r][ev] < minRatio ) {
      lvr = r; minRatio = B[r] / A[r][ev];
    }
      }
    }
    if (lvr < 0) return true;      // unbounded
    int lv = var_r[lvr];       // leaving variable
    var_r[lvr] = ev;
    double p = A[lvr][ev];
    for (int c=0; c<n; c++) A[lvr][c] /= p; B[lvr] /= p;
    for (int r=0; r<m; r++) {
      if ( r != lvr && fabs (A[r][ev]) > EPS ) {
    double p2 = A[r][ev];
    for (int c=0; c<n; c++) A[r][c] -= p2 * A[lvr][c];
    B[r] -= p2 * B[lvr];
      }
    }
    if ( fabs (C[ev]) > EPS ) {
      double p2 = C[ev];
      for (int c=0; c<n; c++) C[c] -= p2 * A[lvr][c];
      obj -= p2 * B[lvr];
    }
 if (DEBUG) {
      for (int c=0; c<n; c++) cout << C[c] << "\t"; cout << obj << endl;
      for (int r=0; r<m; r++) {
    for (int c=0; c<n; c++) cout << A[r][c] << "\t";
    cout << B[r] << endl;
      }
      cout << endl;
    }
  }
  return false;
}


void Milp::LuSolver ( vector <vector <double> > & A, // matrix A
         vector <double>& B,           // b
         vector <double>& X            // x
         )
{
  int n = A.size();
  if (X.size() != n)
    X.resize (n);
vector <vector <double> > L ( n, vector<double> (n) );
  vector <vector <double> > U ( n, vector<double> (n) );
  for (int i=0; i<n; i++)
    L[i][i] = 1.0; // diagonals of L are 1's
  copy ( A[0].begin(), A[0].end(), U[0].begin() );
   for (int k=0; k<n-1; k++) {
    for (int i=k+1; i<n; i++) {
      double t = A[i][k];
      for (int j=0; j<k; j++)
    t -= ( L[i][j] * U[j][k] );
      L[i][k] = t / U[k][k];
    }
    for (int j=k+1; j<n; j++) {
      double t = A[k+1][j];
      for (int i=0; i<k+1; i++)
    t -= ( L[k+1][i] * U[i][j] );
      U[k+1][j] = t;
    }
  }
  for (int k=0; k<n; k++) {
    X[k] = B[k];
    for (int j=0; j<k; j++)
      X[k] -= ( X[j] * L[k][j] );
  }
  for (int k=n-1; k>=0; k--) {
    for (int j=k+1; j<n; j++)
      X[k] -= ( X[j] * U[k][j] );
    X[k] /= U[k][k];
  }
}


int Milp::preprocess ( vector <vector <double> > & A,     // constraint matrix
         vector <double>& B,               // right hand side
         vector <double>& X                // unknowns
         )
{
  int m = A.size ();                 // # of constraints
  int n = A[0].size ();              // # of variables
  vector <bool> is_redundant (m, false);  // flags for redundant constraint
  for (int r=0; r<m; r++) {
    bool allZero = true;
    for (int c=0; c<n; c++)
      if (fabs(A[r][c]) > EPS) { allZero = false; break; }
    if (allZero) {
      if (fabs(B[r]) > EPS) return -1;
      else is_redundant[r] = true;
    }
  }
  for (int i=0; i<m; i++) if (!is_redundant[i]) {
    for (int j=i+1; j<m; j++) if (!is_redundant[j]) {
      int c;
      double ratio = 0.0;
            for (c=0; c<n; c++) {
    if ( fabs(A[i][c]) < EPS && fabs(A[j][c]) < EPS )       // both are 0
      continue;
    else if ( fabs(A[i][c]) < EPS && fabs(A[j][c]) > EPS || // one is 0
         fabs(A[i][c]) > EPS && fabs(A[j][c]) < EPS )
      break;
    else {                                            // both are nonzero
      if ( fabs(ratio) < EPS )
        ratio = A[i][c] / A[j][c];
      else {
        if ( fabs (A[i][c]/A[j][c] - ratio) > EPS )
          break;
      }
    }
      }  if (c == n) {
    if ( fabs(B[i]) < EPS && fabs(B[j]) < EPS ||
         fabs(B[j]) > EPS && fabs (B[i]/B[j] - ratio) < EPS )
      is_redundant[j] = true;
    else return -1;              // inconsistency detected
      }
    }
  }
  int r;
  for (int c=0; c<n; c++)
     r = column_identity (A, c);
if(r==-1)
count==0;
else
count==1;

  int redundancies = count;
  if (redundancies > 0) {
    int ir = 0;
    for (int i=0; i<m; i++) {
      if (!is_redundant[i]) {
    if (ir < i) {                // overiding
      copy (A[i].begin(), A[i].end(), A[ir].begin());
      B[ir] = B[i];
    }
    ir++;
      }
    }
    for (int i=0; i<redundancies; i++) {
      A.erase (A.end()-1);
      B.erase (B.end()-1);
    }
  }
 m -= redundancies;
  if (m >= n) {
    vector <vector <double> > A0 (n, vector<double> (n));
    vector <double> B0 (n);
    for (int r=0; r<n; r++) {
      copy (A[r].begin(), A[r].end(), A0[r].begin());
      B0[r] = B[r];
    }

    LuSolver (A0, B0, X);
    bool nonNegative = true;
    for (int c=0; c<n; c++)
      if (X[c] < 0) { nonNegative = false; break; }
    if (!nonNegative)
      return -1;
    bool consistent = true;
    for (int r=n; r<m; r++) {
      double lhs = 0.0;
      for (int c=0; c<n; c++)
    lhs += A[r][c] * X[c];
      if ( fabs (lhs - B[r]) > EPS ) { // unsatisfied constraint c
    consistent = false;
    break;
      }
    }
    return (consistent ? -2 : -1);
  }
        return redundancies;
}


int Milp::simplex ( const vector <vector <double> > & A,
          const vector <double>& B,
          const vector <double>& C,
          vector <double>& X,
          double & obj
          )
{
  int m = A.size();
  int n = A[0].size();
  if (!m || m != B.size() || n != C.size()) {
     exit(1);
  }
  if (X.size() != n) X.resize(n);
  fill (X.begin(), X.end(), 0);
  vector <vector <double> > A0 ( m, vector<double>(n) );
  vector <double> B0 (m);
  for (int r=0; r<m; r++)
    copy (A[r].begin(), A[r].end(), A0[r].begin() );
  copy ( B.begin(), B.end(), B0.begin() );
  int ret_val = preprocess (A0, B0, X);
  int redundancies;
  if (ret_val == -1)
    return -1;
  else if (ret_val == -2)
    return 1;
  else
    redundancies = ret_val;

  m = A0.size ();
 vector <bool> is_basic (n, false);
  vector <int> var_r (m, -1);
  int num_var_basic = 0;
  for (int c=0; c<n; c++) {
    int r = column_identity (A, c);

    if (r >= 0 && var_r[r] < 0) {
      is_basic[c] = true;
      var_r[r] = c;
      num_var_basic++;
    }
  }
  vector <vector <double> > A2 ( m, vector<double>(n) );
  vector <double> B2 (m);
  vector <double> C2 (n);
  for (int r=0; r<m; r++)
    copy ( A0[r].begin(), A0[r].end(), A2[r].begin() );
  copy ( B0.begin(), B0.end(), B2.begin() );
  for (int c=0; c<n; c++)
    C2[c] = -C[c];
  obj = 0;
  if (num_var_basic < m) {
    int n1 = n;
    vector <vector <double> > A1 (m, vector<double>(n) );
    vector <double> B1 (m);
    vector <double> C1 (n, 0);       // phase I new objective vector
 for (int r=0; r<m; r++)
      copy ( A0[r].begin(), A0[r].end(), A1[r].begin() );
    copy ( B0.begin(), B0.end(), B1.begin() );
    for (int i=0; i<m; i++) {
      if (var_r[i] < 0) {
    for (int r=0; r<m; r++) {
      if (r == i) A1[r].push_back (1);
      else A1[r].push_back (0);
    }
    C1.push_back (1);
    var_r[i] = n1;
    n1++;
      }
    }
     C1.resize (n1, 1);                //  objective vector size adjusted
canonicalize (A1, B1, C1, var_r, obj);
    bool unbounded = pivoting (A1, B1, C1, var_r, obj);  // pivoting
    if (unbounded) {
      exit (1);
    }    bool feasible = (fabs(obj) < EPS) ? true : false;
    if (!feasible) return 0;
    for (int r=0; r<m; r++) {
      for (int c=0; c<n; c++)
    A2[r][c] = A1[r][c];
      B2[r] = B1[r];
    }
  }
  canonicalize (A2, B2, C2, var_r, obj);
  bool unbounded = pivoting (A2, B2, C2, var_r, obj);
 for (int r=0; r<m; r++)              // r.h.s. is the basic solution
    X[var_r[r]] = B2[r];
  return ( unbounded ? -1 : 1 );
}


void Milp::update_inputs(){
  n = atoi(ui->variable_input->text().toStdString().c_str());
  m = atoi(ui->constraints_input->text().toStdString().c_str());
  ui->variable_label->setText(ui->variable_input->text());
  ui->constraints_label->setText(ui->constraints_input->text());
}

void Milp::reset_inputs(){
  n = m = 0;
  ui->variable_label->setText(ui->variable_input->text());
  ui->constraints_label->setText(ui->constraints_input->text());
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
    ui->table_input->setColumnCount(n+1);
    for(int i=0;i<n;i++){
      string str = "var" + to_string(i+1);
      QString s = QString::fromUtf8(str.c_str());
      ui->table_input->setItem(0, i, new QTableWidgetItem(s));
    }
    ui->table_input->setItem(0, n, new QTableWidgetItem("constant"));
    for(int i=1;i<m+1;i++){
      for(int j=0;j<n+1;j++){
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
    for(int j=0;j<n+1;j++){
        if(j==n){
          if(ui->table_input->item(i,j)!=NULL){
            b[i-1] = (double)atoi(ui->table_input->item(i,j)->text().toStdString().c_str());
            //cout << A[i-1][j] << endl;
          }else{
            b[i-1] = 0;
          }
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
    ui->table_optimize->setColumnCount(n);
    for(int i=0;i<n;i++){
      string str = "var" + to_string(i+1);
      QString s = QString::fromUtf8(str.c_str());
      ui->table_optimize->setItem(0, i, new QTableWidgetItem(s));
    }
    for(int i=0;i<n;i++){
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
  for(int i=0;i<n;i++){
    if(i==n-1){
        if(ui->table_optimize->item(1,i)!=NULL){
          v = (double)atoi(ui->table_optimize->item(1,i)->text().toStdString().c_str());
          //cout << A[i-1][j] << endl;
        }else{
          v = 0;
        }
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
  //print_constraint_equations();
  vector <vector <double> > Ar (m, vector<double>(n));
  vector <double> Br (m), Cr(n);

  for(int i =0;i<m;i++){
      for(int j =0;j < n;j++){
        Ar[i][j] = A[i][j];
        //cout << Ar[i][j] << " ";
      }
      //cout << endl;
      Br[i] = b[i];
  }
  for(int i=0;i<n;i++){
      Cr[i] = c[i];
//      cout << Cr[i] << " ";
  }
//  cout << endl;
//  for(int i =0;i<m;i++){
//      cout << Br[i] << " ";
//  }
//  cout << endl;
  unbounded = simplex(Ar, Br, Cr, X, obj);
  //print_result_coefficients();
  display_result_coefficients();
}

void Milp::print_result_coefficients(){

}

void Milp::display_result_coefficients(){
  if (unbounded==-1)
  {
    ui->result_label->setText("The result is: Linear program is unbounded");
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
      string str = to_string(X[i]);
      QString s = QString::fromUtf8(str.c_str());
      ui->table_result->setItem(1, i, new QTableWidgetItem(s));
    }
    string str = to_string(obj);
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


