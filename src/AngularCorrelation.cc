// $Id: AngularCorrelation.cc,v 1.1.1.1 2007/01/13 19:44:02 jordan Exp $

// The functions to evaluate angular correlation coefficients.

#include "AngularCorrelation.hh"
#include "AngularCorrelation_coeff.icc"

#include <iostream>
#include <fstream>
#include <cmath>
using std::cos;
using std::pow;
using std::cout;
using std::endl;


// constructor
Angular_Correlation::Angular_Correlation(bool Verbose_in)
  : Verbose(Verbose_in) {
  if (Verbose) {
    cout << "Default Angular_Correlation constructor is being called." << endl;
  }
}

Angular_Correlation::Angular_Correlation(float io, float i1, float i2, int l1, int l2, float del1, float del2, bool Verbose_in)
  : Verbose(Verbose_in) {
  if (Verbose) {
    cout << "Angular_Correlation object is being constructed." << endl;
  }

  ReInit(io, i1, i2, l1, l2, del1, del2);
}

// destructor
Angular_Correlation::~Angular_Correlation()
{ }

// setup member variables
bool Angular_Correlation::ReInit(float io, float i1, float i2, int l1, int l2, float del1, float del2) {
  I = io; I_1 = i1; I_2 = i2; L_1 = l1; L_2 = l2; Del_1 = del1; Del_2 = del2;

  parameters_in_range = true;

  // check the range of the index variables.
  // if the condition is not satisfied, abort the run.
  if (I < 1.0 || I_1 < 0.0 || I_2 < 0.0 || L_1 < 1 || L_2 < 1 ||
      I > 4.5 || I_1 > 7.5 || I_2 > 7.5 || L_1 > 4 || L_2 > 4) {
     parameters_in_range = false;

     if (true) {
      cout << " Angular_Correlation::ReInit Warning!!! " <<
                " : Check I, I_1, I_2, L_1, and L_2 for" <<
                " Angular Correlation Evaluation " << endl;
     if ((I < 1.0) || (I > 4.5)) {
       cout << "Intermediate angular momentum is out-of-range." << endl;
       cout << "I = " << I << endl;
       cout << "Valid range is 1.0 to 4.5." << endl;
     } else if ((I_1 < 0.0) || (I_1 > 7.5)) {
       cout << "Initial state angular momentum is out-of-range." << endl;
       cout << "I_1 = " << I_1 << endl;
       cout << "Valid range is 0 to 7.5." << endl;
     } else if ((I_2 < 0.0) || (I_2 > 7.5)) {
       cout << "Final state angular momentum is out-of-range." << endl;
       cout << "I_2 = " << I_2 << endl;
       cout << "Valid range is 0 to 7.5." << endl;
     } else if ((L_1 < 1) || (L_1 > 4)) {
       cout << "Initial gamma angular momentum is out-of-range." << endl;
       cout << "L_1 = " << L_1 << endl;
       cout << "Valid range is 1 to 4." << endl;
     } else if ((L_2 < 1) || (L_2 > 4)) {
       cout << "Final gamma angular momentum is out-of-range." << endl;
       cout << "L_2 = " << L_2 << endl;
       cout << "Valid range is 1 to 4." << endl;
     }
    }
  }


  if (parameters_in_range) {
    CoeffCalc();

    // Initialize tables of integrated P2 and P4 Legendre polynomials.

    const float dx = 2.0/(NUM_ANGCOR_ENTRIES - 1);

    for (int i = 0; i < NUM_ANGCOR_ENTRIES; i++) {
      float x = -1.0 + dx*i;
      x_table[i] = x;
      P2_integral_table[i] = P2_integral_func(x);
      P4_integral_table[i] = P4_integral_func(x);
    }
  }

  return parameters_in_range;
}


// calculate A21, A22, etc.
void Angular_Correlation::CoeffCalc() {
  int   J = int(2*I   + 0.0001);
  int J_1 = int(2*I_1 + 0.0001);
  int J_2 = int(2*I_2 + 0.0001);

  A21 = (F2[II2[J][J_1]][LL2[L_1][L_1]] + 2.0*Del_1*
         F2[II2[J][J_1]][LL2[L_1][L_1+1]] + Del_1*Del_1*
         F2[II2[J][J_1]][LL2[L_1+1][L_1+1]])/(1.0+Del_1*Del_1);

  A22 = (F2[II2[J][J_2]][LL2[L_2][L_2]] + 2.0*Del_2*
         F2[II2[J][J_2]][LL2[L_2][L_2+1]] + Del_2*Del_2*
         F2[II2[J][J_2]][LL2[L_2+1][L_2+1]])/(1.0+Del_2*Del_2);

  A41 = (F4[II4[J][J_1]][LL4[L_1][L_1]] + 2.0*Del_1*
         F4[II4[J][J_1]][LL4[L_1][L_1+1]] + Del_1*Del_1*
         F4[II4[J][J_1]][LL4[L_1+1][L_1+1]])/(1.0+Del_1*Del_1);

  A42 = (F4[II4[J][J_2]][LL4[L_2][L_2]] + 2.0*Del_2*
         F4[II4[J][J_2]][LL4[L_2][L_2+1]] + Del_2*Del_2*
         F4[II4[J][J_2]][LL4[L_2+1][L_2+1]])/(1.0+Del_2*Del_2);

  if (Verbose) {
    cout << " In Angular_Correlation::CoeffCalc()." << endl;

    cout << " I   = " << I   << endl;
    cout << " I_1 = " << I_1 << endl;
    cout << " I_2 = " << I_2 << endl;

    cout << " J   = " << J   << endl;
    cout << " J_1 = " << J_1 << endl;
    cout << " J_2 = " << J_2 << endl;

    cout << " A21 = " << A21 << endl;
    cout << " A22 = " << A22 << endl;
    cout << " A41 = " << A41 << endl;
    cout << " A42 = " << A42 << endl;
  }
}

// Sample from angular correlation distribution
float Angular_Correlation::Sample(float y_rnd) {
  // Input: Random variable, y_rnd, drawn from uniform
  // distribution on [0, 1].
  // Output: Value of cos(theta) sampled from
  // angular correlation distribution.

  int indx = Locate(y_rnd);

  if (Verbose) {
    cout << "In Angular_Correlation::Sample." << endl;
    cout << " y_rnd = " << y_rnd << endl;
    cout << " indx  = " << indx << endl;
    cout << " x_table[0]: " << x_table[0] << endl;
    cout << " x_table[last]: " << x_table[NUM_ANGCOR_ENTRIES-1]
      << endl;
    cout << " Table_func(0): " << Table_func(0) << endl;
    cout << " A21 = " << A21 << endl;
    cout << " A22 = " << A22 << endl;
    cout << " A41 = " << A41 << endl;
    cout << " A42 = " << A42 << endl;
    cout << " P2_integral_table[0]: "
      << P2_integral_table[0] << endl;
    cout << " P4_integral_table[0]: "
      << P4_integral_table[0] << endl;
  }

  float Xtab[2];
  float Sum[2];

  Xtab[0] = x_table[indx];
  Xtab[1] = x_table[indx+1];
  Sum[0]  = Table_func(indx);
  Sum[1]  = Table_func(indx+1);

  float dX_dSum = (Xtab[1] - Xtab[0])/(Sum[1] - Sum[0]);

  float x_interp = Xtab[0] + (y_rnd - Sum[0])*dX_dSum;

  if (Verbose) {
    cout << "  x_interp = " << x_interp << endl;
  }

  return x_interp;
}

// evaluate the angular correlation.
float Angular_Correlation::Evaluation(float theta) {
  float P2;
  float P4;

  P2 = (3.0*(pow(cos(theta), 2)) - 1.0)/2.0;
  P4 = (35.0*(pow(cos(theta), 4)) - 30.0*(pow(cos(theta), 2)) + 3.0)/8.0;

  return (1.0 + A21*A22*P2 + A41*A42*P4);
}


// Integral of Legendre function P2(x) with respect to x
float Angular_Correlation::P2_integral_func(float x) {
  // P2(x) = (1/2)(3*x^2 - 1)
  // --> Integral{P2(x) dx} = (1/2)(x^3 - x)

  float f = 0.5 * (pow(x, 3) - x);

  return f;
}

// Integral of Legendre function P4(x) with respect to x
float Angular_Correlation::P4_integral_func(float x) {
  // P4(x) = (1/8)(35*x^4 - 30*x^2 + 3)
  // --> Integral{P4(x) dx} = (1/8)(7*x^5 - 10*x^3 + 3*x)

  float f = 0.125 * (7.0*pow(x, 5) - 10.0*pow(x, 3) + 3.0*x);

  return f;
}


float Angular_Correlation::Table_func(int i) {
  float f = 0.5*(1.0 + x_table[i]
          + A21*A22*P2_integral_table[i]
          + A41*A42*P4_integral_table[i]);

  return f;
}

int Angular_Correlation::Locate(float y) {
  // Adapted from Numerical Recipes routine LOCATE for searching of
  // monotonically-ordered table.

  const int nmax = NUM_ANGCOR_ENTRIES;

  int JL = 0;
  int JU = nmax + 1;

  while ((JU-JL) > 1) {
    int JM = (JU + JL)/2;

    int i = JM-1;
    float f = Table_func(i);

    if (y > f)
      JL = JM;
    else
      JU = JM;
  }

  return JL-1;
}
