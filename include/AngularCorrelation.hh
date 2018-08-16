#ifndef AngCorr_h
#define AngCorr_h 1

const int NUM_ANGCOR_ENTRIES = 100;

// The header file to evaluate the angular correlation.

class Angular_Correlation {
 private: // declare index variables.
  float I;
  float I_1;
  float I_2;
  int   L_1;
  int   L_2;
  float Del_1;
  float Del_2;

  float A21;
  float A22;
  float A41;
  float A42;

  // declare arrays for the coefficient tables.
  static int   II2[10][16];
  static int   II4[10][16];
  static int   LL2[6][6];
  static int   LL4[6][6];
  static float F2[63][10];
  static float F4[49][8];

  // declare arrays for the integrated-Legendre polynomial tables
  static float P2_integral_table[NUM_ANGCOR_ENTRIES];
  static float P4_integral_table[NUM_ANGCOR_ENTRIES];
  static float x_table[NUM_ANGCOR_ENTRIES];

  bool Verbose;
  bool parameters_in_range;

  float P2_integral_func(float x);
  float P4_integral_func(float x);

  // calculate A21, A22, A41, A42
  void CoeffCalc();

  // locate entry in lookup table
  int Locate(float y);

  // calculate table-based value of normalized sampling-distribution
  float Table_func(int i);

 public:
  // set the constructor function.
  Angular_Correlation(bool Verbose_in = false);

  Angular_Correlation(float io, float i1, float i2, int l1, int l2, float del1, float del2, bool Verbose_in = false);

  // set the destructor function.
  ~Angular_Correlation();

  // re-initialize member variables
  bool ReInit(float io, float i1, float i2, int l1, int l2, float del1, float del2);

  bool ValidParameters() {return parameters_in_range;}

  // evaluate the angular  correlation.
  float Evaluation(float);

  // sample from angular correlation distribution
  float Sample(float x_rnd);

  // functions to retrieve the private variable values.
  float GetI()     const {return I;     }
  float GetI_1()   const {return I_1;   }
  float GetI_2()   const {return I_2;   }
  int   GetL_1()   const {return L_1;   }
  int   GetL_2()   const {return L_2;   }
  float GetDel_1() const {return Del_1; }
  float GetDel_2() const {return Del_2; }

  // functions to set the private variable values.
  void SetI(float io)       {I     = io;   }
  void SetI_1(float i1)     {I_1   = i1;   }
  void SetI_2(float i2)     {I_2   = i2;   }
  void SetL_1(int l1)       {L_1   = l1;   }
  void SetL_2(int l2)       {L_2   = l2;   }
  void SetDel_1(float del1) {Del_1 = del1; }
  void SetDel_2(float del2) {Del_2 = del2; }
};

#endif
