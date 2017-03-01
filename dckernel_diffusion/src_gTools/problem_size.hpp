#ifndef PROBLEM_SIZE_H
#define PROBLEM_SIZE_H

#include <string>

//-------------------------------------------------------------------------------
//
// Problem size and global parameters
//
//-------------------------------------------------------------------------------

  int ADM_NSYS     = 32;
  int ADM_MAXFNAME = 1024;
  int ADM_LOG_FID  = 6;

  //--- Identifier of triangle element (i-axis-side or j-axis side)
  int TI  = 1;
  int TJ  = 2;

  //--- Identifier of line element (i-axis-side, ij-axis side, or j-axis side)
  int AI  = 1;
  int AIJ = 2;
  int AJ  = 3;

  //--- Identifier of 1 variable
  int K0  = 1;

  int ADM_nxyz = 3; // dimension of the spacial vector

  //--- region
  static const int ADM_lall      = 1;     // number of regular region per process
  int ADM_lall_pl   = 2;     // number of pole    region per process

  //--- horizontal grid
int ADM_gall      = 16900; // number of horizontal grid per regular region
int ADM_iall      = 130; // number of i points per regular region
int ADM_jall      = 130;//16900; // number of j points per regular region
int ADM_gall_1d   = 13;//130;   // number of horizontal grid (1D)
int ADM_gmin      = 2;     // start index of 1D horizontal grid
int ADM_gmax      = 129;   // end   index of 1D horizontal grid

  int ADM_gall_pl   = 6;     // number of horizontal grid for pole region
  int ADM_gslf_pl   = 1;     // index for pole point
  int ADM_gmin_pl   = 2;     // start index of grid around the pole point
  int ADM_gmax_pl   = 6;     // end   index of grid around the pole point

  //--- vertical grid
int ADM_vlayer    = 40;    // number of vertical layer
int ADM_kall      = 42;    // number of vertical grid
int ADM_kmin      = 2;     // start index of vertical grid
int ADM_kmax      = 41;    // end   index of vertical grid

  // NOTE: to run with a pole region and pentagon
  // set the following values to TRUE
  bool ADM_have_pl     = false; // this ID manages pole region?
  bool ADM_have_sgp[1] = {false}; // region have singlar point?
    //bool ADM_have_pl     = true; // this ID manages pole region?
  //bool ADM_have_sgp[1] = {true}; // region have singlar point?

  //--- mod_grd
  int GRD_XDIR = 1;
  int GRD_YDIR = 2;
  int GRD_ZDIR = 3;

  std::string vgrid_fname ="./vgrid40_600m_24km.dat";

  //--- mod_gmtr
  int GMTR_P_nmax_var = 10;

  int P_AREA  = 1;
  int P_RAREA = 2;
  int GMTR_P_IX    = 3;
  int GMTR_P_IY    = 4;
  int GMTR_P_IZ    = 5;
  int GMTR_P_JX    = 6;
  int GMTR_P_JY    = 7;
  int GMTR_P_JZ    = 8;
  int GMTR_P_LAT   = 9;
  int GMTR_P_LON   = 10;

  int GMTR_T_nmax_var = 7;

  int T_AREA  = 1;
  int T_RAREA = 2;
  int W1    = 3;
  int W2    = 4;
  int W3    = 5;
  int GMTR_T_LAT   = 6;
  int GMTR_T_LON   = 7;

  int GMTR_A_nmax_var    = 12;
  int GMTR_A_nmax_var_pl = 18;

  int HNX  = 1;
  int HNY  = 2;
  int HNZ  = 3;
  int HTX  = 4;
  int HTY  = 5;
  int HTZ  = 6;
  int TNX  = 7;
  int TNY  = 8;
  int TNZ  = 9;
  int TTX  = 10;
  int TTY  = 11;
  int TTZ  = 12;

  int TN2X = 13;
  int TN2Y = 14;
  int TN2Z = 15;
  int TT2X = 16;
  int TT2Y = 17;
  int TT2Z = 18;

  int SET_iteration = 1;
  int SET_prc_me    = 1;

  #endif
