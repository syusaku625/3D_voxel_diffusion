#ifndef _SHAPE_FUNCTION_H_
#define _SHAPE_FUNCTION_H_

#include<vector>
#include"allocation.hpp"

class ShapeFunction2D{
    public:
        static void C2D4_N(std::vector<double> &N,const double &g1,const double &g2)
        {
            N[0] = 2.5e-1 * (1e+0-g1) * (1e+0-g2);
            N[1] = 2.5e-1 * (1e+0+g1) * (1e+0-g2);
            N[2] = 2.5e-1 * (1e+0+g1) * (1e+0+g2);
            N[3] = 2.5e-1 * (1e+0-g1) * (1e+0+g2);
        }
        static void C2D4_dNdr(std::vector<std::vector<double>> &dNdr,const double &g1,const double &g2)
        {
          dNdr[0][0] = -2.5e-1 * (1e+0-g2);
          dNdr[0][1] = -2.5e-1 * (1e+0-g1);
          dNdr[1][0] =  2.5e-1 * (1e+0-g2);
          dNdr[1][1] = -2.5e-1 * (1e+0+g1);
          dNdr[2][0] =  2.5e-1 * (1e+0+g2);
          dNdr[2][1] =  2.5e-1 * (1e+0+g1);
          dNdr[3][0] = -2.5e-1 * (1e+0+g2);
          dNdr[3][1] =  2.5e-1 * (1e+0-g1);
        }
};

class ShapeFunction3D{
 public:
  static void C3D4_N(ARRAY1D<double> &N,const double &L1,const double &L2,const double &L3,const double &L4)
  {
    N(0)=L1;
    N(1)=L2;
    N(2)=L3;
    N(3)=L4;
  }
  static void C3D8_N(ARRAY1D<double> &N,const double &g1,const double &g2,const double &g3)
  {
  N(0)= 1.25e-1 * (1e0-g1) * (1e0-g2) * (1e0-g3);
  N(1) = 1.25e-1 * (1e0+g1) * (1e0-g2) * (1e0-g3);
  N(2) = 1.25e-1 * (1e0+g1) * (1e0+g2) * (1e0-g3);
  N(3) = 1.25e-1 * (1e0-g1) * (1e0+g2) * (1e0-g3);
  N(4) = 1.25e-1 * (1e0-g1) * (1e0-g2) * (1e0+g3);
  N(5) = 1.25e-1 * (1e0+g1) * (1e0-g2) * (1e0+g3);
  N(6) = 1.25e-1 * (1e0+g1) * (1e0+g2) * (1e0+g3);
  N(7) = 1.25e-1 * (1e0-g1) * (1e0+g2) * (1e0+g3);
  }
  static void C3D4_dNdr(ARRAY2D<double> &dNdr,const double &L1,const double &L2,const double &L3,const double &L4)
  {
    dNdr(0,0)=-1e0; dNdr(0,1)=-1e0; dNdr(0,2)=-1e0;
    dNdr(1,0)= 1e0; dNdr(1,1)= 0e0; dNdr(1,2)= 0e0;
    dNdr(2,0)= 0e0; dNdr(2,1)= 1e0; dNdr(2,2)= 0e0;
    dNdr(3,0)= 0e0; dNdr(3,1)= 0e0; dNdr(3,2)= 1e0;
  }
  static void C3D8_dNdr(ARRAY2D<double> &dNdr,const double &g1,const double &g2,const double &g3)
  {
  dNdr(0,0) = -1.25e-1 * (1e0-g2) * (1e0-g3);
  dNdr(0,1) = -1.25e-1 * (1e0-g1) * (1e0-g3);
  dNdr(0,2) = -1.25e-1 * (1e0-g1) * (1e0-g2);
  dNdr(1,0) =  1.25e-1 * (1e0-g2) * (1e0-g3);
  dNdr(1,1) = -1.25e-1 * (1e0+g1) * (1e0-g3);
  dNdr(1,2) = -1.25e-1 * (1e0+g1) * (1e0-g2);
  dNdr(2,0) =  1.25e-1 * (1e0+g2) * (1e0-g3);
  dNdr(2,1) =  1.25e-1 * (1e0+g1) * (1e0-g3);
  dNdr(2,2) = -1.25e-1 * (1e0+g1) * (1e0+g2);
  dNdr(3,0) = -1.25e-1 * (1e0+g2) * (1e0-g3);
  dNdr(3,1) =  1.25e-1 * (1e0-g1) * (1e0-g3);
  dNdr(3,2) = -1.25e-1 * (1e0-g1) * (1e0+g2);
  dNdr(4,0) = -1.25e-1 * (1e0-g2) * (1e0+g3);
  dNdr(4,1) = -1.25e-1 * (1e0-g1) * (1e0+g3);
  dNdr(4,2) =  1.25e-1 * (1e0-g1) * (1e0-g2);
  dNdr(5,0) =  1.25e-1 * (1e0-g2) * (1e0+g3);
  dNdr(5,1) = -1.25e-1 * (1e0+g1) * (1e0+g3);
  dNdr(5,2) =  1.25e-1 * (1e0+g1) * (1e0-g2);
  dNdr(6,0) =  1.25e-1 * (1e0+g2) * (1e0+g3);
  dNdr(6,1) =  1.25e-1 * (1e0+g1) * (1e0+g3);
  dNdr(6,2) =  1.25e-1 * (1e0+g1) * (1e0+g2);
  dNdr(7,0) = -1.25e-1 * (1e0+g2) * (1e0+g3);
  dNdr(7,1) =  1.25e-1 * (1e0-g1) * (1e0+g3);
  dNdr(7,2) =  1.25e-1 * (1e0-g1) * (1e0+g2);
  }

};

#endif