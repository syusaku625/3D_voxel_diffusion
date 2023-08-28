#include "three_dimension.hpp"
#include"shapefunction.hpp"
#include"gauss.hpp"

using namespace std;

void ThreeDimensionalDiffusion::transform_point_data_to_cell_data(ARRAY1D<double> &element_C, ARRAY1D<double> &C)
{
  //Gauss g;
  //for(int i=0; i<numOfElm; i++){
  //  for(int j=0; j<element[0].size(); j++){
  //    ARRAY1D<double> N(element[0].size()); 
  //    ARRAY2D<double> dNdr(element[0].size(),3);
  //    vector<vector<double>> dxdr(3, vector<double>(3)), drdx(3, vector<double>(3));
  //    ShapeFunction3D::C3D4_N(N,0.5,0.5,0.5);
  //    ShapeFunction3D::C3D8_dNdr(dNdr,g.point[0][j],g.point[1][j],g.point[2][j]);
  //    calc_dxdr(i, node, element, dxdr, dNdr);
  //    double detJ = calcDeterminant_3x3(dxdr);
  //    for(int j=0; j<element[0].size(); j++){
  //      element_C(i) += N(j) * C(element[i][j]) * detJ;  
  //    }
  //  }
  //} 
}