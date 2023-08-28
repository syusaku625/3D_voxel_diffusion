#include"three_dimension.hpp"
#include"shapefunction.hpp"

using namespace std;

void ThreeDimensionalDiffusion::set_gauss_point()
{
    gauss.resize(1);
    gauss[0].resize(4);
    
    gauss[0][0] = 2.5e-1;
    gauss[0][1] = 2.5e-1;
    gauss[0][2] = 2.5e-1;
    gauss[0][3] = 2.5e-1;   
}

void ThreeDimensionalDiffusion::matrix_initialize()
{
    mass_centerlization.allocate(numOfNode);
    for(int i=0; i<numOfNode; i++){
      mass_centerlization(i) = 0e0;
    }

    C.allocate(numOfNode);
    for(int i=0; i<numOfNode; i++){
      C(i) = 0e0;
    }

    K.resize(numOfElm);
    MASS.resize(numOfElm);
    for(int i=0; i<numOfElm; i++){
      K[i].allocate(4,4);
      MASS[i].allocate(4,4);
    }

    for(int i=0; i<numOfElm; i++){
      for(int j=0; j<4; j++){
        for(int k=0; k<4; k++){
          K[i](j,k) = 0e0;
          MASS[i](j,k) = 0e0;
        }
      }
    }
}

void ThreeDimensionalDiffusion::calc_dxdr(int ic, std::vector<std::vector<double>> &node, std::vector<std::vector<int>> &element, std::vector<std::vector<double>> &dxdr, ARRAY2D<double> &dNdr)
{
  for(int k=0;k<3;k++){
    for(int l=0;l<3;l++){
      dxdr[k][l] = 0e0;
      for(int p=0;p<4;p++){
        dxdr[k][l] += dNdr(p,l) * node[element[ic][p]][k];
      }
    }
  }
}

void ThreeDimensionalDiffusion::calc_inverse_matrix_3x3(std::vector<std::vector<double>> &dxdr, std::vector<std::vector<double>> &drdx, double det)
{

  drdx[0][0] = dxdr[1][1]*dxdr[2][2] - dxdr[1][2]*dxdr[2][1];
  drdx[0][1] = dxdr[0][2]*dxdr[2][1] - dxdr[0][1]*dxdr[2][2];
  drdx[0][2] = dxdr[0][1]*dxdr[1][2] - dxdr[0][2]*dxdr[1][1];
  drdx[1][0] = dxdr[1][2]*dxdr[2][0] - dxdr[1][0]*dxdr[2][2];
  drdx[1][1] = dxdr[0][0]*dxdr[2][2] - dxdr[0][2]*dxdr[2][0];
  drdx[1][2] = dxdr[0][2]*dxdr[1][0] - dxdr[0][0]*dxdr[1][2];
  drdx[2][0] = dxdr[1][0]*dxdr[2][1] - dxdr[1][1]*dxdr[2][0];
  drdx[2][1] = dxdr[0][1]*dxdr[2][0] - dxdr[0][0]*dxdr[2][1];
  drdx[2][2] = dxdr[0][0]*dxdr[1][1] - dxdr[0][1]*dxdr[1][0];

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) drdx[i][j] = drdx[i][j] / det;
  }
}

void ThreeDimensionalDiffusion::calc_dNdx(ARRAY2D<double> &dNdx, ARRAY2D<double> &dNdr, vector<vector<double>> &drdx)
{
  for(int k=0; k<4; k++){
    for(int l=0; l<3; l++){
      dNdx(k,l) = 0.0;
      for(int p=0; p<3; p++){
        dNdx(k,l) += dNdr(k,p) * drdx[p][l];
      }
    }
  }
}

double ThreeDimensionalDiffusion::calcDeterminant_3x3(vector<vector<double>> &a)
{
  double det  = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[2][0] * a[0][1] * a[1][2]
              - a[2][0] * a[1][1] * a[0][2] - a[1][0] * a[0][1] * a[2][2] - a[0][0] * a[2][1] * a[1][2];
  return det;
}

void ThreeDimensionalDiffusion::boundary_setting(vector<double> R, double time_t)
{
  for(int i=0; i<numOfElm; i++){
    if(phi(i)>0.000001){
      if(time_t>=60.0 && time_t<=120.0){ 
        for(int j=0; j<element[i].size(); j++){
          C(element[i][j]) = ((1.0/60.0)*time_t-1.0-dt*R[element[i][j]]);
        }
      }
      if(time_t>120.0) {
        for(int j=0; j<element[i].size(); j++){
          C(element[i][j]) = ((-1.0/480.0)*time_t+(5.0/4.0)-dt*R[element[i][j]]);
        }
      }
    }
  }
}

void ThreeDimensionalDiffusion::time_step(vector<double> &R, double time_t)
{
  vector<ARRAY2D<double>> LHS(numOfElm);
  vector<ARRAY2D<double>> tmp_RHS(numOfElm);
  vector<ARRAY1D<double>> RHS(numOfElm);

  for(int i=0;i<numOfElm;i++){
    int numOfNodeInElm = element[i].size();
    LHS[i].allocate(numOfNodeInElm,numOfNodeInElm);
    tmp_RHS[i].allocate(numOfNodeInElm,numOfNodeInElm);
    RHS[i].allocate(numOfNodeInElm);
  }

  #pragma omp parallel for
  for(int i=0;i<numOfElm;i++){
    for(int j=0;j<element[i].size();j++){
      RHS[i](j) = 0e0;
      for(int k=0;k<element[i].size();k++){
        LHS[i](j,k) = 0e0;
        tmp_RHS[i](j,k) = 0e0;
      }
    }
  }

  cout << "preprocess" << endl;

  #pragma omp parallel for
  for(int i=0; i<numOfElm; i++){
    //cout << i << endl;
    ARRAY1D<double> N(4); 
    ARRAY2D<double> dNdr(4,3);
    vector<vector<double>> dxdr(3, vector<double>(3)), drdx(3, vector<double>(3));
    ARRAY2D<double> dNdx(4,3);
    ShapeFunction3D::C3D4_N(N,gauss[0][0],gauss[0][1],gauss[0][2],gauss[0][3]);
    ShapeFunction3D::C3D4_dNdr(dNdr,gauss[0][0],gauss[0][1],gauss[0][2],gauss[0][3]);
    
    calc_dxdr(i, node, element, dxdr, dNdr);
    double detJ = calcDeterminant_3x3(dxdr);

    //ofs << detJ << endl;
    calc_inverse_matrix_3x3(dxdr, drdx, detJ);
    calc_dNdx(dNdx, dNdr, drdx);

    for(int k=0; k<4; k++){
      for(int l=0; l<4; l++){
        for(int q=0; q<3; q++){
          K[i](k,l) += diffusion_coefficient * (dNdx(k,q)*dNdx(l,q)) * phi(i) * detJ * 1e0/6e0;
        }
        MASS[i](k,l) = N(k) * N(l) * phi(i) * detJ *1e0/6e0;
      }
    }

    cout << "calclated matrix" << endl;
    
    for(int k=0;k<element[i].size();k++){
      for(int l=0;l<element[i].size();l++){
        LHS[i](k,l) += (MASS[i](k,l)/dt + 0.5 * K[i](k,l)) * detJ * 1e0/6e0;
        tmp_RHS[i](k,l) += (MASS[i](k,l) / dt - 0.5 * K[i](k,l)) * detJ * 1e0/6e0;
      }
    }

    cout << "input LHS RHS" << endl;

    for(int k=0; k<element[i].size(); k++){
      for(int l=0; l<element[i].size(); l++){
        RHS[i](k) += tmp_RHS[i](k,l) * C(l);
      }
    }
  }

  pardiso.set_CSR_value1D(LHS,element,numOfNode,numOfElm,inb);
  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++) pardiso.b[i] = 0e0;

  cout << "pardiso set" << endl;
  ARRAY2D<int> bd_ip(1,1);
  bd_ip(0,0) = 0;

  pardiso.set_CSR_dirichlet_boundary_condition1D(1,bd_ip);
  pardiso.b[0] = 10;

  //for(int i=0;i<numOfElm;i++){
  //  for(int q=0;q<element[i].size();q++) pardiso.b[element[i][q]] += RHS[i](q);
  //}

  cout << "solve" << endl;
  pardiso.main(numOfNode,8);
  cout << "end solve" << endl;
  for(int i=0;i<numOfNode;i++){
    C(i) = pardiso.x[i];
    if(C(i)<0.0){
      C(i) = 0e0;
    }
  }
}

void ThreeDimensionalDiffusion::calc_K_and_M_matrix()
{
  //int numOfInElm = 4;
  //
  //#pragma omp parallel for
  //for(int i=0; i<numOfElm; i++){
  //  cout << i << endl;
  //  ARRAY1D<double> N(4); 
  //  ARRAY2D<double> dNdr(4,3);
  //  vector<vector<double>> dxdr(3, vector<double>(3)), drdx(3, vector<double>(3));
  //  ARRAY2D<double> dNdx(4,3);
  //  ShapeFunction3D::C3D4_N(N,gauss[0][0],gauss[0][1],gauss[0][2],gauss[0][3]);
  //  ShapeFunction3D::C3D4_dNdr(dNdr,gauss[0][0],gauss[0][1],gauss[0][2],gauss[0][3]);
  //  calc_dxdr(i, node, element, dxdr, dNdr);
  //  double detJ = calcDeterminant_3x3(dxdr);
  //  calc_inverse_matrix_3x3(dxdr, drdx, detJ);
  //  calc_dNdx(dNdx, dNdr, drdx);
  //  for(int k=0; k<4; k++){
  //    for(int l=0; l<4; l++){
  //      for(int q=0; q<3; q++){
  //        K[i](k,l) += diffusion_coefficient * (dNdx(k,q)*dNdx(l,q)) * phi(i) * detJ * 1e0/6e0;
  //      }
  //      MASS[i](k,l) = N(k) * N(l) * phi(i) * detJ *1e0/6e0;
  //    }
  //  }
  //}
  //PARDISO.set_CSR_value1D(LHS,element,numOfNode,numOfElm,inb);
  //
  //#pragma omp parallel for
  //for(int i=0;i<numOfNode;i++) PARDISO.b[i] = 0e0;
//
  //for(int i=0;i<numOfElm;i++){
  //  for(int q=0;q<element[i].size();q++) PARDISO.b[element[i].node[q]] += RHS[i](q);
  //}
//
  //PARDISO.main(numOfNode,8);
  //for(int i=0;i<numOfNode;i++) C(i) = PARDISO.x[i];
}