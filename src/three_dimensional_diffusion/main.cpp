#include"three_dimension.hpp"
#include"shapefunction.hpp"
#include"gauss.hpp"

using namespace std;

//void calc_matrix(ThreeDimensionalDiffusion &Fluid, ThreeDimensionalDiffusion &Solid)
//{
//  int numOfInElm = 4;
//  vector<ARRAY2D<double>> Fluid_K(Fluid.numOfElm);
//  vector<ARRAY2D<double>> Solid_K(Fluid.numOfElm);
//  for(int i=0; i<Fluid.numOfElm; i++){
//    Fluid_K[i].allocate(4,4);
//    Solid_K[i].allocate(4,4);
//  }
//
//  #pragma omp parallel for
//  for(int i=0; i<Fluid.numOfElm; i++){
//    cout << i << endl;
//    ARRAY1D<double> N(4); 
//    ARRAY2D<double> dNdr(4,3);
//    vector<vector<double>> dxdr(3, vector<double>(3)), drdx(3, vector<double>(3));
//    ARRAY2D<double> dNdx(4,3);
//    ShapeFunction3D::C3D4_N(N,Fluid.gauss[0][0],Fluid.gauss[0][1],Fluid.gauss[0][2],Fluid.gauss[0][3]);
//    ShapeFunction3D::C3D4_dNdr(dNdr,Fluid.gauss[0][0],Fluid.gauss[0][1],Fluid.gauss[0][2],Fluid.gauss[0][3]);
//    Fluid.calc_dxdr(i, Fluid.node, Fluid.element, dxdr, dNdr);
//    Fluid.calc_inverse_matrix_3x3(dxdr, drdx);
//    double detJ = Fluid.calcDeterminant_3x3(dxdr);
//    Fluid.calc_dNdx(dNdx, dNdr, drdx);
//    for(int k=0; k<4; k++){
//      for(int l=0; l<4; l++){
//        for(int q=0; q<3; q++){
//            Fluid_K[i](k,l) += Fluid.diffusion_coefficient * (dNdx(k,q)*dNdx(l,q)) * Fluid.phi(i) * detJ * 1e0/6e0;
//            Solid_K[i](k,l) += Solid.diffusion_coefficient * (dNdx(k,q)*dNdx(l,q)) * Solid.phi(i) * detJ * 1e0/6e0;
//        }
//        Fluid.mass_centerlization(Fluid.element[i][k]) += N(k) * N(l) * Fluid.phi(i) * detJ * 1e0/6e0;
//        Solid.mass_centerlization(Fluid.element[i][k]) += N(k) * N(l) * Solid.phi(i) * detJ * 1e0/6e0;
//      }
//    }
//  }
//
//  Fluid.set_CSR_value1D(Fluid_K,Fluid.element,Fluid.numOfNode,Fluid.numOfElm,Fluid.inb);
//  Solid.set_CSR_value1D(Solid_K,Solid.element,Solid.numOfNode,Solid.numOfElm,Solid.inb);
//}

void calc_mass_matrix_each_phase(ThreeDimensionalDiffusion &Fluid, vector<vector<double>> &mass)
{
  Gauss g;
  mass.resize(Fluid.numOfElm);
  int numOfNodeInElm = Fluid.element[0].size();
  #pragma omp parallel for
  for(int i=0; i<Fluid.numOfElm; i++){
    mass[i].resize(numOfNodeInElm);
    for(int j=0; j<numOfNodeInElm; j++){
      ARRAY1D<double> N(numOfNodeInElm); 
      ARRAY2D<double> dNdr(numOfNodeInElm,3);
      vector<vector<double>> dxdr(3, vector<double>(3)), drdx(3, vector<double>(3));
      ARRAY2D<double> dNdx(numOfNodeInElm, 3);
      ShapeFunction3D::C3D8_N(N,g.point[0][j],g.point[1][j],g.point[2][j]);
      ShapeFunction3D::C3D8_dNdr(dNdr,g.point[0][j],g.point[1][j],g.point[2][j]);
      Fluid.calc_dxdr(i, Fluid.node, Fluid.element, dxdr, dNdr);
      double detJ = Fluid.calcDeterminant_3x3(dxdr);
      Fluid.calc_inverse_matrix_3x3(dxdr, drdx, detJ);
      Fluid.calc_dNdx(dNdx, dNdr, drdx);
      for(int k=0; k<numOfNodeInElm; k++){
        for(int l=0; l<numOfNodeInElm; l++){
          mass[i][k] += N(k) * N(l) * detJ;
        }
      }
    }
  }
}

ARRAY1D<double> mat_matrix_vec(vector<double> &M, const vector<double> &C_value)
{
  ARRAY1D<double> ans(C_value.size());
  for(int i=0; i<C_value.size(); i++){
    ans(i) = M[i] * C_value[i];
  }
  return ans;
}

int main(int argc, char *argv[])
{
    if(argc!=3){
        cout << "please input tp file & omp_threads" << endl;
        exit(1);
    }
    string input_file = argv[1];
    omp_set_num_threads(stoi(argv[2]));

    ThreeDimensionalDiffusion Fluid("F"), Solid("S"), Vessel("V");

    int ierror;
    if ((ierror = Fluid.tp.read(input_file)) != TP_NO_ERROR){
        printf("\tError at reading '%s' file\n", input_file.c_str());
        exit(1);
    }
    if ((ierror = Solid.tp.read(input_file)) != TP_NO_ERROR){
        printf("\tError at reading '%s' file\n", input_file.c_str());
        exit(1);
    }
    if ((ierror = Vessel.tp.read(input_file)) != TP_NO_ERROR){
        printf("\tError at reading '%s' file\n", input_file.c_str());
        exit(1);
    }

    Fluid.input_geometry(); Solid.input_geometry(); Vessel.input_geometry();
    
    Fluid.input_coupling_coefficient(); Solid.input_coupling_coefficient(); Vessel.input_coupling_coefficient();
    
    Fluid.input_phi(); Solid.input_phi(); Vessel.input_phi();
    
    Fluid.input_time(); Solid.input_time(); Vessel.input_time();
    
    Fluid.input_diffusion_coefficient(); Solid.input_diffusion_coefficient();
    Fluid.pardiso.initialize(Fluid.numOfNode); 
    Solid.pardiso.initialize(Solid.numOfNode); 
    Fluid.pardiso.CSR_initialize(Fluid.inb,Fluid.numOfNode,1); 
    Solid.pardiso.CSR_initialize(Solid.inb,Solid.numOfNode,1); 


    Vessel.export_vtu("vessel_test.vtu");
    Fluid.export_vtu("fluid_test.vtu");
    Solid.export_vtu("solid_test.vtu");

    //fem
    Fluid.set_gauss_point(); Solid.set_gauss_point(); Vessel.set_gauss_point();
    Vessel.matrix_initialize(); Fluid.matrix_initialize(); Solid.matrix_initialize();

    //csr initialize
    Fluid.CSR_initialize(Fluid.inb,Fluid.numOfNode,1); 
    Solid.CSR_initialize(Solid.inb,Solid.numOfNode,1); 

    int loop = Fluid.computational_time/Fluid.dt;

    //vector<vector<double>> mass_fluid, mass_solid, mass_vessel;
    //calc_mass_matrix_each_phase(Fluid, mass_fluid);
    //calc_mass_matrix_each_phase(Solid, mass_solid);
    //calc_mass_matrix_each_phase(Vessel, mass_vessel);

    //ARRAY1D<double> C_sum(Vessel.numOfNode);

    cout << "main loop" << endl;

    vector<double> Fluid_R(Fluid.numOfNode, 0.0);
    vector<double> Solid_R(Solid.numOfNode, 0.0);
    vector<double> Vessel_R(Vessel.numOfNode, 0.0);

    for(int i=0; i<loop; i++){
      cout << "time: "  << i*Fluid.dt << endl;  
      if(i*Fluid.dt<60.0) continue;
      ARRAY1D<double> element_C_vessel(Vessel.numOfElm), element_C_Fluid(Fluid.numOfElm), element_C_Solid(Solid.numOfElm);
      //Vessel.transform_point_data_to_cell_data(element_C_vessel, Vessel.C);
      //Fluid.transform_point_data_to_cell_data(element_C_Fluid, Fluid.C);
      //Solid.transform_point_data_to_cell_data(element_C_Solid, Solid.C);
      //vector<vector<double>> all_Q_vc, all_Q_vi, all_Q_ci;

      for(int j=0; j<Fluid.numOfElm; j++){
        //vector<double> Q_vc(Fluid.element[0].size()), Q_vi(Fluid.element[0].size()), Q_ci(Fluid.element[0].size());
        //vector<double> C_vessel(Fluid.element[0].size(), 0.0), C_CSF(Fluid.element[0].size(), 0.0), C_ISF(Fluid.element[0].size(), 0.0);
        //for(int k=0; k<Fluid.element[0].size(); k++){
        //  C_vessel[k] = Vessel.C(Vessel.element[j][k]); 
        //  C_CSF[k] = Fluid.C(Fluid.element[j][k]); 
        //  C_ISF[k] = Solid.C(Solid.element[j][k]);
        //}
        ////vessel to fluid
        //if(element_C_vessel(j)>=element_C_Fluid(j)){
        //  ARRAY1D<double> MvCv = mat_matrix_vec(mass_vessel[j],C_vessel);
        //  ARRAY1D<double> McCc = mat_matrix_vec(mass_fluid[j],C_CSF);
        //  for(int l=0; l<Fluid.element[0].size(); l++){
        //    Q_vc[l] = Fluid.coupling_coefficient_vc * sqrt(Vessel.phi(j)) * (MvCv(l) - McCc(l));
        //  }
        //  all_Q_vc.push_back(Q_vc);
        //}
        ////csf to vessel
        //else if(element_C_vessel(j)<element_C_Fluid(j)){
        //  ARRAY1D<double> MvCv = mat_matrix_vec(mass_vessel[j],C_vessel);
        //  ARRAY1D<double> McCc = mat_matrix_vec(mass_fluid[j],C_CSF);
        //  for(int l=0; l<Fluid.element[0].size(); l++){
        //    Q_vc[l] = Fluid.coupling_coefficient_cv * sqrt(Vessel.phi(j)) * (MvCv(l) - McCc(l));
        //  }
        //  all_Q_vc.push_back(Q_vc);
        //}
        //else{
        //  cout << "error csf to vessel" << endl;
        //  cout << element_C_vessel(j) << " " << element_C_Fluid(j) << endl;
        //  exit(1);
        //}
        //vessel to solid
        //if(element_C_vessel(j)>=element_C_Solid(j)){
        //  ARRAY1D<double> MvCv = mat_matrix_vec(mass_vessel[j],C_vessel);
        //  ARRAY1D<double> MiCi = mat_matrix_vec(mass_solid[j],C_ISF);
        //  for(int l=0; l<4; l++){
        //    Q_vi[l] = Solid.coupling_coefficient_vi*sqrt(Vessel.phi(j))*(MvCv(l)-MiCi(l));
        //  }
        //  all_Q_vi[j]=(Q_vi);
        //}
        ////solid to vessel
        //else if(element_C_vessel(j)<element_C_Solid(j)){
        //  ARRAY1D<double> MvCv = mat_matrix_vec(mass_vessel[j],C_vessel);
        //  ARRAY1D<double> MiCi = mat_matrix_vec(mass_solid[j],C_ISF);
        //  for(int l=0; l<4; l++){
        //    Q_vi[l] = Solid.coupling_coefficient_iv*sqrt(Vessel.phi(j))*(MvCv(l)-MiCi(l));
        //  }
        //  all_Q_vi[j]=(Q_vi);
        //}
        //else{
        //  cout << "error solid to vessel" << endl;
        //  cout << element_C_vessel(j) << " " << element_C_Solid(j) << endl;
        //  exit(1);
        //}
        //fluid to solid
        //if(element_C_Fluid(j)>=element_C_Solid(j)){
        //  ARRAY1D<double> McCc = mat_matrix_vec(mass_fluid[j],C_CSF);
        //  ARRAY1D<double> MiCi = mat_matrix_vec(mass_solid[j],C_ISF);
        //  for(int l=0; l<4; l++){
        //    Q_ci[l] = Fluid.coupling_coefficient_ci*sqrt(Fluid.phi(j))*(McCc(l)-MiCi(l));
        //  }
        //  all_Q_ci[j]=(Q_ci);
        //}
        ////solid to fluid
        //else if(element_C_Fluid(j)<element_C_Solid(j)){
        //  ARRAY1D<double> McCc = mat_matrix_vec(mass_fluid[j],C_CSF);
        //  ARRAY1D<double> MiCi = mat_matrix_vec(mass_solid[j],C_ISF);
        //  for(int l=0; l<4; l++){
        //    Q_ci[l] = Fluid.coupling_coefficient_ic*sqrt(Fluid.phi(j))*(McCc(l)-MiCi(l));
        //  }
        //  all_Q_ci[j]=(Q_ci);
        //}
        //else{
        //  cout << "error solid to fluid" << endl;
        //  cout << element_C_Fluid(j) << " " << element_C_Solid(j) << endl;
        //  exit(1);
        //}
      }
//
//
      //for(int j=0; j<Fluid.numOfElm; j++){
      //  for(int k=0; k<4; k++){
      //    Fluid_R[Fluid.element[j][k]] += (-all_Q_vc[j][k]+all_Q_ci[j][k]);
      //    Solid_R[Solid.element[j][k]] += (-all_Q_vi[j][k]-all_Q_ci[j][k]);
      //    Vessel_R[Vessel.element[j][k]] += (all_Q_vc[j][k]+all_Q_vi[j][k]);
      //  }
      //}

      //vector<double> Fluid_R_n_1(Fluid.numOfNode, 0.0);
      //vector<double> Vessel_R_n_1(Vessel.numOfNode, 0.0);
      //vector<double> Solid_R_n_1(Solid.numOfNode, 0.0);

      //for(int j=0; j<Fluid.numOfNode; j++){
      //  Fluid_R_n_1[j] = Fluid_R[j];
      //  Vessel_R_n_1[j] = Vessel_R[j];
      //  //Solid_R_n_1[j] = Solid_R[j];
      //  Fluid_R[j] = 0e0;
      //  Vessel_R[j] = 0e0;
      //  //Solid_R[j] = 0e0;
      //}
//
      //for(int j=0; j<Fluid.numOfElm; j++){
      //  for(int k=0; k<Fluid.element[0].size(); k++){
      //    Fluid_R[Fluid.element[j][k]] += (-all_Q_vc[j][k]);
      //    Vessel_R[Vessel.element[j][k]] += (all_Q_vc[j][k]);
      //  }
      //}

      //Vessel.boundary_setting(Vessel_R, Vessel.dt*i);
      cout << "check" << endl;
      Fluid.time_step(Fluid_R, Vessel.dt*i);
      //Solid.time_step(Solid_R, Vessel.dt*i);
      //cout << i << endl;
      if(i%1==0){
        string export_file = "export_vessel_" + to_string((i)) + ".vtu"; 
        Vessel.export_vtu(export_file);
        export_file = "export_fluid_" + to_string(i) + ".vtu"; 
        Fluid.export_vtu(export_file);
        export_file = "export_solid_" + to_string(i) + ".vtu"; 
        //Solid.export_vtu(export_file);
      }
      //if(i==1) exit(1);
    }
}