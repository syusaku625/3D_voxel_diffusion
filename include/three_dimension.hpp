#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include<vector>
#include<sys/stat.h>
#include"TextParser.h"
#include<string>
#include<sstream>
#include<fstream>
#include<iostream>
#include<omp.h>
#include<set>
#include"mkl.h"
#include<stdio.h>
#include<stdlib.h>
#include"allocation.hpp"
#include"pardiso_solver.h"

class ThreeDimensionalDiffusion{
    private:
        std::string material_judge;
    
    public:
        TextParser tp;
        int numOfNode, numOfElm;
        double computational_time, dt;
        double diffusion_coefficient;
        double coupling_coefficient_vc, coupling_coefficient_cv, coupling_coefficient_ci, coupling_coefficient_ic, coupling_coefficient_vi, coupling_coefficient_iv;
        std::vector<std::vector<double>> node;
        std::vector<std::vector<int>> element;
        std::vector<int> elementType;
        ARRAY1D<double> phi;
        ARRAY1D<double> C;
        std::vector<ARRAY2D<double>> K;
        std::vector<ARRAY2D<double>> MASS;
        std::vector<std::vector<int>> ieb, inb;
        int nnz;
        std::vector<int> ptr, index;
        std::vector<double> value;
        PARDISO_solver pardiso;
    
        std::vector<std::vector<double>> gauss;
        ARRAY1D<double> mass_centerlization;

        ThreeDimensionalDiffusion(std::string mat)
        {
            material_judge = mat;
        }
    
        void input_geometry();
        void input_coupling_coefficient();
        void read_node(std::string filename);  
        void read_element(std::string filename);
        void read_elementType(std::string filename);
        void input_phi();
        void input_time();
        void input_diffusion_coefficient();
        void calc_adjacent_elements();
        void calc_adjacent_nodes();
        void CSR_initialize(const std::vector<std::vector<int>> &inb,const int &numOfNode,const int &dim);
        void CSR_ptr_initialize(const std::vector<std::vector<int>> &inb,const int &numOfNode,const int &dim);
        void CSR_index_initialize(const std::vector<std::vector<int>> &inb,const int &numOfNode,const int &dim);
        void set_CSR_value1D(std::vector<ARRAY2D<double>> &K,const std::vector<std::vector<int>> &element,const int &numOfNode,
                               const int &numOfElm,const std::vector<std::vector<int>> &inb);
        //finite_element_module
        void set_gauss_point();
        void matrix_initialize();
        void calc_dxdr(int ic, std::vector<std::vector<double>> &node, std::vector<std::vector<int>> &element, std::vector<std::vector<double>> &dxdr, ARRAY2D<double> &dNdr);
        void calc_inverse_matrix_3x3(std::vector<std::vector<double>> &dxdr, std::vector<std::vector<double>> &drdx, double det);
        void calc_dNdx(ARRAY2D<double> &dNdx, ARRAY2D<double> &dNdr, std::vector<std::vector<double>> &drdx);
        double calcDeterminant_3x3(std::vector<std::vector<double>> &a);
        void boundary_setting(std::vector<double> R, double time_t);
        void time_step(std::vector<double> &R, double time_t);
        void export_vtu(const std::string &file);
        void calc_K_and_M_matrix();

        //O17_transport
        void transform_point_data_to_cell_data(ARRAY1D<double> &element_C, ARRAY1D<double> &C);
};

#endif