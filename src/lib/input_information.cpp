#include "three_dimension.hpp"

using namespace std;

void ThreeDimensionalDiffusion::read_node(string filename)
{
    ifstream ifs(filename);
    if(!ifs){
        cout << "error reading " + filename << endl;
        exit(1);
    }
    string str;
    while(getline(ifs, str)){
        istringstream ss(str);
        vector<double> tmp_node;
        for(int i=0; i<3; i++){
            getline(ss, str, ' ');
            tmp_node.push_back(stod(str)*1e-3);
        }
        node.push_back(tmp_node);
    }
    numOfNode = node.size();
}

void ThreeDimensionalDiffusion::read_element(string filename)
{
    ifstream ifs(filename);
    if(!ifs){
        cout << "error reading " + filename << endl;
        exit(1);
    }
    string str;
    while(getline(ifs, str)){
        istringstream ss(str);
        vector<int> tmp_element;
        for(int i=0; i<5; i++){
            getline(ss, str, ' ');
            if(i==4){
              elementType.push_back(stoi(str));
              continue;
            }
            tmp_element.push_back(stod(str));
        }
        element.push_back(tmp_element);
    }
    numOfElm = element.size();
}

void ThreeDimensionalDiffusion::read_elementType(string filename)
{
    ifstream ifs(filename);
    if(!ifs){
        cout << "error reading " + filename << endl;
        exit(1);
    }
    string str;
    while(getline(ifs,str)){
        elementType.push_back(stoi(str));
    }
}

void ThreeDimensionalDiffusion::input_geometry()
{
    string base_dir, node_file, element_file, elementType_file;
    string str,base_label,label;
    base_label = "/Geometry";
    label = base_label + "/base_dir";
    if ( !tp.getInspectedValue(label,base_dir)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    label = base_label + "/node_file";
    if ( !tp.getInspectedValue(label,node_file)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    label = base_label + "/element_file";
    if ( !tp.getInspectedValue(label,element_file)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    label = base_label + "/elementType";
    if ( !tp.getInspectedValue(label,elementType_file)){
      cout << label << " is not set" << endl;
      exit(0);
    }

    node_file = base_dir + "/" + node_file;
    element_file = base_dir + "/" + element_file;
    elementType_file = base_dir + "/" + elementType_file;

    read_node(node_file);
    read_element(element_file);

    ieb.resize(numOfNode,vector<int>(0));
    inb.resize(numOfNode,vector<int>(0));

    calc_adjacent_elements();
    calc_adjacent_nodes();
}

void ThreeDimensionalDiffusion::input_phi()
{
    phi.allocate(element.size());
    for(int i=0; i<elementType.size(); i++){
        string base_label, label;
        double vessel_phi_tmp, solid_phi_tmp, fluid_phi_tmp;
        if(elementType[i]==4) base_label = "/WhiteMatter";
        else if(elementType[i]==3) base_label = "/GreyMatter";
        else if(elementType[i]==2) base_label = "/CSF";

        label = base_label + "/vessel";
        if ( !tp.getInspectedValue(label,vessel_phi_tmp)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/solid";
        if ( !tp.getInspectedValue(label,solid_phi_tmp)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/fluid";
        if ( !tp.getInspectedValue(label,fluid_phi_tmp)){
          cout << label << " is not set" << endl;
          exit(0);
        }

        if(material_judge=="F"){
            phi(i) = fluid_phi_tmp;
        }
        else if(material_judge=="S"){
            phi(i) = solid_phi_tmp;
        }
        else if(material_judge=="V"){
            phi(i) = vessel_phi_tmp;
        }
    }
}

void ThreeDimensionalDiffusion::input_time()
{
    string base_label, label;
    double dt_tmp, time_tmp;
    base_label = "/Condition";
    label = base_label + "/dt";
    if ( !tp.getInspectedValue(label,dt_tmp)){
      cout << label << " is not set" << endl;
      exit(0);
    }

    label = base_label + "/computational_time";
    if ( !tp.getInspectedValue(label,time_tmp)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    dt = dt_tmp;
    computational_time = time_tmp;
}

void ThreeDimensionalDiffusion::input_diffusion_coefficient()
{
  string base_label, label;
  base_label = "/Condition";
  if(material_judge == "F"){
    label = base_label + "/fluid_diffusion_coefficient";
    if ( !tp.getInspectedValue(label,diffusion_coefficient)){
      cout << label << " is not set" << endl;
      exit(0);
    }
  }
  else if(material_judge == "S"){
    label = base_label + "/solid_diffusion_coefficient";
    if ( !tp.getInspectedValue(label,diffusion_coefficient)){
      cout << label << " is not set" << endl;
      exit(0);
    }
  }
}

void ThreeDimensionalDiffusion::input_coupling_coefficient()
{
  string base_label, label;
  base_label = "/Coupling_coefficient";
  label = base_label + "/coupling_coefficient_vc";
  if ( !tp.getInspectedValue(label,coupling_coefficient_vc)){
    cout << label << " is not set" << endl;
    exit(1);
  }
  label = base_label + "/coupling_coefficient_cv";
  if ( !tp.getInspectedValue(label,coupling_coefficient_cv)){
    cout << label << " is not set" << endl;
    exit(1);
  }
  label = base_label + "/coupling_coefficient_ci";
  if ( !tp.getInspectedValue(label,coupling_coefficient_ci)){
    cout << label << " is not set" << endl;
    exit(1);
  }
  label = base_label + "/coupling_coefficient_ic";
  if ( !tp.getInspectedValue(label,coupling_coefficient_ic)){
    cout << label << " is not set" << endl;
    exit(1);
  }
  label = base_label + "/coupling_coefficient_vi";
  if ( !tp.getInspectedValue(label,coupling_coefficient_vi)){
    cout << label << " is not set" << endl;
    exit(1);
  }
  label = base_label + "/coupling_coefficient_iv";
  if ( !tp.getInspectedValue(label,coupling_coefficient_iv)){
    cout << label << " is not set" << endl;
    exit(1);
  }
}