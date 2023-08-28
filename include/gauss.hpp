#ifndef _GAUSS_H
#define _GAUSS_H

class Gauss{
public:
    Gauss(){
        point[0][0] = -0.577350269; point[1][0] = -0.577350269; point[2][0]=-0.577350269; weight[0]=1.0;
        point[0][1] = 0.577350269; point[1][1] = -0.577350269; point[2][1]=-0.577350269; weight[1]=1.0;
        point[0][2] = 0.577350269; point[1][2] = 0.577350269; point[2][2]=-0.577350269; weight[2]=1.0;
        point[0][3] = -0.577350269; point[1][3] = 0.577350269; point[2][3]=-0.577350269; weight[3]=1.0;
        point[0][4] = -0.577350269; point[1][4] = -0.577350269; point[2][4]=0.577350269; weight[4]=1.0;
        point[0][5] = 0.577350269; point[1][5] = -0.577350269; point[2][5]=0.577350269; weight[5]=1.0;
        point[0][6] = 0.577350269; point[1][6] = 0.577350269; point[2][6]=0.577350269; weight[6]=1.0;
        point[0][7] = -0.577350269; point[1][7] = 0.577350269; point[2][7]=0.577350269; weight[7]=1.0;
    }
    double point[3][8],weight[8];
};

#endif