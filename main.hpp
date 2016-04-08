#include <iostream>
#include <math.h>
#include <cstdlib>
#include <fstream>
using namespace std;

// 0 - x coordinate, 1 - y coordinate
double **nodePos;
int **patchNode;
int **edgeInfo;
int **edgeLink;
int **edgeLinkS;
int *flag_bd;
int **coefMat;

int  **patch_edge;
double *thetaEdge;
double *thetaIn;

double **E_field_pos;
double **H_field_pos;
double *area;
