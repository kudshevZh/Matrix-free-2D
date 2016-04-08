#include"main.hpp"
void InstArry(int N_nodes, int N_patch){
//Coordinates of vertices:
// 1 - x coordinate, 2- y coordinate
    nodePos = new double*[N_nodes];
    for(int i = 1;i<=N_nodes;i++){
        nodePos[i] = new double [2];
        nodePos[i][1] = nodePos[i][2] = 0;
    }

// array of faces with three corresponding verticies
    patchNode = new int*[N_patch];
    for(int i = 1;i<=N_patch;i++){
        patchNode[i] = new int [3];
        patchNode[i][1] = patchNode[i][2] = patchNode[i][3] = 0;
    }
// edges infromation array:
// 1- corresponding patch number,
// 2- edges numeration (edge1: 2-3, edge2:3-1, edge3: 1-2)
// 3,4 - edge vertices
    int maxN = 3*N_patch;
    edgeInfo = new int *[maxN];
    edgeLink = new int *[maxN];

    for(int i = 1;i<=maxN;i++){
        edgeInfo[i] = new int [4];
        edgeLink[i] = new int [2];
        edgeInfo[i][1] = edgeInfo[i][2] = edgeInfo[i][3] = edgeInfo[i][4] = 0;
        edgeLink[i][1] = edgeLink[i][2] = 0;
    }

    patch_edge = new int *[N_patch];
    for(int i = 1;i<=N_patch;i++){
		patch_edge[i] = new int [3];
        patch_edge[i][1] = patch_edge[i][2] = patch_edge[i][3] = 0.0;
	}

	coefMat = new int *[N_patch];
    area = new double [N_patch];
    for(int i = 1;i<=N_patch;i++){
        coefMat[i] = new int[8+10+1];
        area[i] = 0.0;
    }
}

void InstArryAdd(int N_patch, int N_edge, int &Nb){

//1 - begin vertix, 2- end vertix
    edgeLinkS = new int *[N_edge];
    for(int i=1;i<=N_edge;i++){
        edgeLinkS[i] = new int[2];
        edgeLinkS[i][1] = edgeLinkS[i][2]  = 0;
    }
//Sort out edgeLink according to N_edge number
    int m=1;
    for(int i=1;i<=N_edge;i++){
        for(int j=m;j<=3*N_patch;j++){
            if((edgeLink[j][1]!=0)&&(edgeLink[j][2]!=0)){
                edgeLinkS[i][1] = edgeLink[j][1];
                edgeLinkS[i][2] = edgeLink[j][2];
                m=j+1;
                break;
            }
        }
    }
// Array of theta angle of the edge relativaly to x axis
    thetaEdge = new double [N_edge];
    thetaIn = new double [N_patch+1];

    for(int i=1;i<=N_patch;i++){
        for (int j=1;j<=3*N_patch;j++){
            if((edgeInfo[j][1]==i)&&(edgeInfo[j][2]==1)){
                thetaIn[i] = atan2((nodePos[edgeInfo[j][4]][2]- nodePos[edgeInfo[j][3]][2]),(nodePos[edgeInfo[j][4]][1]- nodePos[edgeInfo[j][3]][1]));
                if(thetaIn[i]<0){thetaIn[i] = 2*3.14-fabs(thetaIn[i]);}
            }
        }
    }

    for(int i=1;i<=N_edge;i++){
        thetaEdge[i] = atan2((nodePos[edgeLinkS[i][2]][2]- nodePos[edgeLinkS[i][1]][2]),(nodePos[edgeLinkS[i][2]][1]- nodePos[edgeLinkS[i][1]][1]));
        if(thetaEdge[i]<0){thetaEdge[i] = 2*3.14-fabs(thetaEdge[i]);}

        if(nodePos[edgeLink[i][2]][1]==nodePos[edgeLink[i][1]][1]){
            if(nodePos[edgeLink[i][2]][2]>nodePos[edgeLink[i][1]][2]){
                thetaEdge[i] = 1.57;  //90 degree
            }else if(nodePos[edgeLink[i][2]][2]<nodePos[edgeLink[i][1]][2]){
                thetaEdge[i] = 3*1.57;  //90 degree
            }

        }
        if(nodePos[edgeLink[i][2]][2]==nodePos[edgeLink[i][1]][2]){
            if(nodePos[edgeLink[i][2]][1]>nodePos[edgeLink[i][1]][1]){
                thetaEdge[i] = 0;  //90 degree
            }else if(nodePos[edgeLink[i][2]][1]<nodePos[edgeLink[i][1]][1]){
                thetaEdge[i] = 2*1.57;  //90 degree
            }
        }
    }

    int **flag_tmp;
    flag_tmp = new int *[N_edge];
    flag_bd = new int [N_edge];
    for (int i=1;i<=N_edge;i++){
        flag_tmp[i] = new int [2];
        flag_tmp[i][1] = flag_tmp[i][2] = 0;
    }

    for(int i =1; i<=N_patch;i++){
        for(int j=1;j<=3;j++){
            if(flag_tmp[patch_edge[i][j]][1] == 0){
                flag_tmp[patch_edge[i][j]][1] = i;
            } else{
                flag_tmp[patch_edge[i][j]][2] = i;
            }
        }
    }

    for(int i=1;i<=N_edge;i++){
        if(flag_tmp[i][2]==0){
            flag_bd[i] = 1;
            Nb++;
        }
   }


//Installation of electric field position array
    E_field_pos = new double *[2*N_edge+2*N_patch];
    for(int i=1;i<=2*N_edge+2*N_patch;i++){
        E_field_pos[i] = new double [2];
        E_field_pos[i][1] = E_field_pos[i][2] = 0.0;
    }

    H_field_pos = new double *[4*N_edge+4*N_patch];
    for(int i=1;i<=4*N_edge+4*N_patch;i++){
        H_field_pos[i] = new double [2];
        H_field_pos[i][1] = H_field_pos[i][2] = 0.0;
    }

}
int sizeOfFile (char* nameOfFile, int colNum){
   	ifstream idata;
	idata.open(nameOfFile, std::ios_base::in);
	int j=0;
    double a[colNum];
	while(!idata.eof()){
        for(int k=0;k<colNum;k++){
            idata>>a[k];
        }
		j++;
	}
	idata.close();
    return j-1;
}

void readNodePos(){
   	ifstream idata;
	idata.open("nodepos.dat", std::ios_base::in);
	int j=1;
    while(!idata.eof()){
        idata>>nodePos[j][1];
        idata>>nodePos[j][2];
        j++;
	}
	idata.close();
}

void readpatchNode(){
   	ifstream idata;
	idata.open("patch_node.dat", std::ios_base::in);
	int j=1;
	while(!idata.eof()){
        idata>>patchNode[j][1];
        idata>>patchNode[j][2];
        idata>>patchNode[j][3];
        j++;
	}
	idata.close();
}

void writeGrid(char* nameOfFileNodePos, char* nameOfFilePatch, int N_nodes, int N_patch){
    ofstream outNodePos, outPatch;
	outNodePos.open(nameOfFileNodePos);
	outPatch.open(nameOfFilePatch);

	for (int i=1;i<=N_nodes;i++){
       outNodePos<<" "<<nodePos[i][1]<<" "<<nodePos[i][2]<<endl;
	}

	for (int i=1;i<=N_patch;i++){
       outPatch<<" "<<patchNode[i][1]<<" "<<patchNode[i][2]<<" "<<patchNode[i][3]<<endl;
	}
}

void writeFiledPos(int N_edge,int N_patch){
    ofstream outE_field, outH_field;


    outE_field.open("E_field.dat");
    outH_field.open("H_field.dat");

    for(int i=1;i<=2*N_edge+2*N_patch;i++){
        outE_field<< E_field_pos[i][1]<<" "<<E_field_pos[i][2]<<endl;
    }

    for(int i=1;i<=4*N_edge+4*N_patch;i++){
        outH_field<< H_field_pos[i][1]<<" "<<H_field_pos[i][2]<<endl;
    }
}
int next(int j){
    if (j<3){
        return j+1;
    }else{
        return 1;
    }
 }

void CalcEFieldPos(int N_edge, int N_patch){
    for(int i=1; i<=N_edge;i++){
        E_field_pos[2*i-1][1] = nodePos[edgeLinkS[i][1]][1] + (nodePos[edgeLinkS[i][2]][1]-nodePos[edgeLinkS[i][1]][1])/3.0;
        E_field_pos[2*i-1][2] = nodePos[edgeLinkS[i][1]][2] + (nodePos[edgeLinkS[i][2]][2]-nodePos[edgeLinkS[i][1]][2])/3.0;

        E_field_pos[2*i][1] = nodePos[edgeLinkS[i][1]][1] + (nodePos[edgeLinkS[i][2]][1]-nodePos[edgeLinkS[i][1]][1])*2.0/3.0;
        E_field_pos[2*i][2] = nodePos[edgeLinkS[i][1]][2] + (nodePos[edgeLinkS[i][2]][2]-nodePos[edgeLinkS[i][1]][2])*2.0/3.0;

    }
    for(int i=1; i<=N_patch;i++){
        E_field_pos[2*N_edge+2*i-1][1] = (nodePos[patchNode[i][1]][1] + nodePos[patchNode[i][2]][1] + nodePos[patchNode[i][3]][1])/3.0;
        E_field_pos[2*N_edge+2*i][1] = E_field_pos[2*N_edge+2*i-1][1];

        E_field_pos[2*N_edge+2*i][2] = (nodePos[patchNode[i][1]][2] + nodePos[patchNode[i][2]][2] + nodePos[patchNode[i][3]][2])/3.0;
        E_field_pos[2*N_edge+2*i-1][2] = E_field_pos[2*N_edge+2*i][2];
    }
}

void CalcHfieldPos(int N_edge, int N_patch, double d){
    double COS = 0.0;
    double SIN = 0.0;

    for(int i=1; i<=N_edge;i++){

        COS = d*cos(thetaEdge[i]+1.57);
        SIN = d*sin(thetaEdge[i]+1.57);

        H_field_pos[4*i-3][1] = E_field_pos[2*i-1][1] + COS;
        H_field_pos[4*i-2][1] = E_field_pos[2*i-1][1] - COS;
        H_field_pos[4*i-1][1] = E_field_pos[2*i][1] + COS;
        H_field_pos[4*i][1] = E_field_pos[2*i][1] - COS;

        H_field_pos[4*i-3][2] = E_field_pos[2*i-1][2] + SIN;
        H_field_pos[4*i-2][2] = E_field_pos[2*i-1][2] - SIN;
        H_field_pos[4*i-1][2] = E_field_pos[2*i][2] + SIN;
        H_field_pos[4*i][2] = E_field_pos[2*i][2] - SIN;
    }

    for(int i=1; i<=N_patch;i++){
        H_field_pos[4*N_edge+4*i-3][1] = E_field_pos[2*N_edge+2*i-1][1] + d*cos(thetaIn[i]+1.57);
        H_field_pos[4*N_edge+4*i-2][1] = E_field_pos[2*N_edge+2*i-1][1] - d*cos(thetaIn[i]+1.57);
        H_field_pos[4*N_edge+4*i-1][1] = E_field_pos[2*N_edge+2*i][1] + d*cos(thetaIn[i]);
        H_field_pos[4*N_edge+4*i][1] = E_field_pos[2*N_edge+2*i][1] - d*cos(thetaIn[i]);


        H_field_pos[4*N_edge+4*i-3][2] = E_field_pos[2*N_edge+2*i-1][2] + d*sin(thetaIn[i]+1.57);
        H_field_pos[4*N_edge+4*i-2][2] = E_field_pos[2*N_edge+2*i-1][2] - d*sin(thetaIn[i]+1.57);
        H_field_pos[4*N_edge+4*i-1][2] = E_field_pos[2*N_edge+2*i][2] + d*sin(thetaIn[i]);
        H_field_pos[4*N_edge+4*i][2] = E_field_pos[2*N_edge+2*i][2] - d*sin(thetaIn[i]);
    }
}

void sortArry(int N_patch){
    int tempArry;

    for (int i=1;i<=3*N_patch;i++){
        for (int j=i+1;j<=3*N_patch;j++){
            if(edgeInfo[j][3]<edgeInfo[i][3]){
                for(int k=1;k<=4;k++){
                    tempArry = edgeInfo[j][k];
                    edgeInfo[j][k] = edgeInfo[i][k];
                    edgeInfo[i][k] = tempArry;
                }
            }
        }
    }

    for (int i=1;i<=3*N_patch;i++){
        for (int j=i+1;j<=3*N_patch;j++){
            if(edgeInfo[j][3]==edgeInfo[i][3]){
                if(edgeInfo[j][4]<=edgeInfo[i][4]){
                    for(int k=1;k<=4;k++){
                        tempArry = edgeInfo[j][k];
                        edgeInfo[j][k] = edgeInfo[i][k];
                        edgeInfo[i][k] = tempArry;
                    }
                }
            }
        }
    }
 }

bool PointInTriangle(int patchNumber, int HfieldNum){
        double a = 0.0;
        double b = 0.0;
        double c = 0.0;

        double x1 = nodePos[edgeLinkS[patch_edge[patchNumber][2]][2]][1];
        double y1 = nodePos[edgeLinkS[patch_edge[patchNumber][2]][2]][2];

        double x2 = nodePos[edgeLinkS[patch_edge[patchNumber][1]][1]][1];
        double y2 = nodePos[edgeLinkS[patch_edge[patchNumber][1]][1]][2];

        double x3 = nodePos[edgeLinkS[patch_edge[patchNumber][1]][2]][1];
        double y3 = nodePos[edgeLinkS[patch_edge[patchNumber][1]][2]][2];

        a = ((y2-y3)*(H_field_pos[HfieldNum][1]-x3)+(x3-x2)*(H_field_pos[HfieldNum][2]-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3));
        b = ((y3-y1)*(H_field_pos[HfieldNum][1]-x3)+(x1-x3)*(H_field_pos[HfieldNum][2]-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3));
        c = 1.0-a-b;

        if((a>0)&&(b>0)&&(c>0) ) {return true;} else { return false;}
}

int main(){
    int N_nodes = sizeOfFile("nodepos.dat",2);
    int N_patch = sizeOfFile("patch_node.dat",3);

    InstArry(N_nodes,N_patch);
    readNodePos();
    readpatchNode();

    int curnum = 0;
    for(int i=1;i<=N_patch;i++){
        for(int j=1;j<=3;j++){
            curnum = (i-1)*3+j;
            edgeInfo[curnum][1] = i;
            edgeInfo[curnum][2] = j;
            if (patchNode[i][next(j)]<patchNode[i][next(next(j))]){
                edgeInfo[curnum][3] = patchNode[i][next(j)];
                edgeInfo[curnum][4] = patchNode[i][next(next(j))];
            }else{
                edgeInfo[curnum][3] = patchNode[i][next(next(j))];
                edgeInfo[curnum][4] = patchNode[i][next(j)];
            }
        }
    }


    sortArry(N_patch);

    edgeLink[1][1] = edgeInfo[1][3];
    edgeLink[1][2] = edgeInfo[1][4];
    patch_edge[edgeInfo[1][1]][edgeInfo[1][2]] = 1;

    int N_edge = 1;

    for(int i=2;i<=3*N_patch;i++){
        if((edgeInfo[i][3]!=edgeInfo[i-1][3])||(edgeInfo[i][4]!=edgeInfo[i-1][4])){
            N_edge++;
            edgeLink[i][1] = edgeInfo[i][3];
            edgeLink[i][2] = edgeInfo[i][4];
            patch_edge[edgeInfo[i][1]][edgeInfo[i][2]] = N_edge;
        }else{
            patch_edge[edgeInfo[i][1]][edgeInfo[i][2]] = N_edge;
        }
    }
    int Nb = 0;
    InstArryAdd(N_patch,N_edge,Nb);

    CalcEFieldPos(N_edge,N_patch);
    CalcHfieldPos(N_edge,N_patch, 0.01);

//Output all the data (test)
    writeGrid("pos_nodes.dat","vertiic.dat",N_nodes,N_patch);
    writeFiledPos(N_edge,N_patch);

//form a general matrix of coefs
    for(int i = 1;i<=N_patch;i++){
 // all electric field global numbers corresponding to givern patch
        for(int j=1;j<=3;j++){
            coefMat[i][2*j-1] = 2*patch_edge[i][j]-1;
            coefMat[i][2*j] = 2*patch_edge[i][j];
        }
        coefMat[i][7] = 2*N_edge+2*i-1;
        coefMat[i][8] = 2*N_edge+2*i;

        for(int j=1;j<=3;j++){
            if (PointInTriangle(i, 4*patch_edge[i][j]-3)== true){
                coefMat[i][2*(j+4)-1] = 4*patch_edge[i][j]-3;
            } else{
                coefMat[i][2*(j+4)-1] = 4*patch_edge[i][j]-2;
            }

            if (PointInTriangle(i, 4*patch_edge[i][j]-1)== true){
                coefMat[i][2*(j+4)] = 4*patch_edge[i][j]-1;
            } else{
                coefMat[i][2*(j+4)] = 4*patch_edge[i][j];
            }
        }
        coefMat[i][15] = 4*N_edge+4*i-3;
        coefMat[i][16] = 4*N_edge+4*i-2;

        coefMat[i][17] = 4*N_edge+4*i-1;
        coefMat[i][18] = 4*N_edge+4*i;
// area of the triangle

        double x1 = nodePos[edgeLinkS[patch_edge[i][2]][2]][1];
        double y1 = nodePos[edgeLinkS[patch_edge[i][2]][2]][2];

        double x2 = nodePos[edgeLinkS[patch_edge[i][3]][2]][1];
		double y2 = nodePos[edgeLinkS[patch_edge[i][3]][2]][2];

		double x3 = nodePos[edgeLinkS[patch_edge[i][1]][2]][1];
		double y3 = nodePos[edgeLinkS[patch_edge[i][1]][2]][2];

		area[i] = fabs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
       //area[i] =x1;
    }

	ofstream out1;
	out1.open("Matrix.dat");
	for(int i = 1;i<=N_patch;i++){
		out1<<i<<" "<<coefMat[i][1]<<" "<<coefMat[i][2]<<" "<<coefMat[i][3]<<" "<<coefMat[i][4]<<" "<<coefMat[i][5]<<" "<<coefMat[i][6]<<" "<<coefMat[i][7]<<" "<<coefMat[i][8]<<" "<<coefMat[i][9]<<" "<<coefMat[i][10]<<" "<<coefMat[i][11]<<" "<<coefMat[i][12]<<" "<<coefMat[i][13]<<" "<<coefMat[i][14]<<" "<<coefMat[i][15]<<" "<<coefMat[i][16]<<" "<<coefMat[i][17]<<" "<<coefMat[i][18]<<" "<<area[i]<<endl;
		//out1<<i<<" "<<nodePos[edgeLinkS[patch_edge[i][2]][2]][2]<<" "<<nodePos[edgeLinkS[patch_edge[i][1]][1]][2]<<" "<<nodePos[edgeLinkS[patch_edge[i][1]][2]][2]<<endl;
	}

    return 0;
}
