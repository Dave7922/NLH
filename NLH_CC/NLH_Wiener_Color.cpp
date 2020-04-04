#include <math.h>
#include "mex.h"
#include <stdio.h>

// 声明输入变量
double *N1,*N2,*N3,*Ns,*Nstep,*dismap,*sigma;
float *gamma;
double *im_o1,*im_o2,*im_o3,*im_p1,*im_p2,*im_p3,*im_noise;
int    M,N;
// 声明输出变量
float *imr1,*imr2,*imr3;


void FindMin(float* OrgArray, int OrgNum, int* ResIndexArray, int ResNum)
{
    int i,j,k;
    float* ResArray = new float[ResNum];
    for ( i = 0; i < ResNum; i++ )  *(ResArray+i) = 10000000.0;

    for ( i = 0; i < OrgNum; i++ )
        {
            if ( *(OrgArray+i) >= *(ResArray+ResNum-1) )   continue;
            else if ( *(OrgArray+i) <= *ResArray )
                {
                    for ( j = ResNum-1; j > 0; j-- )
                        {
                            *(ResArray+j) = *(ResArray+j-1);
                            *(ResIndexArray+j) = *(ResIndexArray+j-1);
                        }
                   *ResArray = *(OrgArray+i);
                   *ResIndexArray = i;
                   continue;
                }
            else
                {
                    for ( j = 0; j < ResNum-1; j++ )
                        {
                            if ( *(ResArray+j) < *(OrgArray+i) && *(OrgArray+i) < *(ResArray+j+1) )
                                {
                                    for ( k = ResNum-1; k > j+1; k-- )
                                        {
                                            *(ResArray+k) = *(ResArray+k-1);
                                            *(ResIndexArray+k) = *(ResIndexArray+k-1);
                                        }
                                    *(ResArray+j+1) = *(OrgArray+i);
                                    *(ResIndexArray+j+1) = i;
                                    break;
                                }
                        }
                }
        }
    delete [] ResArray;
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //参数个数检查
    if ( nrhs != 15 ) 
		mexErrMsgTxt("15 input arguments required."); 
	if ( nlhs != 3 )
		mexErrMsgTxt("3 output arguments required.");
    // 定义输入变量
    N1       = mxGetPr(prhs[0]);
    N2       = mxGetPr(prhs[1]);
    N3       = mxGetPr(prhs[2]);
    Ns       = mxGetPr(prhs[3]);
    Nstep    = mxGetPr(prhs[4]);
    im_o1     = mxGetPr(prhs[5]);
    im_o2     = mxGetPr(prhs[6]);
    im_o3     = mxGetPr(prhs[7]);
    M        = mxGetM(prhs[5]);
    N        = mxGetN(prhs[5]);
    gamma    = (float*)mxGetPr(prhs[8]);
    sigma    = mxGetPr(prhs[9]);
    dismap     = mxGetPr(prhs[10]);
    im_p1     = mxGetPr(prhs[11]);
    im_p2     = mxGetPr(prhs[12]);
    im_p3     = mxGetPr(prhs[13]);
    im_noise = mxGetPr(prhs[14]);

  // 定义输出变量
    plhs[0] = mxCreateNumericMatrix(M,N,mxGetClassID(prhs[8]),mxREAL);
    imr1     = (float*)mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericMatrix(M,N,mxGetClassID(prhs[8]),mxREAL);
    imr2     = (float*)mxGetPr(plhs[1]);
    plhs[2] = mxCreateNumericMatrix(M,N,mxGetClassID(prhs[8]),mxREAL);
    imr3     = (float*)mxGetPr(plhs[2]);
    
    int N12 = int(*N1);
    int N22 = int(*N2);
    int N32 = int(*N3);
    int Ns1 = int(*Ns);
    int Nstep1  = int(*Nstep);
    
        
    int i,j,k,l,r;
    float weight,weight1,weight2;
   // 中间矩阵空间定义
    // im_ps，M行M列
    float** im_o1s;
    im_o1s = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_o1s[i] = new float[N+2*Ns1]; 
    // im_ps，M行M列
    float** im_o2s;
    im_o2s = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_o2s[i] = new float[N+2*Ns1]; 
    // im_ps，M行M列
    float** im_o3s;
    im_o3s = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_o3s[i] = new float[N+2*Ns1];
    // im_ps，M行M列
    float** im_p1s;
    im_p1s = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_p1s[i] = new float[N+2*Ns1];
    
    // im_ps，M行M列
    float** im_p2s;
    im_p2s = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_p2s[i] = new float[N+2*Ns1];
    // im_ps，M行M列
    float** im_p3s;
    im_p3s = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_p3s[i] = new float[N+2*Ns1];
   
    
    // im_noises，M行M列
    float** im_noises;
    im_noises = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_noises[i] = new float[N+2*Ns1];
    
      // im_ps，M行M列
    float** dp_ps;
    dp_ps = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   dp_ps[i] = new float[N+2*Ns1];
    
    
    // im_r，M行M列
    float** im_r1;
    im_r1 = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_r1[i] = new float[N+2*Ns1];
     // im_r，M行M列
    float** im_r2;
    im_r2 = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_r2[i] = new float[N+2*Ns1];
     // im_r，M行M列
    float** im_r3;
    im_r3 = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_r3[i] = new float[N+2*Ns1];
     // W，M行M列
    float** W;
    W = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   W[i] = new float[N+2*Ns1]; 
   
     // W，M行M列
    float** W1;
    W1 = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   W1[i] = new float[N+2*Ns1]; 
  
     // W，M行M列
    float** W2;
    W2 = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   W2[i] = new float[N+2*Ns1]; 
    
    // dis1，Ns*Ns-1行3列
    float** dis1;
    dis1 = new float*[(Ns1)*(Ns1)-1];
    for ( i = 0; i < (Ns1)*(Ns1)-1; i++ )   dis1[i] = new float[3];
           
               
    
	// dis1Array
	float* dis1Array = new float[(Ns1)*(Ns1)-1];
	
    
	// IndexArray
    int* IndexArrays = new int[N32];
    
     //d3
    int* dis5 = new int[N32];
    
     // Cood，N2行3列
   int** Cood;
   Cood = new int*[N22];
   for ( i = 0; i < N22; i++ )   Cood[i] = new int[3];
    
   
   
   
    float** groups;
    groups = new float*[N32];
    for ( i = 0; i < N32; i++ )       groups[i] = new float[N22];
    float** groups1;
    groups1 = new float*[N32];
    for ( i = 0; i < N32; i++ )       groups1[i] = new float[N22];
    float** groups2;
    groups2 = new float*[N32];
    for ( i = 0; i < N32; i++ )       groups2[i] = new float[N22];
    float** groupss;
    groupss = new float*[N32];
    for ( i = 0; i < N32; i++ )       groupss[i] = new float[N22];             
     float** groupss1;
    groupss1 = new float*[N32];
    for ( i = 0; i < N32; i++ )       groupss1[i] = new float[N22];  
     float** groupss2;
    groupss2 = new float*[N32];
    for ( i = 0; i < N32; i++ )       groupss2[i] = new float[N22];  
      
    float** groupss_t;
    groupss_t = new float*[N32];
    for ( i = 0; i < N32; i++ )   groupss_t[i] = new float[N22];
     float** groupss_t1;
    groupss_t1 = new float*[N32];
    for ( i = 0; i < N32; i++ )   groupss_t1[i] = new float[N22];
     float** groupss_t2;
    groupss_t2 = new float*[N32];
    for ( i = 0; i < N32; i++ )   groupss_t2[i] = new float[N22];
    float** groups_t;
    groups_t = new float*[N32];
    for ( i = 0; i < N32; i++ )   groups_t[i] = new float[N22];
    
     float** groups_t1;
    groups_t1 = new float*[N32];
    for ( i = 0; i < N32; i++ )   groups_t1[i] = new float[N22];
      float** groups_t2;
    groups_t2 = new float*[N32];
    for ( i = 0; i < N32; i++ )   groups_t2[i] = new float[N22];  
    
     // group，N1N12N22
    float*** group;
    group = new float**[N12];
    for ( i = 0; i < N12; i++ )
        {
            group[i] = new float*[N12];
            for ( j = 0; j < N12; j++ )
                group[i][j] = new float[N22];
        }
    // group，N1N12N22
    float*** group1;
    group1 = new float**[N12];
    for ( i = 0; i < N12; i++ )
        {
            group1[i] = new float*[N12];
            for ( j = 0; j < N12; j++ )
                group1[i][j] = new float[N22];
        }
   // group，N1N12N22
    float*** group2;
    group2 = new float**[N12];
    for ( i = 0; i < N12; i++ )
        {
            group2[i] = new float*[N12];
            for ( j = 0; j < N12; j++ )
                group2[i][j] = new float[N22];
        }
    
    
    float** groupb;
    groupb = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   groupb[i] = new float[N22];   
    
    float** groupb1;
    groupb1 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   groupb1[i] = new float[N22]; 
    
     float** groupb2;
    groupb2 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   groupb2[i] = new float[N22];   
    
    // im_n，N1行N1列
    float** im_n;
    im_n = new float*[N12];
    for ( i = 0; i < N12; i++ )   im_n[i] = new float[N12];
  
    // im_n，N1行N1列
    float** im_n1;
    im_n1 = new float*[N12];
    for ( i = 0; i < N12; i++ )   im_n1[i] = new float[N12];
    
    // im_n，N1行N1列
    float** im_n2;
    im_n2 = new float*[N12];
    for ( i = 0; i < N12; i++ )   im_n2[i] = new float[N12];       
    
    
     float** DD;
    DD = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   DD[i] = new float[N22];
    float** DD1;
    DD1 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   DD1[i] = new float[N22];
    float** DD2;
    DD2 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   DD2[i] = new float[N22];
    float** WW;
    WW = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   WW[i] = new float[N22];  
    float** WW1;
    WW1 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   WW1[i] = new float[N22];  
    float** WW2;
    WW2 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   WW2[i] = new float[N22];  
    
    float** groupsh;
    groupsh = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )       groupsh[i] = new float[N22];
    float** groupsh1;
    groupsh1 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )       groupsh1[i] = new float[N22];
    float** groupsh2;
    groupsh2 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )       groupsh2[i] = new float[N22];
    float** groupsh_t;
    groupsh_t = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   groupsh_t[i] = new float[N22];    
    float** groupsh_t1;
    groupsh_t1 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   groupsh_t1[i] = new float[N22]; 
    float** groupsh_t2;
    groupsh_t2 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   groupsh_t2[i] = new float[N22]; 
      
        
    
    float** groupsw;
    groupsw = new float*[N32];
    for ( i = 0; i < N32; i++ )       groupsw[i] = new float[N22];
    
    float** groupsw1;
    groupsw1 = new float*[N32];
    for ( i = 0; i < N32; i++ )       groupsw1[i] = new float[N22];
    
    float** groupsw2;
    groupsw2 = new float*[N32];
    for ( i = 0; i < N32; i++ )       groupsw2[i] = new float[N22];
    
    
    float** groupssw;
    groupssw = new float*[N32];
    for ( i = 0; i < N32; i++ )       groupssw[i] = new float[N22];
    
    float** groupssw1;
    groupssw1 = new float*[N32];
    for ( i = 0; i < N32; i++ )       groupssw1[i] = new float[N22];
    
    float** groupssw2;
    groupssw2 = new float*[N32];
    for ( i = 0; i < N32; i++ )       groupssw2[i] = new float[N22];
    
    
     float** groupss_tw;
    groupss_tw = new float*[N32];
    for ( i = 0; i < N32; i++ )   groupss_tw[i] = new float[N22];
    
     float** groupss_tw1;
    groupss_tw1 = new float*[N32];
    for ( i = 0; i < N32; i++ )   groupss_tw1[i] = new float[N22];
     float** groupss_tw2;
    groupss_tw2 = new float*[N32];
    for ( i = 0; i < N32; i++ )   groupss_tw2[i] = new float[N22];
    float** groups_tw;
    groups_tw = new float*[N32];
    for ( i = 0; i < N32; i++ )   groups_tw[i] = new float[N22]; 
    float** groups_tw1;
    groups_tw1 = new float*[N32];
    for ( i = 0; i < N32; i++ )   groups_tw1[i] = new float[N22]; 
      float** groups_tw2;
    groups_tw2 = new float*[N32];
    for ( i = 0; i < N32; i++ )   groups_tw2[i] = new float[N22]; 
     // dis2，N2-1行3列
    float** dis2;
    dis2 = new float*[N22-1];
    for ( i = 0; i < N22-1; i++ )   dis2[i] = new float[3];
    
    // IndexArray
    int* IndexArray = new int[N22-1];
    
    float** W_wiener;
    W_wiener = new float*[N32];
    for ( i = 0; i < N32; i++ )       W_wiener[i] = new float[N22]; 
    float** W_wiener1;
    W_wiener1 = new float*[N32];
    for ( i = 0; i < N32; i++ )       W_wiener1[i] = new float[N22];
    float** W_wiener2;
    W_wiener2 = new float*[N32];
    for ( i = 0; i < N32; i++ )       W_wiener2[i] = new float[N22];
    
       
    // group，N1N12N22
    float*** groupw;
    groupw = new float**[N12];
    for ( i = 0; i < N12; i++ )
        {
            groupw[i] = new float*[N12];
            for ( j = 0; j < N12; j++ )
                groupw[i][j] = new float[N22];
        }
   // group，N1N12N22
    float*** groupw1;
    groupw1 = new float**[N12];
    for ( i = 0; i < N12; i++ )
        {
            groupw1[i] = new float*[N12];
            for ( j = 0; j < N12; j++ )
                groupw1[i][j] = new float[N22];
        }
    // group，N1N12N22
    float*** groupw2;
    groupw2 = new float**[N12];
    for ( i = 0; i < N12; i++ )
        {
            groupw2[i] = new float*[N12];
            for ( j = 0; j < N12; j++ )
                groupw2[i][j] = new float[N22];
        }
      
   
    float** groupbw;
    groupbw = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   groupbw[i] = new float[N22];   
    
     float** groupbw1;
    groupbw1 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   groupbw1[i] = new float[N22]; 
     float** groupbw2;
    groupbw2 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   groupbw2[i] = new float[N22]; 
    float** groupshw;
    groupshw = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )       groupshw[i] = new float[N22];
     float** groupshw1;
    groupshw1 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )       groupshw1[i] = new float[N22];
     float** groupshw2;
    groupshw2 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )       groupshw2[i] = new float[N22];
    
    float** groupsh_tw;
    groupsh_tw = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   groupsh_tw[i] = new float[N22];  
     float** groupsh_tw1;
    groupsh_tw1 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   groupsh_tw1[i] = new float[N22];  
     float** groupsh_tw2;
    groupsh_tw2 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   groupsh_tw2[i] = new float[N22];  
    
     // im_nn，N1行N1列
    float** im_nn;
    im_nn = new float*[N12];
    for ( i = 0; i < N12; i++ )   im_nn[i] = new float[N12];
    
    // im_n，N1行N1列
    float** im_nw;
    im_nw = new float*[N12];
    for ( i = 0; i < N12; i++ )   im_nw[i] = new float[N12];
    // im_n，N1行N1列
    float** im_nw1;
    im_nw1 = new float*[N12];
    for ( i = 0; i < N12; i++ )   im_nw1[i] = new float[N12];
    // im_n，N1行N1列
    float** im_nw2;
    im_nw2 = new float*[N12];
    for ( i = 0; i < N12; i++ )   im_nw2[i] = new float[N12];
    
    float** dis3;
    dis3 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   dis3[i] = new float[N12*N12];
    // dis4，N1*N1行2列
    float* dis4 = new float[N12*N12];
    // dis1Array
	float* dis1Arrays = new float[N12*N12];
    
   
    
     
    // x,y
    int x[409600],y[409600];
    // 输出变量初始化
	for ( i = 0; i < M*N; i++ )
		{
			*(imr1+i) = 0.0;
            *(imr2+i) = 0.0;
            *(imr3+i) = 0.0;
		}
            
    
   for ( i = 0; i < M+2*Ns1; i++ )
          for ( j = 0; j < N+2*Ns1; j++ )
              {
              im_p1s[i][j] = 0.0;
              im_p2s[i][j] = 0.0;
              im_p3s[i][j] = 0.0;
              im_o1s[i][j] = 0.0;
              im_o2s[i][j] = 0.0;
              im_o3s[i][j] = 0.0;
              im_r1[i][j]  = 0.0;
              im_r2[i][j]  = 0.0;
              im_r3[i][j]  = 0.0;
              dp_ps[i][j] = 0.0;
              W[i][j]     = 0.0;
              W1[i][j]     = 0.0;
              W2[i][j]     = 0.0;
              im_noises[i][j] = 0.0;
              }
   
       for ( i = Ns1; i < M+Ns1; i++ )
          for ( j = Ns1; j < N+Ns1; j++ )
              {
              im_p1s[i][j] = *(im_p1+(j-Ns1)*M+(i-Ns1));
              im_p2s[i][j] = *(im_p2+(j-Ns1)*M+(i-Ns1));
              im_p3s[i][j] = *(im_p3+(j-Ns1)*M+(i-Ns1));
              im_o1s[i][j] = *(im_o1+(j-Ns1)*M+(i-Ns1));
              im_o2s[i][j] = *(im_o2+(j-Ns1)*M+(i-Ns1));
              im_o3s[i][j] = *(im_o3+(j-Ns1)*M+(i-Ns1));
              dp_ps[i][j]  = *(dismap+(j-Ns1)*M+(i-Ns1));
              im_noises[i][j] = *(im_noise+(j-Ns1)*M+(i-Ns1));
              }
 
  // 周期延拓
    for ( i = 0; i < Ns1; i++ )
          for ( j = 0; j < N+2*Ns1; j++ )
              {
              im_p1s[i][j] = im_p1s[i+Ns1][j];
              im_p1s[i+M+Ns1][j] = im_p1s[i+M][j];
              im_p2s[i][j] = im_p2s[i+Ns1][j];
              im_p2s[i+M+Ns1][j] = im_p2s[i+M][j];
              im_p3s[i][j] = im_p3s[i+Ns1][j];
              im_p3s[i+M+Ns1][j] = im_p3s[i+M][j];
              im_o1s[i][j] = im_o1s[i+Ns1][j];
              im_o1s[i+M+Ns1][j] = im_o1s[i+M][j];
              im_o2s[i][j] = im_o2s[i+Ns1][j];
              im_o2s[i+M+Ns1][j] = im_o2s[i+M][j];
              im_o3s[i][j] = im_o3s[i+Ns1][j];
              im_o3s[i+M+Ns1][j] = im_o3s[i+M][j];
              dp_ps[i][j] = dp_ps[i+Ns1][j];
              dp_ps[i+M+Ns1][j] = dp_ps[i+M][j];
              im_noises[i][j] = im_noises[i+Ns1][j];
              im_noises[i+M+Ns1][j] = im_noises[i+M][j];
              }
    for ( i = 0; i < M+2*Ns1; i++ )
          for ( j = 0; j < Ns1; j++ )
              {
              im_p1s[i][j] = im_p1s[i][j+Ns1];
              im_p1s[i][j+N+Ns1] = im_p1s[i][j+N];
              im_p2s[i][j] = im_p2s[i][j+Ns1];
              im_p2s[i][j+N+Ns1] = im_p2s[i][j+N];
              im_p3s[i][j] = im_p3s[i][j+Ns1];
              im_p3s[i][j+N+Ns1] = im_p3s[i][j+N];
              im_o1s[i][j] = im_o1s[i][j+Ns1];
              im_o1s[i][j+N+Ns1] = im_o1s[i][j+N];
              im_o2s[i][j] = im_o2s[i][j+Ns1];
              im_o2s[i][j+N+Ns1] = im_o2s[i][j+N];
              im_o3s[i][j] = im_o3s[i][j+Ns1];
              im_o3s[i][j+N+Ns1] = im_o3s[i][j+N];
              dp_ps[i][j] = dp_ps[i][j+Ns1];
              dp_ps[i][j+N+Ns1] = dp_ps[i][j+N];
              im_noises[i][j] = im_noises[i][j+Ns1];
              im_noises[i][j+N+Ns1] = im_noises[i][j+N];
              }
   

    
                         float sqrt_N2_1 =  1.0/sqrt(32.0);
                         float sqrt_N2   =  1.0/sqrt(16.0);
                         float sqrt_N2_2  = 1.0/sqrt(8.0);  
                         float sqrt_N2_4  = 1.0/sqrt(4.0);
                         float sqrt_N2_8  = 1.0/sqrt(2.0);
                         float sqrt_N2_3  = 1.0/sqrt(64.0);
                         float sqrt_N2_5  = 1.0/sqrt(128.0);
                         float THR        = (*sigma)*(*sigma);
                         
                         
                      
                         float weightb,weightb1,weightb2, meand;
                         int NN, N_num;
                         
                         float dthr = (*sigma*255)*(*sigma*255)*0.000174;
                        
                         float group_tm = 0.0;
     
                         
	// 算法主体
    for ( int a = Ns1; a <= M+Ns1; a += Nstep1 )
        for ( int b = Ns1; b <= N+Ns1; b += Nstep1 )
            {
  ////////////////////////////////////////////////////////////          	
      float nnn = 0.0;
            
           for ( i = a-1; i < a+N12-1; i++ )
                    for ( j = b-1; j < b+N12-1; j++ )
                    {
                        if ( im_r1[i][j] != 0.0 )
                          nnn += 1.0;
                    }
            
              if ( nnn == N12*N12)
                  continue;
      
   ////////////////////////////////////////////////////////////   
     weightb = 0.0;           
     weightb1 = 0.0;
     weightb2 = 0.0;
     meand = 0.0;

     N_num = 0.0;
    
     NN = 8.0;
    
    
    for ( i = a-1; i < a + NN-1; i++ )
                    for ( j = b-1; j < b + NN-1; j++ )
                               meand += dp_ps[i][j];
  
                               meand /= 64.00000000;

         if ( meand <= dthr*1.1 )
           {
           N22 = 64;
           N_num = 64;
           }
         else if ( meand > dthr*1.1 && meand <= dthr*1.5 )
           {
           N22 = 64;
           N_num = 32;
           }
         else
           {
            N22 = 16;
            N_num = 8;
           }   
 
 ///////////////////////////////////////////////////////////               
                
                for ( i = 0; i < N22; i++ )
                    for ( j = 0; j < 3; j++ )
                        Cood[i][j] = 0.0;  
                for ( i = 0; i < N12; i++ )
                    for ( j = 0; j < N12; j++ )
                        for ( k = 0; k < N22; k++ )
                               {
                            group[i][j][k] = 0.0;
                            groupw[i][j][k] = 0.0;
                             group1[i][j][k] = 0.0;
                            groupw1[i][j][k] = 0.0;
                             group2[i][j][k] = 0.0;
                            groupw2[i][j][k] = 0.0;
                             }
                for ( i = a; i <= a+N12-1; i++ )
                    for ( j = b; j <= b+N12-1; j++ )
                        {
                            im_n[i-a][j-b] = im_o1s[i-1][j-1];
                            group[i-a][j-b][0] = im_n[i-a][j-b];
                            im_nw[i-a][j-b] = im_p1s[i-1][j-1];
                            groupw[i-a][j-b][0] = im_nw[i-a][j-b];
                            im_n1[i-a][j-b] = im_o2s[i-1][j-1];
                            group1[i-a][j-b][0] = im_n1[i-a][j-b];
                            im_nw1[i-a][j-b] = im_p2s[i-1][j-1];
                            groupw1[i-a][j-b][0] = im_nw1[i-a][j-b]; 
                            im_n2[i-a][j-b] = im_o3s[i-1][j-1];
                            group2[i-a][j-b][0] = im_n2[i-a][j-b];
                            im_nw2[i-a][j-b] = im_p3s[i-1][j-1];
                            groupw2[i-a][j-b][0] = im_nw2[i-a][j-b];
                            im_nn[i-a][j-b] = im_noises[i-1][j-1];
                         }
                                             
                Cood[0][0] = 0.0;
                Cood[0][1] = a-1;
                Cood[0][2] = b-1;
                int Q = 0;
                for ( k = 1; k <= (Ns1-1)/2; k++ )   // 半径
                    {
                        for ( i = -k; i <= k-1; i++ )
                            {
                                x[Q] = a-k-1;
                                y[Q] = b+i-1;
                                Q++;
                            }
                        for ( i = -k; i <= k-1; i++ )
                            {
                                x[Q] = a+i-1;
                                y[Q] = b+k-1;
                                Q++;
                            }
                        for ( i = -k; i <= k-1; i++ )
                            {
                                x[Q] = a+k-1;
                                y[Q] = b-i-1;
                                Q++;
                            }
                        for ( i = -k; i <= k-1; i++ )
                            {
                                x[Q] = a-i-1;
                                y[Q] = b-k-1;
                                Q++;
                            }
                    }
               
                 int xk,yk; 
                for ( k = 0; k < Q; k++ )
                    {
                     float  sum  = 0.0,minu = 0.0; 
                          xk = x[k];
                          yk = y[k];
                       for ( i = xk; i < xk+N12; i++ )
                            for ( j = yk; j < yk+N12; j++ )
                          {
                              minu = im_noises[i][j]-im_nn[i-xk][j-yk];
                                 sum += minu*minu;
                            }
                        
                        dis1[k][0] = xk;
                        dis1[k][1] = yk;
                        dis1[k][2] = sum;
                    }
                
                
                for ( i = 0; i < (Ns1)*(Ns1)-1; i++ )   dis1Array[i] = dis1[i][2];
                FindMin(dis1Array,(Ns1)*(Ns1)-1,IndexArray,N22-1);
               
                
                for ( i = 0; i < N22-1; i++ )
                    {
                        dis2[i][0] = dis1[IndexArray[i]][0];
                        dis2[i][1] = dis1[IndexArray[i]][1];
                        dis2[i][2] = dis1[IndexArray[i]][2];
                    }
                
       
                int x,y;                
                for ( k = 0; k < N22-1; k++ )
                    {    
                         x = dis2[k][0];
						 y = dis2[k][1];
                        for ( i = x; i <= x+N12-1; i++ )
                             for ( j = y; j <= y+N12-1; j++ )
                                {
							     group[i-x][j-y][k+1] = im_o1s[i][j];
                                 groupw[i-x][j-y][k+1] = im_p1s[i][j];
                                 group1[i-x][j-y][k+1] = im_o2s[i][j];
                                 groupw1[i-x][j-y][k+1] = im_p2s[i][j];
                                 group2[i-x][j-y][k+1] = im_o3s[i][j];
                                 groupw2[i-x][j-y][k+1] = im_p3s[i][j];
                                }
                        Cood[k+1][0] = k+1;
                        Cood[k+1][1] = x;
                        Cood[k+1][2] = y;
                    }
                
 
                
 //将每组中的所有块按列扫描给groupb         
                            
                      for ( k = 0; k < N12; k++ )
                           for ( j = 0; j < N12; j++ )
                                for ( i = 0; i < N22; i++ )
                                    {
                                      groupb[k*N12+j][i] = group[k][j][i];
                                      groupbw[k*N12+j][i] = groupw[k][j][i];
                                      groupb1[k*N12+j][i] = group1[k][j][i];
                                      groupbw1[k*N12+j][i] = groupw1[k][j][i];
                                      groupb2[k*N12+j][i] = group2[k][j][i];
                                      groupbw2[k*N12+j][i] = groupw2[k][j][i];
                                     }
                               
                             
        for (k = 0; k < N12*N12; k++)
               for ( i = 0; i < N12*N12; i++) 
                     dis3[i][k] = 0.0;  
                
       
               for ( k = 0; k < N12*N12; k++)
                      for ( i = 0; i < N12*N12; i++)
                            { 
                            float sum = 0.0,minu = 0.0;
                              for ( j = 0; j < N22; j++)
                            {
                                minu = groupbw[k][j]-groupbw[i][j];
                                sum += minu*minu;
                              }
                             dis3[k][i] = sum;
                           }
     //      printf("%d,%d,%d,%d\n",dis3[0][0],dis3[0][1],dis3[5][9],dis3[10][11]); 
             
   
       // Haar正变换       
     for ( j = 0; j < N12*N12; j++ )
         for ( i = 0; i < N22; i++ )
        {
           groupsh[j][i] = 0.0;
           groupshw[j][i] = 0.0;
           groupsh_t[j][i] = 0.0;
           groupsh_tw[j][i] = 0.0;
           groupsh1[j][i] = 0.0;
           groupshw1[j][i] = 0.0;
           groupsh_t1[j][i] = 0.0;
           groupsh_tw1[j][i] = 0.0;
           groupsh2[j][i] = 0.0;
           groupshw2[j][i] = 0.0;
           groupsh_t2[j][i] = 0.0;
           groupsh_tw2[j][i] = 0.0;
          } 
   
    
     for ( j = 0; j < N12*N12; j++ )
         for ( i = 0; i < N22; i++ )
        {
           groupsh[j][i] = groupb[j][i];
           groupshw[j][i] = groupbw[j][i];
           groupsh1[j][i] = groupb1[j][i];
           groupshw1[j][i] = groupbw1[j][i];
           groupsh2[j][i] = groupb2[j][i];
           groupshw2[j][i] = groupbw2[j][i];
         }
                
               
                
                          if (N22 == 16.0)
                             {
                          for ( j = 0; j < N12*N12; j++ )
                                
                                         {  
                          for ( i = 0; i < N22; i++ )
                             
                             groupsh_t[j][0] += groupsh[j][i];
                             groupsh_t[j][0] = sqrt_N2*groupsh_t[j][0];
                                                         
                         for ( i = 0; i < 8; i++ )

                             groupsh_t[j][1] += groupsh[j][i];
                             
                             
                          for ( i = 8; i < 16; i++ )

                              groupsh_t[j][1] -= groupsh[j][i];
                              groupsh_t[j][1] *= sqrt_N2;
                          for ( i = 0; i < 4; i++ )

                              groupsh_t[j][2] += groupsh[j][i];
                             
                          for ( i = 4; i < 8; i++ )

                              groupsh_t[j][2] -= groupsh[j][i];
                             groupsh_t[j][2] *= sqrt_N2_2;
                          for ( i = 8; i < 12; i++ )

                              groupsh_t[j][3] += groupsh[j][i];
                          for ( i = 12; i < 16; i++ )

                              groupsh_t[j][3] -= groupsh[j][i];
                              groupsh_t[j][3] *= sqrt_N2_2;
                          
                         for ( i = 0; i < 2; i++ )

                              groupsh_t[j][4] += groupsh[j][i];
                          for ( i = 2; i < 4; i++ )

                              groupsh_t[j][4] -= groupsh[j][i];
                              groupsh_t[j][4] *= sqrt_N2_4;
                              
                          for ( i = 4; i < 6; i++ )

                              groupsh_t[j][5] += groupsh[j][i];
                          for ( i = 6; i < 8; i++ )

                              groupsh_t[j][5] -= groupsh[j][i];
                              groupsh_t[j][5] *= sqrt_N2_4;
                           for ( i = 8; i < 10; i++ )

                              groupsh_t[j][6] += groupsh[j][i]; 
                           for ( i = 10; i < 12; i++ )

                              groupsh_t[j][6] -= groupsh[j][i]; 
                              groupsh_t[j][6] *= sqrt_N2_4;
                           for ( i = 12; i < 14; i++ )

                              groupsh_t[j][7] += groupsh[j][i]; 
                           for ( i = 14; i < 16; i++ ) 
 
                              groupsh_t[j][7] -= groupsh[j][i];  
                              groupsh_t[j][7] *= sqrt_N2_4;
                          for ( l = 0; l < 8; l++ )
                              groupsh_t[j][l+8] = sqrt_N2_8*(groupsh[j][2*l]-groupsh[j][2*l+1]);
                              
                                     }
                              }
                
                
                 if (N22 == 16.0)
                             {
                          for ( j = 0; j < N12*N12; j++ )
                                
                                         {  
                          for ( i = 0; i < N22; i++ )
                             
                             groupsh_t1[j][0] += groupsh1[j][i];
                             groupsh_t1[j][0] = sqrt_N2*groupsh_t1[j][0];
                                                         
                         for ( i = 0; i < 8; i++ )

                             groupsh_t1[j][1] += groupsh1[j][i];
                             
                             
                          for ( i = 8; i < 16; i++ )

                              groupsh_t1[j][1] -= groupsh1[j][i];
                              groupsh_t1[j][1] *= sqrt_N2;
                          for ( i = 0; i < 4; i++ )

                              groupsh_t1[j][2] += groupsh1[j][i];
                             
                          for ( i = 4; i < 8; i++ )

                              groupsh_t1[j][2] -= groupsh1[j][i];
                             groupsh_t1[j][2] *= sqrt_N2_2;
                          for ( i = 8; i < 12; i++ )

                              groupsh_t1[j][3] += groupsh1[j][i];
                          for ( i = 12; i < 16; i++ )

                              groupsh_t1[j][3] -= groupsh1[j][i];
                              groupsh_t1[j][3] *= sqrt_N2_2;
                          
                         for ( i = 0; i < 2; i++ )

                              groupsh_t1[j][4] += groupsh1[j][i];
                          for ( i = 2; i < 4; i++ )

                              groupsh_t1[j][4] -= groupsh1[j][i];
                              groupsh_t1[j][4] *= sqrt_N2_4;
                              
                          for ( i = 4; i < 6; i++ )

                              groupsh_t1[j][5] += groupsh1[j][i];
                          for ( i = 6; i < 8; i++ )

                              groupsh_t1[j][5] -= groupsh1[j][i];
                              groupsh_t1[j][5] *= sqrt_N2_4;
                           for ( i = 8; i < 10; i++ )

                              groupsh_t1[j][6] += groupsh1[j][i]; 
                           for ( i = 10; i < 12; i++ )

                              groupsh_t1[j][6] -= groupsh1[j][i]; 
                              groupsh_t1[j][6] *= sqrt_N2_4;
                           for ( i = 12; i < 14; i++ )

                              groupsh_t1[j][7] += groupsh1[j][i]; 
                           for ( i = 14; i < 16; i++ ) 
 
                              groupsh_t1[j][7] -= groupsh1[j][i];  
                              groupsh_t1[j][7] *= sqrt_N2_4;
                          for ( l = 0; l < 8; l++ )
                              groupsh_t1[j][l+8] = sqrt_N2_8*(groupsh1[j][2*l]-groupsh1[j][2*l+1]);
                              
                                     }
                              }
                
                         if (N22 == 16.0)
                             {
                          for ( j = 0; j < N12*N12; j++ )
                                
                                         {  
                          for ( i = 0; i < N22; i++ )
                             
                             groupsh_t2[j][0] += groupsh2[j][i];
                             groupsh_t2[j][0] = sqrt_N2*groupsh_t2[j][0];
                                                         
                         for ( i = 0; i < 8; i++ )

                             groupsh_t2[j][1] += groupsh2[j][i];
                             
                             
                          for ( i = 8; i < 16; i++ )

                              groupsh_t2[j][1] -= groupsh2[j][i];
                              groupsh_t2[j][1] *= sqrt_N2;
                          for ( i = 0; i < 4; i++ )

                              groupsh_t2[j][2] += groupsh2[j][i];
                             
                          for ( i = 4; i < 8; i++ )

                              groupsh_t2[j][2] -= groupsh2[j][i];
                             groupsh_t2[j][2] *= sqrt_N2_2;
                          for ( i = 8; i < 12; i++ )

                              groupsh_t2[j][3] += groupsh2[j][i];
                          for ( i = 12; i < 16; i++ )

                              groupsh_t2[j][3] -= groupsh2[j][i];
                              groupsh_t2[j][3] *= sqrt_N2_2;
                          
                         for ( i = 0; i < 2; i++ )

                              groupsh_t2[j][4] += groupsh2[j][i];
                          for ( i = 2; i < 4; i++ )

                              groupsh_t2[j][4] -= groupsh2[j][i];
                              groupsh_t2[j][4] *= sqrt_N2_4;
                              
                          for ( i = 4; i < 6; i++ )

                              groupsh_t2[j][5] += groupsh2[j][i];
                          for ( i = 6; i < 8; i++ )

                              groupsh_t2[j][5] -= groupsh2[j][i];
                              groupsh_t2[j][5] *= sqrt_N2_4;
                           for ( i = 8; i < 10; i++ )

                              groupsh_t2[j][6] += groupsh2[j][i]; 
                           for ( i = 10; i < 12; i++ )

                              groupsh_t2[j][6] -= groupsh2[j][i]; 
                              groupsh_t2[j][6] *= sqrt_N2_4;
                           for ( i = 12; i < 14; i++ )

                              groupsh_t2[j][7] += groupsh2[j][i]; 
                           for ( i = 14; i < 16; i++ ) 
 
                              groupsh_t2[j][7] -= groupsh2[j][i];  
                              groupsh_t2[j][7] *= sqrt_N2_4;
                          for ( l = 0; l < 8; l++ )
                              groupsh_t2[j][l+8] = sqrt_N2_8*(groupsh2[j][2*l]-groupsh2[j][2*l+1]);
                              
                                     }
                              }
                
                
                
                        if (N22 == 64.0)
                             {
                          for ( j = 0; j < N12*N12; j++ )
                              
                                   {   
                         for ( i = 0; i < N22; i++ )
                            groupsh_t[j][0] += (sqrt_N2_3)*groupsh[j][i];
                          for ( i = 0; i < 32; i++ )
                            groupsh_t[j][1] += (sqrt_N2_3)*groupsh[j][i];
                          for ( i = 32; i < N22; i++ )
                             groupsh_t[j][1] -= (sqrt_N2_3)*groupsh[j][i];
                          for ( i = 0; i < 16; i++ )
                             groupsh_t[j][2] += (sqrt_N2_1)*groupsh[j][i];
                          for ( i = 16; i < 32; i++ )
                             groupsh_t[j][2] -= (sqrt_N2_1)*groupsh[j][i];
                          for ( i = 32; i < 48; i++ )
                             groupsh_t[j][3] += (sqrt_N2_1)*groupsh[j][i];
                          for ( i = 48; i < N22; i++ )
                             groupsh_t[j][3] -= (sqrt_N2_1)*groupsh[j][i];
                          for ( i = 0; i < 8; i++ )
                             groupsh_t[j][4] += (sqrt_N2)*groupsh[j][i];
                          for ( i = 8; i < 16; i++ )
                             groupsh_t[j][4] -= (sqrt_N2)*groupsh[j][i];
                          for ( i = 16; i < 24; i++ )
                             groupsh_t[j][5] += (sqrt_N2)*groupsh[j][i];
                          for ( i = 24; i < 32; i++ )
                             groupsh_t[j][5] -= (sqrt_N2)*groupsh[j][i];
                          for ( i = 32; i < 40; i++ )
                             groupsh_t[j][6] += (sqrt_N2)*groupsh[j][i]; 
                          for ( i = 40; i < 48; i++ )
                             groupsh_t[j][6] -= (sqrt_N2)*groupsh[j][i]; 
                          for ( i = 48; i < 56; i++ )
                             groupsh_t[j][7] += (sqrt_N2)*groupsh[j][i]; 
                          for ( i = 56; i < N22; i++ )
                             groupsh_t[j][7] -= (sqrt_N2)*groupsh[j][i];
                          for ( i = 0; i < 4; i++ )
                             groupsh_t[j][8] += (sqrt_N2_2)*groupsh[j][i];
                          for ( i = 4; i < 8; i++ )
                             groupsh_t[j][8] -= (sqrt_N2_2)*groupsh[j][i];
                          for ( i = 8; i < 12; i++ )
                             groupsh_t[j][9] += (sqrt_N2_2)*groupsh[j][i];
                          for ( i = 12; i < 16; i++ )
                             groupsh_t[j][9] -= (sqrt_N2_2)*groupsh[j][i];
                          for ( i = 16; i < 20; i++ )
                             groupsh_t[j][10] += (sqrt_N2_2)*groupsh[j][i]; 
                          for ( i = 20; i < 24; i++ )
                             groupsh_t[j][10] -= (sqrt_N2_2)*groupsh[j][i]; 
                          for ( i = 24; i < 28; i++ )
                             groupsh_t[j][11] += (sqrt_N2_2)*groupsh[j][i]; 
                          for ( i = 28; i < 32; i++ )
                             groupsh_t[j][11] -= (sqrt_N2_2)*groupsh[j][i];
                          for ( i = 32; i < 36; i++ )
                             groupsh_t[j][12] += (sqrt_N2_2)*groupsh[j][i];
                          for ( i = 36; i < 40; i++ )
                             groupsh_t[j][12] -= (sqrt_N2_2)*groupsh[j][i];
                          for ( i = 40; i < 44; i++ )
                             groupsh_t[j][13] += (sqrt_N2_2)*groupsh[j][i];
                          for ( i = 44; i < 48; i++ )
                             groupsh_t[j][13] -= (sqrt_N2_2)*groupsh[j][i];
                          for ( i = 48; i < 52; i++ )
                             groupsh_t[j][14] += (sqrt_N2_2)*groupsh[j][i]; 
                          for ( i = 52; i < 56; i++ )
                             groupsh_t[j][14] -= (sqrt_N2_2)*groupsh[j][i]; 
                          for ( i = 56; i < 60; i++ )
                             groupsh_t[j][15] += (sqrt_N2_2)*groupsh[j][i]; 
                          for ( i = 60; i < 64; i++ )
                             groupsh_t[j][15] -= (sqrt_N2_2)*groupsh[j][i];
                          for ( i = 0; i < 2; i++ )
                             groupsh_t[j][16] += (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 2; i < 4; i++ )
                             groupsh_t[j][16] -= (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 4; i < 6; i++ )
                             groupsh_t[j][17] += (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 6; i < 8; i++ )
                             groupsh_t[j][17] -= (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 8; i < 10; i++ )
                             groupsh_t[j][18] += (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 10; i < 12; i++ )
                             groupsh_t[j][18] -= (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 12; i < 14; i++ )
                             groupsh_t[j][19] += (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 14; i < 16; i++ )
                             groupsh_t[j][19] -= (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 16; i < 18; i++ )
                             groupsh_t[j][20] += (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 18; i < 20; i++ )
                             groupsh_t[j][20] -= (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 20; i < 22; i++ )
                             groupsh_t[j][21] += (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 22; i < 24; i++ )
                             groupsh_t[j][21] -= (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 24; i < 26; i++ )
                             groupsh_t[j][22] += (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 26; i < 28; i++ )
                             groupsh_t[j][22] -= (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 28; i < 30; i++ )
                             groupsh_t[j][23] += (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 30; i < 32; i++ )
                             groupsh_t[j][23] -= (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 32; i < 34; i++ )
                             groupsh_t[j][24] += (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 34; i < 36; i++ )
                             groupsh_t[j][24] -= (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 36; i < 38; i++ )
                             groupsh_t[j][25] += (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 38; i < 40; i++ )
                             groupsh_t[j][25] -= (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 40; i < 42; i++ )
                             groupsh_t[j][26] += (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 42; i < 44; i++ )
                             groupsh_t[j][26] -= (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 44; i < 46; i++ )
                             groupsh_t[j][27] += (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 46; i < 48; i++ )
                             groupsh_t[j][27] -= (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 48; i < 50; i++ )
                             groupsh_t[j][28] += (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 50; i < 52; i++ )
                             groupsh_t[j][28] -= (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 52; i < 54; i++ )
                             groupsh_t[j][29] += (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 54; i < 56; i++ )
                             groupsh_t[j][29] -= (sqrt_N2_4)*groupsh[j][i];
                           for ( i = 56; i < 58; i++ )
                             groupsh_t[j][30] += (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 58; i < 60; i++ )
                             groupsh_t[j][30] -= (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 60; i < 62; i++ )
                             groupsh_t[j][31] += (sqrt_N2_4)*groupsh[j][i];
                          for ( i = 62; i < 64; i++ )
                             groupsh_t[j][31] -= (sqrt_N2_4)*groupsh[j][i];
                          for ( l = 0; l < 32; l++ )
                             groupsh_t[j][l+32] = (sqrt_N2_8)*groupsh[j][2*l]-(sqrt_N2_8)*groupsh[j][2*l+1]; 
                                     }
                             } 
                       if (N22 == 64.0)
                             {
                          for ( j = 0; j < N12*N12; j++ )
                              
                                   {   
                         for ( i = 0; i < N22; i++ )
                            groupsh_t1[j][0] += (sqrt_N2_3)*groupsh1[j][i];
                          for ( i = 0; i < 32; i++ )
                            groupsh_t1[j][1] += (sqrt_N2_3)*groupsh1[j][i];
                          for ( i = 32; i < N22; i++ )
                             groupsh_t1[j][1] -= (sqrt_N2_3)*groupsh1[j][i];
                          for ( i = 0; i < N22/4; i++ )
                             groupsh_t1[j][2] += (sqrt_N2_1)*groupsh1[j][i];
                          for ( i = 16; i < 32; i++ )
                             groupsh_t1[j][2] -= (sqrt_N2_1)*groupsh1[j][i];
                          for ( i = 32; i < 48; i++ )
                             groupsh_t1[j][3] += (sqrt_N2_1)*groupsh1[j][i];
                          for ( i = 48; i < N22; i++ )
                             groupsh_t1[j][3] -= (sqrt_N2_1)*groupsh1[j][i];
                          for ( i = 0; i < 8; i++ )
                             groupsh_t1[j][4] += (sqrt_N2)*groupsh1[j][i];
                          for ( i = 8; i < 16; i++ )
                             groupsh_t1[j][4] -= (sqrt_N2)*groupsh1[j][i];
                          for ( i = 16; i < 24; i++ )
                             groupsh_t1[j][5] += (sqrt_N2)*groupsh1[j][i];
                          for ( i = 24; i < 32; i++ )
                             groupsh_t1[j][5] -= (sqrt_N2)*groupsh1[j][i];
                          for ( i = 32; i < 40; i++ )
                             groupsh_t1[j][6] += (sqrt_N2)*groupsh1[j][i]; 
                          for ( i = 40; i < 48; i++ )
                             groupsh_t1[j][6] -= (sqrt_N2)*groupsh1[j][i]; 
                          for ( i = 48; i < 56; i++ )
                             groupsh_t1[j][7] += (sqrt_N2)*groupsh1[j][i]; 
                          for ( i = 56; i < N22; i++ )
                             groupsh_t1[j][7] -= (sqrt_N2)*groupsh1[j][i];
                          for ( i = 0; i < 4; i++ )
                             groupsh_t1[j][8] += (sqrt_N2_2)*groupsh1[j][i];
                          for ( i = 4; i < 8; i++ )
                             groupsh_t1[j][8] -= (sqrt_N2_2)*groupsh1[j][i];
                          for ( i = 8; i < 12; i++ )
                             groupsh_t1[j][9] += (sqrt_N2_2)*groupsh1[j][i];
                          for ( i = 12; i < 16; i++ )
                             groupsh_t1[j][9] -= (sqrt_N2_2)*groupsh1[j][i];
                          for ( i = 16; i < 20; i++ )
                             groupsh_t1[j][10] += (sqrt_N2_2)*groupsh1[j][i]; 
                          for ( i = 20; i < 24; i++ )
                             groupsh_t1[j][10] -= (sqrt_N2_2)*groupsh1[j][i]; 
                          for ( i = 24; i < 28; i++ )
                             groupsh_t1[j][11] += (sqrt_N2_2)*groupsh1[j][i]; 
                          for ( i = 28; i < 32; i++ )
                             groupsh_t1[j][11] -= (sqrt_N2_2)*groupsh1[j][i];
                          for ( i = 32; i < 36; i++ )
                             groupsh_t1[j][12] += (sqrt_N2_2)*groupsh1[j][i];
                          for ( i = 36; i < 40; i++ )
                             groupsh_t1[j][12] -= (sqrt_N2_2)*groupsh1[j][i];
                          for ( i = 40; i < 44; i++ )
                             groupsh_t1[j][13] += (sqrt_N2_2)*groupsh1[j][i];
                          for ( i = 44; i < 48; i++ )
                             groupsh_t1[j][13] -= (sqrt_N2_2)*groupsh1[j][i];
                          for ( i = 48; i < 52; i++ )
                             groupsh_t1[j][14] += (sqrt_N2_2)*groupsh1[j][i]; 
                          for ( i = 52; i < 56; i++ )
                             groupsh_t1[j][14] -= (sqrt_N2_2)*groupsh1[j][i]; 
                          for ( i = 56; i < 60; i++ )
                             groupsh_t1[j][15] += (sqrt_N2_2)*groupsh1[j][i]; 
                          for ( i = 60; i < 64; i++ )
                             groupsh_t1[j][15] -= (sqrt_N2_2)*groupsh1[j][i];
                          for ( i = 0; i < 2; i++ )
                             groupsh_t1[j][16] += (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 2; i < 4; i++ )
                             groupsh_t1[j][16] -= (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 4; i < 6; i++ )
                             groupsh_t1[j][17] += (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 6; i < 8; i++ )
                             groupsh_t1[j][17] -= (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 8; i < 10; i++ )
                             groupsh_t1[j][18] += (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 10; i < 12; i++ )
                             groupsh_t1[j][18] -= (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 12; i < 14; i++ )
                             groupsh_t1[j][19] += (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 14; i < 16; i++ )
                             groupsh_t1[j][19] -= (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 16; i < 18; i++ )
                             groupsh_t1[j][20] += (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 18; i < 20; i++ )
                             groupsh_t1[j][20] -= (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 20; i < 22; i++ )
                             groupsh_t1[j][21] += (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 22; i < 24; i++ )
                             groupsh_t1[j][21] -= (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 24; i < 26; i++ )
                             groupsh_t1[j][22] += (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 26; i < 28; i++ )
                             groupsh_t1[j][22] -= (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 28; i < 30; i++ )
                             groupsh_t1[j][23] += (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 30; i < 32; i++ )
                             groupsh_t1[j][23] -= (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 32; i < 34; i++ )
                             groupsh_t1[j][24] += (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 34; i < 36; i++ )
                             groupsh_t1[j][24] -= (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 36; i < 38; i++ )
                             groupsh_t1[j][25] += (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 38; i < 40; i++ )
                             groupsh_t1[j][25] -= (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 40; i < 42; i++ )
                             groupsh_t1[j][26] += (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 42; i < 44; i++ )
                             groupsh_t1[j][26] -= (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 44; i < 46; i++ )
                             groupsh_t1[j][27] += (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 46; i < 48; i++ )
                             groupsh_t1[j][27] -= (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 48; i < 50; i++ )
                             groupsh_t1[j][28] += (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 50; i < 52; i++ )
                             groupsh_t1[j][28] -= (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 52; i < 54; i++ )
                             groupsh_t1[j][29] += (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 54; i < 56; i++ )
                             groupsh_t1[j][29] -= (sqrt_N2_4)*groupsh1[j][i];
                           for ( i = 56; i < 58; i++ )
                             groupsh_t1[j][30] += (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 58; i < 60; i++ )
                             groupsh_t1[j][30] -= (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 60; i < 62; i++ )
                             groupsh_t1[j][31] += (sqrt_N2_4)*groupsh1[j][i];
                          for ( i = 62; i < 64; i++ )
                             groupsh_t1[j][31] -= (sqrt_N2_4)*groupsh1[j][i];
                          for ( l = 0; l < 32; l++ )
                             groupsh_t1[j][l+32] = (sqrt_N2_8)*groupsh1[j][2*l]-(sqrt_N2_8)*groupsh1[j][2*l+1]; 
                                     }
                             } 
          
                       if (N22 == 64.0)
                             {
                          for ( j = 0; j < N12*N12; j++ )
                              
                                   {   
                         for ( i = 0; i < N22; i++ )
                            groupsh_t2[j][0] += (sqrt_N2_3)*groupsh2[j][i];
                          for ( i = 0; i < 32; i++ )
                            groupsh_t2[j][1] += (sqrt_N2_3)*groupsh2[j][i];
                          for ( i = 32; i < N22; i++ )
                             groupsh_t2[j][1] -= (sqrt_N2_3)*groupsh2[j][i];
                          for ( i = 0; i < N22/4; i++ )
                             groupsh_t2[j][2] += (sqrt_N2_1)*groupsh2[j][i];
                          for ( i = 16; i < 32; i++ )
                             groupsh_t2[j][2] -= (sqrt_N2_1)*groupsh2[j][i];
                          for ( i = 32; i < 48; i++ )
                             groupsh_t2[j][3] += (sqrt_N2_1)*groupsh2[j][i];
                          for ( i = 48; i < N22; i++ )
                             groupsh_t2[j][3] -= (sqrt_N2_1)*groupsh2[j][i];
                          for ( i = 0; i < 8; i++ )
                             groupsh_t2[j][4] += (sqrt_N2)*groupsh2[j][i];
                          for ( i = 8; i < 16; i++ )
                             groupsh_t2[j][4] -= (sqrt_N2)*groupsh2[j][i];
                          for ( i = 16; i < 24; i++ )
                             groupsh_t2[j][5] += (sqrt_N2)*groupsh2[j][i];
                          for ( i = 24; i < 32; i++ )
                             groupsh_t2[j][5] -= (sqrt_N2)*groupsh2[j][i];
                          for ( i = 32; i < 40; i++ )
                             groupsh_t2[j][6] += (sqrt_N2)*groupsh2[j][i]; 
                          for ( i = 40; i < 48; i++ )
                             groupsh_t2[j][6] -= (sqrt_N2)*groupsh2[j][i]; 
                          for ( i = 48; i < 56; i++ )
                             groupsh_t2[j][7] += (sqrt_N2)*groupsh2[j][i]; 
                          for ( i = 56; i < N22; i++ )
                             groupsh_t2[j][7] -= (sqrt_N2)*groupsh2[j][i];
                          for ( i = 0; i < 4; i++ )
                             groupsh_t2[j][8] += (sqrt_N2_2)*groupsh2[j][i];
                          for ( i = 4; i < 8; i++ )
                             groupsh_t2[j][8] -= (sqrt_N2_2)*groupsh2[j][i];
                          for ( i = 8; i < 12; i++ )
                             groupsh_t2[j][9] += (sqrt_N2_2)*groupsh2[j][i];
                          for ( i = 12; i < 16; i++ )
                             groupsh_t2[j][9] -= (sqrt_N2_2)*groupsh2[j][i];
                          for ( i = 16; i < 20; i++ )
                             groupsh_t2[j][10] += (sqrt_N2_2)*groupsh2[j][i]; 
                          for ( i = 20; i < 24; i++ )
                             groupsh_t2[j][10] -= (sqrt_N2_2)*groupsh2[j][i]; 
                          for ( i = 24; i < 28; i++ )
                             groupsh_t2[j][11] += (sqrt_N2_2)*groupsh2[j][i]; 
                          for ( i = 28; i < 32; i++ )
                             groupsh_t2[j][11] -= (sqrt_N2_2)*groupsh2[j][i];
                          for ( i = 32; i < 36; i++ )
                             groupsh_t2[j][12] += (sqrt_N2_2)*groupsh2[j][i];
                          for ( i = 36; i < 40; i++ )
                             groupsh_t2[j][12] -= (sqrt_N2_2)*groupsh2[j][i];
                          for ( i = 40; i < 44; i++ )
                             groupsh_t2[j][13] += (sqrt_N2_2)*groupsh2[j][i];
                          for ( i = 44; i < 48; i++ )
                             groupsh_t2[j][13] -= (sqrt_N2_2)*groupsh2[j][i];
                          for ( i = 48; i < 52; i++ )
                             groupsh_t2[j][14] += (sqrt_N2_2)*groupsh2[j][i]; 
                          for ( i = 52; i < 56; i++ )
                             groupsh_t2[j][14] -= (sqrt_N2_2)*groupsh2[j][i]; 
                          for ( i = 56; i < 60; i++ )
                             groupsh_t2[j][15] += (sqrt_N2_2)*groupsh2[j][i]; 
                          for ( i = 60; i < 64; i++ )
                             groupsh_t2[j][15] -= (sqrt_N2_2)*groupsh2[j][i];
                          for ( i = 0; i < 2; i++ )
                             groupsh_t2[j][16] += (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 2; i < 4; i++ )
                             groupsh_t2[j][16] -= (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 4; i < 6; i++ )
                             groupsh_t2[j][17] += (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 6; i < 8; i++ )
                             groupsh_t2[j][17] -= (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 8; i < 10; i++ )
                             groupsh_t2[j][18] += (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 10; i < 12; i++ )
                             groupsh_t2[j][18] -= (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 12; i < 14; i++ )
                             groupsh_t2[j][19] += (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 14; i < 16; i++ )
                             groupsh_t2[j][19] -= (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 16; i < 18; i++ )
                             groupsh_t2[j][20] += (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 18; i < 20; i++ )
                             groupsh_t2[j][20] -= (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 20; i < 22; i++ )
                             groupsh_t2[j][21] += (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 22; i < 24; i++ )
                             groupsh_t2[j][21] -= (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 24; i < 26; i++ )
                             groupsh_t2[j][22] += (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 26; i < 28; i++ )
                             groupsh_t2[j][22] -= (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 28; i < 30; i++ )
                             groupsh_t2[j][23] += (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 30; i < 32; i++ )
                             groupsh_t2[j][23] -= (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 32; i < 34; i++ )
                             groupsh_t2[j][24] += (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 34; i < 36; i++ )
                             groupsh_t2[j][24] -= (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 36; i < 38; i++ )
                             groupsh_t2[j][25] += (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 38; i < 40; i++ )
                             groupsh_t2[j][25] -= (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 40; i < 42; i++ )
                             groupsh_t2[j][26] += (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 42; i < 44; i++ )
                             groupsh_t2[j][26] -= (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 44; i < 46; i++ )
                             groupsh_t2[j][27] += (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 46; i < 48; i++ )
                             groupsh_t2[j][27] -= (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 48; i < 50; i++ )
                             groupsh_t2[j][28] += (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 50; i < 52; i++ )
                             groupsh_t2[j][28] -= (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 52; i < 54; i++ )
                             groupsh_t2[j][29] += (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 54; i < 56; i++ )
                             groupsh_t2[j][29] -= (sqrt_N2_4)*groupsh2[j][i];
                           for ( i = 56; i < 58; i++ )
                             groupsh_t2[j][30] += (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 58; i < 60; i++ )
                             groupsh_t2[j][30] -= (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 60; i < 62; i++ )
                             groupsh_t2[j][31] += (sqrt_N2_4)*groupsh2[j][i];
                          for ( i = 62; i < 64; i++ )
                             groupsh_t2[j][31] -= (sqrt_N2_4)*groupsh2[j][i];
                          for ( l = 0; l < 32; l++ )
                             groupsh_t2[j][l+32] = (sqrt_N2_8)*groupsh2[j][2*l]-(sqrt_N2_8)*groupsh2[j][2*l+1]; 
                                     }
                             } 
                
                
               // Haar正变换for wiener  
                
                                
                          if (N22 == 16.0)
                             {
                          for ( j = 0; j < N12*N12; j++ )
                                
                                         {  
                          for ( i = 0; i < N22; i++ )
                             
                             groupsh_tw[j][0] += groupshw[j][i];
                             groupsh_tw[j][0] = sqrt_N2*groupsh_tw[j][0];
                                                         
                         for ( i = 0; i < 8; i++ )

                             groupsh_tw[j][1] += groupshw[j][i];
                             
                             
                          for ( i = 8; i < 16; i++ )

                              groupsh_tw[j][1] -= groupshw[j][i];
                              groupsh_tw[j][1] *= sqrt_N2;
                          for ( i = 0; i < 4; i++ )

                              groupsh_tw[j][2] += groupshw[j][i];
                             
                          for ( i = 4; i < 8; i++ )

                              groupsh_tw[j][2] -= groupshw[j][i];
                             groupsh_tw[j][2] *= sqrt_N2_2;
                          for ( i = 8; i < 12; i++ )

                              groupsh_tw[j][3] += groupshw[j][i];
                          for ( i = 12; i < 16; i++ )

                              groupsh_tw[j][3] -= groupshw[j][i];
                              groupsh_tw[j][3] *= sqrt_N2_2;
                          
                         for ( i = 0; i < 2; i++ )

                              groupsh_tw[j][4] += groupshw[j][i];
                          for ( i = 2; i < 4; i++ )

                              groupsh_tw[j][4] -= groupshw[j][i];
                              groupsh_tw[j][4] *= sqrt_N2_4;
                              
                          for ( i = 4; i < 6; i++ )

                              groupsh_tw[j][5] += groupshw[j][i];
                          for ( i = 6; i < 8; i++ )

                              groupsh_tw[j][5] -= groupshw[j][i];
                              groupsh_tw[j][5] *= sqrt_N2_4;
                           for ( i = 8; i < 10; i++ )

                              groupsh_tw[j][6] += groupshw[j][i]; 
                           for ( i = 10; i < 12; i++ )

                              groupsh_tw[j][6] -= groupshw[j][i]; 
                              groupsh_tw[j][6] *= sqrt_N2_4;
                           for ( i = 12; i < 14; i++ )

                              groupsh_tw[j][7] += groupshw[j][i]; 
                           for ( i = 14; i < 16; i++ ) 
 
                              groupsh_tw[j][7] -= groupshw[j][i];  
                              groupsh_tw[j][7] *= sqrt_N2_4;
                          for ( l = 0; l < 8; l++ )
                              groupsh_tw[j][l+8] = sqrt_N2_8*(groupshw[j][2*l]-groupshw[j][2*l+1]);
                              
                                     }
                              }
                
                           if (N22 == 16.0)
                             {
                          for ( j = 0; j < N12*N12; j++ )
                                
                                         {  
                          for ( i = 0; i < N22; i++ )
                             
                             groupsh_tw1[j][0] += groupshw1[j][i];
                             groupsh_tw1[j][0] = sqrt_N2*groupsh_tw1[j][0];
                                                         
                         for ( i = 0; i < 8; i++ )

                             groupsh_tw1[j][1] += groupshw1[j][i];
                             
                             
                          for ( i = 8; i < 16; i++ )

                              groupsh_tw1[j][1] -= groupshw1[j][i];
                              groupsh_tw1[j][1] *= sqrt_N2;
                          for ( i = 0; i < 4; i++ )

                              groupsh_tw1[j][2] += groupshw1[j][i];
                             
                          for ( i = 4; i < 8; i++ )

                              groupsh_tw1[j][2] -= groupshw1[j][i];
                             groupsh_tw1[j][2] *= sqrt_N2_2;
                          for ( i = 8; i < 12; i++ )

                              groupsh_tw1[j][3] += groupshw1[j][i];
                          for ( i = 12; i < 16; i++ )

                              groupsh_tw1[j][3] -= groupshw1[j][i];
                              groupsh_tw1[j][3] *= sqrt_N2_2;
                          
                         for ( i = 0; i < 2; i++ )

                              groupsh_tw1[j][4] += groupshw1[j][i];
                          for ( i = 2; i < 4; i++ )

                              groupsh_tw1[j][4] -= groupshw1[j][i];
                              groupsh_tw1[j][4] *= sqrt_N2_4;
                              
                          for ( i = 4; i < 6; i++ )

                              groupsh_tw1[j][5] += groupshw1[j][i];
                          for ( i = 6; i < 8; i++ )

                              groupsh_tw1[j][5] -= groupshw1[j][i];
                              groupsh_tw1[j][5] *= sqrt_N2_4;
                           for ( i = 8; i < 10; i++ )

                              groupsh_tw1[j][6] += groupshw1[j][i]; 
                           for ( i = 10; i < 12; i++ )

                              groupsh_tw1[j][6] -= groupshw1[j][i]; 
                              groupsh_tw1[j][6] *= sqrt_N2_4;
                           for ( i = 12; i < 14; i++ )

                              groupsh_tw1[j][7] += groupshw1[j][i]; 
                           for ( i = 14; i < 16; i++ ) 
 
                              groupsh_tw1[j][7] -= groupshw1[j][i];  
                              groupsh_tw1[j][7] *= sqrt_N2_4;
                          for ( l = 0; l < 8; l++ )
                              groupsh_tw1[j][l+8] = sqrt_N2_8*(groupshw1[j][2*l]-groupshw1[j][2*l+1]);
                              
                                     }
                              }
                
                          if (N22 == 16.0)
                             {
                          for ( j = 0; j < N12*N12; j++ )
                                
                                         {  
                          for ( i = 0; i < N22; i++ )
                             
                             groupsh_tw2[j][0] += groupshw2[j][i];
                             groupsh_tw2[j][0] = sqrt_N2*groupsh_tw2[j][0];
                                                         
                         for ( i = 0; i < 8; i++ )

                             groupsh_tw2[j][1] += groupshw2[j][i];
                             
                             
                          for ( i = 8; i < 16; i++ )

                              groupsh_tw2[j][1] -= groupshw2[j][i];
                              groupsh_tw2[j][1] *= sqrt_N2;
                          for ( i = 0; i < 4; i++ )

                              groupsh_tw2[j][2] += groupshw2[j][i];
                             
                          for ( i = 4; i < 8; i++ )

                              groupsh_tw2[j][2] -= groupshw2[j][i];
                             groupsh_tw2[j][2] *= sqrt_N2_2;
                          for ( i = 8; i < 12; i++ )

                              groupsh_tw2[j][3] += groupshw2[j][i];
                          for ( i = 12; i < 16; i++ )

                              groupsh_tw2[j][3] -= groupshw2[j][i];
                              groupsh_tw2[j][3] *= sqrt_N2_2;
                          
                         for ( i = 0; i < 2; i++ )

                              groupsh_tw2[j][4] += groupshw2[j][i];
                          for ( i = 2; i < 4; i++ )

                              groupsh_tw2[j][4] -= groupshw2[j][i];
                              groupsh_tw2[j][4] *= sqrt_N2_4;
                              
                          for ( i = 4; i < 6; i++ )

                              groupsh_tw2[j][5] += groupshw2[j][i];
                          for ( i = 6; i < 8; i++ )

                              groupsh_tw2[j][5] -= groupshw2[j][i];
                              groupsh_tw2[j][5] *= sqrt_N2_4;
                           for ( i = 8; i < 10; i++ )

                              groupsh_tw2[j][6] += groupshw2[j][i]; 
                           for ( i = 10; i < 12; i++ )

                              groupsh_tw2[j][6] -= groupshw2[j][i]; 
                              groupsh_tw2[j][6] *= sqrt_N2_4;
                           for ( i = 12; i < 14; i++ )

                              groupsh_tw2[j][7] += groupshw2[j][i]; 
                           for ( i = 14; i < 16; i++ ) 
 
                              groupsh_tw2[j][7] -= groupshw2[j][i];  
                              groupsh_tw2[j][7] *= sqrt_N2_4;
                          for ( l = 0; l < 8; l++ )
                              groupsh_tw2[j][l+8] = sqrt_N2_8*(groupshw2[j][2*l]-groupshw2[j][2*l+1]);
                              
                                     }
                              }
                
                         
                      
                
                          if (N22 == 64.0)
                             {
                          for ( j = 0; j < N12*N12; j++ )
                              
                                   {   
                          for ( i = 0; i < N22; i++ )
                            groupsh_tw[j][0] += (sqrt_N2_3)*groupshw[j][i];
                          for ( i = 0; i < 32; i++ )
                            groupsh_tw[j][1] += (sqrt_N2_3)*groupshw[j][i];
                          for ( i = 32; i < N22; i++ )
                             groupsh_tw[j][1] -= (sqrt_N2_3)*groupshw[j][i];
                          for ( i = 0; i < N22/4; i++ )
                             groupsh_tw[j][2] += (sqrt_N2_1)*groupshw[j][i];
                          for ( i = 16; i < 32; i++ )
                             groupsh_tw[j][2] -= (sqrt_N2_1)*groupshw[j][i];
                          for ( i = 32; i < 48; i++ )
                             groupsh_tw[j][3] += (sqrt_N2_1)*groupshw[j][i];
                          for ( i = 48; i < N22; i++ )
                             groupsh_tw[j][3] -= (sqrt_N2_1)*groupshw[j][i];
                          for ( i = 0; i < 8; i++ )
                             groupsh_tw[j][4] += (sqrt_N2)*groupshw[j][i];
                          for ( i = 8; i < 16; i++ )
                             groupsh_tw[j][4] -= (sqrt_N2)*groupshw[j][i];
                          for ( i = 16; i < 24; i++ )
                             groupsh_tw[j][5] += (sqrt_N2)*groupshw[j][i];
                          for ( i = 24; i < 32; i++ )
                             groupsh_tw[j][5] -= (sqrt_N2)*groupshw[j][i];
                          for ( i = 32; i < 40; i++ )
                             groupsh_tw[j][6] += (sqrt_N2)*groupshw[j][i]; 
                          for ( i = 40; i < 48; i++ )
                             groupsh_tw[j][6] -= (sqrt_N2)*groupshw[j][i]; 
                          for ( i = 48; i < 56; i++ )
                             groupsh_tw[j][7] += (sqrt_N2)*groupshw[j][i]; 
                          for ( i = 56; i < N22; i++ )
                             groupsh_tw[j][7] -= (sqrt_N2)*groupshw[j][i];
                          for ( i = 0; i < 4; i++ )
                             groupsh_tw[j][8] += (sqrt_N2_2)*groupshw[j][i];
                          for ( i = 4; i < 8; i++ )
                             groupsh_tw[j][8] -= (sqrt_N2_2)*groupshw[j][i];
                          for ( i = 8; i < 12; i++ )
                             groupsh_tw[j][9] += (sqrt_N2_2)*groupshw[j][i];
                          for ( i = 12; i < 16; i++ )
                             groupsh_tw[j][9] -= (sqrt_N2_2)*groupshw[j][i];
                          for ( i = 16; i < 20; i++ )
                             groupsh_tw[j][10] += (sqrt_N2_2)*groupshw[j][i]; 
                          for ( i = 20; i < 24; i++ )
                             groupsh_tw[j][10] -= (sqrt_N2_2)*groupshw[j][i]; 
                          for ( i = 24; i < 28; i++ )
                             groupsh_tw[j][11] += (sqrt_N2_2)*groupshw[j][i]; 
                          for ( i = 28; i < 32; i++ )
                             groupsh_tw[j][11] -= (sqrt_N2_2)*groupshw[j][i];
                          for ( i = 32; i < 36; i++ )
                             groupsh_tw[j][12] += (sqrt_N2_2)*groupshw[j][i];
                          for ( i = 36; i < 40; i++ )
                             groupsh_tw[j][12] -= (sqrt_N2_2)*groupshw[j][i];
                          for ( i = 40; i < 44; i++ )
                             groupsh_tw[j][13] += (sqrt_N2_2)*groupshw[j][i];
                          for ( i = 44; i < 48; i++ )
                             groupsh_tw[j][13] -= (sqrt_N2_2)*groupshw[j][i];
                          for ( i = 48; i < 52; i++ )
                             groupsh_tw[j][14] += (sqrt_N2_2)*groupshw[j][i]; 
                          for ( i = 52; i < 56; i++ )
                             groupsh_tw[j][14] -= (sqrt_N2_2)*groupshw[j][i]; 
                          for ( i = 56; i < 60; i++ )
                             groupsh_tw[j][15] += (sqrt_N2_2)*groupshw[j][i]; 
                          for ( i = 60; i < 64; i++ )
                             groupsh_tw[j][15] -= (sqrt_N2_2)*groupshw[j][i];
                          for ( i = 0; i < 2; i++ )
                             groupsh_tw[j][16] += (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 2; i < 4; i++ )
                             groupsh_tw[j][16] -= (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 4; i < 6; i++ )
                             groupsh_tw[j][17] += (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 6; i < 8; i++ )
                             groupsh_tw[j][17] -= (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 8; i < 10; i++ )
                             groupsh_tw[j][18] += (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 10; i < 12; i++ )
                             groupsh_tw[j][18] -= (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 12; i < 14; i++ )
                             groupsh_tw[j][19] += (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 14; i < 16; i++ )
                             groupsh_tw[j][19] -= (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 16; i < 18; i++ )
                             groupsh_tw[j][20] += (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 18; i < 20; i++ )
                             groupsh_tw[j][20] -= (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 20; i < 22; i++ )
                             groupsh_tw[j][21] += (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 22; i < 24; i++ )
                             groupsh_tw[j][21] -= (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 24; i < 26; i++ )
                             groupsh_tw[j][22] += (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 26; i < 28; i++ )
                             groupsh_tw[j][22] -= (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 28; i < 30; i++ )
                             groupsh_tw[j][23] += (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 30; i < 32; i++ )
                             groupsh_tw[j][23] -= (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 32; i < 34; i++ )
                             groupsh_tw[j][24] += (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 34; i < 36; i++ )
                             groupsh_tw[j][24] -= (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 36; i < 38; i++ )
                             groupsh_tw[j][25] += (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 38; i < 40; i++ )
                             groupsh_tw[j][25] -= (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 40; i < 42; i++ )
                             groupsh_tw[j][26] += (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 42; i < 44; i++ )
                             groupsh_tw[j][26] -= (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 44; i < 46; i++ )
                             groupsh_tw[j][27] += (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 46; i < 48; i++ )
                             groupsh_tw[j][27] -= (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 48; i < 50; i++ )
                             groupsh_tw[j][28] += (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 50; i < 52; i++ )
                             groupsh_tw[j][28] -= (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 52; i < 54; i++ )
                             groupsh_tw[j][29] += (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 54; i < 56; i++ )
                             groupsh_tw[j][29] -= (sqrt_N2_4)*groupshw[j][i];
                           for ( i = 56; i < 58; i++ )
                             groupsh_tw[j][30] += (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 58; i < 60; i++ )
                             groupsh_tw[j][30] -= (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 60; i < 62; i++ )
                             groupsh_tw[j][31] += (sqrt_N2_4)*groupshw[j][i];
                          for ( i = 62; i < 64; i++ )
                             groupsh_tw[j][31] -= (sqrt_N2_4)*groupshw[j][i];
                          for ( l = 0; l < 32; l++ )
                             groupsh_tw[j][l+32] = (sqrt_N2_8)*groupshw[j][2*l]-(sqrt_N2_8)*groupshw[j][2*l+1]; 
                                     }
                             }               
                                   
                           if (N22 == 64.0)
                             {
                          for ( j = 0; j < N12*N12; j++ )
                              
                                   {   
                          for ( i = 0; i < N22; i++ )
                            groupsh_tw1[j][0] += (sqrt_N2_3)*groupshw1[j][i];
                          for ( i = 0; i < 32; i++ )
                            groupsh_tw1[j][1] += (sqrt_N2_3)*groupshw1[j][i];
                          for ( i = 32; i < N22; i++ )
                             groupsh_tw1[j][1] -= (sqrt_N2_3)*groupshw1[j][i];
                          for ( i = 0; i < N22/4; i++ )
                             groupsh_tw1[j][2] += (sqrt_N2_1)*groupshw1[j][i];
                          for ( i = 16; i < 32; i++ )
                             groupsh_tw1[j][2] -= (sqrt_N2_1)*groupshw1[j][i];
                          for ( i = 32; i < 48; i++ )
                             groupsh_tw1[j][3] += (sqrt_N2_1)*groupshw1[j][i];
                          for ( i = 48; i < N22; i++ )
                             groupsh_tw1[j][3] -= (sqrt_N2_1)*groupshw1[j][i];
                          for ( i = 0; i < 8; i++ )
                             groupsh_tw1[j][4] += (sqrt_N2)*groupshw1[j][i];
                          for ( i = 8; i < 16; i++ )
                             groupsh_tw1[j][4] -= (sqrt_N2)*groupshw1[j][i];
                          for ( i = 16; i < 24; i++ )
                             groupsh_tw1[j][5] += (sqrt_N2)*groupshw1[j][i];
                          for ( i = 24; i < 32; i++ )
                             groupsh_tw1[j][5] -= (sqrt_N2)*groupshw1[j][i];
                          for ( i = 32; i < 40; i++ )
                             groupsh_tw1[j][6] += (sqrt_N2)*groupshw1[j][i]; 
                          for ( i = 40; i < 48; i++ )
                             groupsh_tw1[j][6] -= (sqrt_N2)*groupshw1[j][i]; 
                          for ( i = 48; i < 56; i++ )
                             groupsh_tw1[j][7] += (sqrt_N2)*groupshw1[j][i]; 
                          for ( i = 56; i < N22; i++ )
                             groupsh_tw1[j][7] -= (sqrt_N2)*groupshw1[j][i];
                          for ( i = 0; i < 4; i++ )
                             groupsh_tw1[j][8] += (sqrt_N2_2)*groupshw1[j][i];
                          for ( i = 4; i < 8; i++ )
                             groupsh_tw1[j][8] -= (sqrt_N2_2)*groupshw1[j][i];
                          for ( i = 8; i < 12; i++ )
                             groupsh_tw1[j][9] += (sqrt_N2_2)*groupshw1[j][i];
                          for ( i = 12; i < 16; i++ )
                             groupsh_tw1[j][9] -= (sqrt_N2_2)*groupshw1[j][i];
                          for ( i = 16; i < 20; i++ )
                             groupsh_tw1[j][10] += (sqrt_N2_2)*groupshw1[j][i]; 
                          for ( i = 20; i < 24; i++ )
                             groupsh_tw1[j][10] -= (sqrt_N2_2)*groupshw1[j][i]; 
                          for ( i = 24; i < 28; i++ )
                             groupsh_tw1[j][11] += (sqrt_N2_2)*groupshw1[j][i]; 
                          for ( i = 28; i < 32; i++ )
                             groupsh_tw1[j][11] -= (sqrt_N2_2)*groupshw1[j][i];
                          for ( i = 32; i < 36; i++ )
                             groupsh_tw1[j][12] += (sqrt_N2_2)*groupshw1[j][i];
                          for ( i = 36; i < 40; i++ )
                             groupsh_tw1[j][12] -= (sqrt_N2_2)*groupshw1[j][i];
                          for ( i = 40; i < 44; i++ )
                             groupsh_tw1[j][13] += (sqrt_N2_2)*groupshw1[j][i];
                          for ( i = 44; i < 48; i++ )
                             groupsh_tw1[j][13] -= (sqrt_N2_2)*groupshw1[j][i];
                          for ( i = 48; i < 52; i++ )
                             groupsh_tw1[j][14] += (sqrt_N2_2)*groupshw1[j][i]; 
                          for ( i = 52; i < 56; i++ )
                             groupsh_tw1[j][14] -= (sqrt_N2_2)*groupshw1[j][i]; 
                          for ( i = 56; i < 60; i++ )
                             groupsh_tw1[j][15] += (sqrt_N2_2)*groupshw1[j][i]; 
                          for ( i = 60; i < 64; i++ )
                             groupsh_tw1[j][15] -= (sqrt_N2_2)*groupshw1[j][i];
                          for ( i = 0; i < 2; i++ )
                             groupsh_tw1[j][16] += (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 2; i < 4; i++ )
                             groupsh_tw1[j][16] -= (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 4; i < 6; i++ )
                             groupsh_tw1[j][17] += (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 6; i < 8; i++ )
                             groupsh_tw1[j][17] -= (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 8; i < 10; i++ )
                             groupsh_tw1[j][18] += (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 10; i < 12; i++ )
                             groupsh_tw1[j][18] -= (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 12; i < 14; i++ )
                             groupsh_tw1[j][19] += (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 14; i < 16; i++ )
                             groupsh_tw1[j][19] -= (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 16; i < 18; i++ )
                             groupsh_tw1[j][20] += (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 18; i < 20; i++ )
                             groupsh_tw1[j][20] -= (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 20; i < 22; i++ )
                             groupsh_tw1[j][21] += (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 22; i < 24; i++ )
                             groupsh_tw1[j][21] -= (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 24; i < 26; i++ )
                             groupsh_tw1[j][22] += (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 26; i < 28; i++ )
                             groupsh_tw1[j][22] -= (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 28; i < 30; i++ )
                             groupsh_tw1[j][23] += (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 30; i < 32; i++ )
                             groupsh_tw1[j][23] -= (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 32; i < 34; i++ )
                             groupsh_tw1[j][24] += (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 34; i < 36; i++ )
                             groupsh_tw1[j][24] -= (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 36; i < 38; i++ )
                             groupsh_tw1[j][25] += (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 38; i < 40; i++ )
                             groupsh_tw1[j][25] -= (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 40; i < 42; i++ )
                             groupsh_tw1[j][26] += (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 42; i < 44; i++ )
                             groupsh_tw1[j][26] -= (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 44; i < 46; i++ )
                             groupsh_tw1[j][27] += (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 46; i < 48; i++ )
                             groupsh_tw1[j][27] -= (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 48; i < 50; i++ )
                             groupsh_tw1[j][28] += (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 50; i < 52; i++ )
                             groupsh_tw1[j][28] -= (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 52; i < 54; i++ )
                             groupsh_tw1[j][29] += (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 54; i < 56; i++ )
                             groupsh_tw1[j][29] -= (sqrt_N2_4)*groupshw1[j][i];
                           for ( i = 56; i < 58; i++ )
                             groupsh_tw1[j][30] += (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 58; i < 60; i++ )
                             groupsh_tw1[j][30] -= (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 60; i < 62; i++ )
                             groupsh_tw1[j][31] += (sqrt_N2_4)*groupshw1[j][i];
                          for ( i = 62; i < 64; i++ )
                             groupsh_tw1[j][31] -= (sqrt_N2_4)*groupshw1[j][i];
                          for ( l = 0; l < 32; l++ )
                             groupsh_tw1[j][l+32] = (sqrt_N2_8)*groupshw1[j][2*l]-(sqrt_N2_8)*groupshw1[j][2*l+1]; 
                                     }
                             }                   
                
                          if (N22 == 64.0)
                             {
                          for ( j = 0; j < N12*N12; j++ )
                              
                                   {   
                          for ( i = 0; i < N22; i++ )
                            groupsh_tw2[j][0] += (sqrt_N2_3)*groupshw2[j][i];
                          for ( i = 0; i < 32; i++ )
                            groupsh_tw2[j][1] += (sqrt_N2_3)*groupshw2[j][i];
                          for ( i = 32; i < N22; i++ )
                             groupsh_tw2[j][1] -= (sqrt_N2_3)*groupshw2[j][i];
                          for ( i = 0; i < N22/4; i++ )
                             groupsh_tw2[j][2] += (sqrt_N2_1)*groupshw2[j][i];
                          for ( i = 16; i < 32; i++ )
                             groupsh_tw2[j][2] -= (sqrt_N2_1)*groupshw2[j][i];
                          for ( i = 32; i < 48; i++ )
                             groupsh_tw2[j][3] += (sqrt_N2_1)*groupshw2[j][i];
                          for ( i = 48; i < N22; i++ )
                             groupsh_tw2[j][3] -= (sqrt_N2_1)*groupshw2[j][i];
                          for ( i = 0; i < 8; i++ )
                             groupsh_tw2[j][4] += (sqrt_N2)*groupshw2[j][i];
                          for ( i = 8; i < 16; i++ )
                             groupsh_tw2[j][4] -= (sqrt_N2)*groupshw2[j][i];
                          for ( i = 16; i < 24; i++ )
                             groupsh_tw2[j][5] += (sqrt_N2)*groupshw2[j][i];
                          for ( i = 24; i < 32; i++ )
                             groupsh_tw2[j][5] -= (sqrt_N2)*groupshw2[j][i];
                          for ( i = 32; i < 40; i++ )
                             groupsh_tw2[j][6] += (sqrt_N2)*groupshw2[j][i]; 
                          for ( i = 40; i < 48; i++ )
                             groupsh_tw2[j][6] -= (sqrt_N2)*groupshw2[j][i]; 
                          for ( i = 48; i < 56; i++ )
                             groupsh_tw2[j][7] += (sqrt_N2)*groupshw2[j][i]; 
                          for ( i = 56; i < N22; i++ )
                             groupsh_tw2[j][7] -= (sqrt_N2)*groupshw2[j][i];
                          for ( i = 0; i < 4; i++ )
                             groupsh_tw2[j][8] += (sqrt_N2_2)*groupshw2[j][i];
                          for ( i = 4; i < 8; i++ )
                             groupsh_tw2[j][8] -= (sqrt_N2_2)*groupshw2[j][i];
                          for ( i = 8; i < 12; i++ )
                             groupsh_tw2[j][9] += (sqrt_N2_2)*groupshw2[j][i];
                          for ( i = 12; i < 16; i++ )
                             groupsh_tw2[j][9] -= (sqrt_N2_2)*groupshw2[j][i];
                          for ( i = 16; i < 20; i++ )
                             groupsh_tw2[j][10] += (sqrt_N2_2)*groupshw2[j][i]; 
                          for ( i = 20; i < 24; i++ )
                             groupsh_tw2[j][10] -= (sqrt_N2_2)*groupshw2[j][i]; 
                          for ( i = 24; i < 28; i++ )
                             groupsh_tw2[j][11] += (sqrt_N2_2)*groupshw2[j][i]; 
                          for ( i = 28; i < 32; i++ )
                             groupsh_tw2[j][11] -= (sqrt_N2_2)*groupshw2[j][i];
                          for ( i = 32; i < 36; i++ )
                             groupsh_tw2[j][12] += (sqrt_N2_2)*groupshw2[j][i];
                          for ( i = 36; i < 40; i++ )
                             groupsh_tw2[j][12] -= (sqrt_N2_2)*groupshw2[j][i];
                          for ( i = 40; i < 44; i++ )
                             groupsh_tw2[j][13] += (sqrt_N2_2)*groupshw2[j][i];
                          for ( i = 44; i < 48; i++ )
                             groupsh_tw2[j][13] -= (sqrt_N2_2)*groupshw2[j][i];
                          for ( i = 48; i < 52; i++ )
                             groupsh_tw2[j][14] += (sqrt_N2_2)*groupshw2[j][i]; 
                          for ( i = 52; i < 56; i++ )
                             groupsh_tw2[j][14] -= (sqrt_N2_2)*groupshw2[j][i]; 
                          for ( i = 56; i < 60; i++ )
                             groupsh_tw2[j][15] += (sqrt_N2_2)*groupshw2[j][i]; 
                          for ( i = 60; i < 64; i++ )
                             groupsh_tw2[j][15] -= (sqrt_N2_2)*groupshw2[j][i];
                          for ( i = 0; i < 2; i++ )
                             groupsh_tw2[j][16] += (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 2; i < 4; i++ )
                             groupsh_tw2[j][16] -= (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 4; i < 6; i++ )
                             groupsh_tw2[j][17] += (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 6; i < 8; i++ )
                             groupsh_tw2[j][17] -= (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 8; i < 10; i++ )
                             groupsh_tw2[j][18] += (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 10; i < 12; i++ )
                             groupsh_tw2[j][18] -= (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 12; i < 14; i++ )
                             groupsh_tw2[j][19] += (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 14; i < 16; i++ )
                             groupsh_tw2[j][19] -= (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 16; i < 18; i++ )
                             groupsh_tw2[j][20] += (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 18; i < 20; i++ )
                             groupsh_tw2[j][20] -= (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 20; i < 22; i++ )
                             groupsh_tw2[j][21] += (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 22; i < 24; i++ )
                             groupsh_tw2[j][21] -= (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 24; i < 26; i++ )
                             groupsh_tw2[j][22] += (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 26; i < 28; i++ )
                             groupsh_tw2[j][22] -= (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 28; i < 30; i++ )
                             groupsh_tw2[j][23] += (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 30; i < 32; i++ )
                             groupsh_tw2[j][23] -= (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 32; i < 34; i++ )
                             groupsh_tw2[j][24] += (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 34; i < 36; i++ )
                             groupsh_tw2[j][24] -= (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 36; i < 38; i++ )
                             groupsh_tw2[j][25] += (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 38; i < 40; i++ )
                             groupsh_tw2[j][25] -= (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 40; i < 42; i++ )
                             groupsh_tw2[j][26] += (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 42; i < 44; i++ )
                             groupsh_tw2[j][26] -= (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 44; i < 46; i++ )
                             groupsh_tw2[j][27] += (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 46; i < 48; i++ )
                             groupsh_tw2[j][27] -= (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 48; i < 50; i++ )
                             groupsh_tw2[j][28] += (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 50; i < 52; i++ )
                             groupsh_tw2[j][28] -= (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 52; i < 54; i++ )
                             groupsh_tw2[j][29] += (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 54; i < 56; i++ )
                             groupsh_tw2[j][29] -= (sqrt_N2_4)*groupshw2[j][i];
                           for ( i = 56; i < 58; i++ )
                             groupsh_tw2[j][30] += (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 58; i < 60; i++ )
                             groupsh_tw2[j][30] -= (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 60; i < 62; i++ )
                             groupsh_tw2[j][31] += (sqrt_N2_4)*groupshw2[j][i];
                          for ( i = 62; i < 64; i++ )
                             groupsh_tw2[j][31] -= (sqrt_N2_4)*groupshw2[j][i];
                          for ( l = 0; l < 32; l++ )
                             groupsh_tw2[j][l+32] = (sqrt_N2_8)*groupshw2[j][2*l]-(sqrt_N2_8)*groupshw2[j][2*l+1]; 
                                     }
                             }                   
     for ( i = 0; i < N12*N12; i++ )
          for ( j = 0; j < N22; j++ ) 
             {     
                 DD[i][j] = 0.0;
                 WW[i][j] = 0.00000001;
                 DD1[i][j] = 0.00000001;
                 WW1[i][j] = 0.00000001;
                 DD2[i][j] = 0.00000001;
                 WW2[i][j] = 0.00000001;
             } 
                
    
              
 // 开始每一组内的行匹配                 
 for (k = 0; k < N12*N12; k++)    //开始循环
                {  
            
        if ( DD[k][0] != 0.0 )              
        continue;              
    
     for ( i = 0; i < N32; i++ )
          for ( j = 0; j < N22; j++ ) 
              { 
                  groups[i][j] = 0.0;
                  groups_t[i][j] = 0.0;
                  groupsw[i][j] = 0.0;
                  groups_tw[i][j] = 0.0;
                  groups1[i][j] = 0.0;
                  groups_t1[i][j] = 0.0;
                  groupsw1[i][j] = 0.0;
                  groups_tw1[i][j] = 0.0;
                  groups2[i][j] = 0.0;
                  groups_t2[i][j] = 0.0;
                  groupsw2[i][j] = 0.0;
                  groups_tw2[i][j] = 0.0;
              }  
           
     
           for (i = 0; i < N12*N12; i++)
                    dis4[i] = 0.0;
     
           for (i = 0; i < N12*N12; i++)
                    dis4[i] = dis3[i][k];
         
              for (i = 0; i < N32; i++)
                dis5[i] = 0.0;
     
                             
                   dis4[k] = 10000.0;
               
         for ( i = 0; i < N12*N12; i++ )   dis1Arrays[i] = dis4[i];   
             FindMin(dis1Arrays,N12*N12,IndexArrays,N32);
  
            
           dis5[0] = k;
             
           for ( i = 0; i < N32-1; i++ )
                dis5[i+1] = IndexArrays[i];  
           
     
            for ( i = 0; i < N32; i++)
                  for ( j = 0; j < N22; j++)
                       {
                           groups[i][j] = groupsh_t[dis5[i]][j];
                           groupsw[i][j] = groupsh_tw[dis5[i]][j];
                           groups1[i][j] = groupsh_t1[dis5[i]][j];
                           groupsw1[i][j] = groupsh_tw1[dis5[i]][j];
                           groups2[i][j] = groupsh_t2[dis5[i]][j];
                           groupsw2[i][j] = groupsh_tw2[dis5[i]][j];
                       }
                      
            
            
          // 纵向 Haar 正变换
             
                  
                          for ( j = 0; j < N22; j++ )
                                      {   
                          for ( i = 0; i < N32; i++ )
                              groups_t[0][j] += groups[i][j];
                              groups_t[0][j] *= (sqrt_N2_2);
                          for ( i = 0; i < 4; i++ )
                             groups_t[1][j] += groups[i][j];
                          for ( i = 4; i < N32; i++ )
                              groups_t[1][j] -= groups[i][j];
                              groups_t[1][j] *= (sqrt_N2_2);
                          for ( i = 0; i < 2; i++ )
                              groups_t[2][j] += groups[i][j];
                          for ( i = 2; i < 4; i++ )
                              groups_t[2][j] -= groups[i][j];
                              groups_t[2][j] *= (sqrt_N2_4);
                          for ( i = 4; i < 6; i++ )
                              groups_t[3][j] += groups[i][j];
                          for ( i = 6; i < N32; i++ )
                              groups_t[3][j] -= groups[i][j];
                              groups_t[3][j] *= (sqrt_N2_4);
                          for ( l = 0; l < 4; l++ )
                              groups_t[l+4][j] = (sqrt_N2_8)*(groups[2*l][j]-groups[2*l+1][j]);
                              
                                     }
                             
                          for ( j = 0; j < N22; j++ )
                                      {   
                          for ( i = 0; i < N32; i++ )
                              groups_t1[0][j] += groups1[i][j];
                              groups_t1[0][j] *= (sqrt_N2_2);
                          for ( i = 0; i < 4; i++ )
                             groups_t1[1][j] += groups1[i][j];
                          for ( i = 4; i < N32; i++ )
                              groups_t1[1][j] -= groups1[i][j];
                              groups_t1[1][j] *= (sqrt_N2_2);
                          for ( i = 0; i < 2; i++ )
                              groups_t1[2][j] += groups1[i][j];
                          for ( i = 2; i < 4; i++ )
                              groups_t1[2][j] -= groups1[i][j];
                              groups_t1[2][j] *= (sqrt_N2_4);
                          for ( i = 4; i < 6; i++ )
                              groups_t1[3][j] += groups1[i][j];
                          for ( i = 6; i < N32; i++ )
                              groups_t1[3][j] -= groups1[i][j];
                              groups_t1[3][j] *= (sqrt_N2_4);
                          for ( l = 0; l < 4; l++ )
                              groups_t1[l+4][j] = (sqrt_N2_8)*(groups1[2*l][j]-groups1[2*l+1][j]);
                              
                                     }
                            
                          for ( j = 0; j < N22; j++ )
                                      {   
                          for ( i = 0; i < N32; i++ )
                              groups_t2[0][j] += groups2[i][j];
                              groups_t2[0][j] *= (sqrt_N2_2);
                          for ( i = 0; i < 4; i++ )
                             groups_t2[1][j] += groups2[i][j];
                          for ( i = 4; i < N32; i++ )
                              groups_t2[1][j] -= groups2[i][j];
                              groups_t2[1][j] *= (sqrt_N2_2);
                          for ( i = 0; i < 2; i++ )
                              groups_t2[2][j] += groups2[i][j];
                          for ( i = 2; i < 4; i++ )
                              groups_t2[2][j] -= groups2[i][j];
                              groups_t2[2][j] *= (sqrt_N2_4);
                          for ( i = 4; i < 6; i++ )
                              groups_t2[3][j] += groups2[i][j];
                          for ( i = 6; i < N32; i++ )
                              groups_t2[3][j] -= groups2[i][j];
                              groups_t2[3][j] *= (sqrt_N2_4);
                          for ( l = 0; l < 4; l++ )
                              groups_t2[l+4][j] = (sqrt_N2_8)*(groups2[2*l][j]-groups2[2*l+1][j]);
                              
                                     }
                         
                        
             
                
                           // 纵向 Haar 正变换for wiener
             
                 
             
                   
                          for ( j = 0; j < N22; j++ )
                                      {   
                          for ( i = 0; i < N32; i++ )
                              groups_tw[0][j] += groupsw[i][j];
                              groups_tw[0][j] *= (sqrt_N2_2);
                          for ( i = 0; i < 4; i++ )
                              groups_tw[1][j] += groupsw[i][j];
                          for ( i = 4; i < N32; i++ )
                              groups_tw[1][j] -= groupsw[i][j];
                              groups_tw[1][j] *= (sqrt_N2_2);
                          for ( i = 0; i < 2; i++ )
                              groups_tw[2][j] += groupsw[i][j];
                          for ( i = 2; i < 4; i++ )
                              groups_tw[2][j] -= groupsw[i][j];
                              groups_tw[2][j] *= (sqrt_N2_4);
                          for ( i = 4; i < 6; i++ )
                              groups_tw[3][j] += groupsw[i][j];
                          for ( i = 6; i < N32; i++ )
                              groups_tw[3][j] -= groupsw[i][j];
                              groups_tw[3][j] *= (sqrt_N2_4);
                          for ( l = 0; l < 4; l++ )
                              groups_tw[l+4][j] = (sqrt_N2_8)*(groupsw[2*l][j]-groupsw[2*l+1][j]);
                              
                                     }
                            
                          for ( j = 0; j < N22; j++ )
                                      {   
                          for ( i = 0; i < N32; i++ )
                              groups_tw1[0][j] +=  groupsw1[i][j];
                              groups_tw1[0][j] *= (sqrt_N2_2);
                          for ( i = 0; i < 4; i++ )
                             groups_tw1[1][j] +=  groupsw1[i][j];
                          for ( i = 4; i < N32; i++ )
                              groups_tw1[1][j] -=  groupsw1[i][j];
                              groups_tw1[1][j] *= (sqrt_N2_2);
                          for ( i = 0; i < 2; i++ )
                              groups_tw1[2][j] +=  groupsw1[i][j];
                          for ( i = 2; i < 4; i++ )
                              groups_tw1[2][j] -=  groupsw1[i][j];
                              groups_tw1[2][j] *= (sqrt_N2_4);
                          for ( i = 4; i < 6; i++ )
                              groups_tw1[3][j] +=  groupsw1[i][j];
                          for ( i = 6; i < N32; i++ )
                              groups_tw1[3][j] -=  groupsw1[i][j];
                              groups_tw1[3][j] *= (sqrt_N2_4);
                          for ( l = 0; l < 4; l++ )
                              groups_tw1[l+4][j] = (sqrt_N2_8)*(groupsw1[2*l][j]-groupsw1[2*l+1][j]);
                              
                                     }
                          
                          for ( j = 0; j < N22; j++ )
                                      {   
                          for ( i = 0; i < N32; i++ )
                              groups_tw2[0][j] +=  groupsw2[i][j];
                              groups_tw2[0][j] *= (sqrt_N2_2);
                          for ( i = 0; i < 4; i++ )
                              groups_tw2[1][j] +=  groupsw2[i][j];
                          for ( i = 4; i < N32; i++ )
                              groups_tw2[1][j] -=  groupsw2[i][j];
                              groups_tw2[1][j] *= (sqrt_N2_2);
                          for ( i = 0; i < 2; i++ )
                              groups_tw2[2][j] +=  groupsw2[i][j];
                          for ( i = 2; i < 4; i++ )
                              groups_tw2[2][j] -=  groupsw2[i][j];
                              groups_tw2[2][j] *= (sqrt_N2_4);
                          for ( i = 4; i < 6; i++ )
                              groups_tw2[3][j] +=  groupsw2[i][j];
                          for ( i = 6; i < N32; i++ )
                              groups_tw2[3][j] -=  groupsw2[i][j];
                              groups_tw2[3][j] *= (sqrt_N2_4);
                          for ( l = 0; l < 4; l++ )
                              groups_tw2[l+4][j] = (sqrt_N2_8)*(groupsw2[2*l][j]-groupsw2[2*l+1][j]);
                              
                                     }
                         
  /////////////////////////////////////////////////////////////////////////////  Wiener filtering         
              weight = 0.00000000; weight1 = 0.00000000; weight2 = 0.00000000;
             
             
             for ( l = 0; l < N32; l++ )
                    for ( j = 0; j < N22; j++ )  
                               {
                                 group_tm = groups_tw[l][j]*groups_tw[l][j]+0.00000001;
                                 W_wiener[l][j] = group_tm /(group_tm+THR);
                                 weight += W_wiener[l][j]*W_wiener[l][j];
                                 groups_t[l][j]  = W_wiener[l][j] *groups_t[l][j];
                                 
                                group_tm = groups_tw1[l][j]*groups_tw1[l][j]+0.00000001;
                                W_wiener1[l][j] = group_tm/(group_tm+THR);
                                weight1 += W_wiener1[l][j]*W_wiener1[l][j];
                                groups_t1[l][j]  = W_wiener1[l][j] *groups_t1[l][j];
                                
                                group_tm = groups_tw2[l][j]*groups_tw2[l][j]+0.00000001;
                                W_wiener2[l][j] = group_tm/(group_tm+THR);
                                weight2 += W_wiener2[l][j]*W_wiener2[l][j];
                                groups_t2[l][j]  = W_wiener2[l][j] *groups_t2[l][j];
                                }
                  
                
                           weight = 1.00000000/(weight);
         
                           weight1 = 1.00000000/(weight1);
   
                           weight2 = 1.0000000/(weight2);           
                           
                           
//               for ( l = 1; l < N32; l++ )
//                     for ( j = 0; j < 1; j++ )  
//                                {
//                                 
//                                  groups_t[l][j]  = groups_t[l][j]-0.05;
//                                  
//                                 
//                                 groups_t1[l][j]  = 0.0 +groups_t1[l][j]-0.05;
//                                 
//                                
//                                 groups_t2[l][j]  = 0.0 +groups_t2[l][j]-0.05;
//                                 }             
                           
                                                           
               weightb += weight;
               weightb1 += weight1;
               weightb2 += weight2;
   ///////////////////////////////////////////////////////////////////////////////////////      
              // Haar小块逆变换    
                      
               
                            for ( j = 0; j < N22; j++ )
                             
                                {
                        groups[0][j] = (sqrt_N2_2)*(groups_t[0][j]+groups_t[1][j])+(sqrt_N2_4)*groups_t[2][j]+(sqrt_N2_8)*groups_t[4][j];
                        groups[1][j] = (sqrt_N2_2)*(groups_t[0][j]+groups_t[1][j])+(sqrt_N2_4)*groups_t[2][j]-(sqrt_N2_8)*groups_t[4][j];               
                        groups[2][j] = (sqrt_N2_2)*(groups_t[0][j]+groups_t[1][j])-(sqrt_N2_4)*groups_t[2][j]+(sqrt_N2_8)*groups_t[5][j];
                        groups[3][j] = (sqrt_N2_2)*(groups_t[0][j]+groups_t[1][j])-(sqrt_N2_4)*groups_t[2][j]-(sqrt_N2_8)*groups_t[5][j]; 
                        groups[4][j] = (sqrt_N2_2)*(groups_t[0][j]-groups_t[1][j])+(sqrt_N2_4)*groups_t[3][j]+(sqrt_N2_8)*groups_t[6][j];
                        groups[5][j] = (sqrt_N2_2)*(groups_t[0][j]-groups_t[1][j])+(sqrt_N2_4)*groups_t[3][j]-(sqrt_N2_8)*groups_t[6][j];               
                        groups[6][j] = (sqrt_N2_2)*(groups_t[0][j]-groups_t[1][j])-(sqrt_N2_4)*groups_t[3][j]+(sqrt_N2_8)*groups_t[7][j];
                        groups[7][j] = (sqrt_N2_2)*(groups_t[0][j]-groups_t[1][j])-(sqrt_N2_4)*groups_t[3][j]-(sqrt_N2_8)*groups_t[7][j]; 
                        groups1[0][j] = (sqrt_N2_2)*(groups_t1[0][j]+groups_t1[1][j])+(sqrt_N2_4)*groups_t1[2][j]+(sqrt_N2_8)*groups_t1[4][j];
                        groups1[1][j] = (sqrt_N2_2)*(groups_t1[0][j]+groups_t1[1][j])+(sqrt_N2_4)*groups_t1[2][j]-(sqrt_N2_8)*groups_t1[4][j];               
                        groups1[2][j] = (sqrt_N2_2)*(groups_t1[0][j]+groups_t1[1][j])-(sqrt_N2_4)*groups_t1[2][j]+(sqrt_N2_8)*groups_t1[5][j];
                        groups1[3][j] = (sqrt_N2_2)*(groups_t1[0][j]+groups_t1[1][j])-(sqrt_N2_4)*groups_t1[2][j]-(sqrt_N2_8)*groups_t1[5][j]; 
                        groups1[4][j] = (sqrt_N2_2)*(groups_t1[0][j]-groups_t1[1][j])+(sqrt_N2_4)*groups_t1[3][j]+(sqrt_N2_8)*groups_t1[6][j];
                        groups1[5][j] = (sqrt_N2_2)*(groups_t1[0][j]-groups_t1[1][j])+(sqrt_N2_4)*groups_t1[3][j]-(sqrt_N2_8)*groups_t1[6][j];               
                        groups1[6][j] = (sqrt_N2_2)*(groups_t1[0][j]-groups_t1[1][j])-(sqrt_N2_4)*groups_t1[3][j]+(sqrt_N2_8)*groups_t1[7][j];
                        groups1[7][j] = (sqrt_N2_2)*(groups_t1[0][j]-groups_t1[1][j])-(sqrt_N2_4)*groups_t1[3][j]-(sqrt_N2_8)*groups_t1[7][j]; 
                        groups2[0][j] = (sqrt_N2_2)*(groups_t2[0][j]+groups_t2[1][j])+(sqrt_N2_4)*groups_t2[2][j]+(sqrt_N2_8)*groups_t2[4][j];
                        groups2[1][j] = (sqrt_N2_2)*(groups_t2[0][j]+groups_t2[1][j])+(sqrt_N2_4)*groups_t2[2][j]-(sqrt_N2_8)*groups_t2[4][j];               
                        groups2[2][j] = (sqrt_N2_2)*(groups_t2[0][j]+groups_t2[1][j])-(sqrt_N2_4)*groups_t2[2][j]+(sqrt_N2_8)*groups_t2[5][j];
                        groups2[3][j] = (sqrt_N2_2)*(groups_t2[0][j]+groups_t2[1][j])-(sqrt_N2_4)*groups_t2[2][j]-(sqrt_N2_8)*groups_t2[5][j]; 
                        groups2[4][j] = (sqrt_N2_2)*(groups_t2[0][j]-groups_t2[1][j])+(sqrt_N2_4)*groups_t2[3][j]+(sqrt_N2_8)*groups_t2[6][j];
                        groups2[5][j] = (sqrt_N2_2)*(groups_t2[0][j]-groups_t2[1][j])+(sqrt_N2_4)*groups_t2[3][j]-(sqrt_N2_8)*groups_t2[6][j];               
                        groups2[6][j] = (sqrt_N2_2)*(groups_t2[0][j]-groups_t2[1][j])-(sqrt_N2_4)*groups_t2[3][j]+(sqrt_N2_8)*groups_t2[7][j];
                        groups2[7][j] = (sqrt_N2_2)*(groups_t2[0][j]-groups_t2[1][j])-(sqrt_N2_4)*groups_t2[3][j]-(sqrt_N2_8)*groups_t2[7][j]; 
                             }
                       
               
                    
                  
               // 小块逆变换结束
                 
             
               
               for ( i = 0; i < N32; i++ )
                    for ( j = 0; j < N22; j++ )
                    {
                       DD[dis5[i]][j] += groups[i][j]*weight;
                       WW[dis5[i]][j] += weight;
                       DD1[dis5[i]][j] += groups1[i][j]*weight1;
                       WW1[dis5[i]][j] += weight1;
                       DD2[dis5[i]][j] += groups2[i][j]*weight2;
                       WW2[dis5[i]][j] += weight2;
                    }
               
              
         }
     // 大块逆变换开始   
                
                   if (N22 == 16.0)
                           {
                        
                            for ( j = 0; j < N12*N12; j++ )
                            {
                        groupsh[j][0]  = sqrt_N2*(DD[j][0]+DD[j][1])+sqrt_N2_2*DD[j][2]+sqrt_N2_4*DD[j][4]+sqrt_N2_8*DD[j][8];
                        groupsh[j][1]  = sqrt_N2*(DD[j][0]+DD[j][1])+sqrt_N2_2*DD[j][2]+sqrt_N2_4*DD[j][4]-sqrt_N2_8*DD[j][8];
                        groupsh[j][2]  = sqrt_N2*(DD[j][0]+DD[j][1])+sqrt_N2_2*DD[j][2]-sqrt_N2_4*DD[j][4]+sqrt_N2_8*DD[j][9];
                        groupsh[j][3]  = sqrt_N2*(DD[j][0]+DD[j][1])+sqrt_N2_2*DD[j][2]-sqrt_N2_4*DD[j][4]-sqrt_N2_8*DD[j][9];
                        groupsh[j][4]  = sqrt_N2*(DD[j][0]+DD[j][1])-sqrt_N2_2*DD[j][2]+sqrt_N2_4*DD[j][5]+sqrt_N2_8*DD[j][10];
                        groupsh[j][5]  = sqrt_N2*(DD[j][0]+DD[j][1])-sqrt_N2_2*DD[j][2]+sqrt_N2_4*DD[j][5]-sqrt_N2_8*DD[j][10];
                        groupsh[j][6]  = sqrt_N2*(DD[j][0]+DD[j][1])-sqrt_N2_2*DD[j][2]-sqrt_N2_4*DD[j][5]+sqrt_N2_8*DD[j][11];
                        groupsh[j][7]  = sqrt_N2*(DD[j][0]+DD[j][1])-sqrt_N2_2*DD[j][2]-sqrt_N2_4*DD[j][5]-sqrt_N2_8*DD[j][11];
                        groupsh[j][8]  = sqrt_N2*(DD[j][0]-DD[j][1])+sqrt_N2_2*DD[j][3]+sqrt_N2_4*DD[j][6]+sqrt_N2_8*DD[j][12];
                        groupsh[j][9]  = sqrt_N2*(DD[j][0]-DD[j][1])+sqrt_N2_2*DD[j][3]+sqrt_N2_4*DD[j][6]-sqrt_N2_8*DD[j][12];
                        groupsh[j][10] = sqrt_N2*(DD[j][0]-DD[j][1])+sqrt_N2_2*DD[j][3]-sqrt_N2_4*DD[j][6]+sqrt_N2_8*DD[j][13];
                        groupsh[j][11] = sqrt_N2*(DD[j][0]-DD[j][1])+sqrt_N2_2*DD[j][3]-sqrt_N2_4*DD[j][6]-sqrt_N2_8*DD[j][13];
                        groupsh[j][12] = sqrt_N2*(DD[j][0]-DD[j][1])-sqrt_N2_2*DD[j][3]+sqrt_N2_4*DD[j][7]+sqrt_N2_8*DD[j][14];
                        groupsh[j][13] = sqrt_N2*(DD[j][0]-DD[j][1])-sqrt_N2_2*DD[j][3]+sqrt_N2_4*DD[j][7]-sqrt_N2_8*DD[j][14];
                        groupsh[j][14] = sqrt_N2*(DD[j][0]-DD[j][1])-sqrt_N2_2*DD[j][3]-sqrt_N2_4*DD[j][7]+sqrt_N2_8*DD[j][15];
                        groupsh[j][15] = sqrt_N2*(DD[j][0]-DD[j][1])-sqrt_N2_2*DD[j][3]-sqrt_N2_4*DD[j][7]-sqrt_N2_8*DD[j][15];
                        groupsh1[j][0]  = sqrt_N2*(DD1[j][0]+DD1[j][1])+sqrt_N2_2*DD1[j][2]+sqrt_N2_4*DD1[j][4]+sqrt_N2_8*DD1[j][8];
                        groupsh1[j][1]  = sqrt_N2*(DD1[j][0]+DD1[j][1])+sqrt_N2_2*DD1[j][2]+sqrt_N2_4*DD1[j][4]-sqrt_N2_8*DD1[j][8];
                        groupsh1[j][2]  = sqrt_N2*(DD1[j][0]+DD1[j][1])+sqrt_N2_2*DD1[j][2]-sqrt_N2_4*DD1[j][4]+sqrt_N2_8*DD1[j][9];
                        groupsh1[j][3]  = sqrt_N2*(DD1[j][0]+DD1[j][1])+sqrt_N2_2*DD1[j][2]-sqrt_N2_4*DD1[j][4]-sqrt_N2_8*DD1[j][9];
                        groupsh1[j][4]  = sqrt_N2*(DD1[j][0]+DD1[j][1])-sqrt_N2_2*DD1[j][2]+sqrt_N2_4*DD1[j][5]+sqrt_N2_8*DD1[j][10];
                        groupsh1[j][5]  = sqrt_N2*(DD1[j][0]+DD1[j][1])-sqrt_N2_2*DD1[j][2]+sqrt_N2_4*DD1[j][5]-sqrt_N2_8*DD1[j][10];
                        groupsh1[j][6]  = sqrt_N2*(DD1[j][0]+DD1[j][1])-sqrt_N2_2*DD1[j][2]-sqrt_N2_4*DD1[j][5]+sqrt_N2_8*DD1[j][11];
                        groupsh1[j][7]  = sqrt_N2*(DD1[j][0]+DD1[j][1])-sqrt_N2_2*DD1[j][2]-sqrt_N2_4*DD1[j][5]-sqrt_N2_8*DD1[j][11];
                        groupsh1[j][8]  = sqrt_N2*(DD1[j][0]-DD1[j][1])+sqrt_N2_2*DD1[j][3]+sqrt_N2_4*DD1[j][6]+sqrt_N2_8*DD1[j][12];
                        groupsh1[j][9]  = sqrt_N2*(DD1[j][0]-DD1[j][1])+sqrt_N2_2*DD1[j][3]+sqrt_N2_4*DD1[j][6]-sqrt_N2_8*DD1[j][12];
                        groupsh1[j][10] = sqrt_N2*(DD1[j][0]-DD1[j][1])+sqrt_N2_2*DD1[j][3]-sqrt_N2_4*DD1[j][6]+sqrt_N2_8*DD1[j][13];
                        groupsh1[j][11] = sqrt_N2*(DD1[j][0]-DD1[j][1])+sqrt_N2_2*DD1[j][3]-sqrt_N2_4*DD1[j][6]-sqrt_N2_8*DD1[j][13];
                        groupsh1[j][12] = sqrt_N2*(DD1[j][0]-DD1[j][1])-sqrt_N2_2*DD1[j][3]+sqrt_N2_4*DD1[j][7]+sqrt_N2_8*DD1[j][14];
                        groupsh1[j][13] = sqrt_N2*(DD1[j][0]-DD1[j][1])-sqrt_N2_2*DD1[j][3]+sqrt_N2_4*DD1[j][7]-sqrt_N2_8*DD1[j][14];
                        groupsh1[j][14] = sqrt_N2*(DD1[j][0]-DD1[j][1])-sqrt_N2_2*DD1[j][3]-sqrt_N2_4*DD1[j][7]+sqrt_N2_8*DD1[j][15];
                        groupsh1[j][15] = sqrt_N2*(DD1[j][0]-DD1[j][1])-sqrt_N2_2*DD1[j][3]-sqrt_N2_4*DD1[j][7]-sqrt_N2_8*DD1[j][15];
                        groupsh2[j][0]  = sqrt_N2*(DD2[j][0]+DD2[j][1])+sqrt_N2_2*DD2[j][2]+sqrt_N2_4*DD2[j][4]+sqrt_N2_8*DD2[j][8];
                        groupsh2[j][1]  = sqrt_N2*(DD2[j][0]+DD2[j][1])+sqrt_N2_2*DD2[j][2]+sqrt_N2_4*DD2[j][4]-sqrt_N2_8*DD2[j][8];
                        groupsh2[j][2]  = sqrt_N2*(DD2[j][0]+DD2[j][1])+sqrt_N2_2*DD2[j][2]-sqrt_N2_4*DD2[j][4]+sqrt_N2_8*DD2[j][9];
                        groupsh2[j][3]  = sqrt_N2*(DD2[j][0]+DD2[j][1])+sqrt_N2_2*DD2[j][2]-sqrt_N2_4*DD2[j][4]-sqrt_N2_8*DD2[j][9];
                        groupsh2[j][4]  = sqrt_N2*(DD2[j][0]+DD2[j][1])-sqrt_N2_2*DD2[j][2]+sqrt_N2_4*DD2[j][5]+sqrt_N2_8*DD2[j][10];
                        groupsh2[j][5]  = sqrt_N2*(DD2[j][0]+DD2[j][1])-sqrt_N2_2*DD2[j][2]+sqrt_N2_4*DD2[j][5]-sqrt_N2_8*DD2[j][10];
                        groupsh2[j][6]  = sqrt_N2*(DD2[j][0]+DD2[j][1])-sqrt_N2_2*DD2[j][2]-sqrt_N2_4*DD2[j][5]+sqrt_N2_8*DD2[j][11];
                        groupsh2[j][7]  = sqrt_N2*(DD2[j][0]+DD2[j][1])-sqrt_N2_2*DD2[j][2]-sqrt_N2_4*DD2[j][5]-sqrt_N2_8*DD2[j][11];
                        groupsh2[j][8]  = sqrt_N2*(DD2[j][0]-DD2[j][1])+sqrt_N2_2*DD2[j][3]+sqrt_N2_4*DD2[j][6]+sqrt_N2_8*DD2[j][12];
                        groupsh2[j][9]  = sqrt_N2*(DD2[j][0]-DD2[j][1])+sqrt_N2_2*DD2[j][3]+sqrt_N2_4*DD2[j][6]-sqrt_N2_8*DD2[j][12];
                        groupsh2[j][10] = sqrt_N2*(DD2[j][0]-DD2[j][1])+sqrt_N2_2*DD2[j][3]-sqrt_N2_4*DD2[j][6]+sqrt_N2_8*DD2[j][13];
                        groupsh2[j][11] = sqrt_N2*(DD2[j][0]-DD2[j][1])+sqrt_N2_2*DD2[j][3]-sqrt_N2_4*DD2[j][6]-sqrt_N2_8*DD2[j][13];
                        groupsh2[j][12] = sqrt_N2*(DD2[j][0]-DD2[j][1])-sqrt_N2_2*DD2[j][3]+sqrt_N2_4*DD2[j][7]+sqrt_N2_8*DD2[j][14];
                        groupsh2[j][13] = sqrt_N2*(DD2[j][0]-DD2[j][1])-sqrt_N2_2*DD2[j][3]+sqrt_N2_4*DD2[j][7]-sqrt_N2_8*DD2[j][14];
                        groupsh2[j][14] = sqrt_N2*(DD2[j][0]-DD2[j][1])-sqrt_N2_2*DD2[j][3]-sqrt_N2_4*DD2[j][7]+sqrt_N2_8*DD2[j][15];
                        groupsh2[j][15] = sqrt_N2*(DD2[j][0]-DD2[j][1])-sqrt_N2_2*DD2[j][3]-sqrt_N2_4*DD2[j][7]-sqrt_N2_8*DD2[j][15];
                                }
                            }
                
                  
                
                    if (N22 == 64.0)
                          {
                            for ( j = 0; j < N12*N12; j++ )
                         
                                {
                        groupsh[j][0]  = (sqrt_N2_3)*(DD[j][0]+DD[j][1])+(sqrt_N2_1)*DD[j][2]+(sqrt_N2)*DD[j][4]+(sqrt_N2_2)*DD[j][8]+(sqrt_N2_4)*DD[j][16]+(sqrt_N2_8)*DD[j][32];
                        groupsh[j][1]  = (sqrt_N2_3)*(DD[j][0]+DD[j][1])+(sqrt_N2_1)*DD[j][2]+(sqrt_N2)*DD[j][4]+(sqrt_N2_2)*DD[j][8]+(sqrt_N2_4)*DD[j][16]-(sqrt_N2_8)*DD[j][32];
                        groupsh[j][2]  = (sqrt_N2_3)*(DD[j][0]+DD[j][1])+(sqrt_N2_1)*DD[j][2]+(sqrt_N2)*DD[j][4]+(sqrt_N2_2)*DD[j][8]-(sqrt_N2_4)*DD[j][16]+(sqrt_N2_8)*DD[j][33];
                        groupsh[j][3]  = (sqrt_N2_3)*(DD[j][0]+DD[j][1])+(sqrt_N2_1)*DD[j][2]+(sqrt_N2)*DD[j][4]+(sqrt_N2_2)*DD[j][8]-(sqrt_N2_4)*DD[j][16]-(sqrt_N2_8)*DD[j][33];
                        groupsh[j][4]  = (sqrt_N2_3)*(DD[j][0]+DD[j][1])+(sqrt_N2_1)*DD[j][2]+(sqrt_N2)*DD[j][4]-(sqrt_N2_2)*DD[j][8]+(sqrt_N2_4)*DD[j][17]+(sqrt_N2_8)*DD[j][34];
                        groupsh[j][5]  = (sqrt_N2_3)*(DD[j][0]+DD[j][1])+(sqrt_N2_1)*DD[j][2]+(sqrt_N2)*DD[j][4]-(sqrt_N2_2)*DD[j][8]+(sqrt_N2_4)*DD[j][17]-(sqrt_N2_8)*DD[j][34];
                        groupsh[j][6]  = (sqrt_N2_3)*(DD[j][0]+DD[j][1])+(sqrt_N2_1)*DD[j][2]+(sqrt_N2)*DD[j][4]-(sqrt_N2_2)*DD[j][8]-(sqrt_N2_4)*DD[j][17]+(sqrt_N2_8)*DD[j][35];
                        groupsh[j][7]  = (sqrt_N2_3)*(DD[j][0]+DD[j][1])+(sqrt_N2_1)*DD[j][2]+(sqrt_N2)*DD[j][4]-(sqrt_N2_2)*DD[j][8]-(sqrt_N2_4)*DD[j][17]-(sqrt_N2_8)*DD[j][35];
                        groupsh[j][8]  = (sqrt_N2_3)*(DD[j][0]+DD[j][1])+(sqrt_N2_1)*DD[j][2]-(sqrt_N2)*DD[j][4]+(sqrt_N2_2)*DD[j][9]+(sqrt_N2_4)*DD[j][18]+(sqrt_N2_8)*DD[j][36];
                        groupsh[j][9]  = (sqrt_N2_3)*(DD[j][0]+DD[j][1])+(sqrt_N2_1)*DD[j][2]-(sqrt_N2)*DD[j][4]+(sqrt_N2_2)*DD[j][9]+(sqrt_N2_4)*DD[j][18]-(sqrt_N2_8)*DD[j][36];
                        groupsh[j][10] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])+(sqrt_N2_1)*DD[j][2]-(sqrt_N2)*DD[j][4]+(sqrt_N2_2)*DD[j][9]-(sqrt_N2_4)*DD[j][18]+(sqrt_N2_8)*DD[j][37];
                        groupsh[j][11] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])+(sqrt_N2_1)*DD[j][2]-(sqrt_N2)*DD[j][4]+(sqrt_N2_2)*DD[j][9]-(sqrt_N2_4)*DD[j][18]-(sqrt_N2_8)*DD[j][37];
                        groupsh[j][12] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])+(sqrt_N2_1)*DD[j][2]-(sqrt_N2)*DD[j][4]-(sqrt_N2_2)*DD[j][9]+(sqrt_N2_4)*DD[j][19]+(sqrt_N2_8)*DD[j][38];
                        groupsh[j][13] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])+(sqrt_N2_1)*DD[j][2]-(sqrt_N2)*DD[j][4]-(sqrt_N2_2)*DD[j][9]+(sqrt_N2_4)*DD[j][19]-(sqrt_N2_8)*DD[j][38];
                        groupsh[j][14] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])+(sqrt_N2_1)*DD[j][2]-(sqrt_N2)*DD[j][4]-(sqrt_N2_2)*DD[j][9]-(sqrt_N2_4)*DD[j][19]+(sqrt_N2_8)*DD[j][39];
                        groupsh[j][15] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])+(sqrt_N2_1)*DD[j][2]-(sqrt_N2)*DD[j][4]-(sqrt_N2_2)*DD[j][9]-(sqrt_N2_4)*DD[j][19]-(sqrt_N2_8)*DD[j][39];
                        groupsh[j][16] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])-(sqrt_N2_1)*DD[j][2]+(sqrt_N2)*DD[j][5]+(sqrt_N2_2)*DD[j][10]+(sqrt_N2_4)*DD[j][20]+(sqrt_N2_8)*DD[j][40];
                        groupsh[j][17] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])-(sqrt_N2_1)*DD[j][2]+(sqrt_N2)*DD[j][5]+(sqrt_N2_2)*DD[j][10]+(sqrt_N2_4)*DD[j][20]-(sqrt_N2_8)*DD[j][40];
                        groupsh[j][18] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])-(sqrt_N2_1)*DD[j][2]+(sqrt_N2)*DD[j][5]+(sqrt_N2_2)*DD[j][10]-(sqrt_N2_4)*DD[j][20]+(sqrt_N2_8)*DD[j][41];
                        groupsh[j][19] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])-(sqrt_N2_1)*DD[j][2]+(sqrt_N2)*DD[j][5]+(sqrt_N2_2)*DD[j][10]-(sqrt_N2_4)*DD[j][20]-(sqrt_N2_8)*DD[j][41];
                        groupsh[j][20] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])-(sqrt_N2_1)*DD[j][2]+(sqrt_N2)*DD[j][5]-(sqrt_N2_2)*DD[j][10]+(sqrt_N2_4)*DD[j][21]+(sqrt_N2_8)*DD[j][42];
                        groupsh[j][21] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])-(sqrt_N2_1)*DD[j][2]+(sqrt_N2)*DD[j][5]-(sqrt_N2_2)*DD[j][10]+(sqrt_N2_4)*DD[j][21]-(sqrt_N2_8)*DD[j][42];
                        groupsh[j][22] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])-(sqrt_N2_1)*DD[j][2]+(sqrt_N2)*DD[j][5]-(sqrt_N2_2)*DD[j][10]-(sqrt_N2_4)*DD[j][21]+(sqrt_N2_8)*DD[j][43];
                        groupsh[j][23] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])-(sqrt_N2_1)*DD[j][2]+(sqrt_N2)*DD[j][5]-(sqrt_N2_2)*DD[j][10]-(sqrt_N2_4)*DD[j][21]-(sqrt_N2_8)*DD[j][43];
                        groupsh[j][24] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])-(sqrt_N2_1)*DD[j][2]-(sqrt_N2)*DD[j][5]+(sqrt_N2_2)*DD[j][11]+(sqrt_N2_4)*DD[j][22]+(sqrt_N2_8)*DD[j][44];
                        groupsh[j][25] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])-(sqrt_N2_1)*DD[j][2]-(sqrt_N2)*DD[j][5]+(sqrt_N2_2)*DD[j][11]+(sqrt_N2_4)*DD[j][22]-(sqrt_N2_8)*DD[j][44];
                        groupsh[j][26] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])-(sqrt_N2_1)*DD[j][2]-(sqrt_N2)*DD[j][5]+(sqrt_N2_2)*DD[j][11]-(sqrt_N2_4)*DD[j][22]+(sqrt_N2_8)*DD[j][45];
                        groupsh[j][27] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])-(sqrt_N2_1)*DD[j][2]-(sqrt_N2)*DD[j][5]+(sqrt_N2_2)*DD[j][11]-(sqrt_N2_4)*DD[j][22]-(sqrt_N2_8)*DD[j][45];
                        groupsh[j][28] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])-(sqrt_N2_1)*DD[j][2]-(sqrt_N2)*DD[j][5]-(sqrt_N2_2)*DD[j][11]+(sqrt_N2_4)*DD[j][23]+(sqrt_N2_8)*DD[j][46];
                        groupsh[j][29] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])-(sqrt_N2_1)*DD[j][2]-(sqrt_N2)*DD[j][5]-(sqrt_N2_2)*DD[j][11]+(sqrt_N2_4)*DD[j][23]-(sqrt_N2_8)*DD[j][46];
                        groupsh[j][30] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])-(sqrt_N2_1)*DD[j][2]-(sqrt_N2)*DD[j][5]-(sqrt_N2_2)*DD[j][11]-(sqrt_N2_4)*DD[j][23]+(sqrt_N2_8)*DD[j][47];
                        groupsh[j][31] = (sqrt_N2_3)*(DD[j][0]+DD[j][1])-(sqrt_N2_1)*DD[j][2]-(sqrt_N2)*DD[j][5]-(sqrt_N2_2)*DD[j][11]-(sqrt_N2_4)*DD[j][23]-(sqrt_N2_8)*DD[j][47];
                        groupsh[j][32] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])+(sqrt_N2_1)*DD[j][3]+(sqrt_N2)*DD[j][6]+(sqrt_N2_2)*DD[j][12]+(sqrt_N2_4)*DD[j][24]+(sqrt_N2_8)*DD[j][48];
                        groupsh[j][33] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])+(sqrt_N2_1)*DD[j][3]+(sqrt_N2)*DD[j][6]+(sqrt_N2_2)*DD[j][12]+(sqrt_N2_4)*DD[j][24]-(sqrt_N2_8)*DD[j][48];
                        groupsh[j][34] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])+(sqrt_N2_1)*DD[j][3]+(sqrt_N2)*DD[j][6]+(sqrt_N2_2)*DD[j][12]-(sqrt_N2_4)*DD[j][24]+(sqrt_N2_8)*DD[j][49];
                        groupsh[j][35] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])+(sqrt_N2_1)*DD[j][3]+(sqrt_N2)*DD[j][6]+(sqrt_N2_2)*DD[j][12]-(sqrt_N2_4)*DD[j][24]-(sqrt_N2_8)*DD[j][49];
                        groupsh[j][36] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])+(sqrt_N2_1)*DD[j][3]+(sqrt_N2)*DD[j][6]-(sqrt_N2_2)*DD[j][12]+(sqrt_N2_4)*DD[j][25]+(sqrt_N2_8)*DD[j][50];
                        groupsh[j][37] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])+(sqrt_N2_1)*DD[j][3]+(sqrt_N2)*DD[j][6]-(sqrt_N2_2)*DD[j][12]+(sqrt_N2_4)*DD[j][25]-(sqrt_N2_8)*DD[j][50];
                        groupsh[j][38] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])+(sqrt_N2_1)*DD[j][3]+(sqrt_N2)*DD[j][6]-(sqrt_N2_2)*DD[j][12]-(sqrt_N2_4)*DD[j][25]+(sqrt_N2_8)*DD[j][51];
                        groupsh[j][39] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])+(sqrt_N2_1)*DD[j][3]+(sqrt_N2)*DD[j][6]-(sqrt_N2_2)*DD[j][12]-(sqrt_N2_4)*DD[j][25]-(sqrt_N2_8)*DD[j][51];
                        groupsh[j][40] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])+(sqrt_N2_1)*DD[j][3]-(sqrt_N2)*DD[j][6]+(sqrt_N2_2)*DD[j][13]+(sqrt_N2_4)*DD[j][26]+(sqrt_N2_8)*DD[j][52];
                        groupsh[j][41] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])+(sqrt_N2_1)*DD[j][3]-(sqrt_N2)*DD[j][6]+(sqrt_N2_2)*DD[j][13]+(sqrt_N2_4)*DD[j][26]-(sqrt_N2_8)*DD[j][52];
                        groupsh[j][42] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])+(sqrt_N2_1)*DD[j][3]-(sqrt_N2)*DD[j][6]+(sqrt_N2_2)*DD[j][13]-(sqrt_N2_4)*DD[j][26]+(sqrt_N2_8)*DD[j][53];
                        groupsh[j][43] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])+(sqrt_N2_1)*DD[j][3]-(sqrt_N2)*DD[j][6]+(sqrt_N2_2)*DD[j][13]-(sqrt_N2_4)*DD[j][26]-(sqrt_N2_8)*DD[j][53];
                        groupsh[j][44] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])+(sqrt_N2_1)*DD[j][3]-(sqrt_N2)*DD[j][6]-(sqrt_N2_2)*DD[j][13]+(sqrt_N2_4)*DD[j][27]+(sqrt_N2_8)*DD[j][54];
                        groupsh[j][45] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])+(sqrt_N2_1)*DD[j][3]-(sqrt_N2)*DD[j][6]-(sqrt_N2_2)*DD[j][13]+(sqrt_N2_4)*DD[j][27]-(sqrt_N2_8)*DD[j][54];
                        groupsh[j][46] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])+(sqrt_N2_1)*DD[j][3]-(sqrt_N2)*DD[j][6]-(sqrt_N2_2)*DD[j][13]-(sqrt_N2_4)*DD[j][27]+(sqrt_N2_8)*DD[j][55];
                        groupsh[j][47] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])+(sqrt_N2_1)*DD[j][3]-(sqrt_N2)*DD[j][6]-(sqrt_N2_2)*DD[j][13]-(sqrt_N2_4)*DD[j][27]-(sqrt_N2_8)*DD[j][55];
                        groupsh[j][48] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])-(sqrt_N2_1)*DD[j][3]+(sqrt_N2)*DD[j][7]+(sqrt_N2_2)*DD[j][14]+(sqrt_N2_4)*DD[j][28]+(sqrt_N2_8)*DD[j][56];
                        groupsh[j][49] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])-(sqrt_N2_1)*DD[j][3]+(sqrt_N2)*DD[j][7]+(sqrt_N2_2)*DD[j][14]+(sqrt_N2_4)*DD[j][28]-(sqrt_N2_8)*DD[j][56];
                        groupsh[j][50] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])-(sqrt_N2_1)*DD[j][3]+(sqrt_N2)*DD[j][7]+(sqrt_N2_2)*DD[j][14]-(sqrt_N2_4)*DD[j][28]+(sqrt_N2_8)*DD[j][57];
                        groupsh[j][51] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])-(sqrt_N2_1)*DD[j][3]+(sqrt_N2)*DD[j][7]+(sqrt_N2_2)*DD[j][14]-(sqrt_N2_4)*DD[j][28]-(sqrt_N2_8)*DD[j][57];
                        groupsh[j][52] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])-(sqrt_N2_1)*DD[j][3]+(sqrt_N2)*DD[j][7]-(sqrt_N2_2)*DD[j][14]+(sqrt_N2_4)*DD[j][29]+(sqrt_N2_8)*DD[j][58];
                        groupsh[j][53] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])-(sqrt_N2_1)*DD[j][3]+(sqrt_N2)*DD[j][7]-(sqrt_N2_2)*DD[j][14]+(sqrt_N2_4)*DD[j][29]-(sqrt_N2_8)*DD[j][58];
                        groupsh[j][54] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])-(sqrt_N2_1)*DD[j][3]+(sqrt_N2)*DD[j][7]-(sqrt_N2_2)*DD[j][14]-(sqrt_N2_4)*DD[j][29]+(sqrt_N2_8)*DD[j][59];
                        groupsh[j][55] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])-(sqrt_N2_1)*DD[j][3]+(sqrt_N2)*DD[j][7]-(sqrt_N2_2)*DD[j][14]-(sqrt_N2_4)*DD[j][29]-(sqrt_N2_8)*DD[j][59];
                        groupsh[j][56] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])-(sqrt_N2_1)*DD[j][3]-(sqrt_N2)*DD[j][7]+(sqrt_N2_2)*DD[j][15]+(sqrt_N2_4)*DD[j][30]+(sqrt_N2_8)*DD[j][60];
                        groupsh[j][57] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])-(sqrt_N2_1)*DD[j][3]-(sqrt_N2)*DD[j][7]+(sqrt_N2_2)*DD[j][15]+(sqrt_N2_4)*DD[j][30]-(sqrt_N2_8)*DD[j][60];
                        groupsh[j][58] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])-(sqrt_N2_1)*DD[j][3]-(sqrt_N2)*DD[j][7]+(sqrt_N2_2)*DD[j][15]-(sqrt_N2_4)*DD[j][30]+(sqrt_N2_8)*DD[j][61];
                        groupsh[j][59] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])-(sqrt_N2_1)*DD[j][3]-(sqrt_N2)*DD[j][7]+(sqrt_N2_2)*DD[j][15]-(sqrt_N2_4)*DD[j][30]-(sqrt_N2_8)*DD[j][61];
                        groupsh[j][60] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])-(sqrt_N2_1)*DD[j][3]-(sqrt_N2)*DD[j][7]-(sqrt_N2_2)*DD[j][15]+(sqrt_N2_4)*DD[j][31]+(sqrt_N2_8)*DD[j][62];
                        groupsh[j][61] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])-(sqrt_N2_1)*DD[j][3]-(sqrt_N2)*DD[j][7]-(sqrt_N2_2)*DD[j][15]+(sqrt_N2_4)*DD[j][31]-(sqrt_N2_8)*DD[j][62];
                        groupsh[j][62] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])-(sqrt_N2_1)*DD[j][3]-(sqrt_N2)*DD[j][7]-(sqrt_N2_2)*DD[j][15]-(sqrt_N2_4)*DD[j][31]+(sqrt_N2_8)*DD[j][63];
                        groupsh[j][63] = (sqrt_N2_3)*(DD[j][0]-DD[j][1])-(sqrt_N2_1)*DD[j][3]-(sqrt_N2)*DD[j][7]-(sqrt_N2_2)*DD[j][15]-(sqrt_N2_4)*DD[j][31]-(sqrt_N2_8)*DD[j][63];  
                       groupsh1[j][0]  = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])+(sqrt_N2_1)*DD1[j][2]+(sqrt_N2)*DD1[j][4]+(sqrt_N2_2)*DD1[j][8]+(sqrt_N2_4)*DD1[j][16]+(sqrt_N2_8)*DD1[j][32];
                       groupsh1[j][1]  = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])+(sqrt_N2_1)*DD1[j][2]+(sqrt_N2)*DD1[j][4]+(sqrt_N2_2)*DD1[j][8]+(sqrt_N2_4)*DD1[j][16]-(sqrt_N2_8)*DD1[j][32];
                       groupsh1[j][2]  = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])+(sqrt_N2_1)*DD1[j][2]+(sqrt_N2)*DD1[j][4]+(sqrt_N2_2)*DD1[j][8]-(sqrt_N2_4)*DD1[j][16]+(sqrt_N2_8)*DD1[j][33];
                       groupsh1[j][3]  = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])+(sqrt_N2_1)*DD1[j][2]+(sqrt_N2)*DD1[j][4]+(sqrt_N2_2)*DD1[j][8]-(sqrt_N2_4)*DD1[j][16]-(sqrt_N2_8)*DD1[j][33];
                       groupsh1[j][4]  = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])+(sqrt_N2_1)*DD1[j][2]+(sqrt_N2)*DD1[j][4]-(sqrt_N2_2)*DD1[j][8]+(sqrt_N2_4)*DD1[j][17]+(sqrt_N2_8)*DD1[j][34];
                       groupsh1[j][5]  = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])+(sqrt_N2_1)*DD1[j][2]+(sqrt_N2)*DD1[j][4]-(sqrt_N2_2)*DD1[j][8]+(sqrt_N2_4)*DD1[j][17]-(sqrt_N2_8)*DD1[j][34];
                       groupsh1[j][6]  = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])+(sqrt_N2_1)*DD1[j][2]+(sqrt_N2)*DD1[j][4]-(sqrt_N2_2)*DD1[j][8]-(sqrt_N2_4)*DD1[j][17]+(sqrt_N2_8)*DD1[j][35];
                       groupsh1[j][7]  = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])+(sqrt_N2_1)*DD1[j][2]+(sqrt_N2)*DD1[j][4]-(sqrt_N2_2)*DD1[j][8]-(sqrt_N2_4)*DD1[j][17]-(sqrt_N2_8)*DD1[j][35];
                       groupsh1[j][8]  = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])+(sqrt_N2_1)*DD1[j][2]-(sqrt_N2)*DD1[j][4]+(sqrt_N2_2)*DD1[j][9]+(sqrt_N2_4)*DD1[j][18]+(sqrt_N2_8)*DD1[j][36];
                       groupsh1[j][9]  = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])+(sqrt_N2_1)*DD1[j][2]-(sqrt_N2)*DD1[j][4]+(sqrt_N2_2)*DD1[j][9]+(sqrt_N2_4)*DD1[j][18]-(sqrt_N2_8)*DD1[j][36];
                       groupsh1[j][10] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])+(sqrt_N2_1)*DD1[j][2]-(sqrt_N2)*DD1[j][4]+(sqrt_N2_2)*DD1[j][9]-(sqrt_N2_4)*DD1[j][18]+(sqrt_N2_8)*DD1[j][37];
                       groupsh1[j][11] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])+(sqrt_N2_1)*DD1[j][2]-(sqrt_N2)*DD1[j][4]+(sqrt_N2_2)*DD1[j][9]-(sqrt_N2_4)*DD1[j][18]-(sqrt_N2_8)*DD1[j][37];
                       groupsh1[j][12] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])+(sqrt_N2_1)*DD1[j][2]-(sqrt_N2)*DD1[j][4]-(sqrt_N2_2)*DD1[j][9]+(sqrt_N2_4)*DD1[j][19]+(sqrt_N2_8)*DD1[j][38];
                       groupsh1[j][13] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])+(sqrt_N2_1)*DD1[j][2]-(sqrt_N2)*DD1[j][4]-(sqrt_N2_2)*DD1[j][9]+(sqrt_N2_4)*DD1[j][19]-(sqrt_N2_8)*DD1[j][38];
                       groupsh1[j][14] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])+(sqrt_N2_1)*DD1[j][2]-(sqrt_N2)*DD1[j][4]-(sqrt_N2_2)*DD1[j][9]-(sqrt_N2_4)*DD1[j][19]+(sqrt_N2_8)*DD1[j][39];
                       groupsh1[j][15] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])+(sqrt_N2_1)*DD1[j][2]-(sqrt_N2)*DD1[j][4]-(sqrt_N2_2)*DD1[j][9]-(sqrt_N2_4)*DD1[j][19]-(sqrt_N2_8)*DD1[j][39];
                       groupsh1[j][16] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])-(sqrt_N2_1)*DD1[j][2]+(sqrt_N2)*DD1[j][5]+(sqrt_N2_2)*DD1[j][10]+(sqrt_N2_4)*DD1[j][20]+(sqrt_N2_8)*DD1[j][40];
                       groupsh1[j][17] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])-(sqrt_N2_1)*DD1[j][2]+(sqrt_N2)*DD1[j][5]+(sqrt_N2_2)*DD1[j][10]+(sqrt_N2_4)*DD1[j][20]-(sqrt_N2_8)*DD1[j][40];
                       groupsh1[j][18] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])-(sqrt_N2_1)*DD1[j][2]+(sqrt_N2)*DD1[j][5]+(sqrt_N2_2)*DD1[j][10]-(sqrt_N2_4)*DD1[j][20]+(sqrt_N2_8)*DD1[j][41];
                       groupsh1[j][19] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])-(sqrt_N2_1)*DD1[j][2]+(sqrt_N2)*DD1[j][5]+(sqrt_N2_2)*DD1[j][10]-(sqrt_N2_4)*DD1[j][20]-(sqrt_N2_8)*DD1[j][41];
                       groupsh1[j][20] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])-(sqrt_N2_1)*DD1[j][2]+(sqrt_N2)*DD1[j][5]-(sqrt_N2_2)*DD1[j][10]+(sqrt_N2_4)*DD1[j][21]+(sqrt_N2_8)*DD1[j][42];
                       groupsh1[j][21] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])-(sqrt_N2_1)*DD1[j][2]+(sqrt_N2)*DD1[j][5]-(sqrt_N2_2)*DD1[j][10]+(sqrt_N2_4)*DD1[j][21]-(sqrt_N2_8)*DD1[j][42];
                       groupsh1[j][22] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])-(sqrt_N2_1)*DD1[j][2]+(sqrt_N2)*DD1[j][5]-(sqrt_N2_2)*DD1[j][10]-(sqrt_N2_4)*DD1[j][21]+(sqrt_N2_8)*DD1[j][43];
                       groupsh1[j][23] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])-(sqrt_N2_1)*DD1[j][2]+(sqrt_N2)*DD1[j][5]-(sqrt_N2_2)*DD1[j][10]-(sqrt_N2_4)*DD1[j][21]-(sqrt_N2_8)*DD1[j][43];
                       groupsh1[j][24] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])-(sqrt_N2_1)*DD1[j][2]-(sqrt_N2)*DD1[j][5]+(sqrt_N2_2)*DD1[j][11]+(sqrt_N2_4)*DD1[j][22]+(sqrt_N2_8)*DD1[j][44];
                       groupsh1[j][25] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])-(sqrt_N2_1)*DD1[j][2]-(sqrt_N2)*DD1[j][5]+(sqrt_N2_2)*DD1[j][11]+(sqrt_N2_4)*DD1[j][22]-(sqrt_N2_8)*DD1[j][44];
                       groupsh1[j][26] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])-(sqrt_N2_1)*DD1[j][2]-(sqrt_N2)*DD1[j][5]+(sqrt_N2_2)*DD1[j][11]-(sqrt_N2_4)*DD1[j][22]+(sqrt_N2_8)*DD1[j][45];
                       groupsh1[j][27] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])-(sqrt_N2_1)*DD1[j][2]-(sqrt_N2)*DD1[j][5]+(sqrt_N2_2)*DD1[j][11]-(sqrt_N2_4)*DD1[j][22]-(sqrt_N2_8)*DD1[j][45];
                       groupsh1[j][28] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])-(sqrt_N2_1)*DD1[j][2]-(sqrt_N2)*DD1[j][5]-(sqrt_N2_2)*DD1[j][11]+(sqrt_N2_4)*DD1[j][23]+(sqrt_N2_8)*DD1[j][46];
                       groupsh1[j][29] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])-(sqrt_N2_1)*DD1[j][2]-(sqrt_N2)*DD1[j][5]-(sqrt_N2_2)*DD1[j][11]+(sqrt_N2_4)*DD1[j][23]-(sqrt_N2_8)*DD1[j][46];
                       groupsh1[j][30] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])-(sqrt_N2_1)*DD1[j][2]-(sqrt_N2)*DD1[j][5]-(sqrt_N2_2)*DD1[j][11]-(sqrt_N2_4)*DD1[j][23]+(sqrt_N2_8)*DD1[j][47];
                       groupsh1[j][31] = (sqrt_N2_3)*(DD1[j][0]+DD1[j][1])-(sqrt_N2_1)*DD1[j][2]-(sqrt_N2)*DD1[j][5]-(sqrt_N2_2)*DD1[j][11]-(sqrt_N2_4)*DD1[j][23]-(sqrt_N2_8)*DD1[j][47];
                       groupsh1[j][32] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])+(sqrt_N2_1)*DD1[j][3]+(sqrt_N2)*DD1[j][6]+(sqrt_N2_2)*DD1[j][12]+(sqrt_N2_4)*DD1[j][24]+(sqrt_N2_8)*DD1[j][48];
                       groupsh1[j][33] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])+(sqrt_N2_1)*DD1[j][3]+(sqrt_N2)*DD1[j][6]+(sqrt_N2_2)*DD1[j][12]+(sqrt_N2_4)*DD1[j][24]-(sqrt_N2_8)*DD1[j][48];
                       groupsh1[j][34] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])+(sqrt_N2_1)*DD1[j][3]+(sqrt_N2)*DD1[j][6]+(sqrt_N2_2)*DD1[j][12]-(sqrt_N2_4)*DD1[j][24]+(sqrt_N2_8)*DD1[j][49];
                       groupsh1[j][35] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])+(sqrt_N2_1)*DD1[j][3]+(sqrt_N2)*DD1[j][6]+(sqrt_N2_2)*DD1[j][12]-(sqrt_N2_4)*DD1[j][24]-(sqrt_N2_8)*DD1[j][49];
                       groupsh1[j][36] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])+(sqrt_N2_1)*DD1[j][3]+(sqrt_N2)*DD1[j][6]-(sqrt_N2_2)*DD1[j][12]+(sqrt_N2_4)*DD1[j][25]+(sqrt_N2_8)*DD1[j][50];
                       groupsh1[j][37] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])+(sqrt_N2_1)*DD1[j][3]+(sqrt_N2)*DD1[j][6]-(sqrt_N2_2)*DD1[j][12]+(sqrt_N2_4)*DD1[j][25]-(sqrt_N2_8)*DD1[j][50];
                       groupsh1[j][38] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])+(sqrt_N2_1)*DD1[j][3]+(sqrt_N2)*DD1[j][6]-(sqrt_N2_2)*DD1[j][12]-(sqrt_N2_4)*DD1[j][25]+(sqrt_N2_8)*DD1[j][51];
                       groupsh1[j][39] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])+(sqrt_N2_1)*DD1[j][3]+(sqrt_N2)*DD1[j][6]-(sqrt_N2_2)*DD1[j][12]-(sqrt_N2_4)*DD1[j][25]-(sqrt_N2_8)*DD1[j][51];
                       groupsh1[j][40] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])+(sqrt_N2_1)*DD1[j][3]-(sqrt_N2)*DD1[j][6]+(sqrt_N2_2)*DD1[j][13]+(sqrt_N2_4)*DD1[j][26]+(sqrt_N2_8)*DD1[j][52];
                       groupsh1[j][41] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])+(sqrt_N2_1)*DD1[j][3]-(sqrt_N2)*DD1[j][6]+(sqrt_N2_2)*DD1[j][13]+(sqrt_N2_4)*DD1[j][26]-(sqrt_N2_8)*DD1[j][52];
                       groupsh1[j][42] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])+(sqrt_N2_1)*DD1[j][3]-(sqrt_N2)*DD1[j][6]+(sqrt_N2_2)*DD1[j][13]-(sqrt_N2_4)*DD1[j][26]+(sqrt_N2_8)*DD1[j][53];
                       groupsh1[j][43] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])+(sqrt_N2_1)*DD1[j][3]-(sqrt_N2)*DD1[j][6]+(sqrt_N2_2)*DD1[j][13]-(sqrt_N2_4)*DD1[j][26]-(sqrt_N2_8)*DD1[j][53];
                       groupsh1[j][44] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])+(sqrt_N2_1)*DD1[j][3]-(sqrt_N2)*DD1[j][6]-(sqrt_N2_2)*DD1[j][13]+(sqrt_N2_4)*DD1[j][27]+(sqrt_N2_8)*DD1[j][54];
                       groupsh1[j][45] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])+(sqrt_N2_1)*DD1[j][3]-(sqrt_N2)*DD1[j][6]-(sqrt_N2_2)*DD1[j][13]+(sqrt_N2_4)*DD1[j][27]-(sqrt_N2_8)*DD1[j][54];
                       groupsh1[j][46] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])+(sqrt_N2_1)*DD1[j][3]-(sqrt_N2)*DD1[j][6]-(sqrt_N2_2)*DD1[j][13]-(sqrt_N2_4)*DD1[j][27]+(sqrt_N2_8)*DD1[j][55];
                       groupsh1[j][47] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])+(sqrt_N2_1)*DD1[j][3]-(sqrt_N2)*DD1[j][6]-(sqrt_N2_2)*DD1[j][13]-(sqrt_N2_4)*DD1[j][27]-(sqrt_N2_8)*DD1[j][55];
                       groupsh1[j][48] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])-(sqrt_N2_1)*DD1[j][3]+(sqrt_N2)*DD1[j][7]+(sqrt_N2_2)*DD1[j][14]+(sqrt_N2_4)*DD1[j][28]+(sqrt_N2_8)*DD1[j][56];
                       groupsh1[j][49] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])-(sqrt_N2_1)*DD1[j][3]+(sqrt_N2)*DD1[j][7]+(sqrt_N2_2)*DD1[j][14]+(sqrt_N2_4)*DD1[j][28]-(sqrt_N2_8)*DD1[j][56];
                       groupsh1[j][50] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])-(sqrt_N2_1)*DD1[j][3]+(sqrt_N2)*DD1[j][7]+(sqrt_N2_2)*DD1[j][14]-(sqrt_N2_4)*DD1[j][28]+(sqrt_N2_8)*DD1[j][57];
                       groupsh1[j][51] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])-(sqrt_N2_1)*DD1[j][3]+(sqrt_N2)*DD1[j][7]+(sqrt_N2_2)*DD1[j][14]-(sqrt_N2_4)*DD1[j][28]-(sqrt_N2_8)*DD1[j][57];
                       groupsh1[j][52] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])-(sqrt_N2_1)*DD1[j][3]+(sqrt_N2)*DD1[j][7]-(sqrt_N2_2)*DD1[j][14]+(sqrt_N2_4)*DD1[j][29]+(sqrt_N2_8)*DD1[j][58];
                       groupsh1[j][53] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])-(sqrt_N2_1)*DD1[j][3]+(sqrt_N2)*DD1[j][7]-(sqrt_N2_2)*DD1[j][14]+(sqrt_N2_4)*DD1[j][29]-(sqrt_N2_8)*DD1[j][58];
                       groupsh1[j][54] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])-(sqrt_N2_1)*DD1[j][3]+(sqrt_N2)*DD1[j][7]-(sqrt_N2_2)*DD1[j][14]-(sqrt_N2_4)*DD1[j][29]+(sqrt_N2_8)*DD1[j][59];
                       groupsh1[j][55] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])-(sqrt_N2_1)*DD1[j][3]+(sqrt_N2)*DD1[j][7]-(sqrt_N2_2)*DD1[j][14]-(sqrt_N2_4)*DD1[j][29]-(sqrt_N2_8)*DD1[j][59];
                       groupsh1[j][56] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])-(sqrt_N2_1)*DD1[j][3]-(sqrt_N2)*DD1[j][7]+(sqrt_N2_2)*DD1[j][15]+(sqrt_N2_4)*DD1[j][30]+(sqrt_N2_8)*DD1[j][60];
                       groupsh1[j][57] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])-(sqrt_N2_1)*DD1[j][3]-(sqrt_N2)*DD1[j][7]+(sqrt_N2_2)*DD1[j][15]+(sqrt_N2_4)*DD1[j][30]-(sqrt_N2_8)*DD1[j][60];
                       groupsh1[j][58] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])-(sqrt_N2_1)*DD1[j][3]-(sqrt_N2)*DD1[j][7]+(sqrt_N2_2)*DD1[j][15]-(sqrt_N2_4)*DD1[j][30]+(sqrt_N2_8)*DD1[j][61];
                       groupsh1[j][59] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])-(sqrt_N2_1)*DD1[j][3]-(sqrt_N2)*DD1[j][7]+(sqrt_N2_2)*DD1[j][15]-(sqrt_N2_4)*DD1[j][30]-(sqrt_N2_8)*DD1[j][61];
                       groupsh1[j][60] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])-(sqrt_N2_1)*DD1[j][3]-(sqrt_N2)*DD1[j][7]-(sqrt_N2_2)*DD1[j][15]+(sqrt_N2_4)*DD1[j][31]+(sqrt_N2_8)*DD1[j][62];
                       groupsh1[j][61] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])-(sqrt_N2_1)*DD1[j][3]-(sqrt_N2)*DD1[j][7]-(sqrt_N2_2)*DD1[j][15]+(sqrt_N2_4)*DD1[j][31]-(sqrt_N2_8)*DD1[j][62];
                       groupsh1[j][62] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])-(sqrt_N2_1)*DD1[j][3]-(sqrt_N2)*DD1[j][7]-(sqrt_N2_2)*DD1[j][15]-(sqrt_N2_4)*DD1[j][31]+(sqrt_N2_8)*DD1[j][63];
                       groupsh1[j][63] = (sqrt_N2_3)*(DD1[j][0]-DD1[j][1])-(sqrt_N2_1)*DD1[j][3]-(sqrt_N2)*DD1[j][7]-(sqrt_N2_2)*DD1[j][15]-(sqrt_N2_4)*DD1[j][31]-(sqrt_N2_8)*DD1[j][63]; 
                        groupsh2[j][0]  = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])+(sqrt_N2_1)*DD2[j][2]+(sqrt_N2)*DD2[j][4]+(sqrt_N2_2)*DD2[j][8]+(sqrt_N2_4)*DD2[j][16]+(sqrt_N2_8)*DD2[j][32];
                        groupsh2[j][1]  = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])+(sqrt_N2_1)*DD2[j][2]+(sqrt_N2)*DD2[j][4]+(sqrt_N2_2)*DD2[j][8]+(sqrt_N2_4)*DD2[j][16]-(sqrt_N2_8)*DD2[j][32];
                        groupsh2[j][2]  = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])+(sqrt_N2_1)*DD2[j][2]+(sqrt_N2)*DD2[j][4]+(sqrt_N2_2)*DD2[j][8]-(sqrt_N2_4)*DD2[j][16]+(sqrt_N2_8)*DD2[j][33];
                        groupsh2[j][3]  = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])+(sqrt_N2_1)*DD2[j][2]+(sqrt_N2)*DD2[j][4]+(sqrt_N2_2)*DD2[j][8]-(sqrt_N2_4)*DD2[j][16]-(sqrt_N2_8)*DD2[j][33];
                        groupsh2[j][4]  = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])+(sqrt_N2_1)*DD2[j][2]+(sqrt_N2)*DD2[j][4]-(sqrt_N2_2)*DD2[j][8]+(sqrt_N2_4)*DD2[j][17]+(sqrt_N2_8)*DD2[j][34];
                        groupsh2[j][5]  = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])+(sqrt_N2_1)*DD2[j][2]+(sqrt_N2)*DD2[j][4]-(sqrt_N2_2)*DD2[j][8]+(sqrt_N2_4)*DD2[j][17]-(sqrt_N2_8)*DD2[j][34];
                        groupsh2[j][6]  = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])+(sqrt_N2_1)*DD2[j][2]+(sqrt_N2)*DD2[j][4]-(sqrt_N2_2)*DD2[j][8]-(sqrt_N2_4)*DD2[j][17]+(sqrt_N2_8)*DD2[j][35];
                        groupsh2[j][7]  = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])+(sqrt_N2_1)*DD2[j][2]+(sqrt_N2)*DD2[j][4]-(sqrt_N2_2)*DD2[j][8]-(sqrt_N2_4)*DD2[j][17]-(sqrt_N2_8)*DD2[j][35];
                        groupsh2[j][8]  = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])+(sqrt_N2_1)*DD2[j][2]-(sqrt_N2)*DD2[j][4]+(sqrt_N2_2)*DD2[j][9]+(sqrt_N2_4)*DD2[j][18]+(sqrt_N2_8)*DD2[j][36];
                        groupsh2[j][9]  = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])+(sqrt_N2_1)*DD2[j][2]-(sqrt_N2)*DD2[j][4]+(sqrt_N2_2)*DD2[j][9]+(sqrt_N2_4)*DD2[j][18]-(sqrt_N2_8)*DD2[j][36];
                        groupsh2[j][10] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])+(sqrt_N2_1)*DD2[j][2]-(sqrt_N2)*DD2[j][4]+(sqrt_N2_2)*DD2[j][9]-(sqrt_N2_4)*DD2[j][18]+(sqrt_N2_8)*DD2[j][37];
                        groupsh2[j][11] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])+(sqrt_N2_1)*DD2[j][2]-(sqrt_N2)*DD2[j][4]+(sqrt_N2_2)*DD2[j][9]-(sqrt_N2_4)*DD2[j][18]-(sqrt_N2_8)*DD2[j][37];
                        groupsh2[j][12] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])+(sqrt_N2_1)*DD2[j][2]-(sqrt_N2)*DD2[j][4]-(sqrt_N2_2)*DD2[j][9]+(sqrt_N2_4)*DD2[j][19]+(sqrt_N2_8)*DD2[j][38];
                        groupsh2[j][13] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])+(sqrt_N2_1)*DD2[j][2]-(sqrt_N2)*DD2[j][4]-(sqrt_N2_2)*DD2[j][9]+(sqrt_N2_4)*DD2[j][19]-(sqrt_N2_8)*DD2[j][38];
                        groupsh2[j][14] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])+(sqrt_N2_1)*DD2[j][2]-(sqrt_N2)*DD2[j][4]-(sqrt_N2_2)*DD2[j][9]-(sqrt_N2_4)*DD2[j][19]+(sqrt_N2_8)*DD2[j][39];
                        groupsh2[j][15] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])+(sqrt_N2_1)*DD2[j][2]-(sqrt_N2)*DD2[j][4]-(sqrt_N2_2)*DD2[j][9]-(sqrt_N2_4)*DD2[j][19]-(sqrt_N2_8)*DD2[j][39];
                        groupsh2[j][16] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])-(sqrt_N2_1)*DD2[j][2]+(sqrt_N2)*DD2[j][5]+(sqrt_N2_2)*DD2[j][10]+(sqrt_N2_4)*DD2[j][20]+(sqrt_N2_8)*DD2[j][40];
                        groupsh2[j][17] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])-(sqrt_N2_1)*DD2[j][2]+(sqrt_N2)*DD2[j][5]+(sqrt_N2_2)*DD2[j][10]+(sqrt_N2_4)*DD2[j][20]-(sqrt_N2_8)*DD2[j][40];
                        groupsh2[j][18] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])-(sqrt_N2_1)*DD2[j][2]+(sqrt_N2)*DD2[j][5]+(sqrt_N2_2)*DD2[j][10]-(sqrt_N2_4)*DD2[j][20]+(sqrt_N2_8)*DD2[j][41];
                        groupsh2[j][19] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])-(sqrt_N2_1)*DD2[j][2]+(sqrt_N2)*DD2[j][5]+(sqrt_N2_2)*DD2[j][10]-(sqrt_N2_4)*DD2[j][20]-(sqrt_N2_8)*DD2[j][41];
                        groupsh2[j][20] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])-(sqrt_N2_1)*DD2[j][2]+(sqrt_N2)*DD2[j][5]-(sqrt_N2_2)*DD2[j][10]+(sqrt_N2_4)*DD2[j][21]+(sqrt_N2_8)*DD2[j][42];
                        groupsh2[j][21] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])-(sqrt_N2_1)*DD2[j][2]+(sqrt_N2)*DD2[j][5]-(sqrt_N2_2)*DD2[j][10]+(sqrt_N2_4)*DD2[j][21]-(sqrt_N2_8)*DD2[j][42];
                        groupsh2[j][22] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])-(sqrt_N2_1)*DD2[j][2]+(sqrt_N2)*DD2[j][5]-(sqrt_N2_2)*DD2[j][10]-(sqrt_N2_4)*DD2[j][21]+(sqrt_N2_8)*DD2[j][43];
                        groupsh2[j][23] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])-(sqrt_N2_1)*DD2[j][2]+(sqrt_N2)*DD2[j][5]-(sqrt_N2_2)*DD2[j][10]-(sqrt_N2_4)*DD2[j][21]-(sqrt_N2_8)*DD2[j][43];
                        groupsh2[j][24] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])-(sqrt_N2_1)*DD2[j][2]-(sqrt_N2)*DD2[j][5]+(sqrt_N2_2)*DD2[j][11]+(sqrt_N2_4)*DD2[j][22]+(sqrt_N2_8)*DD2[j][44];
                        groupsh2[j][25] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])-(sqrt_N2_1)*DD2[j][2]-(sqrt_N2)*DD2[j][5]+(sqrt_N2_2)*DD2[j][11]+(sqrt_N2_4)*DD2[j][22]-(sqrt_N2_8)*DD2[j][44];
                        groupsh2[j][26] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])-(sqrt_N2_1)*DD2[j][2]-(sqrt_N2)*DD2[j][5]+(sqrt_N2_2)*DD2[j][11]-(sqrt_N2_4)*DD2[j][22]+(sqrt_N2_8)*DD2[j][45];
                        groupsh2[j][27] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])-(sqrt_N2_1)*DD2[j][2]-(sqrt_N2)*DD2[j][5]+(sqrt_N2_2)*DD2[j][11]-(sqrt_N2_4)*DD2[j][22]-(sqrt_N2_8)*DD2[j][45];
                        groupsh2[j][28] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])-(sqrt_N2_1)*DD2[j][2]-(sqrt_N2)*DD2[j][5]-(sqrt_N2_2)*DD2[j][11]+(sqrt_N2_4)*DD2[j][23]+(sqrt_N2_8)*DD2[j][46];
                        groupsh2[j][29] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])-(sqrt_N2_1)*DD2[j][2]-(sqrt_N2)*DD2[j][5]-(sqrt_N2_2)*DD2[j][11]+(sqrt_N2_4)*DD2[j][23]-(sqrt_N2_8)*DD2[j][46];
                        groupsh2[j][30] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])-(sqrt_N2_1)*DD2[j][2]-(sqrt_N2)*DD2[j][5]-(sqrt_N2_2)*DD2[j][11]-(sqrt_N2_4)*DD2[j][23]+(sqrt_N2_8)*DD2[j][47];
                        groupsh2[j][31] = (sqrt_N2_3)*(DD2[j][0]+DD2[j][1])-(sqrt_N2_1)*DD2[j][2]-(sqrt_N2)*DD2[j][5]-(sqrt_N2_2)*DD2[j][11]-(sqrt_N2_4)*DD2[j][23]-(sqrt_N2_8)*DD2[j][47];
                        groupsh2[j][32] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])+(sqrt_N2_1)*DD2[j][3]+(sqrt_N2)*DD2[j][6]+(sqrt_N2_2)*DD2[j][12]+(sqrt_N2_4)*DD2[j][24]+(sqrt_N2_8)*DD2[j][48];
                        groupsh2[j][33] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])+(sqrt_N2_1)*DD2[j][3]+(sqrt_N2)*DD2[j][6]+(sqrt_N2_2)*DD2[j][12]+(sqrt_N2_4)*DD2[j][24]-(sqrt_N2_8)*DD2[j][48];
                        groupsh2[j][34] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])+(sqrt_N2_1)*DD2[j][3]+(sqrt_N2)*DD2[j][6]+(sqrt_N2_2)*DD2[j][12]-(sqrt_N2_4)*DD2[j][24]+(sqrt_N2_8)*DD2[j][49];
                        groupsh2[j][35] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])+(sqrt_N2_1)*DD2[j][3]+(sqrt_N2)*DD2[j][6]+(sqrt_N2_2)*DD2[j][12]-(sqrt_N2_4)*DD2[j][24]-(sqrt_N2_8)*DD2[j][49];
                        groupsh2[j][36] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])+(sqrt_N2_1)*DD2[j][3]+(sqrt_N2)*DD2[j][6]-(sqrt_N2_2)*DD2[j][12]+(sqrt_N2_4)*DD2[j][25]+(sqrt_N2_8)*DD2[j][50];
                        groupsh2[j][37] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])+(sqrt_N2_1)*DD2[j][3]+(sqrt_N2)*DD2[j][6]-(sqrt_N2_2)*DD2[j][12]+(sqrt_N2_4)*DD2[j][25]-(sqrt_N2_8)*DD2[j][50];
                        groupsh2[j][38] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])+(sqrt_N2_1)*DD2[j][3]+(sqrt_N2)*DD2[j][6]-(sqrt_N2_2)*DD2[j][12]-(sqrt_N2_4)*DD2[j][25]+(sqrt_N2_8)*DD2[j][51];
                        groupsh2[j][39] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])+(sqrt_N2_1)*DD2[j][3]+(sqrt_N2)*DD2[j][6]-(sqrt_N2_2)*DD2[j][12]-(sqrt_N2_4)*DD2[j][25]-(sqrt_N2_8)*DD2[j][51];
                        groupsh2[j][40] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])+(sqrt_N2_1)*DD2[j][3]-(sqrt_N2)*DD2[j][6]+(sqrt_N2_2)*DD2[j][13]+(sqrt_N2_4)*DD2[j][26]+(sqrt_N2_8)*DD2[j][52];
                        groupsh2[j][41] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])+(sqrt_N2_1)*DD2[j][3]-(sqrt_N2)*DD2[j][6]+(sqrt_N2_2)*DD2[j][13]+(sqrt_N2_4)*DD2[j][26]-(sqrt_N2_8)*DD2[j][52];
                        groupsh2[j][42] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])+(sqrt_N2_1)*DD2[j][3]-(sqrt_N2)*DD2[j][6]+(sqrt_N2_2)*DD2[j][13]-(sqrt_N2_4)*DD2[j][26]+(sqrt_N2_8)*DD2[j][53];
                        groupsh2[j][43] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])+(sqrt_N2_1)*DD2[j][3]-(sqrt_N2)*DD2[j][6]+(sqrt_N2_2)*DD2[j][13]-(sqrt_N2_4)*DD2[j][26]-(sqrt_N2_8)*DD2[j][53];
                        groupsh2[j][44] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])+(sqrt_N2_1)*DD2[j][3]-(sqrt_N2)*DD2[j][6]-(sqrt_N2_2)*DD2[j][13]+(sqrt_N2_4)*DD2[j][27]+(sqrt_N2_8)*DD2[j][54];
                        groupsh2[j][45] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])+(sqrt_N2_1)*DD2[j][3]-(sqrt_N2)*DD2[j][6]-(sqrt_N2_2)*DD2[j][13]+(sqrt_N2_4)*DD2[j][27]-(sqrt_N2_8)*DD2[j][54];
                        groupsh2[j][46] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])+(sqrt_N2_1)*DD2[j][3]-(sqrt_N2)*DD2[j][6]-(sqrt_N2_2)*DD2[j][13]-(sqrt_N2_4)*DD2[j][27]+(sqrt_N2_8)*DD2[j][55];
                        groupsh2[j][47] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])+(sqrt_N2_1)*DD2[j][3]-(sqrt_N2)*DD2[j][6]-(sqrt_N2_2)*DD2[j][13]-(sqrt_N2_4)*DD2[j][27]-(sqrt_N2_8)*DD2[j][55];
                        groupsh2[j][48] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])-(sqrt_N2_1)*DD2[j][3]+(sqrt_N2)*DD2[j][7]+(sqrt_N2_2)*DD2[j][14]+(sqrt_N2_4)*DD2[j][28]+(sqrt_N2_8)*DD2[j][56];
                        groupsh2[j][49] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])-(sqrt_N2_1)*DD2[j][3]+(sqrt_N2)*DD2[j][7]+(sqrt_N2_2)*DD2[j][14]+(sqrt_N2_4)*DD2[j][28]-(sqrt_N2_8)*DD2[j][56];
                        groupsh2[j][50] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])-(sqrt_N2_1)*DD2[j][3]+(sqrt_N2)*DD2[j][7]+(sqrt_N2_2)*DD2[j][14]-(sqrt_N2_4)*DD2[j][28]+(sqrt_N2_8)*DD2[j][57];
                        groupsh2[j][51] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])-(sqrt_N2_1)*DD2[j][3]+(sqrt_N2)*DD2[j][7]+(sqrt_N2_2)*DD2[j][14]-(sqrt_N2_4)*DD2[j][28]-(sqrt_N2_8)*DD2[j][57];
                        groupsh2[j][52] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])-(sqrt_N2_1)*DD2[j][3]+(sqrt_N2)*DD2[j][7]-(sqrt_N2_2)*DD2[j][14]+(sqrt_N2_4)*DD2[j][29]+(sqrt_N2_8)*DD2[j][58];
                        groupsh2[j][53] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])-(sqrt_N2_1)*DD2[j][3]+(sqrt_N2)*DD2[j][7]-(sqrt_N2_2)*DD2[j][14]+(sqrt_N2_4)*DD2[j][29]-(sqrt_N2_8)*DD2[j][58];
                        groupsh2[j][54] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])-(sqrt_N2_1)*DD2[j][3]+(sqrt_N2)*DD2[j][7]-(sqrt_N2_2)*DD2[j][14]-(sqrt_N2_4)*DD2[j][29]+(sqrt_N2_8)*DD2[j][59];
                        groupsh2[j][55] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])-(sqrt_N2_1)*DD2[j][3]+(sqrt_N2)*DD2[j][7]-(sqrt_N2_2)*DD2[j][14]-(sqrt_N2_4)*DD2[j][29]-(sqrt_N2_8)*DD2[j][59];
                        groupsh2[j][56] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])-(sqrt_N2_1)*DD2[j][3]-(sqrt_N2)*DD2[j][7]+(sqrt_N2_2)*DD2[j][15]+(sqrt_N2_4)*DD2[j][30]+(sqrt_N2_8)*DD2[j][60];
                        groupsh2[j][57] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])-(sqrt_N2_1)*DD2[j][3]-(sqrt_N2)*DD2[j][7]+(sqrt_N2_2)*DD2[j][15]+(sqrt_N2_4)*DD2[j][30]-(sqrt_N2_8)*DD2[j][60];
                        groupsh2[j][58] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])-(sqrt_N2_1)*DD2[j][3]-(sqrt_N2)*DD2[j][7]+(sqrt_N2_2)*DD2[j][15]-(sqrt_N2_4)*DD2[j][30]+(sqrt_N2_8)*DD2[j][61];
                        groupsh2[j][59] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])-(sqrt_N2_1)*DD2[j][3]-(sqrt_N2)*DD2[j][7]+(sqrt_N2_2)*DD2[j][15]-(sqrt_N2_4)*DD2[j][30]-(sqrt_N2_8)*DD2[j][61];
                        groupsh2[j][60] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])-(sqrt_N2_1)*DD2[j][3]-(sqrt_N2)*DD2[j][7]-(sqrt_N2_2)*DD2[j][15]+(sqrt_N2_4)*DD2[j][31]+(sqrt_N2_8)*DD2[j][62];
                        groupsh2[j][61] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])-(sqrt_N2_1)*DD2[j][3]-(sqrt_N2)*DD2[j][7]-(sqrt_N2_2)*DD2[j][15]+(sqrt_N2_4)*DD2[j][31]-(sqrt_N2_8)*DD2[j][62];
                        groupsh2[j][62] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])-(sqrt_N2_1)*DD2[j][3]-(sqrt_N2)*DD2[j][7]-(sqrt_N2_2)*DD2[j][15]-(sqrt_N2_4)*DD2[j][31]+(sqrt_N2_8)*DD2[j][63];
                        groupsh2[j][63] = (sqrt_N2_3)*(DD2[j][0]-DD2[j][1])-(sqrt_N2_1)*DD2[j][3]-(sqrt_N2)*DD2[j][7]-(sqrt_N2_2)*DD2[j][15]-(sqrt_N2_4)*DD2[j][31]-(sqrt_N2_8)*DD2[j][63]; 
                        
                            }
                       }                                   
                   
     
         for ( i = 0; i < N22; i++ )
                 for ( k = 0; k < N12; k++ )
                           for ( j = 0; j < N12; j++ )
                                {
                                   group[k][j][i] = groupsh[k*N12+j][i]/ WW[k*N12+j][i];              
                                   group1[k][j][i] = groupsh1[k*N12+j][i]/ WW1[k*N12+j][i]; 
                                   group2[k][j][i] = groupsh2[k*N12+j][i]/ WW2[k*N12+j][i]; 
                                 }
         
                
         
                for ( k = 0; k < N_num; k++ )
                     for ( i = Cood[k][1]; i <= Cood[k][1]+N12-1; i++ )
                            for ( j = Cood[k][2]; j <= Cood[k][2]+N12-1; j++ )
                                {
                                 im_r1[i][j] += group[i-Cood[k][1]][j-Cood[k][2]][k]*weightb;
                                 im_r2[i][j] += group1[i-Cood[k][1]][j-Cood[k][2]][k]*weightb1;
                                 im_r3[i][j] += group2[i-Cood[k][1]][j-Cood[k][2]][k]*weightb2;
                                 W[i][j]   += weightb;
                                 W1[i][j]   += weightb1;
                                 W2[i][j]   += weightb2;
                                 } 
                
             
    
    ///////////////////////////////////////////////////////////////////////////////////
    
    
    
        }

     for ( i = Ns1; i < M+Ns1; i++ )
          for ( j = Ns1; j < N+Ns1; j++ )
                         {
              *(imr1+(j-Ns1)*M+(i-Ns1)) = im_r1[i][j] / W[i][j];
               *(imr2+(j-Ns1)*M+(i-Ns1)) = im_r2[i][j] / W1[i][j];
                *(imr3+(j-Ns1)*M+(i-Ns1)) = im_r3[i][j] / W2[i][j];
          }
      

    
   // 矩阵空间释放
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_p1s[i];
    delete [] im_p1s;
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_p2s[i];
    delete [] im_p2s;
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_p3s[i];
    delete [] im_p3s;
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_o1s[i];
    delete [] im_o1s;
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_o2s[i];
    delete [] im_o2s;
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_o3s[i];
    delete [] im_o3s;
    
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] W[i];
    delete [] W;
     for ( i = 0; i < M+2*Ns1; i++ )   delete [] W1[i];
    delete [] W1;
     for ( i = 0; i < M+2*Ns1; i++ )   delete [] W2[i];
    delete [] W2;
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_r1[i];
    delete [] im_r1;
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_r2[i];
    delete [] im_r2;
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_r3[i];
    delete [] im_r3;
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] dp_ps[i];
    delete [] dp_ps;
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_noises[i];
    delete [] im_noises;
   
    for ( i = 0; i < N12*N12; i++ )   delete [] dis3[i];
    delete [] dis3; 
    
    for ( i = 0; i < N12; i++ )
        for ( j = 0; j < N12; j++ )
            delete [] group[i][j];
    for ( i = 0; i < N12; i++ )  delete [] group[i];
    delete [] group;
     for ( i = 0; i < N12; i++ )
        for ( j = 0; j < N12; j++ )
            delete [] group1[i][j];
    for ( i = 0; i < N12; i++ )  delete [] group1[i];
    delete [] group1;
    for ( i = 0; i < N12; i++ )
        for ( j = 0; j < N12; j++ )
            delete [] group2[i][j];
    for ( i = 0; i < N12; i++ )  delete [] group2[i];
    delete [] group2;
    for ( i = 0; i < N12*N12; i++ )   delete [] groupb[i];
    delete [] groupb;  
     for ( i = 0; i < N12*N12; i++ )   delete [] groupb1[i];
    delete [] groupb1;  
     for ( i = 0; i < N12*N12; i++ )   delete [] groupb2[i];
    delete [] groupb2; 
    
    for ( i = 0; i < N12; i++ )
        for ( j = 0; j < N12; j++ )
            delete [] groupw[i][j];
    for ( i = 0; i < N12; i++ )  delete [] groupw[i];
    delete [] groupw;
     
    for ( i = 0; i < N12; i++ )
        for ( j = 0; j < N12; j++ )
            delete [] groupw1[i][j];
    for ( i = 0; i < N12; i++ )  delete [] groupw1[i];
    delete [] groupw1;
    for ( i = 0; i < N12; i++ )
        for ( j = 0; j < N12; j++ )
            delete [] groupw2[i][j];
    for ( i = 0; i < N12; i++ )  delete [] groupw2[i];
    delete [] groupw2;
    
    for ( i = 0; i < N12*N12; i++ )   delete [] groupbw[i];
    delete [] groupbw;   
    
     for ( i = 0; i < N12*N12; i++ )   delete [] groupbw1[i];
    delete [] groupbw1; 
     for ( i = 0; i < N12*N12; i++ )   delete [] groupbw2[i];
    delete [] groupbw2; 
     for ( i = 0; i < N12; i++ )   delete [] im_n[i];
    delete [] im_n;
     for ( i = 0; i < N12; i++ )   delete [] im_n1[i];
    delete [] im_n1;
     for ( i = 0; i < N12; i++ )   delete [] im_n2[i];
    delete [] im_n2;
    
    for ( i = 0; i < N12; i++ )   delete [] im_nw[i];
    delete [] im_nw;
    for ( i = 0; i < N12; i++ )   delete [] im_nw1[i];
    delete [] im_nw1;
    for ( i = 0; i < N12; i++ )   delete [] im_nw2[i];
    delete [] im_nw2;
    
   for ( i = 0; i < N12; i++ )   delete [] im_nn[i];
    delete [] im_nn;
    
      delete [] dis4;
     delete [] dis1Arrays;            
    for ( i = 0; i < N12*N12; i++ )   delete [] DD[i];
    delete [] DD;
      for ( i = 0; i < N12*N12; i++ )   delete [] DD1[i];
    delete [] DD1;
      for ( i = 0; i < N12*N12; i++ )   delete [] DD2[i];
    delete [] DD2;
    for ( i = 0; i < N12*N12; i++ )   delete [] WW[i];
    delete [] WW;   
     for ( i = 0; i < N12*N12; i++ )   delete [] WW1[i];
    delete [] WW1;
     for ( i = 0; i < N12*N12; i++ )   delete [] WW2[i];
    delete [] WW2;
     for ( i = 0; i < N22; i++ )       delete [] Cood[i];
    delete [] Cood;
    
    for ( i = 0; i < N32; i++ )   delete [] groups[i];
    delete [] groups;
    for ( i = 0; i < N32; i++ )   delete [] groupss[i];
    delete [] groupss;
    for ( i = 0; i < N32; i++ )   delete [] groups_t[i];
    delete [] groups_t;
    for ( i = 0; i < N32; i++ )   delete [] groupss_t[i];
    delete [] groupss_t;
    for ( i = 0; i < N32; i++ )   delete [] groupsw[i];
    delete [] groupsw;
    for ( i = 0; i < N32; i++ )   delete [] groupssw[i];
    delete [] groupssw;
    for ( i = 0; i < N32; i++ )   delete [] groups_tw[i];
    delete [] groups_tw;
    for ( i = 0; i < N32; i++ )   delete [] groupss_tw[i];
    delete [] groupss_tw;
    for ( i = 0; i < N32; i++ )   delete [] W_wiener[i];
    delete [] W_wiener;
    
    
     for ( i = 0; i < N32; i++ )   delete [] groups1[i];
    delete [] groups1;
    for ( i = 0; i < N32; i++ )   delete [] groupss1[i];
    delete [] groupss1;
    for ( i = 0; i < N32; i++ )   delete [] groups_t1[i];
    delete [] groups_t1;
    for ( i = 0; i < N32; i++ )   delete [] groupss_t1[i];
    delete [] groupss_t1;
    for ( i = 0; i < N32; i++ )   delete [] groupsw1[i];
    delete [] groupsw1;
    for ( i = 0; i < N32; i++ )   delete [] groupssw1[i];
    delete [] groupssw1;
    for ( i = 0; i < N32; i++ )   delete [] groups_tw1[i];
    delete [] groups_tw1;
    for ( i = 0; i < N32; i++ )   delete [] groupss_tw1[i];
    delete [] groupss_tw1;
    for ( i = 0; i < N32; i++ )   delete [] W_wiener1[i];
    delete [] W_wiener1;
   
     for ( i = 0; i < N32; i++ )   delete [] groups2[i];
    delete [] groups2;
    for ( i = 0; i < N32; i++ )   delete [] groupss2[i];
    delete [] groupss2;
    for ( i = 0; i < N32; i++ )   delete [] groups_t2[i];
    delete [] groups_t2;
    for ( i = 0; i < N32; i++ )   delete [] groupss_t2[i];
    delete [] groupss_t2;
    for ( i = 0; i < N32; i++ )   delete [] groupsw2[i];
    delete [] groupsw2;
    for ( i = 0; i < N32; i++ )   delete [] groupssw2[i];
    delete [] groupssw2;
    for ( i = 0; i < N32; i++ )   delete [] groups_tw2[i];
    delete [] groups_tw2;
    for ( i = 0; i < N32; i++ )   delete [] groupss_tw2[i];
    delete [] groupss_tw2;
    for ( i = 0; i < N32; i++ )   delete [] W_wiener2[i];
    delete [] W_wiener2;
    
    for ( i = 0; i < N22-1; i++ )   delete [] dis2[i];
    delete [] dis2;
    delete [] IndexArray;
   
    for ( i = 0; i < N12*N12; i++ )   delete [] groupsh[i];
    delete [] groupsh;
    for ( i = 0; i < N12*N12; i++ )   delete [] groupsh_t[i];
    delete [] groupsh_t;
    for ( i = 0; i < N12*N12; i++ )   delete [] groupshw[i];
    delete [] groupshw;
    for ( i = 0; i < N12*N12; i++ )   delete [] groupsh_tw[i];
    delete [] groupsh_tw;
    
    for ( i = 0; i < N12*N12; i++ )   delete [] groupsh1[i];
    delete [] groupsh1;
    for ( i = 0; i < N12*N12; i++ )   delete [] groupsh_t1[i];
    delete [] groupsh_t1;
    for ( i = 0; i < N12*N12; i++ )   delete [] groupshw1[i];
    delete [] groupshw1;
    for ( i = 0; i < N12*N12; i++ )   delete [] groupsh_tw1[i];
    delete [] groupsh_tw1;
    
     for ( i = 0; i < N12*N12; i++ )   delete [] groupsh2[i];
    delete [] groupsh2;
    for ( i = 0; i < N12*N12; i++ )   delete [] groupsh_t2[i];
    delete [] groupsh_t2;
    for ( i = 0; i < N12*N12; i++ )   delete [] groupshw2[i];
    delete [] groupshw2;
    for ( i = 0; i < N12*N12; i++ )   delete [] groupsh_tw2[i];
    delete [] groupsh_tw2;
    
    for ( i = 0; i < (Ns1)*(Ns1)-1; i++ )   delete [] dis1[i];
    delete [] dis1;
    
        
    delete [] dis5;
    delete [] dis1Array;
   
    
    delete [] IndexArrays;
   
    
}
















