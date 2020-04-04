#include <math.h>
#include "mex.h"
#include <stdio.h>

// 声明输入变量
double *N1,*N2,*N3,*Ns,*Nstep,*sigma,*impp;
float *Thr;
double *im_p1,*im_p2,*im_p3;
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
    if ( nrhs != 11) 
		mexErrMsgTxt("11 input arguments required."); 
	if ( nlhs != 3 )
		mexErrMsgTxt("3 output arguments required.");
    // 定义输入变量
    N1       = mxGetPr(prhs[0]);
    N2       = mxGetPr(prhs[1]);
    N3       = mxGetPr(prhs[2]);
    Ns       = mxGetPr(prhs[3]);
    Nstep    = mxGetPr(prhs[4]);
    im_p1     = mxGetPr(prhs[5]);
    M        = mxGetM(prhs[5]);
    N        = mxGetN(prhs[5]);
    im_p2     = mxGetPr(prhs[6]);
    im_p3     = mxGetPr(prhs[7]);
    Thr      = (float*)mxGetPr(prhs[8]);
    sigma    = mxGetPr(prhs[9]);
    impp     = mxGetPr(prhs[10]);
    
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
    float weight;
   // 中间矩阵空间定义
    // im_ps，M行M列
    float** im_ps;
    im_ps = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_ps[i] = new float[N+2*Ns1]; 
    
    // im_ps，M行M列
    float** im_ps1;
    im_ps1 = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_ps1[i] = new float[N+2*Ns1]; 
    
    // im_ps，M行M列
    float** im_ps2;
    im_ps2 = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_ps2[i] = new float[N+2*Ns1]; 
    
     float** im_pp;
    im_pp = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_pp[i] = new float[N+2*Ns1]; 
    
  
    // im_r，M行M列
    float** im_r;
    im_r = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_r[i] = new float[N+2*Ns1];
    
    // im_r，M行M列
    float** im_r1;
    im_r1 = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_r1[i] = new float[N+2*Ns1];
    
    // im_r，M行M列
    float** im_r2;
    im_r2 = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_r2[i] = new float[N+2*Ns1];
    
     // W，M行M列
    float** W;
    W = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   W[i] = new float[N+2*Ns1]; 
    
       
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
    // dis1，Ns*Ns-1行3列
    float** dis1;
    dis1 = new float*[(Ns1)*(Ns1)-1];
    for ( i = 0; i < (Ns1)*(Ns1)-1; i++ )   dis1[i] = new float[3];
    // dis2，N2-1行3列
    float** dis2;
    dis2 = new float*[N22-1];
    for ( i = 0; i < N22-1; i++ )   dis2[i] = new float[3];
    
     
        
	// dis1Array
	float* dis1Array = new float[(Ns1)*(Ns1)-1];
	// IndexArray
    int* IndexArray = new int[N22-1];
    
	// IndexArray
    int* IndexArrays = new int[N32];
    
     //d3
    int* dis5 = new int[N32];
                             
  
                                  
                               
                           
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
    
    float** dis3;
    dis3 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   dis3[i] = new float[N12*N12];
    // dis4，N1*N1行2列
    float* dis4 = new float[N12*N12];
    // dis1Array
	float* dis1Arrays = new float[N12*N12];
    
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
              im_ps[i][j] = 0.0;
              im_ps1[i][j] = 0.0;
              im_ps2[i][j] = 0.0;
              im_pp[i][j] = 0.0;
              im_r[i][j]  = 0.0;
              im_r1[i][j]  = 0.0;
              im_r2[i][j]  = 0.0;
              W[i][j]     = 0.0;
              }
   
       for ( i = Ns1; i < M+Ns1; i++ )
          for ( j = Ns1; j < N+Ns1; j++ )
              {
              im_ps[i][j]  = *(im_p1+(j-Ns1)*M+(i-Ns1));
              im_ps1[i][j] = *(im_p2+(j-Ns1)*M+(i-Ns1));
              im_ps2[i][j] = *(im_p3+(j-Ns1)*M+(i-Ns1));
              im_pp[i][j]  = *(impp+(j-Ns1)*M+(i-Ns1));
              }

  // 周期延拓
    for ( i = 0; i < Ns1; i++ )
          for ( j = 0; j < N+2*Ns1; j++ )
              {
              im_ps[i][j] = im_ps[i+Ns1][j];
              im_ps[i+M+Ns1][j] = im_ps[i+M][j];
              im_ps1[i][j] = im_ps1[i+Ns1][j];
              im_ps1[i+M+Ns1][j] = im_ps1[i+M][j];
              im_ps2[i][j] = im_ps2[i+Ns1][j];
              im_ps2[i+M+Ns1][j] = im_ps2[i+M][j];
              im_pp[i][j] = im_pp[i+Ns1][j];
              im_pp[i+M+Ns1][j] = im_pp[i+M][j];
              }
    for ( i = 0; i < M+2*Ns1; i++ )
          for ( j = 0; j < Ns1; j++ )
              {
              im_ps[i][j] = im_ps[i][j+Ns1];
              im_ps[i][j+N+Ns1] = im_ps[i][j+N];
              im_ps1[i][j] = im_ps1[i][j+Ns1];
              im_ps1[i][j+N+Ns1] = im_ps1[i][j+N];
              im_ps2[i][j] = im_ps2[i][j+Ns1];
              im_ps2[i][j+N+Ns1] = im_ps2[i][j+N];
              im_pp[i][j] = im_pp[i][j+Ns1];
              im_pp[i][j+N+Ns1] = im_pp[i][j+N];
              }
 
        
                         float sqrt_N2_1 =  1.0/sqrt(32.0);
                         float sqrt_N2   =  1.0/sqrt(16.0);
                         float sqrt_N2_2  = 1.0/sqrt(8.0);  
                         float sqrt_N2_4  = 1.0/sqrt(4.0);
                         float sqrt_N2_8  = 1.0/sqrt(2.0);
                         float sqrt_N2_3  = 1.0/sqrt(64.0);
    
                         float weightb;
                                                 
                         float NThr = *Thr*(*sigma);
             
                                             
 // 算法主体
    for ( int a = Ns1; a <= M+Ns1; a += Nstep1 )
        for ( int b = Ns1; b <= N+Ns1; b += Nstep1 )
            {
  ////////////////////////////////////////////////////////////          	
              
     weightb = 0.00000001;           

      
      ///////////////////////////////////////////////////////////               
                
                for ( i = 0; i < N22; i++ )
                    for ( j = 0; j < 3; j++ )
                        Cood[i][j] = 0.0;  
                for ( i = 0; i < N12; i++ )
                    for ( j = 0; j < N12; j++ )
                        for ( k = 0; k < N22; k++ )
                               {
                                group[i][j][k] = 0.0;
                                group1[i][j][k] = 0.0;
                                group2[i][j][k] = 0.0;
                               }
                for ( i = a; i <= a+N12-1; i++ )
                    for ( j = b; j <= b+N12-1; j++ )
                        {
                            im_n[i-a][j-b] = im_ps[i-1][j-1];
                            group[i-a][j-b][0] = im_n[i-a][j-b];
                            im_n1[i-a][j-b] = im_ps1[i-1][j-1];
                            group1[i-a][j-b][0] = im_n1[i-a][j-b];
                            im_n2[i-a][j-b] = im_ps2[i-1][j-1];
                            group2[i-a][j-b][0] = im_n2[i-a][j-b];
                            im_n[i-a][j-b] = im_pp[i-1][j-1];
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
                              minu = im_pp[i][j]-im_n[i-xk][j-yk];
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
							     group[i-x][j-y][k+1]  = im_ps[i][j];
                                 group1[i-x][j-y][k+1] = im_ps1[i][j];
                                 group2[i-x][j-y][k+1] = im_ps2[i][j];
                        Cood[k+1][0] = k+1;
                        Cood[k+1][1] = x;
                        Cood[k+1][2] = y;
                         }
                    } 
                
                      
    for ( j = 0; j < N12*N12; j++ )
         for ( i = 0; i < N22; i++ )
        {
           groupsh[j][i] = 0.0;
           groupsh1[j][i] = 0.0;
           groupsh2[j][i] = 0.0;
           groupsh_t[j][i] = 0.0;
           groupsh_t1[j][i] = 0.0;
           groupsh_t2[j][i] = 0.0;
         }             
 //将每组中的所有块按列扫描给groupb         
                            
                      for ( k = 0; k < N12; k++ )
                           for ( j = 0; j < N12; j++ )
                                for ( i = 0; i < N22; i++ )
                                   {
                                    groupb[k*N12+j][i] = group[k][j][i];
                                     groupb1[k*N12+j][i] = group1[k][j][i];
                                      groupb2[k*N12+j][i] = group2[k][j][i];
                                }
                
                          
     for ( j = 0; j < N12*N12; j++ )
         for ( i = 0; i < N22; i++ )
                {
                 groupsh[j][i] = groupb[j][i];
                 groupsh1[j][i] = groupb1[j][i];
                 groupsh2[j][i] = groupb2[j][i];
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
                                minu = groupb[k][j]-groupb[i][j];
                                sum += minu*minu;
                              }
                             dis3[k][i] = sum;
                           }
       //    printf("%d,%d,%d,%d\n",dis3[0][0],dis3[0][1],dis3[5][9],dis3[10][11]); 
             
   
       // Haar正变换       
       
      
    ///////////////////////////////////////////////////////////////////////////////////////
              
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
      
    
    // 结构硬阈值化
     
   ////////////////////////////////////////////////////////////////////////  
                     for ( i = 0; i < N12*N12; i++ )
                        for ( j = N22/2; j < N22/2+1; j++ )
                             {
                                   groupsh_t[i][j] *= 0.0;
                                   groupsh_t1[i][j] *= 0.0;
                                   groupsh_t2[i][j] *= 0.0;
                             }
   ////////////////////////////////////////////////////////////////////////        
            
                                   
           
     for ( i = 0; i < N12*N12; i++ )
          for ( j = 0; j < N22; j++ ) 
             {     
                 DD[i][j] = 0.0;
                 DD1[i][j] = 0.00000001;
                 DD2[i][j] = 0.00000001;
                 WW[i][j] = 0.00000001;
                 WW1[i][j] = 0.00000001;
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
                  groups1[i][j] = 0.0;
                  groups_t1[i][j] = 0.0;
                  groups2[i][j] = 0.0;
                  groups_t2[i][j] = 0.0;
              }  
           
     
           for (i = 0; i < N12*N12; i++)
                    dis4[i] = 0.0;
               
           for (i = 0; i < N32; i++)
                dis5[i] = 0.0;
     
          for (i = 0; i < N12*N12; i++)
                    dis4[i] = dis3[i][k];
                    
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
                       groups1[i][j] = groupsh_t1[dis5[i]][j];
                       groups2[i][j] = groupsh_t2[dis5[i]][j];
                      }
            
  
           
          // 纵向 Haar 正变换
             
             
                          for ( j = 0; j < N22; j++ )   
                             
                                    {   
                          for ( i = 0; i < N32; i++ )
                              groups_t[0][j] += sqrt_N2_4*groups[i][j];
                          for ( i = 0; i < 2; i++ )
                              groups_t[1][j] += sqrt_N2_4*groups[i][j];
                          for ( i = 2; i < N32; i++ )
                              groups_t[1][j] -= sqrt_N2_4*groups[i][j];
                          for ( l = 0; l < 2; l++ )
                              groups_t[l+2][j] = sqrt_N2_8*groups[2*l][j]-sqrt_N2_8*groups[2*l+1][j];
                                    }
        
    
                          for ( j = 0; j < N22; j++ )   
                             
                                    {   
                          for ( i = 0; i < N32; i++ )
                              groups_t1[0][j] += sqrt_N2_4*groups1[i][j];
                          for ( i = 0; i < 2; i++ )
                              groups_t1[1][j] += sqrt_N2_4*groups1[i][j];
                          for ( i = 2; i < N32; i++ )
                              groups_t1[1][j] -= sqrt_N2_4*groups1[i][j];
                          for ( l = 0; l < 2; l++ )
                              groups_t1[l+2][j] = sqrt_N2_8*groups1[2*l][j]-sqrt_N2_8*groups1[2*l+1][j];
                                    }
         
                          for ( j = 0; j < N22; j++ )   
                             
                                    {   
                          for ( i = 0; i < N32; i++ )
                              groups_t2[0][j] += sqrt_N2_4*groups2[i][j];
                          for ( i = 0; i < 2; i++ )
                              groups_t2[1][j] += sqrt_N2_4*groups2[i][j];
                          for ( i = 2; i < N32; i++ )
                              groups_t2[1][j] -= sqrt_N2_4*groups2[i][j];
                          for ( l = 0; l < 2; l++ )
                              groups_t2[l+2][j] = sqrt_N2_8*groups2[2*l][j]-sqrt_N2_8*groups2[2*l+1][j];
                                    }
      
         
              // 硬阈值化
             
                     for ( i = N32/2; i < N32; i++ )
                        for ( j = 1; j < N22; j++ )
                             {
                                   groups_t[i][j] = 0.0;
                                    groups_t1[i][j] = 0.0;
                                     groups_t2[i][j] = 0.0;
                             }
       
            
                    for ( i = 0; i < N32; i++ )
                        for ( j = 1; j < N22; j++ )  
                            {
                            ( fabs(groups_t[i][j]) <= NThr)?(groups_t[i][j] = 0.0):(groups_t[i][j] = groups_t[i][j]);
             
                            ( fabs(groups_t1[i][j]) <= NThr)?(groups_t1[i][j] = 0.0):(groups_t1[i][j] = groups_t1[i][j]);
    
                           ( fabs(groups_t2[i][j]) <= NThr)?(groups_t2[i][j] = 0.0):(groups_t2[i][j] = groups_t2[i][j]);
                        
                             }
             
                weightb += 1.0;
           

         
        
                            for ( j = 0; j < N22; j++ )
                              {
                        groups[0][j] = sqrt_N2_4*(groups_t[0][j]+groups_t[1][j])+sqrt_N2_8*groups_t[2][j];
                        groups[1][j] = sqrt_N2_4*(groups_t[0][j]+groups_t[1][j])-sqrt_N2_8*groups_t[2][j];               
                        groups[2][j] = sqrt_N2_4*(groups_t[0][j]-groups_t[1][j])+sqrt_N2_8*groups_t[3][j];
                        groups[3][j] = sqrt_N2_4*(groups_t[0][j]-groups_t[1][j])-sqrt_N2_8*groups_t[3][j];
                        groups1[0][j] = sqrt_N2_4*(groups_t1[0][j]+groups_t1[1][j])+sqrt_N2_8*groups_t1[2][j];
                        groups1[1][j] = sqrt_N2_4*(groups_t1[0][j]+groups_t1[1][j])-sqrt_N2_8*groups_t1[2][j];               
                        groups1[2][j] = sqrt_N2_4*(groups_t1[0][j]-groups_t1[1][j])+sqrt_N2_8*groups_t1[3][j];
                        groups1[3][j] = sqrt_N2_4*(groups_t1[0][j]-groups_t1[1][j])-sqrt_N2_8*groups_t1[3][j];
                        groups2[0][j] = sqrt_N2_4*(groups_t2[0][j]+groups_t2[1][j])+sqrt_N2_8*groups_t2[2][j];
                        groups2[1][j] = sqrt_N2_4*(groups_t2[0][j]+groups_t2[1][j])-sqrt_N2_8*groups_t2[2][j];               
                        groups2[2][j] = sqrt_N2_4*(groups_t2[0][j]-groups_t2[1][j])+sqrt_N2_8*groups_t2[3][j];
                        groups2[3][j] = sqrt_N2_4*(groups_t2[0][j]-groups_t2[1][j])-sqrt_N2_8*groups_t2[3][j];
                                }
 
             
               // 小块逆变换结束
                 
                    
   
               
               for ( i = 0; i < N32; i++ )
                    for ( j = 0; j < N22; j++ )
                    {
                      DD[dis5[i]][j] += groups[i][j];
                      DD1[dis5[i]][j] += groups1[i][j];
                      DD2[dis5[i]][j] += groups2[i][j];
                      WW[dis5[i]][j] += 1.0;
                      WW1[dis5[i]][j] += 1.0;
                      WW2[dis5[i]][j] += 1.0;
                    }
              
         }
                
     // 大块逆变换开始   
                 
           
                        
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
                          
   
                
         for ( i = 0; i < N22; i++ )
                 for ( k = 0; k < N12; k++ )
                           for ( j = 0; j < N12; j++ )
                               {
                                   group[k][j][i] = groupsh[k*N12+j][i]/ WW[k*N12+j][i];    
                                   group1[k][j][i] = groupsh1[k*N12+j][i]/ WW1[k*N12+j][i]; 
                                   group2[k][j][i] = groupsh2[k*N12+j][i]/ WW2[k*N12+j][i]; 
                               }
                for ( k = 0; k < N22; k++ )
                     for ( i = Cood[k][1]; i <= Cood[k][1]+N12-1; i++ )
                            for ( j = Cood[k][2]; j <= Cood[k][2]+N12-1; j++ )
                                {
                                 im_r[i][j] += group[i-Cood[k][1]][j-Cood[k][2]][k]*weightb;
                                 im_r1[i][j] += group1[i-Cood[k][1]][j-Cood[k][2]][k]*weightb;
                                 im_r2[i][j] += group2[i-Cood[k][1]][j-Cood[k][2]][k]*weightb;
                                 W[i][j]   += weightb;
                                 } 
         
    ///////////////////////////////////////////////////////////////////////////////////
    
   
       
      } 
                         
                         
    for ( i = Ns1; i < M+Ns1; i++ )
          for ( j = Ns1; j < N+Ns1; j++ )
              {
              *(imr1+(j-Ns1)*M+(i-Ns1)) = im_r[i][j] / W[i][j];
               *(imr2+(j-Ns1)*M+(i-Ns1)) = im_r1[i][j] / W[i][j];
               *(imr3+(j-Ns1)*M+(i-Ns1)) = im_r2[i][j] / W[i][j];
              }
      

    
   // 矩阵空间释放
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_ps[i];
    delete [] im_ps;
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_ps1[i];
    delete [] im_ps1;
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_ps2[i];
    delete [] im_ps2;
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] W[i];
    delete [] W;
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_r[i];
    delete [] im_r;
     for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_r1[i];
    delete [] im_r1;
     for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_r2[i];
    delete [] im_r2;
 
     for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_pp[i];
    delete [] im_pp;
    for ( i = 0; i < N22; i++ )       delete [] Cood[i];
    delete [] Cood;
    
     for ( i = 0; i < N32; i++ )   delete [] groups[i];
    delete [] groups;
     for ( i = 0; i < N32; i++ )   delete [] groups1[i];
    delete [] groups1;
     for ( i = 0; i < N32; i++ )   delete [] groups2[i];
    delete [] groups2;
    for ( i = 0; i < N32; i++ )   delete [] groupss[i];
    delete [] groupss;
    for ( i = 0; i < N32; i++ )   delete [] groupss1[i];
    delete [] groupss1;
    for ( i = 0; i < N32; i++ )   delete [] groupss2[i];
    delete [] groupss2;
    for ( i = 0; i < N32; i++ )   delete [] groups_t[i];
    delete [] groups_t;
     for ( i = 0; i < N32; i++ )   delete [] groups_t1[i];
    delete [] groups_t1;
     for ( i = 0; i < N32; i++ )   delete [] groups_t2[i];
    delete [] groups_t2;
    for ( i = 0; i < N32; i++ )   delete [] groupss_t[i];
    delete [] groupss_t;
     for ( i = 0; i < N32; i++ )   delete [] groupss_t1[i];
    delete [] groupss_t1;
     for ( i = 0; i < N32; i++ )   delete [] groupss_t2[i];
    delete [] groupss_t2;
    for ( i = 0; i < N22-1; i++ )   delete [] dis2[i];
    delete [] dis2;
   
    for ( i = 0; i < (Ns1)*(Ns1)-1; i++ )   delete [] dis1[i];
    delete [] dis1;
            
    delete [] dis5;
    
    delete [] dis1Array;
    
    delete [] IndexArray;
    
    delete [] IndexArrays;
    
    for ( i = 0; i < N12*N12; i++ )   delete [] groupsh[i];
    delete [] groupsh;
   for ( i = 0; i < N12*N12; i++ )   delete [] groupsh1[i];
    delete [] groupsh1;
    for ( i = 0; i < N12*N12; i++ )   delete [] groupsh2[i];
    delete [] groupsh2;
    for ( i = 0; i < N12*N12; i++ )   delete [] groupsh_t[i];
    delete [] groupsh_t;
    
    for ( i = 0; i < N12*N12; i++ )   delete [] groupsh_t1[i];
    delete [] groupsh_t1;
    for ( i = 0; i < N12*N12; i++ )   delete [] groupsh_t2[i];
    delete [] groupsh_t2;
    
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
    
    for ( i = 0; i < N12*N12; i++ )   delete [] dis3[i];
    delete [] dis3;   
     
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
    
    for ( i = 0; i < N12; i++ )   delete [] im_n[i];
    delete [] im_n;
    
    for ( i = 0; i < N12; i++ )   delete [] im_n1[i];
    delete [] im_n1;
    
    for ( i = 0; i < N12; i++ )   delete [] im_n2[i];
    delete [] im_n2;
    
    delete [] dis4;
    delete [] dis1Arrays;            
        
   
    
}
















