#include<iostream>
#include"ctc_func.h"
using namespace std;
//function used to arrange dynamic memory for the array and the initialising value is 0
/*void array_construct(float **(&array_input),int n,int m)
{
    array_input=new float*[n];
    for(int i=0;i<n;i++)
    {
        array_input[i]=new float[m];
        for(int j=0;j<m;j++)
        {
        array_input[i][j]=0;
        }
    }
}
*/

//function used to calculate the loss of the CTC
int main()
{
int array_n=4;
int array_m=12;
int seqLen=9;
//initialize the seq lable
int *seq;
int seqini[9]={0,1,1,1,3,0,0,3,2};
seq=new int[9];
for(int i=0;i<9;i++)
    seq[i]=seqini[i];
//initialize the params array_n*array_m
float **params;
float paramsini[4][12]={{5.36795548e-01,-1.54817607e-02,-7.30924569e-01,-1.85322592e-01,-1.01410294e+00,-4.47990250e-01,-2.53889235e+00,-1.52181946e-01,1.73950856e+00,-3.83602525e-01,6.65459487e-01,-1.36707976e+00},{-3.40452685e-01,-4.07361442e-01,1.33850135e-01,-5.97414038e-01,1.12956817e-01,1.31840490e+00,-1.17690008e+00,-1.35829361e+00,1.06414922e-01,-8.01179540e-01,-5.15575275e-01,-8.41503645e-01},{ -4.15641273e-01,1.03007141e+00,1.27524293e+00,-9.64824938e-01,-2.97305065e-02,-5.04219114e-01,8.81265989e-01,7.46957608e-01,-7.91738806e-01,-1.05708677e+00,2.72360956e-01,-1.03370458e+00},{4.19049745e-01,-4.67706590e-01,-7.11855937e-01,2.04669577e-04,1.57870889e+00,-1.42527026e+00,-3.01010860e-01,-1.30156229e+00,2.35249053e-01,-4.15645013e-01,-1.62828244e+00,-1.11468475e+00}};

params=new float*[array_n];
for(int i=0;i<array_n;i++)
{
    params[i]=new float[array_m];
    for(int j=0;j<array_m;j++)
    {
    params[i][j]=paramsini[i][j];
    }
}
///grad calculating module>
float **grad;
//grad calculating module<
//testing the function ctcloss
float llforward,llbackward;
int blank=0;
llforward=ctc_loss(grad,params,array_n,array_m,seq,seqLen,blank,false);
cout<<"###########################"<<endl;
cout<<"the llforward:"<<-llforward<<endl;
cout<<"###########################"<<endl;
cout<<"the grad calculated:"<<endl;
for(int i=0;i<array_n;i++)
{
for(int t=0;t<array_m;t++)
{
    cout<<grad[i][t]<<" ";
}
    cout<<endl;
}
//#######decode procedure#######>
cout<<endl;
vector<int> label_seq;
best_path_decode(params,array_n,array_m,label_seq);
for(int i=1;i<label_seq.size();i++)
cout<<label_seq[i]<<" ";
//#######decode procedure#######<
delete_bidimen_array(params,array_n,array_m);
delete [] seq;
return 0;
}
