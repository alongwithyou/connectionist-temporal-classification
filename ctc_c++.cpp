#include<iostream>
#include<cmath>
#include<vector>
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
void sub_max(float **(&params),int array_n,int array_m)
{
	//search the largest value in the column of the array/matrrix
	for(int i=0;i<array_m;i++)
	{
		float colum_max=params[0][i];
		for(int j=0;j<array_n;j++)
		{
			if(params[j][i]>colum_max)
				colum_max=params[j][i];
		}
		for(int j=0;j<array_n;j++)
		{
			params[j][i]=params[j][i]-colum_max;
		}
	}
}
void array_exp(float **(&params),int array_n,int array_m)
{
	for(int i=0;i<array_n;i++)
		for(int j=0;j<array_m;j++)
			params[i][j]=exp(params[i][j]);
}
void array_sum_column(float (**params),float *(&sum_of_column),int array_n,int array_m)
{
	float tempt=0;
	for(int j=0;j<array_m;j++)
	{
		tempt=0;
		for(int i=0;i<array_n;i++)
		{
		tempt+=params[i][j];	
		}
		sum_of_column[j]=tempt;
	}
}
void array_norm(float **(&params),int array_n,int array_m)
{
	for(int j=0;j<array_m;j++)
	{
		float col_sum=0;
		for(int i=0;i<array_n;i++)
			col_sum+=params[i][j];
		for(int i=0;i<array_n;i++)
			params[i][j]=params[i][j]/col_sum;
	}
}

//*******************there may exit algorithm problem*******************************
float array_sum_sincolumn(float **array,int num_of_columns,int t=0,int start=0,int end=0)
{
	float c=0;
	for(int i=start;i<end;i++)
	{
		c+=array[i][t];
	}
	return c;
}
//*********************************************************
void delte_bidimen_array(float **array,int n,int m)
{
	for(int i=0;i<n;i++)
		delete []array[i];
	delete []array;
}
//array_n stands for numPhones,array_m stands for T,2*seqLen+1 for u
float ctc_loss(float **(&grad),float **params,int array_n,int array_m,int *seq,int seqLen,int blank=0,bool is_prob=true)
{
	//1.initializing the value 
	int L=seqLen*2+1;
	int T=array_m;
	
	float **alphas;
	float **beta;
	//2.new the space for alphas and beta
	alphas=new float *[L];
	for(int i=0;i<L;i++)
	{
		alphas[i]=new float[array_m];
		for(int j=0;j<array_m;j++)
		{
			alphas[i][j]=0;
		}
	}

	beta=new float *[L];
	for(int i=0;i<L;i++)
	{
		beta[i]=new float[array_m];
		for(int j=0;j<array_m;j++)
		{
			beta[i][j]=0;	
		}
	}
	//3.transfer the params into prob form if necessary}
	if (is_prob==false)
	{  
		//substracting the max
	    sub_max(params,array_n,array_m);	
		//exp the vals of the array
		array_exp(params,array_n,array_m);
		//normalise the array in the column direction
		array_norm(params,array_n,array_m);
	}
	//4.initialize alphas and the forwardpass
	alphas[0][0]=params[blank][0];
	alphas[1][0]=params[seq[0]][0];
	float c=0;
	//communicate the column num of alphas to sumfunc
	c=array_sum_sincolumn(alphas,L,0,0,L-1);
	alphas[0][0]=params[blank][0]/c;
	alphas[1][0]=params[seq[0]][0]/c;
	static float llforward=log(c);
	//used for test
	int count=0;
	int testnum=0;
	//used for test
	//5.calculating the val of alphas
	for(int t=1;t<T;t++)
	{
		int start=(0>(L-2*(T-t))) ? 0:(L-2*(T-t));
		int end=((2*t+2)>L) ? L:(2*t+2);
		for(int s=start;s<L;s++)
		{
			int l=(s-1)/2;
			if(s%2 == 0)
			{
				if(s==0)
					alphas[s][t]=alphas[s][t-1]*params[blank][t];
				else
					alphas[s][t]=(alphas[s][t-1]+alphas[s-1][t-1])*params[blank][t];
			}
			else if(s==1 || seq[l]==seq[l-1] )
			{
				alphas[s][t]=(alphas[s][t-1]+alphas[s-1][t-1])*params[seq[l]][t];
			}
			else
			{
				alphas[s][t]=(alphas[s][t-1]+alphas[s-1][t-1]+alphas[s-2][t-1])*params[seq[l]][t];
			}
		}
		//used for debug>
		cout<<"col num:"<<count+1<<"***:";
		count++;
		for(int i=0;i<L;i++)
			cout<<alphas[i][t]<<" ";
		//used for debug<
		c=array_sum_sincolumn(alphas,L,t,start,end);
		//used for debug>
		cout<<endl;
		cout<<"######the c:"<<c<<endl;
		//used for debug<
		for(int i=start;i<end;i++)
		{
			alphas[i][t]=alphas[i][t]/c;
		}
		cout<<endl;
     	
		llforward+=log(c);
		//used for debug
		cout<<"llforward:"<<llforward<<endl;
		//used for debug
	}
        //used for debug>
		for(int t=0;t<T;t++)
		{
			for(int i=0;i<L;i++)
				cout<<alphas[i][t]<<" ";
			cout<<endl;
		}
		//used for debug<
	
	//6.initialise the betas and the backward pass
	beta[L-1][T-1]=params[blank][T-1];
	beta[L-2][T-1]=params[seq[seqLen-1]][T-1];
	c=beta[L-1][T-1]+beta[L-2][T-1];
	for(int i=0;i<L;i++)
		beta[i][T-1]=beta[i][T-1]/c;
	static float llbackward=log(c);
	//debug module
	count=T-2;
	//debug module
	//7.calculating the val of the beta
	for(int t=T-2;t>=0;t--)
    {
		int start=(0>=(L-(T-t)*2)) ? 0:(L-(T-t)*2);
		int end=(L<=(2*t+2)) ? L:(2*t+2);
		for(int s=(end-1);s>=0;s--)
		{
			int l=(s-1)/2;
			if(s%2==0)
			{
				if (s==(L-1))
					beta[s][t]=beta[s][t+1]*params[blank][t];
				else
					beta[s][t]=(beta[s][t+1]+beta[s+1][t+1])*params[blank][t];
			}
			else if((s==(L-2))|| (seq[l]==seq[l+1]))
			beta[s][t]=(beta[s][t+1]+beta[s+1][t+1])*params[seq[l]][t];
			else
			beta[s][t]=(beta[s][t+1]+beta[s+1][t+1]+beta[s+2][t+1])*params[seq[l]][t];
		}
		/*used for debug>
		cout<<"col num:"<<count<<"***:";
		count--;
		for(int i=0;i<L;i++)
			cout<<beta[i][t]<<" ";
		cout<<endl;
		*///used for debug<
		c=array_sum_sincolumn(beta,L,t,start,end);
		for(int i=start;i<end;i++)
		{
			beta[i][t]=beta[i][t]/c;
		}
		llbackward+=log(c);
	}
//###############################################################################################
//###calculate the grad#########################################################################

	grad=new float *[array_n];
	for(int i=0;i<array_n;i++)
	{
		grad[i]=new float [T];
		for(int j=0;j<T;j++)
			grad[i][j]=0;
	}
//claculating the alphas*beta>
	float **ab;
	ab=new float *[L];
	for(int i=0;i<L;i++)
	{
		ab[i]=new float [T];
			for(int t=0;t<T;t++)
				ab[i][t]=alphas[i][t]*beta[i][t];
	}
//claculating the alphas*beta<
	for(int s=0;s<L;s++)
	{
		if(s%2==0)
		{
		for(int t=0;t<T;t++)
		{
			grad[blank][t]+=ab[s][t];
			ab[s][t]=ab[s][t]/params[blank][t];
		}
		}
		else
		{
		for(int t=0;t<T;t++)
		{
			grad[seq[(s-1)/2]][t]+=ab[s][t];
			ab[s][t]=ab[s][t]/(params[seq[(s-1)/2]][t]);
		}
		}
	}
	float *absum;
	absum=new float[T];
	array_sum_column(ab,absum,L,T);
//module used to judge the output>
	float llDiff=abs(llforward-llbackward);
	bool absum_0=false;
	for(int t=0;t<T;t++)
		if(absum[t]==0)
		{
			absum_0=true;
		}
//>
//cout<<"the raw grad##############:"<<endl;
//for(int i=0;i<array_n;i++)
//{
//for(int t=0;t<T;t++)
//	cout<<grad[i][t]<<" ";
//cout<<endl;
//}
//<
	if(llDiff>1e-5|| absum_0)
		{
			cout<<"Diff in forward/backward LL:"<<llDiff<<endl;
			cout<<"zero found in absum"<<endl;
			delte_bidimen_array(alphas,L,T);
			delte_bidimen_array(beta,L,T);

			return -llforward;
		}
	else
	{
//module used to judge the output<
		for(int t=0;t<T;t++)
			for(int i=0;i<array_n;i++)
			{
				grad[i][t]=params[i][t]-grad[i][t]/(params[i][t]*absum[i]);
			}
		delte_bidimen_array(alphas,L,T);
		delte_bidimen_array(beta,L,T);
		return -llforward;
	}
}
//here the input of the func is params
void best_path_decode(float **(&params),int array_n,int array_m,vector<int> &seq,int blank=0)
{
	int *seq_decode;
	seq_decode=new int[array_m];
	for(int t=0;t<array_m;t++)
	{
		int max_phone_of_time_t=0;
		float max_val_of_time_t=params[0][t];
		for(int u=0;u<array_n;u++)
		{
		if(params[u][t]>max_val_of_time_t)
		{
			max_val_of_time_t=params[u][t];
			max_phone_of_time_t=u;
		}
		}
		seq_decode[t]=max_phone_of_time_t;
	}
	seq.push_back(0);//push a blank symble and delete it later
	for(int t=0;t<array_m;t++)
	{
		int &tempt=seq.back();
		if((t==0)&&(seq_decode[t]!=0))
			seq.push_back(seq_decode[t]);
        else if(((t>=1)&&(seq_decode[t]!=0))&&(seq_decode[t]!=tempt))
			seq.push_back(seq_decode[t]);
	}
}
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
delte_bidimen_array(params,array_n,array_m);
delete [] seq;
return 0;
}
