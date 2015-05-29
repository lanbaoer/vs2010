/*图的邻接矩阵表示*/
#include<iostream>
#include<string>
#include<fstream>
using namespace std;

//因为程序是C++，而CLAPACK是f2c程序转换的C语言版本，所以在此处用extern关键字调用
extern"C"
{
#include <f2c.h>
#include <blaswrap.h>
#include <clapack.h>
}

const int MaxGraphSize=347;
class Graph_Matrix
{
private:
	double * edge;
	double * edge1;
	int graphsize;
	int edgesize;
public:
	Graph_Matrix();
	~Graph_Matrix();
	int GraphEmpty();
	void output();
	int NumberOfVertices();
	int NumberOfEdges();
	//void InsertVertex(int & v);
	void InsertEdge(int & v1,int & v2);
	//void DeleteVertex(int & v);
	//void DeleteEdge(int & v1,int & v2);
	void DegreeCentral();
	void EigenCentral();
	void KatzCentral();
	void PageRank();
	void Betweenness();
	void Closeness();
	void GroupCentral();
	void Floyd();
	//void Tarjan(int & u,int & v);
};
Graph_Matrix::Graph_Matrix()
{
	edge1=new double[MaxGraphSize*MaxGraphSize];
	edge=new double[MaxGraphSize*MaxGraphSize];
	int i;
	for(i=0;i<MaxGraphSize*MaxGraphSize;i++)
		edge[i]=0;
	string readfile="g:\\0.edges";
	ifstream infile(readfile.c_str());
	if(!infile)
	{
		cerr<<"error:unable to open input file:"
			<<readfile<<endl;
	}
	char buf[300];
	char ch=infile.get();
	while(ch=='#')
	{
		infile.getline(buf,300);
		ch=infile.get();
	}
	int x,y;
	x=ch-'0';
	while(!infile.eof())
	{
		infile>>y;
        edge[(x-1)*MaxGraphSize+y-1]=1;
		edge[(y-1)*MaxGraphSize+x-1]=1;
		infile>>x;
	}
	infile.close();
	graphsize=MaxGraphSize;
}
Graph_Matrix::~Graph_Matrix()
{
	delete[] edge;
}
int Graph_Matrix::GraphEmpty()
{
	if(graphsize<=0)
		return 1;
	else
		return 0;
}
int Graph_Matrix::NumberOfEdges()
{
	return edgesize;
}
int Graph_Matrix::NumberOfVertices()
{
	return graphsize;
}
void Graph_Matrix::InsertEdge(int & v1,int & v2)
{
	edge[(v1-1)*MaxGraphSize+v2-1]=1;
	edge[(v2-1)*MaxGraphSize+v1-1]=1;
	edgesize++;
}

void Graph_Matrix::DegreeCentral()
{
	string writefile="g:\\degree.txt";
	ofstream outfile(writefile.c_str());
	int i,j,sum=0;
	for(i=0;i<MaxGraphSize;i++)
	{
		for(j=0;j<MaxGraphSize;j++)
		{
			sum+=edge[i*MaxGraphSize+j];
		}
		outfile<<i<<":"<<sum<<endl;
		sum=0;
	}
	outfile.close();
	cout<<"degree success"<<endl;
}

void Graph_Matrix::GroupCentral()
{

}



void Graph_Matrix::EigenCentral()
{
	int y;
	for(y=0;y<MaxGraphSize*MaxGraphSize;y++)
		edge1[y]=edge[y];
	string writefile="g:\\engin.txt";
	ofstream outfile(writefile.c_str());
	char jobvl = 'N';
    char jobvr = 'V';
    integer n =MaxGraphSize;
    integer lda = n;
    double * wr=new double[MaxGraphSize];
    double * wi=new double[MaxGraphSize];
	integer ldvl = MaxGraphSize;
    double * vl=new double[ldvl*n];
	integer ldvr = MaxGraphSize;
    double * vr=new double[ldvr*n];
	//出现stack flow的问题，一是更改了edge类型为double,二是参数设置规范，改成ldvl*n，不知道怎么就好了,come on.
    integer info;
	integer lwork=4*n;
	double * work=new double[lwork];
	dgeev_(&jobvl,&jobvr,&n,edge1,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);
    if(info==0){
        int i = 0;
        int j = 0;
		double max=wr[0];
		int k=0;
		for(i=0;i<n;i++)
		{
			if(wr[i]>max)
			{
				max=wr[i];
				k=i;
			}
		}
		outfile<<"eigenvalue:"<<wr[k]<<endl;
		outfile<<"eigenvector:"<<endl;
		int a=k*MaxGraphSize;
		for(j=0;j<ldvr;j++)				
				outfile<<vr[a+j]<<endl;
		//千万注意此处如果不输出最大值，而输出所有值时，i要从0开始。从1开始是为了输出最大值。
        /*for(i=0;i<n;i++){
            outfile<<"eigenvalue:"<<i<<"  ";
            outfile<<wr[i]<<"+"<<wi[i]<<"i"<<endl;
            outfile<<"right eigenvector:";
			int a=i*MaxGraphSize;
            for(j=0;j<ldvr;j++){				
				outfile<<vr[a+j]<<"  ";}
            outfile<<endl;
		}*/
		/*	if(wr[i]>max)
			{
				max=wr[i];
				k=i;
			}
        }
		outfile<<"the max of the eigenvalue is:"<<max<<endl;
		outfile<<"the corresponding eigenvector is:"<<endl;
		for(j=0;j<ldvr;j++)
		{
			int a=k*n;
			outfile<<vr[a+j]<<endl;
		}*/
        cout<<"eigen success"<<endl;
    }
	outfile.close();
}

void Graph_Matrix::KatzCentral()
{
	
	string writefile="g:\\katz.txt";
	ofstream outfile(writefile.c_str());
	char jobvl = 'N';
    char jobvr = 'V';
    integer n =MaxGraphSize;
    integer lda = n;
    double * wr=new double[MaxGraphSize];
    double * wi=new double[MaxGraphSize];
	integer ldvl = MaxGraphSize;
    double * vl=new double[ldvl*n];
	integer ldvr = MaxGraphSize;
    double * vr=new double[ldvr*n];
	int y;
	for(y=0;y<MaxGraphSize*MaxGraphSize;y++)
		edge1[y]=edge[y];
	//出现stack flow的问题，一是更改了edge类型为double,二是参数设置规范，改成ldvl*n，不知道怎么就好了,come on.
    integer info;
	integer lwork=4*n;
	double * work=new double[lwork];
	dgeev_(&jobvl,&jobvr,&n,edge1,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);
	 int i = 0;
        int j = 0;
		double max=wr[0];
		int k=0;
    if(info==0){
        for(i=1;i<n;i++)
		{
		if(wr[i]>max)
			{
				max=wr[i];
				k=i;
			}
        }		
    }
	double p=0.02,q=0.02;
	//cout<<max<<endl;
	for(y=0;y<MaxGraphSize*MaxGraphSize;y++)
		edge1[y]=edge[y];
		//cout<<"please input term p and q:"<<endl;
		//cin>>p>>q;
		for(i=0;i<MaxGraphSize;i++)
		{
			for(j=0;j<MaxGraphSize;j++)
			{
			    int b=i*MaxGraphSize+j;
				if(i==j)
					edge1[b]=1-p*edge1[b];
				else
					edge1[b]=(-1)*p*edge1[b];	
			}
		}
	//output();
	//cout<<endl;
	integer t1 = MaxGraphSize;  
    integer t2 = MaxGraphSize;  
    integer t3 = MaxGraphSize;  
    integer t4[MaxGraphSize];  
    integer t5;       
    dgetrf_(&t1,&t2,edge1,&t3,t4,&t5);  
	//cout<<"what"<<endl;
	//output();
    double *t6 = new double[t1]();  
    //求普通矩阵的逆矩阵  
    dgetri_(&t1,edge1,&t3,t4,t6,&t2,&t5);
	//double * t7=new double[MaxGraphSize];
	double sumofrow=0;
	//output();
	for(i=0;i<MaxGraphSize;i++)
		{
			for(j=0;j<MaxGraphSize;j++)
			{
			    sumofrow+=edge1[j+i*MaxGraphSize];
			}
			sumofrow=sumofrow*q;
			outfile<<sumofrow<<endl;
			//cout<<sumofrow<<endl;
			sumofrow=0;
	}
	cout<<"katz success"<<endl;
	outfile.close();
}
void Graph_Matrix::output()
{
	int i=0,j=0;
	for(i=0;i<MaxGraphSize*MaxGraphSize;i++)
	{
		cout<<edge[i]<<" ";
		j++;
		if(j%MaxGraphSize==0)
			cout<<endl;
	}
}

void Graph_Matrix::PageRank()
{
	int y;
	int i,j;
	double * degree=new double[MaxGraphSize];
	double sumdegree=0;
	for(y=0;y<MaxGraphSize*MaxGraphSize;y++)
		edge1[y]=edge[y];
	for(i=0;i<MaxGraphSize;i++)
	{
		for(j=0;j<MaxGraphSize;j++)
		{
			sumdegree+=edge1[i*MaxGraphSize+j];
		}
					
			degree[i]=1/sumdegree;
			sumdegree=0;
	}
	//output();
	string writefile="g:\\pagerank.txt";
	ofstream outfile(writefile.c_str());
	double p=0.3,q=0.3;
	int b;
	double c;
		//cout<<"please input term p and q:"<<endl;
		//cin>>p>>q;
		for(i=0;i<MaxGraphSize;i++)
		{
			for(j=0;j<MaxGraphSize;j++)
			{
				b=i*MaxGraphSize+j;
			    c=p*edge1[b]*degree[j];
				if(i==j)
					edge1[b]=1-c;
				else
					edge1[b]=(-1)*c;
			}
		}
/*		for(i=0;i<MaxGraphSize;i++)
	{
		for(j=0;j<MaxGraphSize;j++)
		{
			cout<<edge1[i*MaxGraphSize+j]<<" ";
		}
					cout<<endl;
	}*/
	//output();
	//cout<<endl;
	integer t1 = MaxGraphSize;  
    integer t2 = MaxGraphSize;  
    integer t3 = MaxGraphSize;  
    integer t4[MaxGraphSize];  
    integer t5;       
    dgetrf_(&t1,&t2,edge1,&t3,t4,&t5);  
	//cout<<"what"<<endl;
	//output();
    double *t6 = new double[t1]();  
    //求普通矩阵的逆矩阵  
    dgetri_(&t1,edge1,&t3,t4,t6,&t2,&t5);
	//double * t7=new double[MaxGraphSize];
	double sumofrow=0;
	//output();
	
	for(i=0;i<MaxGraphSize;i++)
		{
			for(j=0;j<MaxGraphSize;j++)
			{
				//注意此处必须要这样顺序，虽然有点糊涂，但这样是对的。感觉都是行主序就可以，不要考虑什么列主序。
			    sumofrow+=edge1[j+i*MaxGraphSize];
			}
			sumofrow=sumofrow*q;
			outfile<<sumofrow<<endl;
			//cout<<sumofrow<<endl;
			sumofrow=0;
	}
	cout<<"pagerank success"<<endl;
	outfile.close();
}

/*void Graph_Matrix::Betweenness()
{
	string writefile="g:\\between.txt";
	ofstream outfile(writefile.c_str());
	int x,y;
	for(x=0;x<MaxGraphSize;x++)
	for(y=0;y<MaxGraphSize;y++)
	{
		if(edge[x*MaxGraphSize+y]==0&&x!=y)
			edge1[x*MaxGraphSize+y]=1000;
		else
			edge1[x*MaxGraphSize+y]=1;
	}
	int path[MaxGraphSize][MaxGraphSize];
	int parent[MaxGraphSize][MaxGraphSize];
	int cross[MaxGraphSize][MaxGraphSize][MaxGraphSize];
	int i,j,k,n=MaxGraphSize;
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			for(k=0;k<n;k++)
			{
				if(i!=j&&i!=k&&j!=k)
					cross[i][j][k]=0;
				else
					cross[i][j][k]=1;
			}
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			if(edge1[i*MaxGraphSize+j]<1000)
				path[i][j]=1;
		}
	}
		for(k=0;k<n;k++)
		{ 
			for(i=0;i<n;i++)
				for(j=0;j<n;j++)
				{
					int a,b,c;
					a=i*MaxGraphSize+j;
					b=i*MaxGraphSize+k;
					c=k*MaxGraphSize+j;
					if(edge1[a]>(edge1[b]+edge1[c])&&i!=j&&i!=k&&j!=k)
					{
						edge1[a]=edge1[b]+edge1[c];						
						if(k!=parent[i][j])
						{
						 path[i][j]=1;
						// cross[i][j][k]=1;
						}
						parent[i][j]=k;
					}
					else
						if((edge1[a]==(edge1[b]+edge1[c]))&&i!=j&&i!=k&&j!=k)
						{
							//parent应该做成三维的，用来存储所有的前继结点。
							//i到j之间经过的路径数应该等于i到k和k到j之间路径数相乘的最大值
							path[i][j]++;
							//cross[i][j][k]++;
						}
				}
		}
		for(i=0;i<n;i++)
		{
		for(j=0;j<n;j++)
		{
			cout<<path[i][j]<<" ";
		}
		cout<<endl;
		}
		for(i=0;i<n;i++)
		{
		for(j=0;j<n;j++)
		{
			for(k=0;k<n;k++)
			{
				    if(i!=k&&j!=k&&i!=j&&edge1[i*MaxGraphSize+j]==edge1[i*MaxGraphSize+k]+edge1[k*MaxGraphSize+j])
						cross[i][j][k]=path[i][k]*path[k][j];
					//cout<<cross[i][j][k]<<" ";
			}
			//cout<<endl;
		}
		//cout<<endl;
		}
		double * close=new double[MaxGraphSize];
		double sum=0;
		for(k=0;k<MaxGraphSize;k++)
		{
			for(i=0;i<MaxGraphSize;i++)
			{
				for(j=0;j<MaxGraphSize;j++)
				{
					if(i!=k&&j!=k&&i!=j)
						sum+=double(cross[i][j][k])/double(path[i][j]);
				}
			}
			//不乘以2是因为1到3和3到1均已经被加进去了
			close[k]=sum;
			outfile<<close[k]<<endl;
			//cout<<close[k]<<endl;
			sum=0;
		}
		outfile.close();
		cout<<"between success"<<endl;
}*/
void Graph_Matrix::Betweenness()
{
	string writefile="g:\\between.txt";
	ofstream outfile(writefile.c_str());
	int x,y;
	for(x=0;x<MaxGraphSize;x++)
	for(y=0;y<MaxGraphSize;y++)
	{
		if(edge[x*MaxGraphSize+y]==0&&x!=y)
			edge1[x*MaxGraphSize+y]=1000;
		else
			edge1[x*MaxGraphSize+y]=1;
	}
	int path[MaxGraphSize][MaxGraphSize];
	int parent[MaxGraphSize][MaxGraphSize];
	//int cross[MaxGraphSize][MaxGraphSize][MaxGraphSize];
	int i,j,k,n=MaxGraphSize;
	/*for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			for(k=0;k<n;k++)
			{
				if(i!=j&&i!=k&&j!=k)
					cross[i][j][k]=0;
				else
					cross[i][j][k]=1;
			}*/
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			if(edge1[i*MaxGraphSize+j]<1000)
				path[i][j]=1;
		}
	}
		for(k=0;k<n;k++)
		{ 
			for(i=0;i<n;i++)
				for(j=0;j<n;j++)
				{
					int a,b,c;
					a=i*MaxGraphSize+j;
					b=i*MaxGraphSize+k;
					c=k*MaxGraphSize+j;
					if(edge1[a]>(edge1[b]+edge1[c])&&i!=j&&i!=k&&j!=k)
					{
						edge1[a]=edge1[b]+edge1[c];						
						if(k!=parent[i][j])
						{
						 path[i][j]=1;
						// cross[i][j][k]=1;
						}
						parent[i][j]=k;
					}
					else
						if((edge1[a]==(edge1[b]+edge1[c]))&&i!=j&&i!=k&&j!=k)
						{
							//parent应该做成三维的，用来存储所有的前继结点。
							//i到j之间经过的路径数应该等于i到k和k到j之间路径数相乘的最大值
							path[i][j]++;
							//cross[i][j][k]++;
						}
				}
		}
		/*for(i=0;i<n;i++)
		{
		for(j=0;j<n;j++)
		{
			cout<<path[i][j]<<" ";
		}
		cout<<endl;
		}*/
		/*for(i=0;i<n;i++)
		{
		for(j=0;j<n;j++)
		{
			for(k=0;k<n;k++)
			{
				    if(i!=k&&j!=k&&i!=j&&edge1[i*MaxGraphSize+j]==edge1[i*MaxGraphSize+k]+edge1[k*MaxGraphSize+j])
						cross[i][j][k]=path[i][k]*path[k][j];
					//cout<<cross[i][j][k]<<" ";
			}
			//cout<<endl;
		}
		//cout<<endl;
		}*/
		double * close=new double[MaxGraphSize];
		double sum=0;
		for(k=0;k<MaxGraphSize;k++)
		{
			for(i=0;i<MaxGraphSize;i++)
			{
				for(j=0;j<MaxGraphSize;j++)
				{
					if(i!=k&&j!=k&&i!=j&&edge1[i*MaxGraphSize+j]==edge1[i*MaxGraphSize+k]+edge1[k*MaxGraphSize+j])
						sum+=double(path[i][k]*path[k][j])/double(path[i][j]);
				}
			}
			//不乘以2是因为1到3和3到1均已经被加进去了
			close[k]=sum;
			outfile<<close[k]<<endl;
			//cout<<close[k]<<endl;
			sum=0;
		}
		outfile.close();
		cout<<"between success"<<endl;
}

void Graph_Matrix::Closeness()
{
	string writefile="g:\\closeness.txt";
	ofstream outfile(writefile.c_str());
	Floyd();
	double * closeness=new double[MaxGraphSize];
	int n=MaxGraphSize-1;
	int i,j;
	double k=0;
	for(i=0;i<MaxGraphSize;i++)
	{
		for(j=0;j<MaxGraphSize;j++)
		{
			k+=edge1[i*MaxGraphSize+j];
		}
 		closeness[i]=n/k;
		outfile<<closeness[i]<<endl;
		//cout<<closeness[i]<<endl;
		k=0;
	}
	outfile.close();
	cout<<"close success"<<endl;

}

void Graph_Matrix::Floyd()
{
	int x,y;
	for(x=0;x<MaxGraphSize;x++)
	for(y=0;y<MaxGraphSize;y++)
	{
		if(edge[x*MaxGraphSize+y]==0&&x!=y)
			edge1[x*MaxGraphSize+y]=1000;
		else
			edge1[x*MaxGraphSize+y]=edge[x*MaxGraphSize+y];
	}
	int parent[MaxGraphSize][MaxGraphSize];
	int i,j,k,n=MaxGraphSize;
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
				parent[i][j]=-1;
		}
	}
		for(k=0;k<n;k++)
		{ 
			for(i=0;i<n;i++)
				for(j=0;j<n;j++)
				{
					int a,b,c;
					a=i*MaxGraphSize+j;
					b=i*MaxGraphSize+k;
					c=k*MaxGraphSize+j;
					if(edge1[a]>(edge1[b]+edge1[c]))
					{
						edge1[a]=edge1[b]+edge1[c];
						parent[i][j]=k;
					} 
				}
		}
/*		for(i=0;i<MaxGraphSize;i++)
		{
		for(y=0;y<MaxGraphSize;y++)
			cout<<edge1[i*MaxGraphSize+y]<<" ";
		cout<<endl;
		}*/
} 