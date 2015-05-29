// test1.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include"Graph.h"
#include "time.h"
#include<vector>
#include <process.h>

const int X=100;
const int MIN=10;
const int MAX=20;
void gCentrality(int**a)
{
	srand( (unsigned)time( NULL ) ); 
	int vertex=rand()%X;
	int num=10+rand()%MAX;
	vector<int> temp(num,-1);
	vector<int> remain(X-num,-1);
	cout<<"vertex:"<<vertex<<"num:"<<num<<endl;
	int i,j,k;
	bool flag=false;
	int ver=vertex;
	for(i=0;i<num;i++)
	{
		for(j=0;j<X;j++)
		{
			if(a[ver][j]==1)
			{
				for(k=0;k<temp.size();k++)
					if(j==temp[k])
						continue;
				temp.push_back(j);
			}
			if(temp.size()==num)
			{
				flag=true;
				break;
			}
		}
		if(flag==true)
			break;
		ver=temp[i];
	}

	
	for(i=0;i<X;i++)
	{
		bool flag1=true;
		for(j=0;j<temp.size();j++)
		{
			if(i==temp[j])
			{
				flag1=false;
				break;
			}
		}
		if(flag1)
			remain.push_back(i);
	}

	int sum=0;
	
	for(i=0;i<remain.size();i++)
	{
		bool flag2=false;
		for(j=0;j<temp.size();j++)
		{
			if(a[i][j]==1)
			{
				flag2=true;
				break;
			}
		}
		if(flag2)
			sum++;
	}

	cout<<"group:"<<sum<<endl;



}

void sizetest(string readfile)
{
	ifstream infile(readfile.c_str());
	if(!infile)
	{
		cerr<<"error:unable to open input file:"
			<<readfile<<endl;
		//return 0;
	}
	char buf[300];
	char ch=infile.get();
	while(ch=='#')
	{
		infile.getline(buf,300);
		ch=infile.get();
	}
	int x,sum,min;
	x=ch-'0';
	sum=x;
	min=x;
	while(!infile.eof())
	{
		infile>>x;
		if(x>sum)
			sum=x;
		if(x<min)
			min=x;
	}
	infile.close();
	cout<<sum<<endl;
	cout<<min<<endl;
}
const int size=2;
void Silhouette()
{
	double x[size][size]={{-10,-5},{5,10}};
	int i,j,k,p,q;
	double a[size][size],b[size][size],s[size][size];
	int cnum=size;
	int ncnum=size;
	double temp1,temp2;
	double sum1=0,sum2=0;
	double silhouette=0;
	//compute a and b
	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
				for(p=0;p<size;p++)
				{
					if(p!=j)
					{
					temp1=x[i][j]-x[i][p];
					sum1+=temp1*temp1;
					}
				}
				a[i][j]=(1/(cnum-1))*sum1;
				sum1=0;

				for(k=0;k<size;k++)
					if(k!=i)
					{
						for(q=0;q<size;q++)
						{
							temp2=x[i][j]-x[k][q];
							sum2+=temp2*temp2;
						}
					}
				b[i][j]=((double)1/ncnum)*sum2;
				sum2=0;

				double temp3=b[i][j];
				if(a[i][j]>temp3)
					temp3=a[i][j];
				s[i][j]=(b[i][j]-a[i][j])/temp3;
				silhouette+=((double)1/4)*s[i][j];
		}		
	}
	cout<<"silhouette index is:"<<silhouette<<endl;
}
int main()
{
	/*clock_t start,finish;
	start=clock();
	Graph_Matrix graph;
	//sizetest("g:\\0.edges");
	graph.DegreeCentral();
	graph.EigenCentral();
	graph.KatzCentral();
	graph.PageRank();
	graph.Betweenness();
	graph.Closeness();
	//graph.Floyd();
	//Silhouette();
	finish=clock();
	double duration=(double)(finish-start)/CLOCKS_PER_SEC;
	cout<<"time:"<<duration<<"seconds"<<endl;*/
	Graph_Matrix graph;
	system("pause");
	return 0;
}    
