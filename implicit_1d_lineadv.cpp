#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

//finite-difference implementation of upwind for linear
//advection i.e. a_t+u*a_x=0
//here use implicit methos so we should solve the matrix equation
//use Xushiliang3.1 Gauss xiaoqu method

class gaus
{
    private:
        int n;
        int ghost1;
        double **a,*b;
    public:
        gaus(int nn, int ghost)
        {
            int i;
            n=nn;
            ghost1 = ghost;
            a = new double * [n]; //malloc 2d array
            for(i=0;i<n;i++)
            {
                a[i] = new double [n];
            }
            b = new double [n];
        }
        void input(double bb[], double cc);//input data
        void gauss();//calculate 
        void output(double celldomianup[]);
        ~gaus()
        {
            int i;
            for(i=0;i<n;i++)
            {
                delete[] a[i];
            }
            delete[] a;
            delete[] b;
            cout << "all allocated memory has been released." << endl;
        }
};


void gaus::input(double bb[],double cc)
{
    int i;
    int j;
    for(i=0;i<n;i++)
    {
        b[i] = bb[i];
    }
    for(i=0;i<n;i++)
    for(j=0;j<n;j++)
        {
            a[i][j] = 0.0;
            if(i==j)
            {
                a[i][j] = 1.0 + cc;
                if(i > 0)
                {
                    a[i][i-1] = -cc;
                }
            }
        }
    a[0][n-1] = -cc;

}


void gaus::gauss()
{
	  int *js,l,k,i,j,is;
      double d,t;
      js = new int[n];
      l=1;
      for (k=0; k<=n-2; k++)
      { 
		  d=0.0;
          for (i=k;i<=n-1;i++)
          for (j=k;j<=n-1;j++)
          { 
			  t=fabs(a[i][j]);
              if (t>d) { d=t; js[k]=j; is=i;}
          }
          if (d+1.0==1.0) l=0;
          else
          { if (js[k]!=k)
              for (i=0;i<=n-1;i++)
              { 
                  t=a[i][k]; 
				  a[i][k]=a[i][js[k]]; 
				  a[i][js[k]]=t;
              }
              if (is!=k)
              { 
				  for (j=k;j<=n-1;j++)
                  { 
                      t=a[k][j]; 
					  a[k][j]=a[is][j]; 
					  a[is][j]=t;
                  }
                  t=b[k]; b[k]=b[is]; b[is]=t;
              }
          }
          if (l==0)
          { 
			  delete [] js;
			  cout <<"\n系数矩阵奇异！无解." <<endl;
              return;
          }
          d=a[k][k];
          for (j=k+1;j<=n-1;j++)
              a[k][j]=a[k][j]/d;
          b[k]=b[k]/d;
          for (i=k+1;i<=n-1;i++)
          { 
			  for (j=k+1;j<=n-1;j++)
                  a[i][j]=a[i][j]-a[i][k]*a[k][j];
              b[i]=b[i]-a[i][k]*b[k];
          }
      }
      d=a[n-1][n-1];
      if (fabs(d)+1.0==1.0)
      { 
		  delete [] js;
		  cout <<"\n系数矩阵奇异！无解." <<endl;
          return;
      }
      b[n-1]=b[n-1]/d;
      for (i=n-2;i>=0;i--)
      { 
		  t=0.0;
          for (j=i+1;j<=n-1;j++)
              t=t+a[i][j]*b[j];
          b[i]=b[i]-t;
      }
      js[n-1]=n-1;
      for (k=n-1;k>=0;k--)
        if (js[k]!=k)
        { 
			t=b[k]; b[k]=b[js[k]]; b[js[k]]=t;
		}
    delete [] js;


}


void gaus::output(double celldomianup[])
{
    int i;
    for(i=0;i<=n-1;i++)
    {
        celldomianup[i] = b[i];
    }

}


/////////////////////////////////////////////
int main()
{
    int i;//loop index
    double u = 1.0;// the velocity of physical quantity >0 is go right
    double C = 10.0;//dt<=C*dx/u C is the CFL numbers
    int ncells = 65;//total 65 cells in computational domain
    //here cell in fact means point in the domian not the interval
    int nghost = 1;//each boundary owns 1 //in total 2 ghost cells
    double xmin = 0.0;//left boundary
    double xmax = 1.0;//right boundary
    double dx = (xmax-xmin)/(ncells-1.0);//65 cells = 64 intervals
    //dx the interval between two cells or two points
    double dt = C*dx/u;//time step
    double t_zhouqi = (xmax-xmin)/u;
    double t_tot = 1.0 * t_zhouqi; //one zhou qi
    double t = 0.0;//time
    int total_cells = ncells+2*nghost;
    double cell_pos[total_cells];
    double cell_val[total_cells];
    double cell_val_update[total_cells]={0.0};
    //index of the left ghost cell is from 0 to nghost-1
    //index of the most left cell is nghost
    //index of the most right cell is nghost+ncells-1
    //index of the right ghost cell is from nghost+ncells to total_cells-1
    gaus ct(ncells,nghost);

    //assign the domain i.e. not include the ghost cells
    for(i=0;i<ncells;i++)
    {
        cell_pos[i+nghost] = xmin + dx*i; //position of each cell
    }
    for(i=nghost;i<=nghost+ncells-1;i++)
    { //each cell's physical value
        if(cell_pos[i]<1.0/3.0)
        {
            cell_val[i] = 0.0;
        }
        else if(cell_pos[i]<2.0/3.0)
        {
            cell_val[i] = 1.0;
        }
        else
        {
            cell_val[i] = 0.0;
        }
    }

    //assign the ghost cells //here use periodic boundary conditions
    for(i=0;i<=nghost-1;i++) //left ghost cells
    {
        cell_pos[i] = -1.0; //no need for pos of ghost cells
        cell_val[i] = cell_val[nghost+ncells-1-(nghost-i)];//in case nghost large than 1
    }//use illustrtion to get the expression in []
    for(i=nghost+ncells;i<=total_cells-1;i++)
    {
        cell_pos[i] = -1.0;
        cell_val[i] = cell_val[i-ncells+1];//use illustrtion to get the expression in []
    }
//    for(i=0;i<total_cells;i++) //check for initial profile
//    {
//        cout << cell_pos[i] << '\t' << cell_val[i] <<endl;
//    }


    while(t<t_tot)
    {
        //update the BCs from last timestep
            //assign the ghost cells //here use periodic boundary conditions
        for(i=0;i<=nghost-1;i++) //left ghost cells
        {
            cell_pos[i] = -1.0; //no need for pos of ghost cells
            cell_val[i] = cell_val[nghost+ncells-1-(nghost-i)];//in case nghost large than 1
        }//use illustrtion to get the expression in []
        for(i=nghost+ncells;i<=total_cells-1;i++)
        {
        cell_pos[i] = -1.0;
        cell_val[i] = cell_val[i-ncells+1];//use illustrtion to get the expression in []
        }
        double cell_in[ncells]={0.0};
        double cell_out[ncells]={0.0};
        
        for(i=nghost;i<=nghost+ncells-1;i++)
        {
            cell_in[i-nghost] = cell_val[i];
        }
        ct.input(cell_in,C);
        ct.gauss();
        ct.output(cell_out);
        //update domain using implicit method

        for(i=nghost;i<=nghost+ncells-1;i++)
        {
            cell_val[i] = cell_out[i-nghost];
        }
        t = t + dt;
    }

    ofstream myfile("implicitc10.dat");
    myfile.precision(10);

    for(i=nghost;i<=nghost+ncells-1;i++)
    {
        myfile << cell_pos[i] << '\t' << cell_val[i] << endl;
    }


    myfile.close();
    return 0;
}






































