//ComPhy Using C++ A book written by KNA from NTU of Athens
//example 1 motion of a particle on the plane 

#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

//all in c g s
int main()
{
    const double PI = 3.1415926;
    double x0=6.0;//yuanxin x
    double y0=6.0;//yuanxin y
    double R=1.0;//radius
    double x;
    double y;
    double Vx;//velocity project in x
    double Vy;//velocity project in y
    double t;
    double t0=0.0;//begin time
    double tf=50.0;//final time
    double dt=0.02;//time step
    double theta;//angle theta determine the position of a particle
    double omega=0.1;//the rad velocity 0.1 here is about 6 degeree per second

    cout << "T = " << 2.0 * PI/omega << endl;
    ofstream myfile("circle.dat");
    myfile.precision(10);
    t=t0;
    while(t<=tf)
    {
        theta = omega*(t-t0);
        x = x0+R*cos(theta);
        y = y0+R*sin(theta);
        Vx = -omega*R*sin(theta);
        Vy = omega*R*cos(theta);
        myfile << t << '\t' << x << '\t' << y << '\t' << Vx << '\t' << Vy << endl;
        t = t + dt;
    }

    return 0;
}

//example 2 单摆

