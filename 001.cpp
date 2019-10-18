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
#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

//all in c g s
int main()
{
    const double PI = 3.1415926;
    const double g = 9.81;

    double l=1.0;//length of line
    double x;
    double y;//coordinate
    double Vx;
    double Vy;//velocity
    double t;
    double t0=0.0;//initial time
    double tf=20;//final time
    double dt=0.01;//time step

    double theta;//in rad
    double theta0=0.314;//initial theta
    double dtheta_dt;
    double omega;//rad velocity

    omega = sqrt(g/l);
    cout << "omega is " << omega << endl;
    cout << "T is " << 2.0*PI/omega << endl;
    ofstream myfile("SimplePendulum.dat");
    myfile.precision(10);

    t=t0;
    while(t<=tf)
    {
        theta = theta0*cos(omega*(t-t0));
        dtheta_dt = -omega*theta0*sin(omega*(t-t0));
        x = l*sin(theta);
        y = -l*cos(theta);
        Vx = l*dtheta_dt*cos(theta);
        Vy = l*dtheta_dt*sin(theta);
        myfile << t << '\t' << x <<'\t' << y << '\t' << Vx << '\t' << Vy
        << '\t' << theta << '\t' << dtheta_dt << endl;
        t = t+dt;
    }

    return 0;
}


//example 3 平抛

#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

//all in c g s
int main()
{
    const double PI = 3.1415926;
    const double g = 9.81;

    double x0;
    double y0;
    double x;
    double y;//coordinate
    double Vx;
    double Vy;//velocity
    double t;
    double t0=0.0;//initial time
    double tf=1.4416;//final time
    double dt=0.001;//time step

    double theta=PI/4.0;//in rad
    double V0=10.0;
    double V0x=V0*cos(theta);
    double V0y=V0*sin(theta);

    cout << "V0x = " << V0x << "V0y = " << V0y << endl;
    ofstream myfile("pingpao.dat");
    myfile.precision(10);

    t=t0;
    while(t<=tf)
    {
        x=V0x*t;
        y=V0y*t-0.5*g*t*t;
        Vx=V0x;
        Vy=V0y-g*t;
        myfile << t << '\t' << x <<'\t' << y << '\t' << Vx << '\t' << Vy << endl;
        t = t+dt;
    }

    return 0;
}
