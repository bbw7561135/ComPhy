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

    omega = sqrt(g/l);//角频率 而非角速度
    cout << "omega is " << omega << endl;
    cout << "T is " << 2.0*PI/omega << endl;
    ofstream myfile("SimplePendulum.dat");
    myfile.precision(10);

    t=t0;
    while(t<=tf)
    {
        theta = theta0*cos(omega*(t-t0)); //由单摆推导的公式 dx^2/dt^2=-w^2x的结果给出x=Asin(wt+phi)
        //结合初始条件x0=A=lsin(theta0) 以及sin(theta)~theta/l可推导出theta的表达式 微分后得到角速度的表达式
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


#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

//all in c g s
int main()
{
    const double PI = 3.1415926;
    const double g = 9.81;
//consider the resistance of air F=-mkV
    double k=5.0;
    double x0;
    double y0;
    double x;
    double y;//coordinate
    double Vx;
    double Vy;//velocity
    double t;
    double t0=0.0;//initial time
    double tf=0.91;//final time
    double dt=0.001;//time step

    double theta=PI/4.0;//in rad
    double V0=10.0;
    double V0x=V0*cos(theta);
    double V0y=V0*sin(theta);

    cout << "V0x = " << V0x << "V0y = " << V0y << endl;
    ofstream myfile("pingpaoairresistance.dat");
    myfile.precision(10);

    t=t0;
    while(t<=tf)
    {
        x=V0x/k*(1.0-exp(-k*t)); //积分 -mkVx=mdVx/dt
        y=(1.0/k)*(V0y+(g/k))*(1.0-exp(-k*t))-(g/k)*t;//积分 -mkVy=mg=mdVy/dt 利用积分表
        Vx=V0x*exp(-k*t); //积分 x(t)=jifen(Vx(x)dt) Vx表达式在上面
        Vy=(V0y+(g/k))*exp(-k*t)-(g/k);//同理积分
        myfile << t << '\t' << x <<'\t' << y << '\t' << Vx << '\t' << Vy << endl;
        t = t+dt;
    }

    return 0;
}

//各向异性的简谐振荡
//Fx=-mw1^2x Fy=mw2^2y 相当于回复力 F=-kx k等于mw1^2 mw2^2
//令x=Acos(wt+phi) y=Asin(wt+phi) 利用初始条件 x0=A y0=0 Vx0=0 Vy0=w2A 以及求导两次与F比较可得
//xt yt //然后求导可得vxt和yxt

#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

//all in c g s
int main()
{
    const double PI = 3.1415926;

    double x0;
    double y0;
    double x;
    double y;//coordinate
    double Vx;
    double Vy;//velocity
    double t;
    double t0=0.0;//initial time
    double tf=10.0;//final time
    double dt=0.01;//time step
    double A=1;//zhenfu

    double o1=3.0;//omega1
    double o2=5.0;//omega2
    double T1=2*PI/o1;
    double T2=2*PI/o2;

    ofstream myfile("xyjianxiezhendang.dat");
    myfile.precision(10);

    t=t0;
    while(t<=tf)
    {
        x = A*cos(o1*t);
        y = A*sin(o2*t);
        Vx = -o1*A*sin(o1*t);
        Vy = o2*A*cos(o2*t);
        myfile << t << '\t' << x <<'\t' << y << '\t' << Vx << '\t' << Vy << endl;
        t = t+dt;
    }

    return 0;
}


//圆锥摆
//初始条件 t=0时 x=r y=0
//则t时后 x=rcos（wt） y=rsin（wt） z=-lcos（theta）
//利用Tz和Txy的表达式 相除可得cos（theta）的表达式 进而sin（theta）也有了
//再利用Tz表达式将T表示出来 带入到Txy中可得r表达式
//微分x y z 可得速度表达式
//利用cos表达式可得 w>=genhao(g/l)

#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

//all in c g s
int main()
{
    const double PI = 3.1415926;
    const double g=9.81;

    double l=1.0;
    double r;
    double x0;
    double y0;
    double x;
    double y;//coordinate
    double z;
    double Vx;
    double Vy;//velocity
    double Vz;
    double t;
    double t0=0.0;//initial time
    double tf=10.0;//final time
    double dt=0.01;//time step

    double omega=6.28;
    double cos_theta=g/(omega*omega*l);
    double sin_theta=sqrt(1-cos_theta*cos_theta);
    z = -l*cos_theta;
    Vz = 0.0;
    r = (g*sin_theta)/(omega*omega*cos_theta);

    ofstream myfile("xyzyuanzhuibai.dat");
    myfile.precision(10);

    t=t0;
    while(t<=tf)
    {
        x = r*cos(omega*t);
        y = r*sin(omega*t);
        Vx = -r*sin(omega*t)*omega;
        Vy = r*cos(omega*t)*omega;
        myfile << t << '\t' << x <<'\t' << y << '\t'<< z << '\t' << Vx << '\t' << Vy <<'\t'<< Vz << endl;
        t = t+dt;
    }

    return 0;
}


//粒子在磁场中运动
//初始条件见page95
//推导注意利用初始条件 特别是在积分加速度得到速度表达式时注意运动的方向 theta设为与Y轴的夹角 尽管最后表达式中没有theta
//最后令x0=V0y/w是一种选择 为的是让xy平面的圆周运动中心在坐标原点

#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

//all in c g s
int main()
{
    const double PI = 3.1415926;


    double x0;
    double y0;
    double z0;
    double x;
    double y;//coordinate
    double z;
    double Vx;
    double Vy;//velocity
    double Vz;
    double V0=1.0;
    double V0x;
    double V0y;
    double V0z;
    double t;
    double t0=0.0;//initial time
    double tf=10.0;//final time
    double dt=0.01;//time step

    double omega=6.28;
    double theta=PI/9.0;//20 degree in rad

    V0y=V0*cos(theta);
    V0z=V0*sin(theta);
    x0=-V0y/omega;

    ofstream myfile("qvb.dat");
    myfile.precision(10);

    t=t0;
    Vz=V0z;
    while(t<=tf)
    {
        x = x0*cos(omega*t);
        y = -x0*sin(omega*t);
        z = V0z*t;
        Vx = V0y*sin(omega*t);
        Vy = V0y*cos(omega*t);
        myfile << t << '\t' << x <<'\t' << y << '\t'<< z << '\t' << Vx << '\t' << Vy <<'\t'<< Vz << endl;
        t = t+dt;
    }

    return 0;
}

//1d 自由盒 在0-L之间 速度为v 自由移动 在0和L处势能无限大 速度反弹为-v

#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

//all in c g s
int main()
{
    const double PI = 3.1415926;


    double L=10.0;
    double x0=0;
    double v0=1.0;
    double t0=0.0;
    double tf=100.0;
    double dt=0.01;
    double t;
    double x;
    double v;

    ofstream myfile("box1d.dat");
    myfile.precision(10);

    t=t0;
    x=x0;
    v=v0;
    while(t<=tf)
    {
        myfile << t << '\t' << x <<'\t' << v << endl;
        x = x + v*dt;
        t = t+dt;
        if(x<0.0 || x > L)
        {
            v = -v;
        }
    }
    myfile.close();
    return 0;
}


