#include<cmath>
#include<cstdio>
#include<iostream>
#include<vector>
#include<windows.h>
#include"eggplot.h"

using namespace std;
using namespace eggp;

int main()
{
    // 设置多孔介质初始物性参数
    /*---------------------------------------------------------------------------
    L0 -- 初始厚度    e0 -- 初始孔隙度    e0x -- 初始时刻x处的孔隙度
    Es -- 多孔介质弹性模量    D0 -- 初始纤维当量直径
    k1 -- 粘性项系数    k2 -- 惯性项系数
    ---------------------------------------------------------------------------*/
    double L0 = 1.0e-3; // m
    double e0 = 0.4;
    double k1 = 150.0;
    double k2 = 1.75;
    double D0 = 6.0e-3; // m
    double Es = 1400; //Pa
    // 设置流体介质物性参数
    /*---------------------------------------------------------------------------
    rhoL -- 流体密度    mu -- 流体动力粘度    u0 -- 流体初始表征速度
    v0 -- 流体初始孔隙速度    v0x -- 初始时刻x处的孔隙速度    ut -- t时刻的表征速度
    g -- 重力加速度    DeltaP -- 压差    Pin -- 入口压力    Pout -- 出口压力
    px -- x处的压力
    ---------------------------------------------------------------------------*/
    double rhoL = 1000.0; // kg/m3
    double mu = 1.01e-3; // Pa.s
    double g = 9.8; // m/s2
    double Pin = 101325.0; // Pa
    double Pout = 71325.0; // Pa
    double DeltaP = Pin - Pout; // Pa
    double v0, u0, vx, ut, px, v1, u1;
    // 设置多孔介质变形后物性参数
    /*---------------------------------------------------------------------------
    L1 -- 平衡厚度    e1 -- 平衡孔隙度    e1x -- 平衡时刻x处的孔隙度
    Es -- 多孔介质弹性模量    D1 -- 平衡纤维当量直径
    ---------------------------------------------------------------------------*/
    double L1 = (1 - DeltaP / Es) * L0; // m 胡克定律
    double e1 = 1 - (1 - e0) * L0 / L1;
    double D1 = D0 * L1 / L0; // m

    // 设置原始Ergun方程计算用中间变量
    /*---------------------------------------------------------------------------
    phi0 -- (1+e0)/e0    A0 -- 原始Ergun方程二次项系数
    B0 -- 原始Ergun方程一次项系数    C0 -- 原始Ergun方程常数项
    ---------------------------------------------------------------------------*/
    double phi0 = (1 + e0) / e0;
    double A0 = k2 * phi0 * rhoL / D0 / g;
    double B0 = k1 * phi0 * phi0 * mu / (D0 * D0) / g;
    double C0 = -DeltaP / L0;

    // 设置修正Ergun方程计算用中间变量
    /*---------------------------------------------------------------------------
    phi1 -- (1+e1)/e1    A1 -- 原始Ergun方程二次项系数
    B1 -- 原始Ergun方程一次项系数    C1 -- 原始Ergun方程常数项
    ---------------------------------------------------------------------------*/
    double phi1 = (1 + e1) / e1;
    double A1 = k2 * phi1 * rhoL / D1 / g;
    double B1 = k1 * phi1 * phi1 * mu / (D1 * D1) / g;
    double C1 = -DeltaP / L1;

    //    // 设置弹性修正Ergun方程中间变量
    //    /*---------------------------------------------------------------------------
    //    phi0 -- (1+e0)/e0    A0 -- 弹性修正Ergun方程二次项系数
    //    B0 -- 弹性修正Ergun方程一次项系数    C0 -- 弹性修正Ergun方程常数项
    //    ---------------------------------------------------------------------------*/
    //    double phi1 = (1 + e0) / (e0 - DeltaP / Es);
    //    double A1 = k2 * phi0 * rhoL / D0 / g;
    //    double B1 = k1 * phi0 * phi0 * mu / (D0 * D0) / g / (1- DeltaP / Es);
    //    double C1 = -DeltaP / L0;

    // 计算原始Ergun方程
    v0 = (-B0 + sqrt(B0 * B0 - 4 * A0 * C0)) / (2 * A0);
    // 计算弹性修正Ergun方程
    v1 = (-B1 + sqrt(B1 * B1 - 4 * A1 * C1)) / (2 * A1);
    // 计算线弹性修正Ergun方程相对变化率
    double vChange = (abs(v1) - v0) / v0 * 100;
    std::cout << "v0:  " << v0 << "\t" << "v1:  " << v1 << "\t" << "m/s" << "\t" << vChange << "\t" << "%" << endl;

    // 原Ergun方程画图
    // 计算Y:v0-X:deltaP-[e0]关系
    int dPNum0; // 压差取点数
    int e0Num; // 孔隙率变化取点数
    dPNum0 = 14;
    e0Num = 4;
    std::vector<std::vector<double> >DeltaPs0(e0Num, std::vector<double>(dPNum0));
    std::vector<std::vector<double> >v0s0(e0Num, std::vector<double>(dPNum0));
    std::vector<double>C0s(dPNum0);
    double* phi0s = new double[e0Num];

    eggp::Eggplot curvePlot0(SCREEN | PNG | EPS | PDF);
    curvePlot0.xlabel("{/Symbol D}P(kPa)");
    curvePlot0.ylabel("v_0(m/s)");
    curvePlot0.grid(true);
    curvePlot0.name("R01");
    curvePlot0.title("v_0 changing with {/Symbol D}P in different {/Symbol e}_0 and rigid body assumption");
    curvePlot0.legend({ "{/Symbol e}_0 = 0.3", "{/Symbol e}_0 = 0.5", "{/Symbol e}_0 = 0.7", "{/Symbol e}_0 = 0.9"});
    
    for (int i = 1; i < e0Num + 1; i++)
    {
        double e0s = 0.1 + 0.2 * i;
        phi0s[i-1] = (1 + e0s) / e0s;
        double A0s = k2 * phi0s[i - 1] * rhoL / D0 / g;
        double B0s = k1 * phi0s[i - 1] * phi0s[i - 1] * mu / (D0 * D0) / g;

        for (int j = 1; j < dPNum0 + 1; j++)
        {
            DeltaPs0[i - 1][j - 1] = (Pin - j * 5000.0)/1000.0; //  转换成kPa
            C0s[j - 1] = -DeltaPs0[i - 1][j - 1] * 1000.0/ L0; //  运算中转换成Pa
            v0s0[i - 1][j - 1] = ((-B0s + sqrt(B0s * B0s - 4 * A0s * C0s[j - 1])) / (2 * A0s));
            std::cout << "e0:  " << e0s << "\t" << "DeltaP " << j << ":  " << DeltaPs0[i - 1][j - 1] << "\t" << "v0 " << ":  " 
                      << v0s0[i - 1][j - 1] <<  endl;
        }
    }
    curvePlot0.plot({ DeltaPs0[0], v0s0[0], DeltaPs0[1], v0s0[1], DeltaPs0[2], v0s0[2], DeltaPs0[3], v0s0[3]});
    std::cout << "Fig1\n" << endl;
    curvePlot0.exec();
    DeltaPs0.clear();
    v0s0.clear();
    

    // 计算原Ergun方程与线弹性修正后的速度与压差Y:v-X:deltaP-[v0-v1]-e0关系
    int dPNum1; // 压差取点数
    int e1Num; // 孔隙率变化取点数
    dPNum1 = 14;
    e1Num = 4;
    std::vector<std::vector<double> >DeltaPs1(e1Num, std::vector<double>(dPNum1));
    std::vector<std::vector<double> >v0s1(e1Num, std::vector<double>(dPNum1));
    std::vector<std::vector <double>>v1s(e1Num, std::vector<double>(dPNum1));
    std::vector<double>C1s(dPNum1);
    double* phi1s = new double[e1Num];

    eggp::Eggplot curvePlot1(SCREEN | PNG | EPS | PDF);
    curvePlot1.xlabel("{/Symbol D}P(kPa)");
    curvePlot1.ylabel("v(m/s)");
    curvePlot1.grid(true);
    curvePlot1.name("L02");
    curvePlot1.title("v_0 and v_1 changing with {/Symbol D}P in {/Symbol e}_0=0.3 and linear elasticity assumption" );
    curvePlot1.legend({ "  v_0", "  v_1"});
    for (int i = 1; i < e1Num + 1; i++)
    {
        double e0s = 0.1 + 0.2 * i;
        phi0s[i - 1] = (1 + e0s) / e0s;
        double A0s = k2 * phi0s[i - 1] * rhoL / D0 / g;
        double B0s = k1 * phi0s[i - 1] * phi0s[i - 1] * mu / (D0 * D0) / g;

        for (int j = 1; j < dPNum1 + 1; j++)
        {
            DeltaPs1[i - 1][j - 1] = (Pin - j * 5000.0)/1000.0;

            C0s[j - 1] = (-DeltaPs1[i - 1][j - 1] * 1000.0 / L0);
            v0s1[i - 1][j - 1] = ((-B0s + sqrt(B0s * B0s - 4 * A0s * C0s[j - 1])) / (2 * A0s));

            double L1s = (1- DeltaPs1[i - 1][j - 1] * 1000.0 / Es) * L0; 
            double e1s = 1 - (1 - e0s) * L0 / L1s;
            double D1s =  D0 * L1s / L0;

            phi1s[i - 1] = (1 + e1s) / e1s;
            double A1s = k2 * phi1s[i - 1] * rhoL / D1s / g;
            double B1s = k1 * phi1s[i - 1] * phi1s[i - 1] * mu / (D1s * D1s) / g;
            C1s[j - 1] = (-DeltaPs1[i - 1][j - 1] * 1000.0 / L1s);
            v1s[i - 1][j - 1] = ((-B1s - sqrt(B1s * B1s - 4 * A1s * C1s[j - 1])) / (2 * A1s));
            std::cout << "e0:  " << e0s << "\t" << "DeltaP " << j << ":  " << DeltaPs1[i - 1][j - 1] << "\t" << "v0 " << ":  " 
                      << v0s1[i - 1][j - 1] << "\t" << "v1 " << ":  " << v1s[i - 1][j - 1] << endl;
        }
    }
    curvePlot1.plot({ DeltaPs1[0], v0s1[0], DeltaPs1[0], v1s[0] });
    std::cout << "Fig2\n" << endl;
    curvePlot1.exec();
    curvePlot1.name("L03");
    curvePlot1.title("v_0 and v_1 changing with {/Symbol D}P in {/Symbol e}_0=0.5 and linear elasticity assumption" );
    curvePlot1.plot({ DeltaPs1[1], v0s1[1], DeltaPs1[1], v1s[1] });
    std::cout << "Fig3\n" << endl;
    curvePlot1.exec();
    curvePlot1.name("L04");
    curvePlot1.title("v_0 and v_1 changing with {/Symbol D}P in {/Symbol e}_0=0.7 and linear elasticity assumption" );
    curvePlot1.plot({ DeltaPs1[2], v0s1[2], DeltaPs1[2], v1s[2] });
    std::cout << "Fig4\n" << endl;
    curvePlot1.exec();
    curvePlot1.name("L05");
    curvePlot1.title("v_0 and v_1 changing with {/Symbol D}P in {/Symbol e}_0=0.9 and linear elasticity assumption" );
    curvePlot1.plot({ DeltaPs1[3], v0s1[3], DeltaPs1[3], v1s[3] });
    std::cout << "Fig5\n" << endl;
    curvePlot1.exec();
    DeltaPs1.clear();
    v0s1.clear();
    v1s.clear();
    C1s.clear();

    // 计算原Ergun方程与粘弹性修正Ergun方程速度与压差Y:v-X:deltaP-[v0-v1]-e0关系
    dPNum0 = 14;
    e0Num = 4;
    std::vector<std::vector<double> >DeltaPs2(e0Num, std::vector<double>(dPNum0));
    std::vector<std::vector<double> >v0s2(e0Num, std::vector<double>(dPNum0));
    std::vector<double>C1s1(dPNum0);
    double* phi2s = new double[e0Num];

    eggp::Eggplot curvePlot3(SCREEN | PNG | EPS | PDF);
    curvePlot3.xlabel("{/Symbol D}P(kPa)");
    curvePlot3.ylabel("v_0(m/s)");
    curvePlot3.grid(true);
    curvePlot3.name("V06");
    curvePlot3.title("v_0 changing with {/Symbol D}P in different {/Symbol e}_0 and viscoelasticity assumption");
    curvePlot3.legend({ "{/Symbol e}_0 = 0.3", "{/Symbol e}_0 = 0.5", "{/Symbol e}_0 = 0.7", "{/Symbol e}_0 = 0.9"});
    
    for (int i = 1; i < e0Num + 1; i++)
    {
//        double e0s = 0.1 + 0.2 * i;
//        phi2s[i-1] = (1 + e0s) / e0s;
//        double A0s = k2 * phi2s[i - 1] * rhoL / D0 / g;
//        double B0s = k1 * phi2s[i - 1] * phi2s[i - 1] * mu / (D0 * D0) / g;
//
        for (int j = 1; j < dPNum0 + 1; j++)
        {
//            DeltaPs0[i - 1][j - 1] = (Pin - j * 5000.0)/1000.0; //  转换成kPa
//            C1s1[j - 1] = -DeltaPs0[i - 1][j - 1] * 1000.0/ L0; //  运算中转换成Pa
//            v0s2[i - 1][j - 1] = ((-B0s + sqrt(B0s * B0s - 4 * A0s * C1s1[j - 1])) / (2 * A0s));
//            std::cout << "e0:  " << e0s << "\t" << "DeltaP " << j << ":  " << DeltaPs0[i - 1][j - 1] << "\t" << "v0 " << ":  " 
//                      << v0s2[i - 1][j - 1] <<  endl;
        }
    }
//    curvePlot3.plot({ DeltaPs0[0], v0s2[0], DeltaPs0[1], v0s2[1], DeltaPs0[2], v0s2[2], DeltaPs0[3], v0s2[3]});
    std::cout << "Fig6\n" << endl;
    curvePlot3.exec();
    DeltaPs2.clear();
    v0s2.clear();


    Sleep(4000);
    delete[] phi0s;
    delete[] phi1s;
    delete[] phi2s;

  return 0;
}
