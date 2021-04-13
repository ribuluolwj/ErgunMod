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
    L -- 初始厚度    e0 -- 初始孔隙度    ex -- 平衡时刻x处的孔隙度
    Es -- 多孔介质弹性模量    D0 -- 初始纤维当量直径
    k1 -- 粘性项系数    k2 -- 惯性项系数
    ---------------------------------------------------------------------------*/
    double L0 = 1.0e-3; // m
    double e0 = 0.4;
    double k1 = 150.0;
    double k2 = 1.75;
    double D0 = 6.0e-3; // m
    double Es = 1400; //Pa
//    double Es = 1400; //Pa
    // 设置流体介质物性参数
    /*---------------------------------------------------------------------------
    rhoL -- 流体密度    mu -- 流体动力粘度    u0 -- 流体初始表征速度
    v0 -- 流体初始孔隙速度    vx -- 平衡时刻x处的孔隙速度    ut -- t时刻的表征速度
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
    L1 -- 初始厚度    e1 -- 初始孔隙度    ex -- 平衡时刻x处的孔隙度
    Es -- 多孔介质弹性模量    D1 -- 初始纤维当量直径
    k1 -- 粘性项系数    k2 -- 惯性项系数
    ---------------------------------------------------------------------------*/
    double L1 = (1- DeltaP / Es) * L0; // m
    double e1 = 1 - (1 - e0) * L0 / L1;
    double k1 = 150.0;
    double k2 = 1.75;
    double D1 = D0 * L1 / L0; // m

    // 设置原始Ergun方程中间变量
    /*---------------------------------------------------------------------------
    phi0 -- (1+e0)/e0    A0 -- 原始Ergun方程二次项系数
    B0 -- 原始Ergun方程一次项系数    C0 -- 原始Ergun方程常数项
    ---------------------------------------------------------------------------*/
    double phi0 = (1 + e0) / e0;
    double A0 = k2 * phi0 * rhoL / D0 / g;
    double B0 = k1 * phi0 * phi0 * mu / (D0 * D0) / g;
    double C0 = -DeltaP / L0;
    // 设置修正Ergun方程中间变量
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
    // 计算弹性修正Ergun方程相对变化率
    double vChange = (v1 - v0) / v0 * 100;

    std::cout << "v0:  " << v0 << "\t" << "v1:  " << v1 << "\t" << "m/s" << "\t" << vChange << "\t" << "%" << endl;
    // 计算v0-deltaP关系
    int dPNum; // delta pressure data points number
    int e0Num; // e0 number
    dPNum = 14;
    e0Num = 4;
    std::vector<double>DeltaPs;
    std::vector<double>v0s,v1s;
    std::vector<double>C0s,C1s;
    double* phi0s = new double[e0Num];
    double* phi1s = new double[e0Num];

    eggp::Eggplot curvePlot(SCREEN|PNG);

    curvePlot.xlabel("{/Symbol D}P");
    curvePlot.ylabel("v_0");
    curvePlot.grid(true);
    for (int j = 1; j < e0Num + 1; j++)
    {
        double e0s = 0.1 + 0.2 * j;
        phi0s[j-1] = (1 + e0s) / e0s;
        double A0s = k2 * phi0s[j-1] * rhoL / D0 / g;
        double B0s = k1 * phi0s[j-1] * phi0s[j-1] * mu / (D0 * D0) / g;

        for (int i = 1; i < dPNum + 1; i++)
        {
            DeltaPs.push_back(Pin - i * 5000.0);

            C0s.push_back(-DeltaPs[i - 1] / L);
            v0s.push_back((-B0s + sqrt(B0s * B0s - 4 * A0s * C0s[i - 1])) / (2 * A0s));

            double L1s = (1- DeltaPs[i-1] / Es) * L0; 
            double e1s = 1 - (1 - e0s) * L0 / L1s;

            phi1s[j-1] = (1 + e1s) / e1s;
            double A1s = k2 * phi1s[j - 1] * rhoL / D0 / g;
            double B1s = k1 * phi1s[j - 1] * phi1s[j-1] * mu / (D0 * D0) / g;
            C1s.push_back(-DeltaPs[i - 1] / L1s);
            v1s.push_back((-B1s + sqrt(B1s * B1s - 4 * A1s * C1s[i - 1])) / (2 * A1s));
            std::cout << "DeltaP " << i << ":  " << DeltaPs[i - 1] << "\t" << "v0 " << ":  " 
                      << v0s[i - 1] << "\t" << "v1 " << ":  " << v1s[i - 1] << endl;
        }
    }
    curvePlot.plot({ DeltaPs, v0s, DeltaPs, v1s });
    //curvePlot.legend({ "{/Symbol 1}=1",
    //                   "{/Symbol 1}=4" });
    curvePlot.exec();
    Sleep(60000);
    delete[] phi0s;

  return 0;
}
