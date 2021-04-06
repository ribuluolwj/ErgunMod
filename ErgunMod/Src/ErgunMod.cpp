#include<cmath>
#include<cstdio>
#include<iostream>
#include<windows.h>
#include "MatPlot.h"

using namespace std;
using namespace MatPlot;


int main()
{
  // 设置多孔介质物性参数
  /*---------------------------------------------------------------------------
  L -- 初始厚度    e0 -- 初始孔隙度    ex -- 平衡时刻x处的孔隙度
  Es -- 多孔介质弹性模量    D0 -- 初始纤维当量直径
  k1 -- 粘性项系数    k2 -- 惯性项系数
  ---------------------------------------------------------------------------*/
  double L = 1.0e-3; // m
  double e0 = 0.7;
  double k1 = 150.0;
  double k2 = 1.75;
  double D0 = 6.0e-3; // m
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
  double v0, u0, vx, ut, px;
  // 设置方程中间变量
  /*---------------------------------------------------------------------------
  phi0 -- (1+e0)/e0    A0 -- 原始Ergun方程二次项系数
  B0 -- 原始Ergun方程一次项系数    C0 -- 原始Ergun方程常数项
  ---------------------------------------------------------------------------*/
  double phi0 = (1 + e0) / e0;
  double A0 = k2 * phi0 * rhoL / D0 / g;
  double B0 = k1 * phi0 * phi0 * mu / (D0 * D0) / g;
  double C0 = - DeltaP / L;
  // 计算原始Ergun方程
  v0 = (-B0 + sqrt(B0 * B0 - 4 * A0 * C0)) / (2 * A0);
  std::cout << v0 << "\t" << "m/s" << endl;
  // 计算v0-deltaP关系
  int dPNum; // data points number
  dPNum = 14;
  double *DeltaPs = new double[dPNum];
  double *v0s = new double[dPNum];
  double *C0s = new double[dPNum];
  MatPlotInit();
  show_control();
  for (int i = 1; i < dPNum+1; i++)
  {
    DeltaPs[i-1] = Pin - i * 5000.0;
    C0s[i-1] = - DeltaPs[i-1] / L;
    v0s[i-1] = (-B0 + sqrt(B0 * B0 - 4 * A0 * C0s[i-1])) / (2 * A0);
    std::cout << "DeltaP " << i << ":" << DeltaPs[i - 1] << "\t" << "v0 " << ":" << v0s[i - 1] << endl;
  }
  plot(DeltaPs, v0s, dPNum, 'r');
  e0 = 0.4;
  for (int i = 1; i < dPNum+1; i++)
  {
    DeltaPs[i-1] = Pin - i * 5000.0;
    C0s[i-1] = - DeltaPs[i-1] / L;
    v0s[i-1] = (-B0 + sqrt(B0 * B0 - 4 * A0 * C0s[i-1])) / (2 * A0);
    std::cout << "DeltaP " << i << ":" << DeltaPs[i - 1] << "\t" << "v0 " << ":" << v0s[i - 1] << endl;
  }
  plot(DeltaPs, v0s, dPNum, 'g');

  Sleep(600000);
  MatPlotClose();
  delete[] DeltaPs;
  delete[] v0s;
  delete[] C0s;

  return 0;
}
