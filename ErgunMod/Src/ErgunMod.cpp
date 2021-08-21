#include<cmath>
#include<cstdio>
#include<iostream>
#include<vector>
#include<windows.h>
#include"eggplot.h"

using namespace std;
using namespace eggp;

class PorousMedia
{
public:
    double PMlength;
    double PMporousR;
    double PMelasticM;
    double PMfibreD;
};
class Fluid
{
public:
    double Fdensity;
    double Fviscosity;
};
int main()
{
//    system("chcp 65001 && cls");
    // 设置变量
    /*---------------------------------------------------------------------------
    ---设置多孔介质介质特性参数
    basisWeight -- 纸幅克重         thick -- 纸幅厚度
    porosity -- 孔隙率              elasticMod -- 纸幅弹性模量
    dryContent -- 纸幅干度          effDiameter -- 纸幅有效直径
    ---设置多孔介质颗粒特性参数
    parLength -- 颗粒长度           parDiameter -- 颗粒直径(圆柱)
    cellDensity -- 纤维素密度       lenDensity -- 长度密度
    ---设置多孔介质中流体特性参数
    airViscosity -- 空气动力粘度    airDensity -- 空气密度
    waterDensity -- 水动力粘度
    ---设置输入输出参数
    pressureInlet -- 输入压力       pressureOutlet -- 输出压力
    pressureDiff -- 压差            realVelocity -- 孔隙速度
    specificVelocity -- 表征速度    realVelocityInite -- 初始孔隙速度
    ---设置Ergun方程系数
    k1 -- 粘性项系数                k2 -- 惯性项系数
    ---设置初始变量
    porosityInit -- 初始孔隙率        thickInit -- 初始厚度
    effDiameterInit -- 有效初始颗粒直径   dryContentInit -- 初始干度
    ---------------------------------------------------------------------------*/
    double basisWeight = 50e-3; // kg/m2
    double thick;
    double thickInit; // m
    double thickChange; // %
    double porosity;
    double porosityInit = 0.4; 
    double elasticMod = 10e6; // Pa
    double dryContentInit = 0.054;
    double dryContent = 0.24;
    double dryContentFinal = 0.94;
    double parLength = 0.72e-3; // m
    double lenDensity = 59e-9; // kg/m3
    double cellDensity = 1.5e3; // kg/m3
    double parDiameter; // m
    double parDiameterInit; // m
    double effDiameter;
    double effDiameterInit;
    double airViscosity = 1.01e-3; // Pa.s
    double airDensity = 1.29; // kg/m3
    double waterDensity = 1000.0; // kg/m3
    double pressureInlet = 101325.0; // Pa
    double pressureOutlet = 61325.0; // Pa
    double pressureDiff;
    double specificVelocity;
    double specificVelocityInit;
    double realVelocity;
    double realVelocityInit;
    double k1 = 150.0;
    double k2 = 1.750;
    double pi = 3.1415927;
    double g = 9.8; // m/s2
    double Reynolds;
    double ReynoldsInit;


    thickInit = basisWeight * dryContentFinal / (1.0 - porosityInit) * (1.0 / cellDensity + (1.0 - dryContentInit) / waterDensity / dryContentInit);
    // 设置变量关系公式
    /*---------------------------------------------------------------------------
    ---------------------------------------------------------------------------*/
    pressureDiff = pressureInlet - pressureOutlet;
    thick = (1 - pressureDiff / elasticMod) * thickInit; // m 胡克定律
    parDiameterInit = sqrt(lenDensity * 4.0 / pi * (dryContentInit / cellDensity + \
        (1.0 - dryContentInit) / waterDensity));
    parDiameter = sqrt(lenDensity * 4.0 / pi * (dryContent / cellDensity + \
        (1.0 - dryContent) / waterDensity));
    effDiameterInit = 3.0 * parDiameterInit * parLength / (2.0 * parLength + parDiameterInit);
    effDiameter = 3.0 * parDiameter * parLength / (2.0 * parLength + parDiameter);
    porosity = porosityInit / (1.0 - pressureDiff / elasticMod);
    thickChange = abs(thick - thickInit) / thickInit * 100.0;

    // 设置原始Ergun方程计算用中间变量
    /*---------------------------------------------------------------------------
    phi0 -- (1+porosityInit)/(porosityInit*effDiameterInit)    A0 -- 原始Ergun方程二次项系数
    B0 -- 原始Ergun方程一次项系数    C0 -- 原始Ergun方程常数项
    ---------------------------------------------------------------------------*/
    double phi0 = (1.0 - porosityInit) / (porosityInit * effDiameterInit);
    double A0 = k2 * phi0 * airDensity;
    double B0 = k1 * phi0 * phi0 * airViscosity;
    double C0 = -pressureDiff / thickInit * g;
    // 设置线弹性修正Ergun方程计算用中间变量
    /*---------------------------------------------------------------------------
    phi1 -- (1+porosity)/(porosity*effDiameter)    A1 -- 原始Ergun方程二次项系数
    B1 -- 原始Ergun方程一次项系数    C1 -- 原始Ergun方程常数项
    ---------------------------------------------------------------------------*/
    double phi1 = (1.0 - porosity) / (porosity * effDiameter);
    double A1 = k2 * phi1 * airDensity;
    double B1 = k1 * phi1 * phi1 * airViscosity;
    double C1 = -pressureDiff / thick * g;
    // 计算原始Ergun方程
    realVelocityInit = (-B0 + sqrt(B0 * B0 - 4 * A0 * C0)) / (2 * A0);
    specificVelocityInit = (-B0 + sqrt(B0 * B0 - 4 * A0 * C0)) / (2 * A0) * porosityInit;
    // 计算弹性修正Ergun方程
    realVelocity = (-B1 + sqrt(B1 * B1 - 4 * A1 * C1)) / (2 * A1);
    specificVelocity = (-B1 + sqrt(B1 * B1 - 4 * A1 * C1)) / (2 * A1) * porosity;
    // 计算线弹性修正Ergun方程相对变化率
    double vChange = abs(abs(realVelocity) - realVelocityInit) / realVelocityInit * 100.0;
    double uChange = abs(abs(specificVelocity) - specificVelocityInit) / specificVelocityInit * 100.0;
    // 计算雷诺数
    ReynoldsInit = airDensity * realVelocityInit * effDiameterInit / airViscosity;
    Reynolds = airDensity * realVelocity * effDiameter / airViscosity;
    // 输出已知量
    std::cout << "Init:" << endl;
    std::cout << "PresDiff  :" << pressureDiff << "  Pa" << "\t" << "PaprThik  :" << thickInit << "  m" << "\n" \
              << "FibrDiam  :" << parDiameterInit << "  m" << "\t" << "FibrLeng  :" << parLength << "  m" << "\n" \
              << "EffcDiam  :" << effDiameterInit << "  m" << "\t" << "ElasModl  :" << elasticMod << "  Pa" << "\n" << endl;
    // 输出计算量
    std::cout << "Deform:" << endl;
    std::cout << "PresDiff  :" << pressureDiff << "  Pa" << "\t" << "PaprThik  :" << thick << "  m" << "\n" \
              << "FibrDiam  :" << parDiameter << "  m" << "\t" << "FibrLeng  :" << parLength << "  m" << "\n" \
              << "EffcDiam  :" << effDiameter << "  m" << "\t" << "ElasModl  :" << elasticMod << "  Pa" << "\n" << endl;
    std::cout << "Compute:" << endl;
    std::cout << "VeloBefr  :" << realVelocityInit << "  m/s" << "\t" << "VeloAftr  :" << realVelocity << "  m/s" << endl;
    std::cout << "RealVeloChange  :" << vChange << "  %" << "\t" << "SpecVeloChange  :" << uChange << "  %" << endl;
    std::cout << "ThickInit  :" << thickInit << "\t" << "ThickChange  :" << thickChange << "  %" << endl;
    std::cout << "ReInit  :" << ReynoldsInit << "\t" << "Re  :" << Reynolds << endl;

    // 原Ergun方程画图
    // 计算Y:v0-X:deltaP-[e0]关系
    int dPNum0; // 压差取点数
    int e0Num; // 孔隙率变化取点数
    dPNum0 = 14;
    e0Num = 4;
    std::vector<std::vector<double> >DeltaPs0(e0Num, std::vector<double>(dPNum0));
    std::vector<std::vector<double> >v0s0(e0Num, std::vector<double>(dPNum0));
    std::vector<std::vector<double> >u0s0(e0Num, std::vector<double>(dPNum0));
    std::vector<double>C0s(dPNum0);
    double* phi0s = new double[e0Num];

    eggp::Eggplot curvePlot0(SCREEN | PNG | EPS | PDF);
    curvePlot0.xlabel("{/Symbol D}P(kPa)");
    curvePlot0.ylabel("v_0(m/s)");
    curvePlot0.grid(true);
    curvePlot0.name("Rigid");
    curvePlot0.title("Air velocity changing with pressure difference in different porosity and rigid body assumption");
    curvePlot0.legend({ "{/Symbol e}_0 = 0.3", "{/Symbol e}_0 = 0.5", "{/Symbol e}_0 = 0.7", "{/Symbol e}_0 = 0.9"});
    
    for (int i = 1; i < e0Num + 1; i++)
    {
        double e0s = 0.2 * i;
        phi0s[i-1] = (1 + e0s) / e0s / effDiameterInit;
        double A0s = k2 * phi0s[i - 1] * airDensity;
        double B0s = k1 * phi0s[i - 1] * phi0s[i - 1] * airViscosity;

        for (int j = 1; j < dPNum0 + 1; j++)
        {
            DeltaPs0[i - 1][j - 1] = (pressureInlet - j * 5000.0 - 30000.0)/1000.0; //  转换成kPa
            C0s[j - 1] = -DeltaPs0[i - 1][j - 1] * 1000.0 * g / thickInit; //  运算中转换成Pa
            v0s0[i - 1][j - 1] = ((-B0s + sqrt(B0s * B0s - 4 * A0s * C0s[j - 1])) / (2 * A0s));
            u0s0[i - 1][j - 1] = ((-B0s + sqrt(B0s * B0s - 4 * A0s * C0s[j - 1])) / (2 * A0s)) * e0s;
            std::cout << "e0:  " << e0s << "\t" << "DeltaP " << j << ":  " << DeltaPs0[i - 1][j - 1] << "\t" << "v0 " << ":  " 
                      << v0s0[i - 1][j - 1] << endl;
        }
    }
//    curvePlot0.plot({ DeltaPs0[0], v0s0[0], DeltaPs0[1], v0s0[1], DeltaPs0[2], v0s0[2], DeltaPs0[3], v0s0[3]});
    curvePlot0.plot({ DeltaPs0[0], u0s0[0], DeltaPs0[1], u0s0[1], DeltaPs0[2], u0s0[2], DeltaPs0[3], u0s0[3]});
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
    std::vector<std::vector<double> >u0s1(e1Num, std::vector<double>(dPNum1));
    std::vector<std::vector <double>>v1s(e1Num, std::vector<double>(dPNum1));
    std::vector<std::vector <double>>u1s(e1Num, std::vector<double>(dPNum1));
    std::vector<std::vector <double>>vChanges(e1Num, std::vector<double>(dPNum1));
    std::vector<std::vector <double>>uChanges(e1Num, std::vector<double>(dPNum1));
    std::vector<double>C1s(dPNum1);
    double* phi1s = new double[e1Num];

    eggp::Eggplot curvePlot1(SCREEN | PNG | EPS | PDF);
    curvePlot1.xlabel("{/Symbol D}P(kPa)");
    curvePlot1.ylabel("v(m/s)");
    curvePlot1.grid(true);
    curvePlot1.name("LinearElasity0_2");
    curvePlot1.title("Air velocity changing with pressure difference in porosity=0.3 and linear elasticity assumption" );
    curvePlot1.legend({ "  v_0", "  v_1"});
    for (int i = 1; i < e1Num + 1; i++)
    {
        double e0s =  0.2 * i;
        phi0s[i - 1] = (1 + e0s) / e0s / effDiameterInit;
        double A0s = k2 * phi0s[i - 1] * airDensity;
        double B0s = k1 * phi0s[i - 1] * phi0s[i - 1] * airViscosity;

        for (int j = 1; j < dPNum1 + 1; j++)
        {
            DeltaPs1[i - 1][j - 1] = (pressureInlet - j * 5000.0 - 30000.0)/1000.0;

            C0s[j - 1] = (-DeltaPs1[i - 1][j - 1] * g * 1000.0 / thickInit);
            v0s1[i - 1][j - 1] = ((-B0s + sqrt(B0s * B0s - 4 * A0s * C0s[j - 1])) / (2 * A0s));
            u0s1[i - 1][j - 1] = ((-B0s + sqrt(B0s * B0s - 4 * A0s * C0s[j - 1])) / (2 * A0s)) * e0s;

            double L1s = (1- DeltaPs1[i - 1][j - 1] * 1000.0 / elasticMod) * thickInit; 
            double e1s = 1 - (1 - e0s) * thickInit / L1s;
//            double D1s =  effDiameter;

            phi1s[i - 1] = (1 + e1s) / e1s / effDiameter;
            double A1s = k2 * phi1s[i - 1] * airDensity;
            double B1s = k1 * phi1s[i - 1] * phi1s[i - 1] * airViscosity;
            C1s[j - 1] = (-DeltaPs1[i - 1][j - 1] * g * 1000.0 / L1s);
            v1s[i - 1][j - 1] = ((-B1s + sqrt(B1s * B1s - 4 * A1s * C1s[j - 1])) / (2 * A1s));
            u1s[i - 1][j - 1] = ((-B1s + sqrt(B1s * B1s - 4 * A1s * C1s[j - 1])) / (2 * A1s)) * e1s;
            vChanges[i - 1][j - 1] = abs(v1s[i - 1][j - 1] - v0s1[i - 1][j - 1]) / v0s1[i - 1][j - 1] * 100.0;
            uChanges[i - 1][j - 1] = abs(v1s[i - 1][j - 1] - v0s1[i - 1][j - 1]) / v0s1[i - 1][j - 1] * e1s * 100.0;
            std::cout << "e0:  " << e0s << "\t" << "DeltaP " << j << ":  " << DeltaPs1[i - 1][j - 1] << "\t" \
                      << "v0:  " << v0s1[i - 1][j - 1] << "\t" << "v1:  " << v1s[i - 1][j - 1]  << "\t" \
                      << "vChange:  " << vChanges[i - 1][j - 1] << "\t" << "uChange:" << uChanges[i - 1][j - 1] << endl;
        }
    }
//    std::cout << "Fig2\n" << endl;
//    curvePlot1.plot({ DeltaPs1[0], v0s1[0], DeltaPs1[0], v1s[0] });
    curvePlot1.plot({ DeltaPs1[0], u0s1[0], DeltaPs1[0], u1s[0] });
    curvePlot1.exec();
//    std::cout << "Fig3\n" << endl;
    curvePlot1.name("LinearElasity0_4");
    curvePlot1.title("Velocity changing with pressure difference in porosity=0.5 and linear elasticity assumption" );
//    curvePlot1.plot({ DeltaPs1[1], v0s1[1], DeltaPs1[1], v1s[1] });
    curvePlot1.plot({ DeltaPs1[1], u0s1[1], DeltaPs1[1], u1s[1] });
    curvePlot1.exec();
//    std::cout << "Fig4\n" << endl;
    curvePlot1.name("LinearElasity0_6");
    curvePlot1.title("Velocity changing with pressure difference in porosity=0.7 and linear elasticity assumption" );
//    curvePlot1.plot({ DeltaPs1[2], v0s1[2], DeltaPs1[2], v1s[2] });
    curvePlot1.plot({ DeltaPs1[2], u0s1[2], DeltaPs1[2], u1s[2] });
    curvePlot1.exec();
//   std::cout << "Fig5\n" << endl;
    curvePlot1.name("LinearElasity0_8");
    curvePlot1.title("Velocity changing with pressure difference in porosity=0.9 and linear elasticity assumption" );
//    curvePlot1.plot({ DeltaPs1[3], v0s1[3], DeltaPs1[3], v1s[3] });
    curvePlot1.plot({ DeltaPs1[3], u0s1[3], DeltaPs1[3], u1s[3] });
    curvePlot1.exec();
    DeltaPs1.clear();
    v0s1.clear();
    v1s.clear();
    C1s.clear();

    // 计算原Ergun方程与粘弹性修正Ergun方程速度与压差Y:v-X:deltaP-[v0-v1]-e0关系
//    dPNum0 = 14;
//    e0Num = 4;
//    std::vector<std::vector<double> >DeltaPs2(e0Num, std::vector<double>(dPNum0));
//    std::vector<std::vector<double> >v0s2(e0Num, std::vector<double>(dPNum0));
//    std::vector<double>C1s1(dPNum0);
//    double* phi2s = new double[e0Num];
//
//    eggp::Eggplot curvePlot3(SCREEN | PNG | EPS | PDF);
//    curvePlot3.xlabel("{/Symbol D}P(kPa)");
//    curvePlot3.ylabel("v_0(m/s)");
//    curvePlot3.grid(true);
//    curvePlot3.name("V06");
//    curvePlot3.title("v_0 changing with {/Symbol D}P in different {/Symbol e}_0 and viscoelasticity assumption");
//    curvePlot3.legend({ "{/Symbol e}_0 = 0.3", "{/Symbol e}_0 = 0.5", "{/Symbol e}_0 = 0.7", "{/Symbol e}_0 = 0.9"});
//    
//    for (int i = 1; i < e0Num + 1; i++)
//    {
////        double e0s = 0.1 + 0.2 * i;
////        phi2s[i-1] = (1 + e0s) / e0s;
////        double A0s = k2 * phi2s[i - 1] * rhoL / D0 / g;
////        double B0s = k1 * phi2s[i - 1] * phi2s[i - 1] * mu / (D0 * D0) / g;
////
//        for (int j = 1; j < dPNum0 + 1; j++)
//        {
////            DeltaPs0[i - 1][j - 1] = (Pin - j * 5000.0)/1000.0; //  转换成kPa
////            C1s1[j - 1] = -DeltaPs0[i - 1][j - 1] * 1000.0/ L0; //  运算中转换成Pa
////            v0s2[i - 1][j - 1] = ((-B0s + sqrt(B0s * B0s - 4 * A0s * C1s1[j - 1])) / (2 * A0s));
////            std::cout << "e0:  " << e0s << "\t" << "DeltaP " << j << ":  " << DeltaPs0[i - 1][j - 1] << "\t" << "v0 " << ":  " 
////                      << v0s2[i - 1][j - 1] <<  endl;
//        }
//    }
////    curvePlot3.plot({ DeltaPs0[0], v0s2[0], DeltaPs0[1], v0s2[1], DeltaPs0[2], v0s2[2], DeltaPs0[3], v0s2[3]});
//    std::cout << "Fig6\n" << endl;
//    curvePlot3.exec();
//    DeltaPs2.clear();
//    v0s2.clear();


    Sleep(4000);
    delete[] phi0s;
    delete[] phi1s;
//    delete[] phi2s;

  return 0;
}
