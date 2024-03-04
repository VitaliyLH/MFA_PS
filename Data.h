#ifndef DATA_H
#define DATA_H

#include <QFile>
#include <QVector>
#include <cmath>
#include <iostream>
#include <QStringList>
#include <QTextStream>
#include <QGenericMatrix>


using namespace  std;

struct Solutions
{
    double _Hz;
    double _HzA;
    double _H2;

    Solutions operator +(const Solutions &sol);
};



class Data
{
public:
    Data();
    ~Data();
    void parse(QFile &file, QFile &Gfile, QFile &Sfile);

    void SetHz(Solutions& sol ,double Hz) { sol._Hz=Hz; }
    void SetHzA(Solutions& sol ,double HzA) { sol._HzA=HzA; }
    void SetH2(Solutions& sol ,double H2) { sol._H2=H2; }

    double Hz(Solutions sol) { return sol._Hz; }
    double HzA(Solutions sol) { return sol._HzA; }
    double H2(Solutions sol) { return sol._H2; }

    void RemoveDuplicates();
    void PrintSolutions(double T, double mu, QFile &file);
    void PrintGlobalSolutions(double T, double mu, QFile &file);
    double Norm(double x1, double x2, double y1, double y2);
    QVector<Solutions> SignCO(double &T, double mu, Solutions Min, Solutions Max);
    QVector<Solutions> SignSF(double &T, double mu, Solutions Min, Solutions Max);
    Solutions SolveSF(double &T, double mu, Solutions sol);
    Solutions SolveCO(double &T, double mu, Solutions sol);
    void RealSolveNO(double &T, double mu);
    void RealSolveCO(double &T, double mu);
    void RealSolveSF(double &T, double mu);

    double OmegaSol(double T, double mu, Solutions sol);
    bool FindMUS(double T, double MuMin, double MuMax, double StepMU, QFile &Gfile);
    void FindPS(double T, double muPS, QFile &Gfile, QFile &Sfile);


private:
    double _MuMin;
    double _MuMax;
    double _StepMU;
    double _TMax;
    double _TMin;
    double _D;
    double _V;
    double _tb;
    QVector<double> _MUCO;
    QVector<double> _MUSF;
    QVector<double> *_Data;
    QVector<double> *_FEnergy;
    QVector<Solutions> *_Solutions;
    QVector<Solutions> *_SolutionsCO;
    QVector<Solutions> *_SolutionsSF;
    double _Accuracy;

    double SqP(double &T, Solutions sol)
    {
        return sqrt((Hz(sol)+HzA(sol))*(Hz(sol)+HzA(sol))+H2(sol)*H2(sol))/T;
    }

    double SqM(double &T, Solutions sol)
    {
        return sqrt((Hz(sol)-HzA(sol))*(Hz(sol)-HzA(sol))+H2(sol)*H2(sol))/T;
    }

    double Zc(double &T, Solutions sol)
    {
        return 4*(1+exp(-_D/T)*cosh(SqP(T,sol)))*(1+exp(-_D/T)*cosh(SqM(T,sol)));
    }

    double PhiP(double &T, Solutions sol)
    {
        return sinh(SqP(T,sol))/(SqP(T,sol)*(1+exp(-_D/T)*cosh(SqP(T,sol))));
    }

    double PhiM(double &T, Solutions sol)
    {
        return sinh(SqM(T,sol))/(SqM(T,sol)*(1+exp(-_D/T)*cosh(SqM(T,sol))));
    }

    double ThetaP(double &T, Solutions sol)
    {
        return (cosh(SqP(T,sol))+exp(-_D/T)-sinh(SqP(T,sol))*(1+exp(-_D/T)*cosh(SqP(T,sol)))/SqP(T,sol))
                /
                (SqP(T,sol)*SqP(T,sol)*(1+exp(-_D/T)*cosh(SqP(T,sol)))*(1+exp(-_D/T)*cosh(SqP(T,sol))));
    }

    double ThetaM(double &T, Solutions sol)
    {
        return (cosh(SqM(T,sol))+exp(-_D/T)-sinh(SqM(T,sol))*(1+exp(-_D/T)*cosh(SqM(T,sol)))/SqM(T,sol))
                /
                (SqM(T,sol)*SqM(T,sol)*(1+exp(-_D/T)*cosh(SqM(T,sol)))*(1+exp(-_D/T)*cosh(SqM(T,sol))));
    }

    double PiP(double &T, Solutions sol)
    {
        return (Hz(sol)+HzA(sol))*(Hz(sol)+HzA(sol))*ThetaP(T,sol)/(T*T)+PhiP(T,sol);
    }

    double PiM(double &T, Solutions sol)
    {
        return (Hz(sol)-HzA(sol))*(Hz(sol)-HzA(sol))*ThetaM(T,sol)/(T*T)+PhiM(T,sol);
    }

    double PsiP(double &T, Solutions sol)
    {
        return (Hz(sol)+HzA(sol))*H2(sol)*ThetaP(T,sol)/(T*T);
    }

    double PsiM(double &T, Solutions sol)
    {
        return (Hz(sol)-HzA(sol))*H2(sol)*ThetaM(T,sol)/(T*T);
    }

    double DeltaP(double &T, Solutions sol)
    {
        return H2(sol)*H2(sol)*ThetaP(T,sol)/(T*T)+PhiP(T,sol);
    }

    double DeltaM(double &T, Solutions sol)
    {
        return H2(sol)*H2(sol)*ThetaM(T,sol)/(T*T)+PhiM(T,sol);
    }

    double n(double &T, Solutions sol)
    {
        return 0.5*exp(-_D/T)*((Hz(sol)+HzA(sol))*PhiP(T,sol)+(Hz(sol)-HzA(sol))*PhiM(T,sol))/T;
    }

    double a(double &T, Solutions sol)
    {
        return 0.5*exp(-_D/T)*((Hz(sol)+HzA(sol))*PhiP(T,sol)-(Hz(sol)-HzA(sol))*PhiM(T,sol))/T;
    }

    double B(double &T, Solutions sol)
    {
        return 0.5*exp(-_D/T)*(H2(sol)*PhiP(T,sol)+H2(sol)*PhiM(T,sol))/T;
    }

    double EqN(double &T, double mu, Solutions sol)
    {
        return Hz(sol)-mu+4*_V*n(T,sol);
    }

    double EqCO(double &T, Solutions sol)
    {
        return  HzA(sol)-4*_V*a(T,sol);
    }

    double EqSF(double &T, Solutions sol)
    {
        return  H2(sol)-2*_tb*B(T,sol);
    }

    double Omega(double &T, double &mu, Solutions sol)
    {
        return  -
                0.5*T*log(Zc(T,sol))
                -
                2*_V*a(T,sol)*a(T,sol)
                -
                _tb*(B(T,sol)*B(T,sol))
                +
                H2(sol)*B(T,sol)
                +
                HzA(sol)*a(T,sol)
                -
                (mu-Hz(sol))*(mu-Hz(sol))/(8*_V);
    }



};

#endif // DATA_H
