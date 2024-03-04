#include "Data.h"

Data::Data()
{
    _Data = new QVector<double>;
    _FEnergy = new QVector<double>;
    _Solutions = new QVector<Solutions>;
    _SolutionsCO = new QVector<Solutions>;
    _SolutionsSF = new QVector<Solutions>;
}

Data::~Data()
{
    delete _Data;
    delete _FEnergy;
    delete _Solutions;
    delete _SolutionsCO;
    delete _SolutionsSF;
}

Solutions Solutions::operator +(const Solutions& sol)
{
    Solutions solp;
    solp._Hz=_Hz+sol._Hz;
    solp._HzA=_HzA+sol._HzA;
    solp._H2=_H2+sol._H2;

    return solp;
}

void Data::parse(QFile &file, QFile &Gfile, QFile &Sfile)
{
    QTextStream stream(&file);
    while (!stream.atEnd()) {

        double data;
        stream>>data;
        _Data->append(data);
    }



    _MuMin=_Data->at(0);
    _MuMax=_Data->at(1);
    _TMin=_Data->at(2);
    _TMax=_Data->at(3);
    _D=_Data->at(4);
    _V=_Data->at(5);
    _tb=_Data->at(6);


    _Accuracy= 0.0000000001;

    double StepT=0.01;

    double MuMin=_MuMin;
    double MuMax=_MuMax;
    double StepMu=0.0005; // Be really very super careful with this parameter :-)

    double MUPS=0;

    FindMUS(_TMin,MuMin,MuMax,StepMu,Gfile);

    for (double T=_TMin; T <= _TMax+0.01*StepT; T+=StepT)
    {
        bool findit=false;

        if (MUPS!=0)
        {
            _MuMin=MUPS-0.2;
            _MuMax=MUPS+0.2;
            _StepMU=0.0005;
        }

        for (double mu=_MuMin; mu <= _MuMax+0.01*_StepMU; mu+=_StepMU)
        {
            RealSolveCO(T,mu);

            PrintGlobalSolutions(T, mu, Gfile);

            RealSolveSF(T,mu);

            PrintGlobalSolutions(T, mu, Gfile);
        }

        int Size=(_MUCO.size() >= _MUSF.size() ? _MUSF.size()-1 : _MUCO.size()-1);

        int i=0;

        while(!findit && i<Size)
        {
            if ( fabs(_MUCO.at(i)-_MUSF.at(i)) < 0.0000001 && fabs(_MUCO.at(i+1)-_MUSF.at(i+1)) < 0.0000001 )
            {

                Solutions CO1=_SolutionsCO->at(i);

                double mu1=_MUCO.at(i);

                double OmegaCO1=Omega(T,mu1,CO1);


                Solutions CO2=_SolutionsCO->at(i+1);

                double mu2=_MUCO.at(i+1);

                double OmegaCO2=Omega(T,mu2,CO2);


                Solutions SF1=_SolutionsSF->at(i);

                double OmegaSF1=Omega(T,mu1,SF1);


                Solutions SF2=_SolutionsSF->at(i+1);

                double OmegaSF2=Omega(T,mu2,SF2);


                if ( ((OmegaCO1 > OmegaSF1 ) && (OmegaCO2 < OmegaSF2)) ||  ((OmegaCO1 < OmegaSF1) && (OmegaCO2 > OmegaSF2)) )
                {

                    if (fabs(OmegaCO1-OmegaCO2) > 0.0000001)
                    {

                        double UP=OmegaCO1*(mu2-mu1)/(OmegaCO2-OmegaCO1)-OmegaSF1*(mu2-mu1)/(OmegaSF2-OmegaSF1);

                        double DOWN=(mu2-mu1)/(OmegaCO2-OmegaCO1)-(mu2-mu1)/(OmegaSF2-OmegaSF1);

                        double OmegaPS=UP/DOWN;

                        MUPS=(OmegaPS-OmegaCO1)*(mu2-mu1)/(OmegaCO2-OmegaCO1)+mu1;
                    }
                    else
                    {
                        MUPS=(OmegaCO1-OmegaSF1)*(mu2-mu1)/(OmegaSF2-OmegaSF1)+mu1;
                    }

                    FindPS(T, MUPS, Gfile, Sfile);

                    findit=true;
                }

            }

            i++;
        }


        if (findit==false) T=_TMax+0.01*StepT;

        _MUCO.clear();
        _MUSF.clear();
        _SolutionsCO->clear();
        _SolutionsSF->clear();

    }

}

bool Data::FindMUS(double T, double MuMin, double MuMax, double StepMU, QFile &Gfile)
{
    bool is_findit=false;

    for (double mu=MuMin; mu <= MuMax+0.01*StepMU; mu+=StepMU)
    {
        RealSolveCO(T,mu);

        PrintGlobalSolutions(T, mu, Gfile);

        RealSolveSF(T,mu);

        PrintGlobalSolutions(T, mu, Gfile);
    }

    int Size=(_MUCO.size() > _MUSF.size() ? _MUSF.size()-1 : _MUCO.size()-1);

    for (int i=0; i < Size; i++)
        for (int j=0; j < Size; j++)
        {
            if ( fabs(_MUCO.at(i)-_MUSF.at(j)) < 0.0000001 && fabs(_MUCO.at(i+1)-_MUSF.at(j+1)) < 0.0000001 )
            {

                Solutions CO1=_SolutionsCO->at(i);

                double mu1=_MUCO.at(i);

                double OmegaCO1=Omega(T,mu1,CO1);


                Solutions CO2=_SolutionsCO->at(i+1);

                double mu2=_MUCO.at(i+1);

                double OmegaCO2=Omega(T,mu2,CO2);


                Solutions SF1=_SolutionsSF->at(i);

                double OmegaSF1=Omega(T,mu1,SF1);


                Solutions SF2=_SolutionsSF->at(i+1);

                double OmegaSF2=Omega(T,mu2,SF2);


                if ( ((OmegaCO1 > OmegaSF1) && (OmegaCO2 < OmegaSF2)) ||  ((OmegaCO1 < OmegaSF1) && (OmegaCO2 > OmegaSF2)) )
                {
                    _MuMin=mu1-50*StepMU;
                    _MuMax=mu2+50*StepMU;
                    _StepMU=(_MuMax-_MuMin)/100.0;

                    is_findit=true;
                }

            }

            if (is_findit)
            {
                i=Size-1;
                j=Size-1;
            }
        }

    _MUCO.clear();
    _MUSF.clear();
    _SolutionsCO->clear();
    _SolutionsSF->clear();

    return is_findit;
}

void Data::FindPS(double T, double muPS, QFile &Gfile, QFile &Sfile)
{
    int FACU=12;

    _MUCO.clear();
    _MUSF.clear();
    _SolutionsCO->clear();
    _SolutionsSF->clear();

    RealSolveCO(T,muPS);

    PrintGlobalSolutions(T, muPS, Gfile);

    RealSolveSF(T,muPS);

    PrintGlobalSolutions(T, muPS, Gfile);


    Solutions CO=_SolutionsCO->at(0);

    double nCO=n(T, CO);

    double OmegaCO=Omega(T, muPS, CO);


    Solutions SF=_SolutionsSF->at(0);

    double nSF=n(T, SF);

    double OmegaSF=Omega(T, muPS, SF);


    double FEnergyCO=OmegaCO+muPS*nCO;
    double FEnergySF=OmegaSF+muPS*nSF;

    QTextStream PrintSolutions( &Sfile );

    if (fabs(OmegaCO-OmegaSF) < 0.0000001 )
    {

        PrintSolutions << T
                       << "\t"
                       << qSetRealNumberPrecision(FACU) << muPS
                       << "\t"
                       << qSetRealNumberPrecision(FACU) << nCO
                       << "\t"
                       << qSetRealNumberPrecision(FACU) << nSF
                       << "\t"
                       << qSetRealNumberPrecision(FACU) << OmegaCO
                       << "\t"
                       << qSetRealNumberPrecision(FACU) << OmegaSF
                       << "\n";

        QTextStream PrintGlobal( &Gfile );

        if (nCO < nSF)
        {
            for (double n=0.0001; n<=nSF; n+=0.001)
            {
                if (n>nCO)
                {
                    PrintGlobal << 4
                                << "\t"
                                << T
                                << "\t"
                                << n
                                << "\t"
                                << qSetRealNumberPrecision(FACU) << ((nSF-n)*FEnergyCO)/(nSF-nCO)+((n-nCO)*FEnergySF)/(nSF-nCO)
                                << "\n";

                }
            }
        }


        if (nCO > nSF)
        {
            for (double n=0.0001; n<=nCO; n+=0.001)
            {
                if (n>nSF)
                {
                    PrintGlobal << 4
                                << "\t"
                                << T
                                << "\t"
                                << n
                                << "\t"
                                << qSetRealNumberPrecision(FACU) << ((nSF-n)*FEnergyCO)/(nSF-nCO)+((n-nCO)*FEnergySF)/(nSF-nCO)
                                << "\n";

                }
            }

        }
    }

}



void Data::RealSolveNO(double &T, double mu)
{
    Solutions Min;
    Solutions Max;

    SetHz(Min,0);
    SetHzA(Min,0);
    SetH2(Min,-0.001);

    SetHz(Max,30);
    SetHzA(Max,0);
    SetH2(Max,0.001);

    QVector<Solutions> sign=SignSF(T, mu, Min, Max);

    for (int i=0; i < sign.size(); i++)
    {
        Solutions sol=SolveSF(T, mu, sign.at(i));

        bool is_equal=false;

        if (_Solutions->size()!=0)
        {
            for (int i=0; i<_Solutions->size();i++)
            {
                double Om1=Omega(T, mu, sol);

                double n1=n(T, sol);

                double Om2=Omega(T, mu, _Solutions->at(i));

                double n2=n(T, _Solutions->at(i));

                if (fabs(Om2-Om1) < 0.0001 && fabs(fabs(n2)-fabs(n1)) < 0.0001 ) is_equal=true;
            }
        }

        double bla1=fabs(EqN(T, mu, sol));
        double bla2=fabs(EqSF(T, sol));

        if (fabs(EqN(T, mu, sol)) > _Accuracy || fabs(EqSF(T, sol)) > _Accuracy ) is_equal=true;
        if (isnan(EqN(T, mu, sol)) || isnan(EqSF(T, sol)) || fabs(Hz(sol)) > 20.0 || fabs(H2(sol)) > 3.0 ) is_equal=true;

        if (!is_equal) _Solutions->append(sol);
    }

}


void Data::RealSolveCO(double &T, double mu)
{
    Solutions Min;
    Solutions Max;

    SetHz(Min,0);
    SetHzA(Min,0.1);
    SetH2(Min,0);

    SetHz(Max,4);
    SetHzA(Max,4*_V+0.1);
    SetH2(Max,0);

    QVector<Solutions> sign=SignCO(T, mu, Min, Max);

    for (int i=0; i < sign.size(); i++)
    {
        Solutions sol=SolveCO(T, mu, sign.at(i));

        bool is_equal=false;

        if (_Solutions->size()!=0)
        {
            for (int i=0; i<_Solutions->size();i++)
            {
                double Om1=Omega(T, mu, sol);

                double n1=n(T, sol);

                double Om2=Omega(T, mu, _Solutions->at(i));

                double n2=n(T, _Solutions->at(i));

                if (fabs(Om2-Om1) < 0.0001 && fabs(fabs(n2)-fabs(n1)) < 0.0001 ) is_equal=true;
            }
        }

        if (fabs(EqN(T, mu, sol)) > _Accuracy || fabs(EqCO(T, sol)) > _Accuracy ) is_equal=true;
        if (isnan(EqN(T, mu, sol)) || isnan(EqCO(T, sol)) || fabs(Hz(sol)) > 20.0 || fabs(HzA(sol)) > 4.0 ) is_equal=true;

        if (!is_equal) _Solutions->append(sol);
    }

}



void Data::RealSolveSF(double &T, double mu)
{
    Solutions Min;
    Solutions Max;

    SetHz(Min,0);
    SetHzA(Min,0);
    SetH2(Min,0.1);

    SetHz(Max,4);
    SetHzA(Max,0);
    SetH2(Max,2*_tb+0.1);

    QVector<Solutions> sign=SignSF(T, mu, Min, Max);

    for (int i=0; i < sign.size(); i++)
    {
        Solutions sol=SolveSF(T, mu, sign.at(i));

        bool is_equal=false;

        if (_Solutions->size()!=0)
        {
            for (int i=0; i<_Solutions->size();i++)
            {
                double Om1=Omega(T, mu, sol);

                double n1=n(T, sol);

                double Om2=Omega(T, mu, _Solutions->at(i));

                double n2=n(T, _Solutions->at(i));

                if (fabs(Om2-Om1) < 0.0001 && fabs(fabs(n2)-fabs(n1)) < 0.0001 ) is_equal=true;
            }
        }

        if (fabs(EqN(T, mu, sol)) > _Accuracy || fabs(EqSF(T, sol)) > _Accuracy ) is_equal=true;
        if (isnan(EqN(T, mu, sol)) || isnan(EqSF(T, sol)) || fabs(Hz(sol)) > 20.0 || fabs(H2(sol)) > 3.0 ) is_equal=true;

        if (!is_equal) _Solutions->append(sol);
    }

}


void Data::PrintSolutions(double T, double mu, QFile &file)
{
    int FACU=12;

    QTextStream Print( &file );

    for (int i=0; i<_Solutions->size();i++)
    {

        int Type=0;

        if (HzA(_Solutions->at(i))>0.0001) Type=1;
        if (H2(_Solutions->at(i))>0.0001) Type=2;


        if (Type==1)
        {
            _MUCO.append(mu);
            _SolutionsCO->append(_Solutions->at(i));
        }
        if (Type==2)
        {
            _MUSF.append(mu);
            _SolutionsSF->append(_Solutions->at(i));
        }

        Print << Type
              << "\t"
              << mu
              << "\t"
              << T
              << "\t"
              << qSetRealNumberPrecision(FACU) << Hz(_Solutions->at(i))
              << "\t"
              << qSetRealNumberPrecision(FACU) << HzA(_Solutions->at(i))
              << "\t"
              << qSetRealNumberPrecision(FACU) << H2(_Solutions->at(i))
              << "\t"
              << qSetRealNumberPrecision(FACU) << Omega(T, mu, _Solutions->at(i))
              << "\t"
              << EqN(T, mu, _Solutions->at(i))
              << "\t"
              << EqCO(T, _Solutions->at(i))
              << "\t"
              << EqSF(T, _Solutions->at(i))
              << "\n";

    }



}

void Data::PrintGlobalSolutions(double T, double mu, QFile &file)
{
    QTextStream Print( &file );

    QVector<double> Om;

    for (int i=0; i<_Solutions->size();i++) Om.append(Omega(T, mu, _Solutions->at(i)));


    qSort(Om);

    for (int i=0; i<_Solutions->size();i++)
    {

        if (fabs(Omega(T, mu, _Solutions->at(i))-Om.at(0))<_Accuracy)
        {

            int Type=0;

            if (HzA(_Solutions->at(i))>0.0001) Type=1;
            if (H2(_Solutions->at(i))>0.0001) Type=2;


            if (Type==1)
            {
                _MUCO.append(mu);
                _SolutionsCO->append(_Solutions->at(i));
            }
            if (Type==2)
            {
                _MUSF.append(mu);
                _SolutionsSF->append(_Solutions->at(i));
            }

            //                        Print << Type
            //                              << "\t"
            //                              << mu
            //                              << "\t"
            //                              << T
            //                              << "\t"
            //                              << qSetRealNumberPrecision(10) << Hz(_Solutions->at(i))
            //                              << "\t"
            //                              << qSetRealNumberPrecision(10) << HzA(_Solutions->at(i))
            //                              << "\t"
            //                              << qSetRealNumberPrecision(10) << H2(_Solutions->at(i))
            //                              << "\t"
            //                              << qSetRealNumberPrecision(10) << Omega(T, mu, _Solutions->at(i))
            //                              << "\t"
            //                              << EqN(T, mu, _Solutions->at(i))
            //                              << "\t"
            //                              << EqCO(T, _Solutions->at(i))
            //                              << "\t"
            //                              << EqSF(T, _Solutions->at(i))
            //                              << "\n";
        }
    }

    _Solutions->clear();

}



double Data::Norm(double x1, double x2, double y1, double y2)
{
    return sqrt((x1-y1)*(x1-y1)+(x2-y2)*(x2-y2));
}

Solutions Data::SolveCO(double &T, double mu, Solutions sol)
{
    Solutions point=sol;

    QGenericMatrix<1,2,double> X0;
    QGenericMatrix<1,2,double> Xi;
    QGenericMatrix<1,2,double> F;
    QGenericMatrix<2,2,double> J;
    QGenericMatrix<2,2,double> InvJ;
    QGenericMatrix<2,2,double> Test;


    X0.fill(0);
    Xi.fill(1);
    F.fill(0);
    J.fill(0);
    InvJ.fill(0);
    Test.fill(0);

    int step=0;

    while( Norm(X0(0,0),X0(1,0),Xi(0,0),Xi(1,0)) > _Accuracy )
    {

        X0(0,0)=Hz(point);
        X0(1,0)=HzA(point);

        F(0,0)=EqN(T, mu, point);
        F(1,0)=EqCO(T, point);


        J(0,0)=1+2*_V*exp(-_D/T)*(PiP(T, point)+PiM(T, point))/T;
        J(0,1)=2*_V*exp(-_D/T)*(PiP(T, point)-PiM(T, point))/T;
        J(1,0)=-2*_V*exp(-_D/T)*(PiP(T, point)-PiM(T, point))/T;
        J(1,1)=1-2*_V*exp(-_D/T)*(PiP(T, point)+PiM(T, point))/T;

        double Det=J(0,0)*J(1,1)-J(0,1)*J(1,0);

        InvJ(0,0)=J(1,1)/Det;
        InvJ(0,1)=-J(0,1)/Det;
        InvJ(1,0)=-J(1,0)/Det;
        InvJ(1,1)=J(0,0)/Det;

        Test=J*InvJ;

        if (fabs(Test(0,0)) > 1-0.0001 && fabs(Test(0,0)) < 1+0.0001)
            if (fabs(Test(1,1)) > 1-0.0001 && fabs(Test(1,1)) < 1+0.0001)
                if (fabs(Test(0,1)) < 0.0001)
                    if (fabs(Test(1,0)) < 0.0001)
                        Xi=X0-InvJ*F;

        SetHz(point,Xi(0,0));
        SetHzA(point,Xi(1,0));
        SetH2(point,H2(point));

        step++;

        if (Hz(point) < 0.0 || step>10000)
        {
            X0(0,0)=Xi(0,0);
            X0(1,0)=Xi(1,0);
        }

    }

    SetHz(point,Xi(0,0));
    SetHzA(point,fabs(Xi(1,0)));
    SetH2(point,H2(point));

    return point;

}

Solutions Data::SolveSF(double &T, double mu, Solutions sol)
{
    Solutions point=sol;

    QGenericMatrix<1,2,double> X0;
    QGenericMatrix<1,2,double> Xi;
    QGenericMatrix<1,2,double> F;
    QGenericMatrix<2,2,double> J;
    QGenericMatrix<2,2,double> InvJ;
    QGenericMatrix<2,2,double> Test;


    X0.fill(0);
    Xi.fill(100);
    F.fill(0);
    J.fill(0);
    InvJ.fill(0);
    Test.fill(0);

    int step=0;

    while( Norm(X0(0,0),X0(1,0),Xi(0,0),Xi(1,0)) > _Accuracy )
    {

        X0(0,0)=Hz(point);
        X0(1,0)=H2(point);

        F(0,0)=EqN(T, mu, point);
        F(1,0)=EqSF(T, point);


        J(0,0)=1+2*_V*exp(-_D/T)*(PiP(T, point)+PiM(T, point))/T;
        J(0,1)=2*_V*exp(-_D/T)*(PsiP(T, point)+PsiM(T, point))/T;
        J(1,0)=-_tb*exp(-_D/T)*(PsiP(T, point)+PsiM(T, point))/T;
        J(1,1)=1-_tb*exp(-_D/T)*(DeltaP(T, point)+DeltaM(T, point))/T;

        double Det=J(0,0)*J(1,1)-J(0,1)*J(1,0);

        InvJ(0,0)=J(1,1)/Det;
        InvJ(0,1)=-J(0,1)/Det;
        InvJ(1,0)=-J(1,0)/Det;
        InvJ(1,1)=J(0,0)/Det;

        Test=J*InvJ;

        if (fabs(Test(0,0)) > 1-0.0001 && fabs(Test(0,0)) < 1+0.0001)
            if (fabs(Test(1,1)) > 1-0.0001 && fabs(Test(1,1)) < 1+0.0001)
                if (fabs(Test(0,1)) < 0.0001)
                    if (fabs(Test(1,0)) < 0.0001)
                        Xi=X0-InvJ*F;

        SetHz(point,Xi(0,0));
        SetHzA(point,HzA(point));
        SetH2(point,Xi(1,0));

        step++;

        if (Hz(point) < 0.0 || step>10000)
        {
            X0(0,0)=Xi(0,0);
            X0(1,0)=Xi(1,0);
        }

    }

    SetHz(point,Xi(0,0));
    SetHzA(point,HzA(point));
    SetH2(point,fabs(Xi(1,0)));

    return point;

}


QVector<Solutions> Data::SignCO(double &T, double mu, Solutions Min, Solutions Max)
{
    QVector<Solutions> pts;

    Solutions p1;
    Solutions p2;

    double SStep=10;
    double signconst=0.01;

    double h2=H2(Min);

    double StepHz=(Hz(Max)-Hz(Min))/SStep;
    double StepHzA=(HzA(Max)-HzA(Min))/SStep;

    for (double P1=Hz(Min);P1 <= Hz(Max); P1+=StepHz)
        for (double P2=HzA(Min);P2 <= HzA(Max); P2+=StepHzA)
        {
            bool is_equal=false;

            SetHz(p1,P1);
            SetHzA(p1,P2);
            SetH2(p1,h2);

            if(SqP(T, p1)!=0.0 && SqM(T, p1)!=0.0)
            {
                SetHz(p2,P1+StepHz);
                SetHzA(p2,P2);
                SetH2(p2,h2);

                if (
                        (EqN(T, mu, p1)*EqN(T, mu, p2)<signconst)
                        &&
                        (EqCO(T, p1)*EqCO(T, p2)<signconst)
                        )
                {
                    Solutions point=p1+p2;

                    SetHz(point,0.5*Hz(point));
                    SetHzA(point,0.5*HzA(point));
                    SetH2(point,0.5*H2(point));

                    if (pts.size()!=0)
                    {
                        for (int i=0; i<pts.size();i++)
                        {
                            double n1=n(T, point);

                            double Om1=Omega(T, mu, point);

                            double n2=n(T, pts.at(i));

                            double Om2=Omega(T, mu, pts.at(i));

                            if ( fabs(Om2-Om1) < 0.001 && fabs(n2-n1) < 0.001 ) is_equal=true;
                        }
                    }

                    if (!is_equal) pts.append(point);
                }

                SetHz(p2,P1);
                SetHzA(p2,P2+StepHzA);
                SetH2(p2,h2);

                if (
                        (EqN(T, mu, p1)*EqN(T, mu, p2)<signconst)
                        &&
                        (EqCO(T, p1)*EqCO(T, p2)<signconst)
                        )
                {
                    Solutions point=p1+p2;

                    SetHz(point,0.5*Hz(point));
                    SetHzA(point,0.5*HzA(point));
                    SetH2(point,0.5*H2(point));

                    if (pts.size()!=0)
                    {
                        for (int i=0; i<pts.size();i++)
                        {
                            double n1=n(T, point);

                            double Om1=Omega(T, mu, point);

                            double n2=n(T, pts.at(i));

                            double Om2=Omega(T, mu, pts.at(i));

                            if ( fabs(Om2-Om1) < 0.001 && fabs(n2-n1) < 0.001 ) is_equal=true;
                        }
                    }

                    if (!is_equal) pts.append(point);
                }

                SetHz(p2,P1+StepHz);
                SetHzA(p2,P2+StepHzA);
                SetH2(p2,h2);

                if (
                        (EqN(T, mu, p1)*EqN(T, mu, p2)<signconst)
                        &&
                        (EqCO(T, p1)*EqCO(T, p2)<signconst)
                        )
                {
                    Solutions point=p1+p2;

                    SetHz(point,0.5*Hz(point));
                    SetHzA(point,0.5*HzA(point));
                    SetH2(point,0.5*H2(point));

                    if (pts.size()!=0)
                    {
                        for (int i=0; i<pts.size();i++)
                        {
                            double n1=n(T, point);

                            double Om1=Omega(T, mu, point);

                            double n2=n(T, pts.at(i));

                            double Om2=Omega(T, mu, pts.at(i));

                            if ( fabs(Om2-Om1) < 0.001 && fabs(n2-n1) < 0.001 ) is_equal=true;
                        }
                    }

                    if (!is_equal) pts.append(point);
                }

            }

        }


    return pts;


}



QVector<Solutions> Data::SignSF(double &T, double mu, Solutions Min, Solutions Max)
{
    QVector<Solutions> pts;

    Solutions p1;
    Solutions p2;

    double SStep=10;
    double signconst=0.01;

    double hza=HzA(Min);

    double StepHz=(Hz(Max)-Hz(Min))/SStep;
    double StepH2=(H2(Max)-H2(Min))/SStep;

    for (double P1=Hz(Min);P1 <= Hz(Max); P1+=StepHz)
        for (double P2=H2(Min);P2 <= H2(Max); P2+=StepH2)
        {
            bool is_equal=false;

            SetHz(p1,P1);
            SetHzA(p1,hza);
            SetH2(p1,P2);

            if(SqP(T, p1)!=0.0 && SqM(T, p1)!=0.0)
            {
                SetHz(p2,P1+StepHz);
                SetHzA(p2,hza);
                SetH2(p2,P2);

                if (
                        (EqN(T, mu, p1)*EqN(T, mu, p2)<signconst)
                        &&
                        (EqSF(T, p1)*EqSF(T, p2)<signconst)
                        )
                {
                    Solutions point=p1+p2;

                    SetHz(point,0.5*Hz(point));
                    SetHzA(point,0.5*HzA(point));
                    SetH2(point,0.5*H2(point));

                    if (pts.size()!=0)
                    {
                        for (int i=0; i<pts.size();i++)
                        {
                            double n1=n(T, point);

                            double Om1=Omega(T, mu, point);

                            double n2=n(T, pts.at(i));

                            double Om2=Omega(T, mu, pts.at(i));

                            if ( fabs(Om2-Om1) < 0.001 && fabs(n2-n1) < 0.001 ) is_equal=true;
                        }
                    }

                    if (!is_equal) pts.append(point);
                }

                SetHz(p2,P1);
                SetHzA(p2,hza);
                SetH2(p2,P2+StepH2);

                if (
                        (EqN(T, mu, p1)*EqN(T, mu, p2)<signconst)
                        &&
                        (EqSF(T, p1)*EqSF(T, p2)<signconst)
                        )
                {
                    Solutions point=p1+p2;

                    SetHz(point,0.5*Hz(point));
                    SetHzA(point,0.5*HzA(point));
                    SetH2(point,0.5*H2(point));

                    if (pts.size()!=0)
                    {
                        for (int i=0; i<pts.size();i++)
                        {
                            double n1=n(T, point);

                            double Om1=Omega(T, mu, point);

                            double n2=n(T, pts.at(i));

                            double Om2=Omega(T, mu, pts.at(i));

                            if ( fabs(Om2-Om1) < 0.001 && fabs(n2-n1) < 0.001 ) is_equal=true;
                        }
                    }

                    if (!is_equal) pts.append(point);
                }

                SetHz(p2,P1+StepHz);
                SetHzA(p2,hza);
                SetH2(p2,P2+StepH2);

                if (
                        (EqN(T, mu, p1)*EqN(T, mu, p2)<signconst)
                        &&
                        (EqSF(T, p1)*EqSF(T, p2)<signconst)
                        )
                {
                    Solutions point=p1+p2;

                    SetHz(point,0.5*Hz(point));
                    SetHzA(point,0.5*HzA(point));
                    SetH2(point,0.5*H2(point));

                    if (pts.size()!=0)
                    {
                        for (int i=0; i<pts.size();i++)
                        {
                            double n1=n(T, point);

                            double Om1=Omega(T, mu, point);

                            double n2=n(T, pts.at(i));

                            double Om2=Omega(T, mu, pts.at(i));

                            if ( fabs(Om2-Om1) < 0.001 && fabs(n2-n1) < 0.001 ) is_equal=true;
                        }
                    }

                    if (!is_equal) pts.append(point);
                }

            }

        }


    return pts;


}
