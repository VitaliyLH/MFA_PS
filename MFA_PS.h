#ifndef TC_H
#define TC_H

#include "Data.h"
#include "complex.h"

class Tc
{
public:
    Tc(QFile &file , QFile &Gfile, QFile &Sfile);
    ~Tc();
    void Print();

private:

    Data data;

};

#endif // TC_H
