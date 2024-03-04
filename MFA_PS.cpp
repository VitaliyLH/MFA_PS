#include "MFA_PS.h"

Tc::Tc(QFile &file, QFile &Gfile, QFile &Sfile)
{
   data.parse(file, Gfile, Sfile);
}


Tc::~Tc()
{
  data.~Data();
}
