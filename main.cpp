#include"MFA_PS.h"
#include <QStringBuilder>
#include <iomanip>


int main(int argc, char *argv[])
{
    if (argc != 2) {
        puts("WTF!");
        return 1;
    }

    Tc* T;

        QFile file(argv[1]);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
            return EXIT_FAILURE;

        QString G("Global");
        QString S("Solutions");

        QString other(argv[1]);

        int index = other.indexOf(".dat");

        QString fileG(other);
        QString fileS(other);

        fileG.insert(index,G);
        fileS.insert(index,S);

        QFile Gfile(fileG);
        QFile Sfile(fileS);

        Gfile.open(QIODevice::WriteOnly | QIODevice::Text);
        Sfile.open(QIODevice::WriteOnly | QIODevice::Text);

        T = new Tc(file , Gfile , Sfile);

        file.close();

        Gfile.close();
        Sfile.close();

        return 0;

}
