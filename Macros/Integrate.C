#include <TF2.h>

using namespace std;

void Integrate(){
    TF2* integrand = new TF2("Integrand", "6/TMath::Pi()/TMath::Power((TMath::Power(x, 2) + TMath::Power(y, 2) + 1), 5./2.)*(1-13.43/16.*x)*(1-13.43/172.*y)", 0., 16./13.43, 0., 172./13.43);

    integrand->SetNpx(10000);
    integrand->SetNpy(10000);

    cout << "Integral = " << integrand->Integral(0., 16./13.43, 0., 172./13.43) << endl;

}