using namespace std;

#define NCHANNELS 16

void CalcLightYield(){
    double trapEff[NCHANNELS] = {0.044, 0.044, 0.0344, 0.0344, 0.044, 0.044, 0.042, 0.042, 0.042, 0.044, 0.044, 0.044, 0.0344, 0.0344, 0.044, 0.044};
    // BC12 spect. peak = 435 nm, SCSF7 spect. peak = 450 nm -> Q.E. circ 20%
    double quantEff = 0.20;
    double deltaEnergy[NCHANNELS] = {0.2245, 0.2389, 0.1211, 0.1226, 0.0619, 0.0623, 0.0503, 0.0503, 0.0504, 0.0631, 0.0632, 0.0633, 0.1270, 0.1269, 0.2531, 0.2459};

    double meanCharge[NCHANNELS] = {0., 0.006126, 0.003629, 0.00259, 0.0007382, 0.0007352, 0.0000675, -0.0007171, 0.0003127, 0.0004528, 0.0003826, 0.001944, 0.00257, 0.005411, 0.007866, 0.}; //nC
    double pedestal[NCHANNELS] = {1, -0.00114, 0.00001555, 0.0001628, -0.0009387, -0.0007868, -0.001655, -0.002598, -0.001915, -0.001185, -0.001157, -0.0004987, -0.001136, -0.000788, -0.0002164, 1.}; //nC
    double peCharge[NCHANNELS] = {1.80, 1.63, 1.57, 1.59, 1.60, 1.77, 1.80, 1.77, 1.76, 1.73, 1.69, 1.64, 1.64, 1.59, 1.59, 1.79}; //pC

    double lightyield[NCHANNELS] = {0.};

    for (int i = 1; i < NCHANNELS - 1; i++){
        lightyield[i] = (meanCharge[i] - pedestal[i])*TMath::Power(10., 3)/(peCharge[i]*quantEff*trapEff[i]*deltaEnergy[i]);
        cout << "Channel " << i << " L.Y. = " << lightyield[i] << endl;
    }

}