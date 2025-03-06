using namespace std;

#define NCHANNELS 16

void CalcLightYield(){
    double trapEff[NCHANNELS] = {0.044, 0.044, 0.0344, 0.0344, 0.044, 0.044, 0.042, 0.042, 0.042, 0.044, 0.044, 0.044, 0.0344, 0.0344, 0.044, 0.044};
    // BC12 spect. peak = 435 nm, SCSF7 spect. peak = 450 nm -> Q.E. circ 20%
    double quantEff = 0.20;
    double deltaEnergy[NCHANNELS] = {0.2245, 0.2389, 0.1211, 0.1226, 0.0619, 0.0623, 0.0503, 0.0503, 0.0504, 0.0631, 0.0632, 0.0633, 0.1270, 0.1269, 0.2531, 0.2459};
    double deltaEnergySigma[NCHANNELS] = {1., 0.1327, 0.0813, 0.0844, 0.0515, 0.0524, 0.0446, 0.0450, 0.0448, 0.0539, 0.0527, 0.0534, 0.0901, 0.0901, 0.1578, 1.};
    int nGeant = 183380;

    double meanCharge[NCHANNELS] = {1., 0.006126, 0.003629, 0.00259, 0.0007382, 0.0007352, 0.0000675, -0.0007171, 0.0003127, 0.0004528, 0.0003826, 0.001944, 0.00257, 0.005411, 0.007866, 1.}; //nC
    double meanChargeSigma[NCHANNELS] = {1., 0.007204, 0.003753, 0.003059, 0.002457, 0.002366, 0.002472, 0.002671, 0.003053, 0.002566, 0.002235, 0.003196, 0.004344, 0.006274, 0.006804, 1.};
    int nData = 10000;

    double pedestal[NCHANNELS] = {1, -0.00114, 0.00001555, 0.0001628, -0.0009387, -0.0007868, -0.001655, -0.002598, -0.001915, -0.001185, -0.001157, -0.0004987, -0.001136, -0.000788, -0.0002164, 1.}; //nC

    double peCharge[NCHANNELS] = {1.80, 1.63, 1.57, 1.59, 1.60, 1.77, 1.80, 1.77, 1.76, 1.73, 1.69, 1.64, 1.64, 1.59, 1.59, 1.79}; //pC
    double gainSigma[NCHANNELS] = {0.005, 0.006, 0.006, 0.006, 0.006, 0.005, 0.005, 0.005, 0.004, 0.005, 0.006, 0.004, 0.005, 0.005, 0.004, 0.004}; //x10^7
    double eCharge = 1.6;

    double lightyield[NCHANNELS] = {0.};
    double lightyieldSigma[NCHANNELS] = {0.};

    for (int i = 1; i < NCHANNELS - 1; i++){
        lightyield[i] = (meanCharge[i] - pedestal[i])*TMath::Power(10., 3)/(peCharge[i]*quantEff*trapEff[i]*deltaEnergy[i]);

        lightyieldSigma[i] = lightyield[i]*TMath::Sqrt(TMath::Power(deltaEnergySigma[i]/(deltaEnergy[i]*TMath::Sqrt(nGeant)), 2.) + TMath::Power(meanChargeSigma[i]/((meanCharge[i] - pedestal[i])*TMath::Sqrt(nData)), 2) + TMath::Power(gainSigma[i]*eCharge/peCharge[i], 2));
        cout << "Channel " << i << " L.Y. = " << lightyield[i] << " pm " << lightyieldSigma[i] << endl;
    }

}