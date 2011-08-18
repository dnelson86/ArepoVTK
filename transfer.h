/*
 * integrator.h
 * dnelson
 */
 
#ifndef AREPO_RT_TRANSFER_H
#define AREPO_RT_TRANSFER_H

#include "ArepoRT.h"
#include "spectrum.h"

#define TF_NUM_QUANT 2

class TransferFunction {
public:
    // construction
    TransferFunction();

    // methods
		bool AddGaussian(string &quantity, float *amp, float *midp, float *sigma);
    bool AddTophat(string &quantity, float *amp, float *midp, float *width);
		bool AddPiecewise();
		
    Spectrum Lve(const float &rho, const float &utherm) const;
		
private:
    // data
    int numFuncs;   // number of independent functions
		int clamping;   // return behavior outside of range
		float range[2]; // 
};

TransferFunction *CreateTransferFunction();


#endif
