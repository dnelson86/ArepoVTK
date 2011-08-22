/*
 * integrator.h
 * dnelson
 */
 
#ifndef AREPO_RT_TRANSFER_H
#define AREPO_RT_TRANSFER_H

#include "ArepoRT.h"
#include "spectrum.h"
#include "util.h"

// TODO: dynamic
#define TF_MAX_FUNCS 10

// TODO: bit field approach, such that TransferFunction keeps a record of all the 
//       SphP fields it will need to evaluate Lve(), and returns this upon request
//       to AdvanceRayOneCell() which will use it to form vals[] of the needed elements

#define TF_VAL_DENS   0
#define TF_VAL_UTHERM 1

class TransferFunc1D {
public:
		TransferFunc1D(short int ty, short int vn, vector<float> &params, vector<Spectrum> &spec);
	
		bool InRange(const float *vals);
		Spectrum Lve(const float *vals) const;
		
private:
		short int valNum;   // 1 - density, 2 - utherm
		short int type;     // 1 - constant, 2 - tophat, 3 - gaussian
		bool clamp;         // false - return zero outside range
		                    // true - extrapolate constant values above/below range
		float range[2];     // min,max of val to calculate TF on
		Spectrum le;
		float gaussParam[2]; // mean, sigma
};

class TransferFunction {
public:
    // construction
    TransferFunction();
		~TransferFunction();

    // methods
		bool AddGaussian(int valNum, float midp, float sigma, Spectrum &spec);
    bool AddTophat(int valNum, float min, float max, Spectrum &spec);
		bool AddConstant(int valNum, Spectrum &spec);
		bool AddPiecewise(int valNum, vector<float> &params);
		
		// evaluation
    Spectrum Lve(const float *vals) const;
		
private:
    // data
    short int numFuncs;
		vector<TransferFunc1D *> f_1D;
};



TransferFunction *CreateTransferFunction();


#endif
