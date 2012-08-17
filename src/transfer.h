/*
 * integrator.h
 * dnelson
 */
 
#ifndef AREPO_RT_TRANSFER_H
#define AREPO_RT_TRANSFER_H

//#include "ArepoRT.h"
#include "spectrum.h"

// TODO: bit field approach, such that TransferFunction keeps a record of all the 
//       SphP fields it will need to evaluate Lve(), and returns this upon request
//       to AdvanceRayOneCell() which will use it to form vals[] of the needed elements

#define TF_VAL_DENS        0
#define TF_VAL_UTHERM      1
#define TF_VAL_PRES        2
#define TF_VAL_ENERGY      3

#define TF_VAL_VEL_X       4
#define TF_VAL_VEL_Y       5
#define TF_VAL_VEL_Z       6
#define TF_VAL_VEL_DIV     7
#define TF_VAL_VEL_CURL    8

#define TF_VAL_POTENTIAL   9

#define TF_VAL_METALLICITY 10
#define TF_VAL_NE          11
#define TF_VAL_SFR         12

#define TF_NUM_VALS        13

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

class TransferFunction /*: public VolumeRegion*/ {
public:
    // construction
    TransferFunction(const Spectrum &sig_a);
		~TransferFunction();

    // methods
		bool AddGaussian(int valNum, float midp, float sigma, Spectrum &spec);
    bool AddTophat(int valNum, float min, float max, Spectrum &spec);
		bool AddConstant(int valNum, Spectrum &spec);
		bool AddPiecewise(int valNum, vector<float> &params);
		
		bool AddParseString(string &addTFstr);
		
		// evaluation
    //Spectrum sigma_a(const Point &p, const Vector &, float) const {    }
    //Spectrum sigma_s(const Point &p, const Vector &, float) const {    }
		Spectrum sigma_t() const { return sig_t; }
    Spectrum Lve(const float *vals) const;
    //Spectrum tau(const Ray &r, float stepSize, float offset) const {   }
		
private:
    // data
    short int numFuncs;
		vector<TransferFunc1D *> f_1D;
		
		map<string,int> valNums;
		
		// tau = scatter + abs
		Spectrum sig_a, sig_s, sig_t;
};

#endif
