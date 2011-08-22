/*
 * transfer.cpp
 * dnelson
 */
 
#include "transfer.h"

TransferFunc1D::TransferFunc1D(short int ty, short int vn, 
															 vector<float> &params, vector<Spectrum> &spec)
{
		valNum = vn;
		type   = ty;
		
		// constant (value weighted)
		if (type == 1) {
				if (spec.empty() || spec.size() > 1) {
						IF_DEBUG(cout << "TF1D: Error! Constant type but spec out of bounds." << endl);
						exit(1111);
				}
				
				range[0] = 0.0;
				range[1] = INFINITY;
				le       = spec[0];
		}
		
		// tophat (value weighted)
		if (type == 2) {
				if (spec.empty() || spec.size() > 1 || params.size() != 2) {
						IF_DEBUG(cout << "TF1D: Error! Tophat type but spec/params out of bounds." << endl);
						exit(1112);
				}
				
				range[0] = params[0];
				range[1] = params[1];
				le       = spec[0];
		}
		
		// gaussian (same width for each RGB channel)
		if (type == 3) {
				if (spec.empty() || spec.size() > 1 || params.size() != 2) {
						IF_DEBUG(cout << "TF1D: Error! Gaussian type but spec/params out of bounds." << endl);
						exit(1113);
				}
		
				range[0]      = params[0] - 3.0f * params[1]; //cut at +/- 3s
				range[1]      = params[0] + 3.0f * params[1];
				gaussParam[0] = params[0]; //mean
				gaussParam[1] = params[1]; //sigma
				le            = spec[0];
		}
		
		// unrecognized type / val
		if (type < 0 || type > 3)
				exit(1110);
		if (valNum < 0 || valNum > 1)
				exit(1112);
}

bool TransferFunc1D::InRange(const float *vals)
{
		//IF_DEBUG(cout << "TF1D InRange(" << range[0] << "," << range[1] << ") test = " << vals[valNum] << endl);
		
		if (vals[valNum] < range[0] || vals[valNum] > range[1])
				return false;

		return true;
}

Spectrum TransferFunc1D::Lve(const float *vals) const
{
		float rgb[3];
		le.ToRGB(&rgb[0]);
		
		// constant (value weighted)
		if (type == 1) {
				rgb[0] *= vals[valNum];
				rgb[1] *= vals[valNum];
				rgb[2] *= vals[valNum];
		}
		
		// tophat (value weighted)
		if (type == 2) {
				rgb[0] *= vals[valNum];
				rgb[1] *= vals[valNum];
				rgb[2] *= vals[valNum];
		}
		
		// gaussian
		if (type == 3) {
				//float fwhm = 2.0 * sqrt(2.0 * log(2.0)) * gaussParam[1];
			  float mult   = exp( -1.0 * (vals[valNum] - gaussParam[0])*(vals[valNum] - gaussParam[0]) / 
																  	(2.0f * gaussParam[1]*gaussParam[1]) );
																		
				rgb[0] *= mult;
				rgb[1] *= mult;
				rgb[2] *= mult;
		}
		
		// unrecognized type
		if (type < 0 || type > 3)
				exit(1111);
				
		return Spectrum::FromRGB(rgb);

}

TransferFunction::TransferFunction()
{
		IF_DEBUG(cout << "TransferFunction() constructor." << endl);
		numFuncs = 0;
}

TransferFunction::~TransferFunction()
{
		IF_DEBUG(cout << "TransferFunction() destructor." << endl);
		
		for (int i=0; i < numFuncs; i++) {
				if (f_1D[i]) delete f_1D[i];
		}
}

Spectrum TransferFunction::Lve(const float *vals) const
{
		Spectrum Lve(0.0f);
		
		// consider each independent transfer function
		for (int i=0; i < numFuncs; i++) {
				// exit early if out of range
				if (!f_1D[i]->InRange(vals)) {
						//IF_DEBUG(cout << " Exiting, vals out of bound." << endl);
						continue;
				}
						
				// evaluate TF
				Lve += f_1D[i]->Lve(vals);		
		}

		return Lve;
}

bool TransferFunction::AddConstant(int valNum, Spectrum &sp)
{
		IF_DEBUG(cout << "TF::AddConstant(" << valNum << ",sp) new numFuncs = " << numFuncs+1 << endl);
		
		TransferFunc1D *f;
		vector<float> params;
		vector<Spectrum> spec;
		
		// set constant type and spectrum
		short int type = 1;
		spec.push_back(sp);
		
		// create and store
		f = new TransferFunc1D(type, valNum, params, spec);
		
		f_1D.push_back(f);
		numFuncs++;
		
		return true;
}

bool TransferFunction::AddTophat(int valNum, float min, float max, Spectrum &sp)
{
		IF_DEBUG(cout << "TF::AddTophat(" << valNum << ",sp) range [" << min << "," << max 
									<< "] new numFuncs = " << numFuncs+1 << endl);
		
		TransferFunc1D *f;
		vector<float> params;
		vector<Spectrum> spec;
		
		// set constant type and spectrum
		short int type = 2;
		spec.push_back(sp);
		params.push_back(min);
		params.push_back(max);
		
		// create and store
		f = new TransferFunc1D(type, valNum, params, spec);
		
		f_1D.push_back(f);
		numFuncs++;
		
		return true;
}

bool TransferFunction::AddGaussian(int valNum, float mean, float sigma, Spectrum &sp)
{
		IF_DEBUG(cout << "TF::AddGaussian(" << valNum << ",sp)"
									<< " mean = " << mean << " sigma = " << sigma << " new numFuncs = " << numFuncs+1 << endl);
		
		TransferFunc1D *f;
		vector<float> params;
		vector<Spectrum> spec;
		
		// set constant type and spectrum
		short int type = 3;
		spec.push_back(sp);
		params.push_back(mean);
		params.push_back(sigma);
		
		// create and store
		f = new TransferFunc1D(type, valNum, params, spec);
		
		f_1D.push_back(f);
		numFuncs++;
		
		return true;
}

TransferFunction *CreateTransferFunction()
{
    return new TransferFunction();
}
