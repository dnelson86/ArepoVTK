/*
 * transfer.cpp
 * dnelson
 */
 
#include "transfer.h"
#include "spectrum.h"
#include "util.h"
#include "fileio.h"

TransferFunc1D::TransferFunc1D(short int ty, short int vn, 
															 vector<float> &params, vector<Spectrum> &spec, string ctName)
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
				
				gaussParam[0] = params[0]; //mean
				gaussParam[1] = params[1]; //sigma
				range[0]      = gaussParam[0] - 4.0f * gaussParam[1]; //cut at +/- 4s
				range[1]      = gaussParam[0] + 4.0f * gaussParam[1];
				le            = spec[0];
		}
		
		// constant (discrete color table)
		if (type == 4) {
				if (!spec.empty() || params.size() != 2 || !ctName.length()) {
						IF_DEBUG(cout << "TF1D: Error! ConstantDiscrete type but spec/params out of bounds." << endl);
						exit(1158);
				}
				
				colorTable  = ctName;
				ctMinMax[0] = params[0];
				ctMinMax[1] = params[1];
				range[0] = 0.0;
				range[1] = INFINITY;
				
				colorTableLen = loadDiscreteColorTable(ctName, &colorTableVals);
				CheckReverse();
				
				ctStep = (ctMinMax[1]-ctMinMax[0]) / colorTableLen;
		}
		
		// tophat (discrete color table)
		if (type == 5) {
				if (!spec.empty() || params.size() != 4 || !ctName.length()) {
						IF_DEBUG(cout << "TF1D: Error! TophatDiscrete type but spec/params out of bounds." << endl);
						exit(1159);
				}
				
				colorTable  = ctName;
				ctMinMax[0] = params[0];
				ctMinMax[1] = params[1];
				range[0]    = params[2];
				range[1]    = params[3];
				
				colorTableLen = loadDiscreteColorTable(ctName, &colorTableVals);
				CheckReverse();
				
				ctStep = (ctMinMax[1]-ctMinMax[0]) / colorTableLen;
		}		
		
		// gaussian (discrete color table)
		if (type == 6) {
				if (!spec.empty() || params.size() != 4 || !ctName.length()) {
						IF_DEBUG(cout << "TF1D: Error! GaussianDiscrete type but spec/params out of bounds." << endl);
						exit(1160);
				}
				
				colorTable  = ctName;
				ctMinMax[0] = params[0];
				ctMinMax[1] = params[1];

				gaussParam[0] = params[2]; //mean
				gaussParam[1] = params[3]; //sigma
				range[0]      = gaussParam[0] - 4.0f * gaussParam[1]; //cut at +/- 4s
				range[1]      = gaussParam[0] + 4.0f * gaussParam[1];
				
				colorTableLen = loadDiscreteColorTable(ctName, &colorTableVals);
				CheckReverse();
				
				ctStep = (ctMinMax[1]-ctMinMax[0]) / colorTableLen;
		}

		// unrecognized type / val
		if (type < 0 || type > 6)
				exit(1110);
		if (valNum < 0 || valNum >= TF_NUM_VALS)
				exit(1112);
}

void TransferFunc1D::CheckReverse()
{
		// if max<min, swap them and reverse the color table
		if (ctMinMax[1] < ctMinMax[0]) {
				float swap = ctMinMax[1];
				ctMinMax[1] = ctMinMax[0];
				ctMinMax[0] = swap;
				
				vector<float> oldCtVals = colorTableVals;
				int count = 0;
				
				for (size_t i = colorTableLen-1; i != (size_t)-1; i-- ) {
						// stride four
						colorTableVals[count++] = oldCtVals[i*4+0];
						colorTableVals[count++] = oldCtVals[i*4+1];
						colorTableVals[count++] = oldCtVals[i*4+2];
						colorTableVals[count++] = oldCtVals[i*4+3];
				}
		}
		
		// calculate a max/min of all the RGB values

		float curMax = 0.0, curMin = INFINITY;
		
		for (size_t i = 0; i < colorTableVals.size()/4; i++) {
				if (colorTableVals[i*4+0] < curMin) curMin = colorTableVals[i*4+0];
				if (colorTableVals[i*4+1] < curMin) curMin = colorTableVals[i*4+1];
				if (colorTableVals[i*4+2] < curMin) curMin = colorTableVals[i*4+2];
				
				if (colorTableVals[i*4+0] > curMax) curMax = colorTableVals[i*4+0];
				if (colorTableVals[i*4+1] > curMax) curMax = colorTableVals[i*4+1];
				if (colorTableVals[i*4+2] > curMax) curMax = colorTableVals[i*4+2];
		}
		
		if (curMin < 0.0 || curMax > 255.0) {
				cout << "ERROR: Corrupt (scale) color table, exiting." << endl;
				exit(1163);
		}
		
		// if the maximum is greater than 2.0 then assume the inputs are in [0,255] and rescale to [0,1]
		if (curMax > 2.0) {
				cout << "Warning: Rescaling input color table to [0,1]!" << endl;
				
				for (size_t i = 0; i < colorTableVals.size()/4; i++) {
						colorTableVals[i*4+0] /= 255.0;
						colorTableVals[i*4+1] /= 255.0;
						colorTableVals[i*4+2] /= 255.0;
						// leave alpha alone
				}
		}		
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
		rgb[0] = 0; rgb[1] = 0; rgb[2] = 0;
		
		// constant or tophat (value weighted)
		if (type == 1 || type == 2) {
				le.ToRGB(&rgb[0]);
				rgb[0] *= vals[valNum];
				rgb[1] *= vals[valNum];
				rgb[2] *= vals[valNum];
		}
		
		// gaussian
		if (type == 3) {
				le.ToRGB(&rgb[0]);
				
				//float fwhm = 2.0 * sqrt(2.0 * log(2.0)) * gaussParam[1];
			  float mult   = exp( -1.0 * (vals[valNum] - gaussParam[0])*(vals[valNum] - gaussParam[0]) / 
																  	(2.0f * gaussParam[1]*gaussParam[1]) );
																		
				rgb[0] *= 1.0 * mult;
				rgb[1] *= 1.0 * mult;
				rgb[2] *= 1.0 * mult;
		}
		
		// constant or tophat (discrete color table)
		if (type == 4 || type == 5) {
				// determine lerp indices and interpolant point
				int left  = Clamp((int)floor((vals[valNum] - ctMinMax[0] ) / ctStep),0,colorTableLen-1);
				int right = Clamp(left + 1,0,colorTableLen-1);

				float t = Clamp((vals[valNum] - ctMinMax[0] - left*ctStep) / ctStep,0.0,1.0);
				
				// lerp alpha and each rgb channel separately
				float alpha = Lerp(t,colorTableVals[left*4+3],colorTableVals[right*4+3]);
				
				rgb[0] = 1.0 * alpha * vals[TF_VAL_DENS] * Lerp(t,colorTableVals[left*4+0],colorTableVals[right*4+0]);
				rgb[1] = 1.0 * alpha * vals[TF_VAL_DENS] * Lerp(t,colorTableVals[left*4+1],colorTableVals[right*4+1]);
				rgb[2] = 1.0 * alpha * vals[TF_VAL_DENS] * Lerp(t,colorTableVals[left*4+2],colorTableVals[right*4+2]);
		}
		
		// gaussian (discrete color table)
		if (type == 6) {
				// determine lerp indices and interpolant point
				int left  = Clamp((int)floor((vals[valNum] - ctMinMax[0] ) / ctStep),0,colorTableLen-1);
				int right = Clamp(left + 1,0,colorTableLen-1);

				float t = Clamp((vals[valNum] - ctMinMax[0] - left*ctStep) / ctStep,0.0,1.0);
				
				// lerp alpha and each rgb channel separately
				float alpha = Lerp(t,colorTableVals[left*4+3],colorTableVals[right*4+3]);
		
			  float mult  = exp( -1.0 * (vals[valNum] - gaussParam[0])*(vals[valNum] - gaussParam[0]) / 
																  	(2.0f * gaussParam[1]*gaussParam[1]) );
																		
				rgb[0] = alpha * mult * Lerp(t,colorTableVals[left*4+0],colorTableVals[right*4+0]);
				rgb[1] = alpha * mult * Lerp(t,colorTableVals[left*4+1],colorTableVals[right*4+1]);
				rgb[2] = alpha * mult * Lerp(t,colorTableVals[left*4+2],colorTableVals[right*4+2]);
		}		
				
		return Spectrum::FromRGB(rgb);

}

TransferFunc1D::~TransferFunc1D()
{
		// free any color tables
		// unused right now
}

TransferFunction::TransferFunction(const Spectrum &sig_a)
{
		IF_DEBUG(cout << "TransferFunction() constructor." << endl);
		numFuncs = 0;
		
		// value name -> index number map
		valNums["Density"] = 0;
		valNums["Utherm"]  = 1;
		valNums["Vel_X"]   = 2;
		valNums["Vel_Y"]   = 3;
		valNums["Vel_Z"]   = 4;
		
		// scattering
		sig_s = 0.0f;
		
		// tau/trans
		sig_t = sig_a + sig_s;
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
				if (!f_1D[i]->InRange(vals))
						continue;
						
				// evaluate TF
				Lve += f_1D[i]->Lve(vals);		
		}

		return Lve;
}

bool TransferFunction::InRange(const float *vals) const
{
		bool flag = false;
		
		// consider each independent transfer function
		for (int i=0; i < numFuncs; i++) {
				if (f_1D[i]->InRange(vals))
						flag = true;
		}
		
		return flag;
}

bool TransferFunction::AddConstant(int valNum, Spectrum &sp)
{
		IF_DEBUG(cout << "TF::AddConstant(" << valNum << ",sp) new numFuncs = " << numFuncs+1 << endl);
		
		TransferFunc1D *f;
		vector<float> params;
		vector<Spectrum> spec;
		string ctName; // dummy
		
		// set constant type and spectrum
		short int type = 1;
		spec.push_back(sp);
		
		// create and store
		f = new TransferFunc1D(type, valNum, params, spec, ctName);
		
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
		string ctName; // dummy
		
		// set constant type and spectrum
		short int type = 2;
		spec.push_back(sp);
		params.push_back(min);
		params.push_back(max);
		
		// create and store
		f = new TransferFunc1D(type, valNum, params, spec, ctName);
		
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
		string ctName; // dummy
		
		// set constant type and spectrum
		short int type = 3;
		spec.push_back(sp);
		params.push_back(mean);
		params.push_back(sigma);
		
		// create and store
		f = new TransferFunc1D(type, valNum, params, spec, ctName);
		
		f_1D.push_back(f);
		numFuncs++;
		
		return true;
}

bool TransferFunction::AddConstantDiscrete(int valNum, string ctName, float ctMin, float ctMax)
{
		IF_DEBUG(cout << "TF::AddConstantDiscrete(" << valNum << "," << ctName << "," << ctMin << "," << ctMax 
									<< ") new numFuncs = " << numFuncs+1 << endl);
		
		TransferFunc1D *f;
		vector<float> params;
		vector<Spectrum> spec;
		
		// set type and parameters
		short int type = 4;
		params.push_back(ctMin);
		params.push_back(ctMax);
		
		// create and store
		f = new TransferFunc1D(type, valNum, params, spec, ctName);
		
		f_1D.push_back(f);
		numFuncs++;
		
		return true;
}

bool TransferFunction::AddTophatDiscrete(int valNum, string ctName, float ctMin, float ctMax, float min, float max)
{
		IF_DEBUG(cout << "TF::AddTophatDiscrete(" << valNum << "," << ctName << "," << ctMin << "," << ctMax 
									<< ") new numFuncs = " << numFuncs+1 << endl);
		
		TransferFunc1D *f;
		vector<float> params;
		vector<Spectrum> spec;
		
		// set type and parameters
		short int type = 5;
		params.push_back(ctMin);
		params.push_back(ctMax);
		params.push_back(min);
		params.push_back(max);
		
		// create and store
		f = new TransferFunc1D(type, valNum, params, spec, ctName);
		
		f_1D.push_back(f);
		numFuncs++;

		return true;
}

bool TransferFunction::AddGaussianDiscrete(int valNum, string ctName, float ctMin, float ctMax, float midp, float sigma)
{
		IF_DEBUG(cout << "TF::AddConstantDiscrete(" << valNum << "," << ctName << "," << ctMin << "," << ctMax 
									<< ") new numFuncs = " << numFuncs+1 << endl);
		
		TransferFunc1D *f;
		vector<float> params;
		vector<Spectrum> spec;
		
		// set type and parameters
		short int type = 6;
		params.push_back(ctMin);
		params.push_back(ctMax);
		params.push_back(midp);
		params.push_back(sigma);
		
		// create and store
		f = new TransferFunc1D(type, valNum, params, spec, ctName);
		
		f_1D.push_back(f);
		numFuncs++;

		return true;
}

bool TransferFunction::AddParseString(string &addTFstr)
{
		if (Config.verbose)
				cout << "Add to TF: " << addTFstr << endl;
				
		int valNum;
		Spectrum spec;
		float rgb[3];
		
	  // split string
		vector<string> p;
		string item;
		char delim;
		delim = ' ';
		
		stringstream ss(addTFstr);
		while(getline(ss, item, delim))
				p.push_back(item);
				
		// check string size
		if (p.size() < 3) {
				cout << "ERROR: addTF string too short: " << addTFstr << endl;
				exit(1122);
		}
		
		// check value name
		if (!valNums.count(p[1])) {
				cout << "ERROR: addTF string has invalid value name: " << p[1] << endl;
				exit(1123);
		}
		
		// check that code was compiled with necessary flags to support this value name
		if (p[1] == "Metallicity") {
#ifndef METALS
				cout << "ERROR: TF on value 'Metallicity' requires METALS compile option in Arepo." << endl;
				exit(1124);
#endif
		}
		
		valNum = valNums[p[1]];
		
		// determine type of TF
		if (p[0] == "constant") {
				if (p.size() == 3) {
						spec = Spectrum::FromNamed(p[2]);
				} else if (p.size() == 5) {
						rgb[0] = atof(p[2].c_str());
						rgb[1] = atof(p[3].c_str());
						rgb[2] = atof(p[4].c_str());
						spec = Spectrum::FromRGB(rgb);
				} else {
						cout << "ERROR: addTF constant string bad: " << addTFstr << endl;
						exit(1126);
				}
				
				AddConstant(valNum,spec);
		
		} else if (p[0] == "gaussian") {
				if (p.size() == 5) {
						spec = Spectrum::FromNamed(p[4]);
				} else if (p.size() == 7) {
						rgb[0] = atof(p[4].c_str());
						rgb[1] = atof(p[5].c_str());
						rgb[2] = atof(p[6].c_str());
						spec = Spectrum::FromRGB(rgb);
				} else {
						cout << "ERROR: addTF gaussian string bad: " << addTFstr << endl;
						exit(1127);
				}
				
				float mean  = atof(p[2].c_str());
				float sigma = atof(p[3].c_str());
				
				AddGaussian(valNum,mean,sigma,spec);
		
		} else if (p[0] == "tophat") {
				if (p.size() == 5) {
						spec = Spectrum::FromNamed(p[4]);
				} else if (p.size() == 7) {
						rgb[0] = atof(p[4].c_str());
						rgb[1] = atof(p[5].c_str());
						rgb[2] = atof(p[6].c_str());
						spec = Spectrum::FromRGB(rgb);
				} else {
						cout << "ERROR: addTF tophat string bad: " << addTFstr << endl;
						exit(1128);
				}
				
				float min = atof(p[2].c_str());
				float max = atof(p[3].c_str());
				
				AddTophat(valNum,min,max,spec);
		} else if (p[0] == "constant_table") {
				if (p.size() != 5) {
						cout << "ERROR: addTF constant_table string bad: " << addTFstr << endl;
						exit(1126);
				}
				
				float ctMin = atof(p[3].c_str());
				float ctMax = atof(p[4].c_str());
				
				AddConstantDiscrete(valNum,p[2],ctMin,ctMax);
		} else if (p[0] == "tophat_table") {
				if (p.size() != 7) {
						cout << "ERROR: addTF tophat_table string bad: " << addTFstr << endl;
						exit(1126);
				}
				
				float ctMin = atof(p[3].c_str());
				float ctMax = atof(p[4].c_str());
				float min   = atof(p[5].c_str());
				float max   = atof(p[6].c_str());
				
				AddTophatDiscrete(valNum,p[2],ctMin,ctMax,min,max);
		
		} else if (p[0] == "gaussian_table") {
				if (p.size() != 7) {
						cout << "ERROR: addTF gaussian_table string bad: " << addTFstr << endl;
						exit(1126);
				}
				
				float ctMin = atof(p[3].c_str());
				float ctMax = atof(p[4].c_str());
				float mean  = atof(p[5].c_str());
				float sigma = atof(p[6].c_str());
				
				AddGaussianDiscrete(valNum,p[2],ctMin,ctMax,mean,sigma);
		}
		
		return true;
}
