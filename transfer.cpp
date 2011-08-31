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
				
				range[0]      = params[0] - 4.0f * params[1]; //cut at +/- 4s
				range[1]      = params[0] + 4.0f * params[1];
				gaussParam[0] = params[0]; //mean
				gaussParam[1] = params[1]; //sigma
				le            = spec[0];
		}

		// unrecognized type / val
		if (type < 0 || type > 3)
				exit(1110);
		if (valNum < 0 || valNum >= TF_NUM_VALS)
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

bool TransferFunction::AddParseString(string &addTFstr)
{
		int valNum;
		Spectrum spec;
		float rgb[3];
		
	  // split string
		vector<string> p;
		string item;
		char delim;
		delim = ' ';
		
		IF_DEBUG(cout << "TF str = " << addTFstr << endl);
		
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
						exit(1124);
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
						exit(1124);
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
						exit(1124);
				}
				
				float min = atof(p[2].c_str());
				float max = atof(p[3].c_str());
				
				AddTophat(valNum,min,max,spec);
		}
		
				
		//Spectrum s1 = Spectrum::FromRGB(Config.rgbEmit);
		//Spectrum s2 = Spectrum::FromNamed("green");
		//Spectrum s3 = Spectrum::FromNamed("blue");
		
		// examples:
		//tf->AddConstant(TF_VAL_DENS,s1);
		//tf->AddTophat(TF_VAL_DENS,5.0,10.0,s1);
		//tf->AddGaussian(TF_VAL_DENS,0.1,0.01,s1);
		
		// "3gaussian"
		//tf->AddGaussian(TF_VAL_DENS,1e-3,2e-4,s1);
		//tf->AddGaussian(TF_VAL_DENS,1e-5,2e-6,s2);
		//tf->AddGaussian(TF_VAL_DENS,1e-7,2e-8,s3);
		
		// "rho_velz"
		//tf->AddGaussian(TF_VAL_DENS,1e-5,2e-6,s2);
		//tf->AddGaussian(TF_VAL_VEL_Z,-10.0,0.1,s1);
		
		// "2gauss"
		//tf->AddGaussian(TF_VAL_DENS,2e-3,3e-4,s1);
		//tf->AddGaussian(TF_VAL_DENS,2e-2,3e-3,s2);

}
