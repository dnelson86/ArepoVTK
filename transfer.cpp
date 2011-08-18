/*
 * transfer.cpp
 * dnelson
 */
 
#include "transfer.h"

TransferFunction::TransferFunction()
{
		IF_DEBUG(cout << "TransferFunction() constructor." << endl);
		numFuncs = 0;
		
		//Spectrum Le = Spectrum::FromRGB(Config.rgbEmit);
}

Spectrum TransferFunction::Lve(const float &rho, const float &utherm) const
{
		float rgb[3];

		// testing only
		rgb[0] = rho/100.0;
		rgb[1] = rho/100.0;
		rgb[2] = rho/100.0;
			
		Spectrum Lve = Spectrum::FromRGB(rgb);
		
		return Lve;
}

TransferFunction *CreateTransferFunction()
{
    return new TransferFunction();
}
