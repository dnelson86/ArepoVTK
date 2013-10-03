/*
 * keyframe.cpp
 * dnelson
 */

#include <sstream> 
#include "util.h"
#include "keyframe.h"

FrameManager::FrameManager(vector<string> kfSet)
{
		IF_DEBUG(cout << "FrameManager() constructor." << endl);
		
		// parse keyframe strings
		for(unsigned int i = 0; i < kfSet.size(); i++)
			FrameManager::AddParseString( kfSet[i] );
			
		// init
		curTime  = -1;
		curFrame = -1;
		
		for(int i=0; i < 3; i++)
		  cameraPosition[i] = Config.cameraPosition[i];
}

void FrameManager::AddParseString(string &addTFstr)
{
	// split string
	vector<string> p;
	string item;
	char delim;
	delim = ' ';
	
	stringstream ss(addTFstr);
	while(getline(ss, item, delim))
			p.push_back(item);
			
	// check string size
	if (p.size() != 5) {
			cout << "ERROR: addKF string too short: " << addTFstr << endl;
			exit(1164);
	}
	
	// parse
	float startTime  = atof(p[0].c_str());
	float stopTime   = atof(p[1].c_str());
	float startVal   = 0.0;
	float stopVal    = atof(p[3].c_str());
	string qName     = p[2];
	string intMethod = p[4];
	
	// temporary sanity checks
	if (qName != "cameraX" && qName != "cameraY" && qName != "cameraZ" && qName != "rotXY") {
		cout << "ERROR: addKF unsupported parameter to keyframe: " << qName << endl;
		exit(1165);
	}
	
	if (intMethod != "linear" && intMethod != "quadratic_inout") {
		cout << "ERROR: addKF unsupported interpolation method: " << intMethod << endl;
		exit(1167);
	}
	
	// get starting value (at t=0)
	if( qName == "cameraX" )
		startVal = Config.cameraPosition[0];
	if( qName == "cameraY" )
		startVal = Config.cameraPosition[1];
	if( qName == "cameraZ" )
		startVal = Config.cameraPosition[2];
	if( qName == "rotXY" )
		startVal = 0.0; // angle in radians
	
	// add to keyframe arrays
	start.push_back(startTime);
	stop.push_back(stopTime);
	start_val.push_back(startVal);
	stop_val.push_back(stopVal);
	method.push_back(intMethod);
	quantity.push_back(qName);
	
	cout << "Added KF: [" << addTFstr << "] now have [" << start.size() << "] keyframes." << endl;
}

void FrameManager::Advance()
{
	curFrame++;
	curTime = curFrame * Config.timePerFrame;
	
	// update imageFilename
	if( Config.numFrames > 1 )
	{
		std::ostringstream s;
		s << curFrame;
		const char padChar = '0';
		string num = s.str();
		if( num.size() < 4 ) // output image numbering
		num.insert(0, 4-num.size(), padChar);
		
		size_t found = Config.imageFile.find_last_of("/");
		Config.imageFile = Config.imageFile.substr(0,found) + "/frame_" + num + ".tga";
	}

	// update all quantities for the next frame
	for(unsigned int i=0; i < quantity.size(); i++)
	{
		// is this keyframe not active?
		if( start[i] > curTime || stop[i] < curTime )
			continue;
			
		float duration = (stop[i]-start[i]);
		float fracTime = (curTime-start[i]) / duration;
		float curVal = 0;
			
		if( method[i] == "linear" )
		{
		  curVal = Lerp(fracTime, start_val[i], stop_val[i]);
		}
		else if( method[i] == "quadratic_inout" )
		{
		  float intVal = 0;
		  if( fracTime < 0.5 )
		    intVal = 2.0 * fracTime * fracTime;
		  if( fracTime >= 0.5 )
		    intVal = (-2.0 * fracTime * fracTime) + (4 * fracTime) - 1.0;

		  curVal = start_val[i] + (stop_val[i]-start_val[i]) * intVal;
		}
		
		// what quantity are we updating?
		if( quantity[i] == "cameraX" )
		  cameraPosition[0] = curVal;
		if( quantity[i] == "cameraY" )
		  cameraPosition[1] = curVal;
		if( quantity[i] == "cameraZ" )
		  cameraPosition[2] = curVal;
			
		// "rotXY": do an orbit/rotation in the z-plane as a function of theta (curVal, from 0 to 2pi)
		if( quantity[i] == "rotXY" )
		{
			// currently hardcoded to do spoonHD ending rotation
			float rad = 3.0;
			cameraPosition[0] = -1.0 * rad * cosf( curVal );
			cameraPosition[1] = 0.5 + rad * sinf( curVal );
		}
		
	}
	
	cout << "Next frame to render: " << curFrame << " (time " << curTime << ") cameraXYZ [" 
	     << cameraPosition[0] << " " << cameraPosition[1] << " " << cameraPosition[2] << "]" << endl;
	cout << "Next frame: " << Config.imageFile << endl;
}

Transform FrameManager::SetCamera()
{
	Transform world2camera;
	
	world2camera = LookAt(Point(cameraPosition), 
	                      Point(Config.cameraLookAt), 
												Vector(Config.cameraUp));
	
	return world2camera;
}

FrameManager::~FrameManager()
{
		IF_DEBUG(cout << "FrameManager() destructor." << endl);
		
		
}
