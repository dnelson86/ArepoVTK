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
	{
	  cameraPosition[i] = Config.cameraPosition[i];
	  cameraLookAt[i] = Config.cameraLookAt[i];
	}
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
	if (qName != "cameraX" && qName != "cameraY" && qName != "cameraZ" && 
	    qName != "lookAtX" && qName != "lookAtY" && qName != "lookAtZ" && 
	    qName != "rotXY" && qName != "rotXZ") {
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
	if( qName == "lookAtX" )
		startVal = Config.cameraLookAt[0];
	if( qName == "lookAtY" )
		startVal = Config.cameraLookAt[1];
	if( qName == "lookAtZ" )
		startVal = Config.cameraLookAt[2];
	if( qName == "rotXY" )
		startVal = -1.0; // angle in radians, override when rotation starts
	if( qName == "rotXZ" )
		startVal = -1.0;
	
	// add to keyframe arrays
	start.push_back(startTime);
	stop.push_back(stopTime);
	start_val.push_back(startVal);
	stop_val.push_back(stopVal);
	method.push_back(intMethod);
	quantity.push_back(qName);
	
	cout << "Added KF: [" << addTFstr << "] now have [" << start.size() << "] keyframes." << endl;
}

void FrameManager::Advance(int curFrame)
{
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
	
		// set start_val for rotation (angle offset) if needed
		float rad = 0.0;

		if( quantity[i] == "rotXY" && start_val[i] < 0.0 )
		{
			start_val[i] = asin( (cameraPosition[0] - cameraLookAt[0]) / rad );
			stop_val[i] += start_val[i];	
			
			rad = sqrt( pow(cameraPosition[0]-cameraLookAt[0],2) +
								  pow(cameraPosition[1]-cameraLookAt[1],2) );
		}
		
		if( quantity[i] == "rotXZ" && start_val[i] < 0.0 )
		{
			start_val[i] = 0.0;
			stop_val[i] += start_val[i];
			
			// this is to re-produce Mark's inside subbox0 movie rotation
			rad = sqrt( pow(cameraPosition[0]-cameraLookAt[0],2) +
								  pow(cameraPosition[2]-cameraLookAt[2],2) );
		}

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
		if( quantity[i] == "lookAtX" )
		  cameraLookAt[0] = curVal;
		if( quantity[i] == "lookAtY" )
		  cameraLookAt[1] = curVal;
		if( quantity[i] == "lookAtZ" )
		  cameraLookAt[2] = curVal;
			
		// "rotXY": do an orbit/rotation in the z-plane as a function of theta (curVal, from 0 to 2pi)
		if( quantity[i] == "rotXY" )
		{
			cout << " curVal : " << curVal << " start_val : " << start_val[i] << " sin: " << sin(curVal+start_val[i]) << endl;
			// progress camera along rotation
			cameraPosition[0] = rad * sin( curVal + start_val[i] ) + cameraLookAt[0];
			cameraPosition[1] = rad * cos( curVal + start_val[i] ) + cameraLookAt[1];
		}
		
		// "rotXZ": do a continual orbit/rotation in the y-plane (curVal, from 0 to 2pi)
		if( quantity[i] == "rotXZ" )
		{
			cout << " curVal : " << curVal << " sin: " << sin(curVal) << endl;
			// progress camera along rotation
			cameraPosition[0] = rad * cos( curVal + start_val[i] ) + cameraLookAt[0];
			cameraPosition[2] = rad * sin( curVal + start_val[i] ) + cameraLookAt[2];
		}
		
	}
	
	cout << "Next frame to render: [" << curFrame << "] " << Config.imageFile << " (time " << curTime 
	     << ") cameraXYZ [" << cameraPosition[0] << " " << cameraPosition[1] << " " << cameraPosition[2] 
	     << "] cameraLookAt [" << cameraLookAt[0] << " " << cameraLookAt[1] << " " << cameraLookAt[2] 
	     << "]" << endl;
}

Transform FrameManager::SetCamera()
{
	Transform world2camera;
	
	world2camera = LookAt(Point(cameraPosition), 
	                      Point(cameraLookAt), 
	                      Vector(Config.cameraUp));
	
	return world2camera;
}

FrameManager::~FrameManager()
{
		IF_DEBUG(cout << "FrameManager() destructor." << endl);
		
}
