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
	if (p.size() != 4) {
			cout << "ERROR: addKF string too short: " << addTFstr << endl;
			exit(1164);
	}
	
	// parse
	float startTime  = 0.0;
	float stopTime   = atof(p[0].c_str());
	float startVal   = 0.0;
	float stopVal    = atof(p[2].c_str());
	string qName     = p[1];
	string intMethod = p[3];
	
	// temporary sanity checks
	if (p[1] != "cameraX" && p[1] != "cameraY" && p[1] != "cameraZ") {
		cout << "ERROR: addKF unsupported parameter to keyframe: " << p[1] << endl;
		exit(1165);
	}
	
	if (p[3] != "linear") {
		cout << "ERROR: addKF only supports linear interpolation." << endl;
		exit(1166);
	}
	
	// get starting value (at t=0)
	if( qName == "cameraX" )
		startVal = Config.cameraPosition[0];
	if( qName == "cameraY" )
		startVal = Config.cameraPosition[1];
	if( qName == "cameraZ" )
		startVal = Config.cameraPosition[2];
	
	// add to keyframe arrays
	start.push_back(startTime);
	stop.push_back(stopTime);
	start_val.push_back(startVal);
	stop_val.push_back(stopVal);
	//method.push_back(intMethod);
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
		
		Config.imageFile = "frames_rot/frame_" + num + ".tga";
	}

	// update all quantities for the next frame
	for(unsigned int i=0; i < quantity.size(); i++)
	{
		// is this keyframe not active?
		if( start[i] > curTime || stop[i] < curTime )
			continue;
			
		float duration = (stop[i]-start[i]);
		float change   = (stop_val[i]-start_val[i]);
		float fracTime = (curTime-start[i]) / duration;
		float curVal;
			
		// linear
		curVal = Lerp(fracTime, start_val[i], stop_val[i]);
		
		// easing
		
/*
	myTween.setup(100, 0, ofGetWidth() - 25, Easing::BounceEaseOut);
	setup(duration, start, change, easing func)

num = _ease(_tween.time, _start, _change, _tween.duration);
float Easing::QuadEaseInOut(float t,float b , float c, float d)
{
if ((t/=d/2) < 1) return ((c/2)*(t*t)) + b;
return -c/2 * (((t-2)*(--t)) - 1) + b;
}
*/
		//if ((fracTime/=1.0/2) < 1)
		//	curVal = ((change/2)*(fracTime*fracTime)) + start_val[i];
		//else
		//	curVal = -change/2 * (((fracTime-2)*(fracTime)) - 1) + start_val[i];
	
		// what quantity are we updating?
		if( quantity[i] == "cameraX" )
		  cameraPosition[0] = curVal;
		if( quantity[i] == "cameraY" )
		  cameraPosition[1] = curVal;
		if( quantity[i] == "cameraZ" )
		  cameraPosition[2] = curVal;
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
