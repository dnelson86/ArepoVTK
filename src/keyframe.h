/*
 * keyframe.h
 * dnelson
 */
 
#ifndef AREPO_RT_KEYFRAME_H
#define AREPO_RT_KEYFRAME_H

#include "transform.h"

class FrameManager {
public:
    // construction
    FrameManager(vector<string> kfSet);
		~FrameManager();

		// return values at current frame
		Transform SetCamera();
		
		// frame management
		void Advance();
		
private:
		void AddParseString(string &addTFstr);
		
    // data
		float curTime;
    int curFrame;
		
		vector<float> start;
		vector<float> stop;
		vector<float> start_val;
		vector<float> stop_val;
		vector<string> method;
		vector<string> quantity;
		
		// quantities
		float cameraPosition[3];
};

#endif
