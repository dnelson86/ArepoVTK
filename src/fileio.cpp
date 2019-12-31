/*
 * fileio.cpp
 * dnelson
 */

#include "fileio.h"
#include "fileio_img.h"

// ConfigSet
void ConfigSet::ReadFile(string cfgfile)
{
		cout << "Reading Configuration File [" << cfgfile << "]." << endl;
		
		ifstream is(cfgfile.c_str());
		if (!is) {
				cout << "ConfigSet: ERROR opening configuration file [" << cfgfile << "]." << endl;
				exit(1113);
		}
		
		// prepare to parse
		const string::size_type skip = delim.length();
		string nextLine = "";
		
		static const char whitespace[] = " \n\t\v\r\f";
		
		while (is || nextLine.length() > 0) {
				// read line
				string line = "";
				
				if (nextLine.length() > 0) {
						line     = nextLine;
						nextLine = "";
				} else {
						getline(is,line);
				}
				
				// ignore inline and full line comments
				line = line.substr(0,line.find(comment));
				
				// split into key = value pair
				string::size_type offset = line.find(delim);
				
				if (offset < string::npos) {
						// get key
						string key = line.substr(0,offset);
						line.replace(0,offset+skip,"");
				
						// trim
						key.erase( 0, key.find_first_not_of(whitespace) );
						key.erase( key.find_last_not_of(whitespace) + 1U );
						
						line.erase( 0, line.find_first_not_of(whitespace) );
						line.erase( line.find_last_not_of(whitespace) + 1U );
				
						// keyframe and transfer function handling (multiple entries for each)
						if (key.substr(0,5) == "addTF")
							tfSet.push_back(line);
						if (key.substr(0,5) == "addKF")
							kfSet.push_back(line);
				
						// store key,value (map type, keys unique)
						parsedParams[key] = line;
						//IF_DEBUG(cout << " [" << key << "] = " << line << endl);
				}
		}
		
		// Input/Output
		imageFile     = readValue<string>("imageFile",  "frame.png");
		rawRGBFile    = readValue<string>("rawRGBFile", "frame.raw.txt");
		filename      = readValue<string>("filename", "none");
		paramFilename = readValue<string>("paramFilename");
		writeRGB8bit  = readValue<bool>("writeRGB8bit",  true); // png
		writeRGB16bit = readValue<bool>("writeRGB16bit", false); // png
		writeRGBA8bit = readValue<bool>("writeRGBA8bit", false); // png

		dumpMeshText   = readValue<bool>("dumpMeshText", false);
		dumpMeshBinary = readValue<bool>("dumpMeshBinary", false);
		dumpMeshCells  = readValue<bool>("dumpMeshCells", false);
		
		// General
		nTasks        = readValue<int> ("nTasks",        1);
		nCores        = readValue<int> ("nCores",        0);
		quickRender   = readValue<bool>("quickRender",   false);
		openWindow    = readValue<bool>("openWindow",    false);
		verbose       = readValue<bool>("verbose",       false);

		// Job Sub-Divison
		totNumJobs    = readValue<int>   ("totNumJobs", 0);
		curJobNum     = -1; // read from commandline
		maskFileBase  = readValue<string>("maskFileBase", "");
		maskPadFac    = readValue<float> ("maskPadFac", 0.0f);
		
		jobExpansionFac = readValue<int> ("jobExpansionFac", 1);
		expandedJobNum  = -1; // read from commandline

		// Frame/Camera
		imageXPixels  = readValue<int>  ("imageXPixels", 500);
		imageYPixels  = readValue<int>  ("imageYPixels", 500);
		swScale       = readValue<float>("swScale",      1.0f);
		cameraType    = readValue<string>("cameraType",  "ortho");
		cameraFOV     = readValue<float>("cameraFOV",    0.0f); // degrees
		splitStrArray( readValue<string>("cameraPosition") , &cameraPosition[0]   );
		splitStrArray( readValue<string>("cameraLookAt")   , &cameraLookAt[0] );
		splitStrArray( readValue<string>("cameraUp")       , &cameraUp[0] );
	
		// Data Processing
		readPartType          = readValue<int>("readPartType", 0);
		splitStrArray( readValue<string>("recenterBoxCoords", "-1 -1 -1") , &recenterBoxCoords[0] );
		convertUthermToKelvin = readValue<bool>("convertUthermToKelvin", false);
		takeLogUtherm         = readValue<bool>("takeLogUtherm", false);
		takeLogDens           = readValue<bool>("takeLogDens", false);

		// Animation
		startFrame    = readValue<int>("startFrame",     0);	
		numFrames     = readValue<int>("numFrames",      1);
		timePerFrame  = readValue<float>("timePerFrame", 1.0);
		minScale      = readValue<float>("minScale", -1.0);
		maxScale      = readValue<float>("maxScale", -1.0);
		minAlpha      = readValue<float>("minAlpha", -1.0);
		maxAlpha      = readValue<float>("maxAlpha", -1.0);
	
		// Render
		drawBBox      = readValue<bool>("drawBBox",      true);
		drawTetra     = readValue<bool>("drawTetra",     true);	
		drawVoronoi   = readValue<bool>("drawVoronoi",   false);
		drawSphere    = readValue<bool>("drawSphere",    false);

		projColDens   = readValue<bool>("projColDens",     false); // write raw values
		nTreeNGB      = readValue<int>("nTreeNGB",             0); // disabled by default
		viStepSize    = readValue<float>("viStepSize",      0.0f); // disabled by default
		rayMaxT       = readValue<float>("rayMaxT",         0.0f);
		
		// rgb triplets input		
		splitStrArray( readValue<string>("rgbLine",     "0.1  0.1  0.1")  , &rgbLine[0]    );
		splitStrArray( readValue<string>("rgbTetra",    "0.01 0.01 0.01") , &rgbTetra[0]   );
		splitStrArray( readValue<string>("rgbVoronoi",  "0.0  0.05 0.0")  , &rgbVoronoi[0] );
		splitStrArray( readValue<string>("rgbAbsorb",   "0.0  0.0  0.0")  , &rgbAbsorb[0]  );
		
		// basic validation
		if (!tfSet.size())
			terminate("Config: no TFs specified, going to be a very boring image.");
		//if (projColDens && !totNumJobs)
		//	terminate("Config: ERROR! projColDens only with totNumJobs>0 (custom load).");
		if (totNumJobs < 0 || totNumJobs > 16*16)
			terminate("Config: ERROR! Strange totNumJobs value, should be >=0 and <256");
		if (jobExpansionFac != 1 && jobExpansionFac != 2 && jobExpansionFac != 4)
			terminate("Config: ERROR! Odd choice of jobExpansionFac.");
			
		// data read/TFs
		if (readPartType != 0 && readPartType != 1)
			terminate("Config: ERROR! Unsupported readPartType.");
		if (takeLogDens && (rgbAbsorb[0] > 0.0 || rgbAbsorb[1] > 0.0 || rgbAbsorb[2] > 0.0))
			terminate("Config: WARNING: Will be using log(density) weighting due to nonzero absorption (maybe ok).");
			
		// render setup validation
		if (viStepSize == 0.0 && nTreeNGB)
			terminate("Config: ERROR! Need to specify viStepSize!=0 if nTreeNGB>0.");
#if !defined(NATURAL_NEIGHBOR_IDW) && !defined(NATURAL_NEIGHBOR_SPHKERNEL)
		if (nTreeNGB)
			terminate("Config: ERROR! Must enable IDW or SPHKERNEL for nTreeNGB>0.");
#endif
			
		// camera type mappings
		if (cameraType == "ortho") { cameraType = "orthographic"; }
		if (cameraType == "persp") { cameraType = "perspective"; }
		if (cameraType == "env")   { cameraType = "environmental"; }
		if (cameraType == "fish")  { cameraType = "fisheye"; }
		if (cameraType == "rift")  { cameraType = "oculusrift"; }
			
		// camera validation
		if (cameraType == "fisheye" && imageXPixels != imageYPixels)
			terminate("Config: ERROR! Fisheye camera only supports square images.");
		if (cameraType == "environmental" && imageXPixels != 2*imageYPixels)
			terminate("Config: ERROR! Environment camera requires image width equal to twice the height.");
		if ((cameraType == "fisheye" || cameraType == "environmental") && swScale != 1.0 )
			terminate("Config: ERROR! swScale not used for fisheye or environmental (leave at 1.0).");
		if (cameraType == "environmental" && cameraFOV != 360.0)
			terminate("Config: ERROR! Environment camera requires 360 degree fov.");
		if (cameraType == "orthographic" && cameraFOV != 0.0)
			terminate("Config: ERROR! FOV not used for ortho camera (leave at 0.0).");
		if (cameraType == "perspective" && (cameraFOV <= 0.0 || cameraFOV >= 180.0))
			terminate("Config: ERROR! Perspective camera expects sane FOV.");
			
		// validation not directly related to config file
#if defined(NATURAL_NEIGHBOR_INTERP) && !defined(NATURAL_NEIGHBOR_INNER)
		terminate("NNI without NN_INNER likely gives incorrect weights.");
#endif
}

void ConfigSet::print()
{
		map_i pi;
		
		cout << endl << "CONFIGURATION PARAMETERS USED:" << endl << endl;
		for (pi = parsedParams.begin(); pi != parsedParams.end(); ++pi) {
				cout << " " << setw(22) << pi->first << " " << delim << " " << pi->second << endl;
		}
		cout << endl;

}

template<class T> T ConfigSet::readValue(const string &key) const
{
		// no default value = required key, die if missing
		map_ci pi = parsedParams.find(key);
		
		if( pi == parsedParams.end() )
			terminate(" ERROR: Required parameter [%s] missing from configuration file.",key.c_str());
		
		return string_to_T<T>( pi->second );
}

template<class T> T ConfigSet::readValue(const string &key, const T &defaultValue) const
{
		// default value supplied = optional param in config file
		map_ci pi = parsedParams.find(key);
		
		if( pi == parsedParams.end() )
				return defaultValue;
				
		return string_to_T<T>( pi->second );
}

template<class T> T ConfigSet::string_to_T(const string &str)
{
		T t;
		istringstream ist(str);
		
		ist >> t; // use >> operator to typecast to T
		return t;
}

template<> inline bool ConfigSet::string_to_T<bool>(const string &str)
{
		// special case conversion to bool
		bool ret   = true;
		string strUpper = str;
		
		// uppercase
		for(string::iterator p = strUpper.begin(); p != strUpper.end(); ++p)
				*p = toupper(*p);
			
		if(strUpper == string("FALSE") || strUpper == string("NO") || strUpper == string("0"))
				ret = false;
			
		return ret;
}

template<> inline string ConfigSet::string_to_T<string>(const string &str)
{
	// special case conversion to string, otherwise spaces get lost with >>
	return str;
}

void ConfigSet::splitStrArray(const string str, float *rgb) //size=3
{
		char *t;
	  t = new char [str.size()+1];
		strcpy (t, str.c_str());

		char *pch = strtok((char*)t, " ,");
		for (int i=0; i < 3; i++) {
				rgb[i] = atof(pch);
				pch = strtok(NULL," ,");
		}
		delete[] t;
}

bool ReadFloatFile(const char *filename, vector<float> *values)
{
    FILE *f = fopen(filename, "r");
    if (!f) {
        return false;
    }

    int c;
    bool inNumber = false;
    char curNumber[32];
    int curNumberPos = 0;
    int lineNumber = 1;
    while ((c = getc(f)) != EOF) {
        if (c == '\n') ++lineNumber;
        if (inNumber) {
            if (isdigit(c) || c == '.' || c == 'e' || c == '-' || c == '+')
                curNumber[curNumberPos++] = c;
            else {
                curNumber[curNumberPos++] = '\0';
                values->push_back(atof(curNumber));
                //Assert(curNumberPos < (int)sizeof(curNumber));
                inNumber = false;
                curNumberPos = 0;
            }
        }
        else {
            if (isdigit(c) || c == '.' || c == '-' || c == '+') {
                inNumber = true;
                curNumber[curNumberPos++] = c;
            }
            else if (c == '#') {
                while ((c = getc(f)) != '\n' && c != EOF)
                    ;
                ++lineNumber;
            }
            else if (!isspace(c)) {
								cout << "Unexpected text found at line " << lineNumber << " of float file "
								     << filename << endl;
            }
        }
    }
    fclose(f);
    return true;
}

bool WriteFloatFile(const char *filename, float *values, int nx, int ny)
{
		// c - text
    FILE *f = fopen(filename, "w");
    if (!f) {
        cout << "WriteFloatFile: Unable to open file " << filename << endl;
        return false;
		}
		
		IF_DEBUG(cout << "WriteFloatFile(" << filename << ") nx = " << nx << " ny = " << ny << endl);
		
		int y, x;
		int offset = 0;
		float amp = 0.0f;
		
		fprintf(f, "%d\n%d\n", nx, ny);
		
    for (y = 0; y < ny; ++y) {
        for (x = 0; x < nx; ++x) {
						fprintf(f, "(%f %f %f)", values[3*offset  ],
						                         values[3*offset+1],
																		 values[3*offset+2]);
						offset++;
				}
				fprintf(f, "\n");
		}
		
		fprintf(f, "amp\n");
		
		offset = 0;
		
    for (y = 0; y < ny; ++y) {
        for (x = 0; x < nx; ++x) {
						amp = sqrt( values[3*offset]*values[3*offset] + values[3*offset+1]*values[3*offset+1]
											+ values[3*offset+2]*values[3*offset+2] );
						fprintf(f, "%f ", amp);
						offset++;
				}
				fprintf(f, "\n");
		}
		
    fclose(f);
		
		// c++ - binary
/*		string fname = filename;
		fname += ".bin";
		ofstream f2(fname.c_str(), ofstream::binary);
		
		if (!f2.is_open()) {
        cout << "WriteFloatFile: Unable to open file " << filename << endl;
        return false;
		}
		
		cout << "WriteFloatFile(" << fname << ") nx = " << nx << " ny = " << ny << endl;		
		
		f2.write((const char*)&nx,sizeof(int));
		f2.write((const char*)&ny,sizeof(int));
		f2.write((const char*)&values,sizeof(float)*3*nx*ny);
		f2.close(); */

    return true;
}

int loadDiscreteColorTable(const string &filename, vector<float> *vals)
{
		// set colortables path
		string filepath = "colortables/" + filename + ".tbl";
		
		if (!ReadFloatFile(filepath.c_str(), vals)) {
				cout << "ERROR: Unable to open color table file: '" << filepath << "', exiting." << endl;
				exit(1161);
		}

		int numEntries = (int)(*vals)[0];
		
		if (numEntries != ((int)vals->size()-1)/4) {
				cout << "ERROR: Corrupt color table file: '" << filepath << "', exiting." << endl;
				exit(1162);
		}
			
		// pop the first entry (the count of the number of entries in this color table)
		vals->erase(vals->begin());

		return numEntries;
}

int parseSceneFile(const string &filename, int &nx, int &ny, int &nz, vector<float> *data)
{
		int ni = 0;
		vector<float> vals;
		
		if (!ReadFloatFile(filename.c_str(), &vals))
				return 0;

		// size
		if (vals.size() >= 4) {
			ni = (int)vals[0];
			nx = (int)vals[1];
			ny = (int)vals[2];
			nz = (int)vals[3];
		}

		// errorcheck size
		if (ni != nx*ny*nz) {
			cout << "Error: parseSceneFile ni = " << ni << " but nx*ny*nz = " << nx*ny*nz << endl;
			return 0;
		}
		if (ni != (int)ni || nx != (int)nx || ny != (int)ny || nz != (int)nz) {
			cout << "Error: parseSceneFile non-integer size, ni nx ny nz = " 
					 << ni << " " << nx << " " << ny << " " << nz << endl;
			return 0;
		}
	
#if defined(DEBUG) && 0
		cout << "parseSceneFile volumeGrid: ni = " << ni 
		              << " nx = " << nx << " ny = " << ny << " nz = " << nz << endl;
									
		cout << "parseSceneFile volumeGrid data:";
		for (unsigned int i=4; i < vals.size(); i++)
		  cout << " " << vals[i];
		cout << endl;
#endif
		
		for (unsigned int i=4; i < vals.size(); i++)
			data->push_back(vals[i]);
		
		return ni;
}

// Image Output

void WriteImage(const string &name, float *pixels, float *alpha, int xRes,
                int yRes, int totalXRes, int totalYRes, int xOffset, int yOffset)
{
    if (name.size() >= 5) {
        uint32_t suffixOffset = name.size() - 4;
				
        if (!strcmp(name.c_str() + suffixOffset, ".tga") ||
            !strcmp(name.c_str() + suffixOffset, ".TGA")) {
            WriteImageTGA(name, pixels, alpha, xRes, yRes, totalXRes,
                          totalYRes, xOffset, yOffset);
            return;
        }
				
        if (!strcmp(name.c_str() + suffixOffset, ".png") ||
            !strcmp(name.c_str() + suffixOffset, ".PNG")) {
            WriteImagePNG(name, pixels, alpha, xRes, yRes, totalXRes,
                          totalYRes, xOffset, yOffset);
            return;
        }
    }
    cout << "WriteImage: Can't determine image file type from suffix of filename '" << name << "'!" << endl;
}
