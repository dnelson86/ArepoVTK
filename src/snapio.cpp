/*
 * snapio.cpp
 * dnelson
 */

#include <algorithm> // min_element/max_element
 
#include "snapio.h"
#include "util.h"
#include "geometry.h"
#include "transform.h"
#include "camera.h"

void ArepoSnapshot::read_ic()
{
	unsigned int i;
	
	// locate file(s) and store list of all snapshot files
	vector<string> snapFilenames;
	vector<unsigned int> numFiles;
	string maskFileName;
	
	string f1 = fileBase + ".hdf5";
	string f2 = fileBase + ".0.hdf5";
	
	if( ifstream(f1.c_str()) ) {
		// single file
		snapFilenames.push_back(f1);
		
		// verify with header
		readGroupAttribute( f1, "Header", "NumFilesPerSnapshot", numFiles );
		if( numFiles[0] != 1 ) {
			cout << "Error! Strange numFiles value: [" << numFiles[0] << "] when expecting 1.";
			exit(1172);
		}
	}	else {
		// multiple files, read total number from hdf5 header of first
		readGroupAttribute( f2, "Header", "NumFilesPerSnapshot", numFiles );
		cout << "Found: [" << numFiles[0] << "] files." << endl;
		
		if( numFiles[0] <= 0 || numFiles[0] > 1024 ) {
			cout << "Error! Strange numFiles value: [" << numFiles[0] << "].";
			exit(1171);
		}
		
		for( i=0; i < numFiles[0]; i++ )
		  snapFilenames.push_back( fileBase + toStr(i) + ".hdf5" );

	}
	
	// if running multiple jobs/doing selective load, load maskFile now, or if it does not exist, create it now then exit
	if( Config.totNumJobs >= 1 )
	{	
		// file format: maskFileBase_hash(hashKey).hdf5
		string hashKey = toStr(Config.totNumJobs) + "_" + toStr(Config.maskPadFac) + "_" + fileBase + "_" + 
										 toStr(Config.rayMaxT) + "_" + toStr(Config.cameraFOV) + "_" + toStr(Config.cameraPosition[0]) + 
										 "-" + toStr(Config.cameraPosition[1]) + "-" + toStr(Config.cameraPosition[2]) + "_" + 
										 toStr(Config.swScale) + "-" + toStr(Config.cameraLookAt[0]) + "-" + toStr(Config.cameraLookAt[1]) + "-" + 
										 toStr(Config.cameraLookAt[2]) + "_" + toStr(Config.cameraUp[0]) + "-" + toStr(Config.cameraUp[1]) + 
										 "-" + toStr(Config.cameraUp[2]) + "_" + toStr(float(Config.imageXPixels)/float(Config.imageYPixels));
		
		maskFileName = Config.maskFileBase + "_" + toStr(hash_djb2(hashKey)) + ".hdf5";
	
		if ( !ifstream(maskFileName.c_str()) ) {
			cout << "Mask file [" << maskFileName << "] not found, making it now." << endl;

			ArepoSnapshot::makeNewMaskFile( maskFileName, snapFilenames );
			
			cout << "Mask file written, exiting." << endl << endl;
			exit(0);
		} else {
			cout << "Reading existing mask file [" << maskFileName << "]." << endl;
		}
	}
	
	// selectively load from all snapshot chunks using the maskfile
	ArepoSnapshot::loadAllChunksWithMask( maskFileName, snapFilenames );
	
}

void ArepoSnapshot::loadAllChunksWithMask( string maskFileName, vector<string> snapFilenames )
{
	unsigned int i, j, k;
	int offset = 0;
	
	// existing maskfile: load Header/partCountsTot
	vector<long long> partCountsTot;
	
	readGroupDataset( maskFileName, "Header", "PartCountsTot", -1, partCountsTot );
	
	if( Config.verbose )
		for( i=0; i < partCountsTot.size(); i++ )
			cout << " Job [" << i << "] partCountsTot = " << partCountsTot[i] << endl;
	
	// allocate (P/SphP)
	All.MaxPart = partCountsTot[Config.curJobNum] / (1.0 - 2 * ALLOC_TOLERANCE);
	All.MaxPartSph = partCountsTot[Config.curJobNum] / (1.0 - 2 * ALLOC_TOLERANCE);
	
	allocate_memory();	
	
	// loop over each chunk
	vector<int> jobIndList;
	vector<float> quantity;	
	vector<int> quantity_id;
	
	for( i=0; i < snapFilenames.size(); i++ )
	{
		jobIndList.clear();
		
		// load File_Y/Job_X_IndList where X=curJobNum
		string groupName = "File_" + toStr(i);
		string objName   = "Job_PartCounts";
		
		vector<unsigned int> partCounts;
		
		readGroupDataset( maskFileName, groupName, objName, -1, partCounts );
		
		if( partCounts.size() != (unsigned int) Config.totNumJobs ) {
			cout << "Error: loadAllChunksWithMask: partCounts.size() = " << partCounts.size() << endl;
			exit(1190);
		}
		
		// reading nothing from this file? then skip
		if( partCounts[Config.curJobNum] == 0 ) {
			if( Config.verbose )
				cout << "File [" << i << ": " << snapFilenames[i] << "] empty load, skipping." << endl;
			
			continue;
		}
		
		// reserve enough space for indice list and float quantity buffer
		jobIndList.reserve( partCounts[Config.curJobNum] );
		quantity.reserve( partCounts[Config.curJobNum] );
		quantity_id.reserve( partCounts[Config.curJobNum] );
		
		if( Config.verbose )
			cout << "File [" << i << ": " << snapFilenames[i] << "] loading [" << partCounts[Config.curJobNum] << "] indices." << endl;
		
		// load index list, verify partCounts[curJobNum] equals the size of the loaded job index list
		objName = "Job_" + toStr(Config.curJobNum) + "_IndList";
		
		readGroupDataset( maskFileName, groupName, objName, -1, jobIndList );
		
		if( jobIndList.size() != partCounts[Config.curJobNum] ) {
			cout << "Error: loadAllChunksWithMask: Index size mismatch with expected." << endl;
			exit(1191);
		}
		
		// init type
		for( j=0; j < jobIndList.size(); j++ )
			P[offset + j].Type = 0;
		
		// read fields using hdf5 point selection, then transfer buffer into P/SphP starting at offset
		// (ParticleIDs, Coordinates, Density, Mass, Utherm, Velocity)
		// based on compile flags: Ne, Metallicity, Sfr
		// note these are all for Gas, for meshing+rendering of different types, will have to change this logic
		groupName = "PartType" + toStr(READ_PARTTYPE);
		
		// make 1D serialized index list of hsize_t type since we will use it several times
		vector<hsize_t> coord;
		coord.reserve( jobIndList.size() );
		
		for( j=0; j < jobIndList.size(); j++ )
			coord.push_back(jobIndList[j]);
		
		// Coordinates (X)
		for( k=0; k < 3; k++ )
		{
			readGroupDatasetSelect( snapFilenames[i], groupName, "Coordinates", coord, k, quantity );
			
			if( Config.verbose ) {
				cout << " pos_" << k << ".size() = " << quantity.size();
				cout << " min = " << *min_element(quantity.begin(), quantity.end());
				cout << " max = " << *max_element(quantity.begin(), quantity.end()) << endl;
			}
			
			for( j=0; j < quantity.size(); j++ )
				P[offset + j].Pos[k] = quantity[j];
		}
		
		// ParticleIDs
		readGroupDatasetSelect( snapFilenames[i], groupName, "ParticleIDs", coord, -1, quantity_id );
		
		for( j=0; j < quantity_id.size(); j++ )
			P[offset + j].ID = quantity_id[j];
		
		// Velocities
		for( k=0; k < 3; k++ )
		{
			readGroupDatasetSelect( snapFilenames[i], groupName, "Velocities", coord, k, quantity );
			
			for( j=0; j < quantity.size(); j++ )
				P[offset + j].Vel[k] = quantity[j];
		}
		
		// Mass
		readGroupDatasetSelect( snapFilenames[i], groupName, "Masses", coord, -1, quantity );
		
		for( j=0; j < quantity.size(); j++ )
			P[offset + j].Mass = quantity[j];
			
		// Density
		readGroupDatasetSelect( snapFilenames[i], groupName, "Density", coord, -1, quantity );
		
		for( j=0; j < quantity.size(); j++ )
			SphP[offset + j].Density = quantity[j];
			
		// Utherm
		readGroupDatasetSelect( snapFilenames[i], groupName, "InternalEnergy", coord, -1, quantity );
		
		for( j=0; j < quantity.size(); j++ )
			SphP[offset + j].Utherm = quantity[j];

		// increment global snapshot offset as we move to next chunk
		offset += partCounts[Config.curJobNum];
	
	} // snapFilenames
	
	// handle other things read_ic() sets (global/local total particle counts, massTable, time/cosmo factors)
	All.TotNumGas = partCountsTot[Config.curJobNum];
	All.TotNumPart = partCountsTot[Config.curJobNum];
	
	//for(i = 0; i < NTYPES; i++)
	//	All.MassTable[i] = header.mass[i];
	
	vector<float> snapTime;
	readGroupAttribute( snapFilenames[0], "Header", "Time", snapTime );
	All.Time = All.TimeBegin = snapTime[0];
	set_cosmo_factors_for_current_time();
	
	NumPart = partCountsTot[Config.curJobNum];
	NumGas = partCountsTot[Config.curJobNum];
	
	// verify our total particle counts equal partCountsTot[curJobNum]
	if( offset != partCountsTot[Config.curJobNum] ) {
		cout << "Error: Failed to read the expected number of particles." << endl;
		exit(1192);
	}
	
}

void ArepoSnapshot::calcSubCamera(struct cameraParams &cP, int jobNum)
{
	float screen[4];
	
	Filter *filter       = CreateBoxFilter();
	Film *film           = CreateFilm(filter);
	
	film->CalculateScreenWindow(screen, jobNum);
			
	// calculate x,y bounds (these, as with Z bounds, are all relative to camera position, in code units)
	cP.yMin = screen[2] - Config.maskPadFac;
	cP.yMax = screen[3] + Config.maskPadFac;
	cP.xMin = screen[0] - Config.maskPadFac;
	cP.xMax = screen[1] + Config.maskPadFac;
	
	// also do z-bounds (same for all sub-cameras)
	cP.zMin = 0 - Config.maskPadFac;
	cP.zMax = Config.rayMaxT + Config.maskPadFac; // TODO, should use min(rayMaxT,dist to exiting box)

	if( Config.verbose )
		cout << " zMinMax [" << cP.zMin << " " << cP.zMax << "] yMinMax [" << cP.yMin << " " << cP.yMax 
				 << "] xMinMax [" << cP.xMin << " " << cP.xMax << "]" << endl;	
	
	delete film;
}

void ArepoSnapshot::makeNewMaskFile( string maskFileName, vector<string> snapFilenames )
{
	unsigned int i, j, k;
	
	// initialize sub-cameras and verify parameters
	vector<struct cameraParams> camParams;
	
	for( j=0; j < (unsigned int)Config.totNumJobs; j++ ) {
		cameraParams cP;
		calcSubCamera(cP,j);
		camParams.push_back(cP);
	}
	
	// X,Y,Z for camera referential system
	Vector Z = Vector(Config.cameraLookAt) - Vector(Config.cameraPosition);
	Z /= Z.Length(); // normalize
	
	Vector X = Cross( Z, Vector(Config.cameraUp) );
	X /= X.Length(); //normalize
	
	Vector Y = Cross( X, Z );
	
	if( Config.verbose ) {
		X.print("X: ");
		Y.print("Y: ");
		Z.print("Z: ");	
	}
	
	Vector camPos(Config.cameraPosition);
	
	// get total particle counts in all files
	vector<unsigned long long> numPartTotal_low, numPartTotal_high, numPartTotal;
	unsigned long long numReadTot = 0;
	
	readGroupAttribute( snapFilenames[0], "Header", "NumPart_Total", numPartTotal_low );
	readGroupAttribute( snapFilenames[0], "Header", "NumPart_Total_HighWord", numPartTotal_high );
	
	if( numPartTotal_low.size() != numPartTotal_high.size() ) {
		cout << "Error: Size mismatch between low and highwords for numPartTotal!" << endl;
		exit(1173);
	}
	
	numPartTotal.reserve( numPartTotal_low.size() );
	numPartTotal.assign( numPartTotal.size(), 0 );
	
	for( i=0; i < numPartTotal_low.size(); i++ ) {
		numPartTotal[i] = numPartTotal_low[i] + (((long long) numPartTotal_high[i]) << 32);
		if( Config.verbose )
			cout << " PartType[" << i << "] numPartTotal = " << numPartTotal[i] << endl;
	}
	
	// make HDF5 mask file (empty)
	createNewFile( maskFileName );	

	// hold total number of particles that each job is responsible for
	vector<long long> partCountsTot ( Config.totNumJobs, 0 );
	
	string groupName = "PartType" + toStr(READ_PARTTYPE);
	
	// allocate space for approximate size of Coordinates per chunk
	vector<float> pos_x, pos_y, pos_z;
	int allocNum = round(numPartTotal[READ_PARTTYPE] / snapFilenames.size() * 1.5);
	
	pos_x.reserve( allocNum );
	pos_y.reserve( allocNum );
	pos_z.reserve( allocNum );
	
	int fillDivisor = floor( numPartTotal[READ_PARTTYPE] * snapFilenames.size() / FILL_NUMPART );
	cout << " fillDivisor: " << fillDivisor << endl;
	
	// loop over each snapshot chunk
	for( i=0; i < snapFilenames.size(); i++ )
	{
		vector<int> partCounts ( Config.totNumJobs, 0 );
		int partCountsAllJobs = 0;
		int count_x, count_y, count_z;
		
		// make group in mask file for this chunk
		string groupNameSave = "File_" + toStr(i);
		createNewGroup( maskFileName, groupNameSave );
		
		// load x,y,z coordinates directly into their vectors
		count_x = readGroupDataset( snapFilenames[i], groupName, "Coordinates", 0, pos_x );
		count_y = readGroupDataset( snapFilenames[i], groupName, "Coordinates", 1, pos_y );
		count_z = readGroupDataset( snapFilenames[i], groupName, "Coordinates", 2, pos_z );
	
		if( count_x != count_y || count_x != count_z ) {
			cout << "Error: readCoordinatesAllFiles: Different number of elements read for x,y,z!" << endl;
			exit(1177);
		}
		if( pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() ) {
			cout << "Error: readCoordinatesAllFiles: Different size vectors for x,y,z!" << endl;
			exit(1181);
		}
		
		// allocate the index list (zeroed per job)
		vector<int> jobIndList;
		jobIndList.reserve( allocNum );
		
		// loop over each jobNum
		for( j=0; j < (unsigned int)Config.totNumJobs; j++ )
		{
			jobIndList.clear();
			
			// used per point
			Vector p, v;
			float pc_x, pc_y, pc_z;
			
			// loop over all particles
			for( k=0; k < (unsigned int)count_x; k++ )
			{
				// decide for each if inside sub-camera frustrum+padding
				// the below logic is currently for ORTHOGRAPHIC only!
				p[0] = pos_x[k];
				p[1] = pos_y[k];
				p[2] = pos_z[k];
				
				v = p - camPos; // vector separation between camera and point
					
				// z-coordinate (periodic)
				pc_z = Dot(v, Z);
				
				int fillCheck = k % fillDivisor;
				
				if( (pc_z < camParams[j].zMin || pc_z > camParams[j].zMax) && 
				    (pc_z - All.BoxSize < camParams[j].zMin || pc_z - All.BoxSize > camParams[j].zMax) &&
						(pc_z + All.BoxSize < camParams[j].zMin || pc_z + All.BoxSize > camParams[j].zMax) &&
						(fillCheck != 0) )
					continue;
					
				// y-coordinate (periodic)
				pc_y = Dot(v, Y);
				//aux = pc_z * tang; //tang = (float)tan(ANG2RAD * angle * 0.5) ;
				
				if( (pc_y < camParams[j].yMin || pc_y > camParams[j].yMax) && 
				    (pc_y - All.BoxSize < camParams[j].yMin || pc_y - All.BoxSize > camParams[j].yMax) &&
						(pc_y + All.BoxSize < camParams[j].yMin || pc_y + All.BoxSize > camParams[j].yMax) &&
						(fillCheck != 0) )
					continue;
					
				// x-coordinate (periodic)
				pc_x = Dot(v, X);
				//aux = aux*ratio; // ratio = aspect ratio? probably
				
				if( (pc_x < camParams[j].xMin || pc_x > camParams[j].xMax) && 
				    (pc_x - All.BoxSize < camParams[j].xMin || pc_x - All.BoxSize > camParams[j].xMax) &&
						(pc_x + All.BoxSize < camParams[j].xMin || pc_x + All.BoxSize > camParams[j].xMax) && 
						(fillCheck != 0) )
					continue;
				
				// made it this far, add to keeper list
				jobIndList.push_back(k);
				partCounts[j]++;
			}
			
			// save jobIndList as dataset to mask HDF5 file
			string objName = "Job_" + toStr(j) + "_IndList";
			
			writeGroupDataset( maskFileName, groupNameSave, objName, jobIndList );
			
			partCountsTot[j] += partCounts[j];
			partCountsAllJobs += partCounts[j];
		}

		numReadTot += count_x;
		
		// save partCounts (per job for this file) to dataset in mask HDF5 file
		writeGroupDataset( maskFileName, groupNameSave, "Job_PartCounts", partCounts );
		
		if( Config.verbose )
			cout << "File [" << i << ": " << snapFilenames[i] << "] done, saved [" 
					 << partCountsAllJobs << "] over all jobs." << endl;
	}
	
	// save partCountsTot to hdf5 file
	createNewGroup( maskFileName, "Header" );
	writeGroupDataset( maskFileName, "Header", "PartCountsTot", partCountsTot );
	
	if( numReadTot != numPartTotal[READ_PARTTYPE] ) {
		cout << "Error: Failed to read all coordinates for mask creation!" << endl;
		exit(1178);
	}
	
	if( Config.verbose )
		for( j=0; j < (unsigned int)Config.totNumJobs; j++ )
			cout << " Job [" << j << "] loading [" << partCountsTot[j] << "] in total (" 
					 << setprecision(4) << toStr(float(partCountsTot[j])/float(numReadTot)*100) << "%)" << endl;

}

// ----------- hdf5 helpers ----------------

template<typename T> hid_t ArepoSnapshot::getDataType( void )
{
  if( typeid(T) == typeid(int) )
    return H5T_NATIVE_INT;

  if( typeid(T) == typeid(unsigned) )
    return H5T_NATIVE_UINT;

  if( typeid(T) == typeid(float) )
    return H5T_NATIVE_FLOAT;

  if( typeid(T) == typeid(double) )
    return H5T_NATIVE_DOUBLE;

	if( typeid(T) == typeid(long long) )
		return H5T_NATIVE_LLONG;

	if( typeid(T) == typeid(unsigned long long) )
		return H5T_NATIVE_ULLONG;

	if( typeid(T) == typeid(size_t) )
		return H5T_NATIVE_ULLONG;

  cout << "Error: getDataType: Unknown type";
	exit(1166);
}

void ArepoSnapshot::getDatasetExtent( string fileName,
																			string groupName,
																			string objName,
																			vector<int> &Extent )
{
	// open file and group
  HDF_FileID = H5Fopen( fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
	HDF_GroupID = H5Gopen( HDF_FileID, groupName.c_str() );

  // open dataset
  HDF_DatasetID = H5Dopen( HDF_GroupID, objName.c_str() );

  // dataset did not exist or was empty
  if( HDF_FileID < 0 || HDF_GroupID < 0 || HDF_DatasetID < 0 )
	{
		H5Fclose( HDF_FileID );
    cout << "getDatasetExtent: dataset [" << groupName << "/" << objName 
		     << "] does not exist/is empty" << endl;
		exit(1167);
  }

  // get space associated with dataset
  HDF_DataspaceID = H5Dget_space( HDF_DatasetID );

  int ndims = H5Sget_simple_extent_ndims( HDF_DataspaceID );

  hsize_t *dimsize = new hsize_t[ndims];

  H5Sget_simple_extent_dims( HDF_DataspaceID, dimsize, NULL );

  Extent.clear();
  for(int i=0; i<ndims; ++i )
    Extent.push_back( dimsize[i] );

  delete[] dimsize;

  H5Sclose( HDF_DataspaceID );
  H5Dclose( HDF_DatasetID );
	H5Gclose( HDF_GroupID );
  H5Fclose( HDF_FileID );
}

template<typename T> void ArepoSnapshot::readGroupAttribute( string fileName,
																														 string groupName,
																														 string objName,
																													 	 vector<T> &Data )
{
  HDF_Type = getDataType<T>();
	hsize_t HDF_StorageSize;

  // open file and group
  HDF_FileID = H5Fopen( fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF_GroupID = H5Gopen( HDF_FileID, groupName.c_str() );
	
	// open attribute
  HDF_AttributeID = H5Aopen_name( HDF_GroupID, objName.c_str() );

  if( HDF_FileID < 0 || HDF_GroupID < 0 || HDF_AttributeID < 0 ){
		H5Fclose( HDF_FileID );
    cout << "readGroupAttribute: attr [" << groupName << "/" << objName 
		     << "] does not exist/is empty" << endl;
		exit(1168);
  }
	
	// get size of the file dataspace
	HDF_DataspaceID = H5Aget_space( HDF_AttributeID );
	int ndims = H5Sget_simple_extent_ndims( HDF_DataspaceID );
	
	hsize_t dimsize[ndims];
	H5Sget_simple_extent_dims( HDF_DataspaceID, dimsize, NULL );
	
  HDF_StorageSize = 1;
  for(int i=0; i<ndims; ++i )
    HDF_StorageSize *= dimsize[i];

	// reserve space in Data
	Data.clear();
	Data.reserve( HDF_StorageSize );
	Data.assign( HDF_StorageSize, (T) 0 );
	
	// read
  H5Aread( HDF_AttributeID, HDF_Type, &Data[0] );

  H5Aclose( HDF_AttributeID );
	H5Sclose( HDF_DataspaceID );
  H5Gclose( HDF_GroupID );
  H5Fclose( HDF_FileID );
}

template<typename T> int ArepoSnapshot::readGroupDataset( string fileName,
																					 string groupName,
																					 string objName,
																					 int objVectorIndex,
																					 vector<T> &Data )
{
  HDF_Type = getDataType<T>();
	hsize_t offset[2], stride[2], count[2], block[2], mem_block[1];
	
	if( objVectorIndex < -1 || objVectorIndex > 2 ) {
		cout << "Error: readGroupDataset: Unexpected objVectorIndex [" << objVectorIndex << "]!" << endl;
		exit(1176);
	}

  // open file and group
  HDF_FileID = H5Fopen( fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF_GroupID = H5Gopen( HDF_FileID, groupName.c_str() );
	
	// open dataset
	HDF_DatasetID = H5Dopen( HDF_GroupID, objName.c_str() );
	
  if( HDF_FileID < 0 || HDF_GroupID < 0 || HDF_DatasetID < 0 ){
		H5Fclose( HDF_FileID );
    cout << "readGroupAttribute: attr [" << groupName << "/" << objName 
		     << "] does not exist/is empty" << endl;
		exit(1175);
  }
	
	// get file dataspace size
	HDF_DataspaceID = H5Dget_space( HDF_DatasetID );
	int ndims = H5Sget_simple_extent_ndims( HDF_DataspaceID );
	
  hsize_t dimsize[ndims];
  H5Sget_simple_extent_dims( HDF_DataspaceID, dimsize, NULL );
	
	// if ndims==1, we should have no objVectorIndex
	if( ndims == 1 && objVectorIndex != -1 ) {
		cout << "Error: readGroupDataset: ndims=1 but objVectorIndex = " << objVectorIndex << endl;
		exit(1188);
	}
	if( ndims == 2 && objVectorIndex == -1 ) {
		cout << "Error: readGroupDataset: ndims==2 but objVectorIndex = -1" << endl;
		exit(1189);
	}
	
	// select hyperslab in file
	if( ndims == 2 ) {
		offset[0] = 0;
		offset[1] = objVectorIndex;
		count[0]  = 1;
		count[1]  = 1;
		stride[0] = 1;
		stride[1] = 1;
		block[0]  = dimsize[0];
		block[1]  = 1;
	}
	if( ndims == 1 ) {
		offset[0] = 0;
		count[0]  = 1;
		stride[0] = 1;
		block[0]  = dimsize[0];
	}
	
	H5Sselect_hyperslab( HDF_DataspaceID, H5S_SELECT_SET, offset, stride, count, block );
	
	// make memory space and do read (starting at dataOffset within Data)
	mem_block[0] = dimsize[0];
	
	HDF_MemspaceID = H5Screate_simple( 1, mem_block, NULL );
	
	// reserve space in Data
	Data.reserve( dimsize[0] ); // reallocate only if capacity is insufficient
	Data.resize( dimsize[0] );
	Data.assign( dimsize[0], (T) 0 );	
	
	H5Dread( HDF_DatasetID, HDF_Type, HDF_MemspaceID, HDF_DataspaceID, H5P_DEFAULT, &Data[0] );
	
	H5Sclose( HDF_MemspaceID );
	H5Sclose( HDF_DataspaceID );
	H5Dclose( HDF_DatasetID );
  H5Gclose( HDF_GroupID );
  H5Fclose( HDF_FileID );	
	
	return (int) dimsize[0];
}

template<typename T> int ArepoSnapshot::readGroupDatasetSelect( string fileName,
																																string groupName,
																																string objName,
																																vector<hsize_t> indList,
																																int objVectorIndex,
																																vector<T> &Data )
{
  HDF_Type = getDataType<T>();
	hsize_t mem_block[1];
	unsigned int i;
	
	if( objVectorIndex < -1 || objVectorIndex > 2 ) {
		cout << "Error: readGroupDatasetSelect: Unexpected objVectorIndex [" << objVectorIndex << "]!" << endl;
		exit(1193);
	}

  // open file and group
  HDF_FileID = H5Fopen( fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF_GroupID = H5Gopen( HDF_FileID, groupName.c_str() );
	
	// open dataset
	HDF_DatasetID = H5Dopen( HDF_GroupID, objName.c_str() );
	
  if( HDF_FileID < 0 || HDF_GroupID < 0 || HDF_DatasetID < 0 ){
		H5Fclose( HDF_FileID );
    cout << "readGroupAttribute: attr [" << groupName << "/" << objName << "] does not exist/is empty" << endl;
		exit(1194);
  }
	
	// get file dataspace size
	HDF_DataspaceID = H5Dget_space( HDF_DatasetID );
	int ndims = H5Sget_simple_extent_ndims( HDF_DataspaceID );
	
  hsize_t dimsize[ndims];
  H5Sget_simple_extent_dims( HDF_DataspaceID, dimsize, NULL );
	
	// if ndims==1, we should have no objVectorIndex
	if( ndims == 1 && objVectorIndex != -1 ) {
		cout << "Error: readGroupDatasetSelect: ndims=1 but objVectorIndex = " << objVectorIndex << endl;
		exit(1195);
	}
	if( ndims == 2 && objVectorIndex == -1 ) {
		cout << "Error: readGroupDatasetSelect: ndims==2 but objVectorIndex = -1" << endl;
		exit(1196);
	}
	
	// for 1D, the serialized coordinate list is directly read
	if( ndims == 1 )
	{
		// make selection in filespace
		if( H5Sselect_elements( HDF_DataspaceID, H5S_SELECT_SET, indList.size(), (const hsize_t *)&indList[0] ) < 0 )
			cout << "ERROR: H5Sselect_elements() failed 1D." << endl;
			
		if( H5Sselect_valid( HDF_DataspaceID ) <= 0 )
			cout << "ERROR: H5Sselect_valid() failed 1D." << endl;
	}
	
	// for 2D, make serialized coordinate list
	if( ndims == 2 ) 
	{
		vector<hsize_t> coord;
		coord.reserve( indList.size() * 2 );
		
		for( i=0; i < indList.size(); i++ ) {
			coord.push_back(indList[i]);
			coord.push_back(objVectorIndex);
		}
		
		// make selection in filespace
		if( H5Sselect_elements( HDF_DataspaceID, H5S_SELECT_SET, coord.size()/2, (const hsize_t *)&coord[0] ) < 0 )
			cout << "ERROR: H5Sselect_elements() failed 2D." << endl;
			
		if( H5Sselect_valid( HDF_DataspaceID ) <= 0 )
			cout << "ERROR: H5Sselect_valid() failed 2D." << endl;
	}
	
	// make memory space and do read (starting at dataOffset within Data)
	mem_block[0] = (hsize_t) indList.size();
	
	HDF_MemspaceID = H5Screate_simple( 1, mem_block, NULL );
	
	// reserve space in Data
	Data.reserve( mem_block[0] ); // reallocate only if capacity is insufficient
	Data.resize( mem_block[0] );
	Data.assign( mem_block[0], (T) 0 );	
	
	H5Dread( HDF_DatasetID, HDF_Type, HDF_MemspaceID, HDF_DataspaceID, H5P_DEFAULT, &Data[0] );
	
	H5Sclose( HDF_MemspaceID );
	H5Sclose( HDF_DataspaceID );
	H5Dclose( HDF_DatasetID );
  H5Gclose( HDF_GroupID );
  H5Fclose( HDF_FileID );
	
	if( Data.size() != indList.size() ) {
		cout << "Error: Size mismatch between Data and indList." << endl;
		exit(1197);
	}
	
	return (int) Data.size();
}

template<typename T> int ArepoSnapshot::writeGroupDataset( string fileName,
																													 string groupName,
																													 string objName,
																													 vector<T> &Data )
{
	HDF_Type = getDataType<T>();
	
	// open file and groups (make if necessary)
	HDF_FileID = H5Fopen( fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
	HDF_GroupID = H5Gopen( HDF_FileID, groupName.c_str() );
	
	//if( HDF_GroupID < 0 )
	//	HDF_GroupID = H5Gcreate( HDF_FileID, groupName.c_str(), 0 );
		
	// make spaces
  HDF_Dims = Data.size();
  HDF_DataspaceID = H5Screate_simple(1, &HDF_Dims, NULL);
  HDF_DatasetID = H5Dcreate( HDF_GroupID, objName.c_str(), HDF_Type, HDF_DataspaceID, H5P_DEFAULT );
				
  if( HDF_FileID < 0 || HDF_GroupID < 0 || HDF_DatasetID < 0 ){
		H5Fclose( HDF_FileID );
    cout << "writeGroupDataset: failed to create [" << groupName << "/" << objName 
		     << "] for write" << endl;
		exit(1183);
  }		
				
	// write and close
  herr_t HDF_Err = H5Dwrite( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Data[0] );
						
  H5Dclose( HDF_DatasetID );
  H5Sclose( HDF_DataspaceID );
	H5Gclose( HDF_GroupID );
  H5Fclose( HDF_FileID );

	return (HDF_Err>=0 ? 1 : 0);
} 

void ArepoSnapshot::createNewFile( string fileName )
{
	HDF_FileID = H5Fcreate( fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
	H5Fclose( HDF_FileID );
}

void ArepoSnapshot::createNewGroup( string fileName, string groupName )
{
	HDF_FileID = H5Fopen( fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
	HDF_GroupID = H5Gcreate( HDF_FileID, groupName.c_str(), 0 );
	H5Gclose( HDF_GroupID );
	H5Fclose( HDF_FileID );
}
