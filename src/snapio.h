/*
 * snapio.h
 * dnelson
 */
 
#ifndef AREPO_RT_SNAPIO_H
#define AREPO_RT_SNAPIO_H

#include "ArepoRT.h"
#include "hdf5.h"

// ArepoSnapshot: handle HDF5 I/O related to snapshots
class ArepoSnapshot
{
public:
  // construction
  ArepoSnapshot(string fB) { fileBase = fB; }

  // methods
  void read_ic();
  template<typename T> void recenterBoxCoords(vector<T> &pos_xyz, int index);

  // HDF5 helpers
  template<typename T> hid_t getDataType(void);   
  bool groupExists(string fileName, string groupName, string objName);
  void getDatasetExtent(string fileName, string groupName, string objName, vector<int> &Extent);
  int getDatasetTypeSize(string fileName, string groupName, string objName);
  bool hasGroupAttribute(string fileName, string groupName, string objName);
  
  template<typename T> void readGroupAttribute(string fileName,
    string groupName, string objName, vector<T> &Data);
  template<typename T> int readGroupDataset(string fileName,
    string groupName, string objName, int objVectorIndex, vector<T> &Data);
  template<typename T> int readGroupDatasetSelect(string fileName,
    string groupName, string objName, vector<hsize_t> indList, int objVectorIndex, vector<T> &Data);                            
  template<typename T> int writeGroupDataset(string fileName,
    string groupName, string objName, vector<T> &Data, int flag2d = 0, int dim0 = 0, int dim1 = 0);  

  void createNewFile(string fileName);
  void createNewGroup(string fileName, string groupName);
  
private:
  // methods  
  unsigned long hash_djb2( string str )
  {
    const char *cstr = str.c_str();
    unsigned long hash = 5381;
    int c;
    
    while ((c = *cstr++))
      hash = ((hash << 5) + hash) + c;
    return hash;
  }
  
  void makeNewMaskFile( string maskFileName, vector<string> snapFilenames );
  void loadAllChunksWithMask( string maskFileName, vector<string> snapFilenames );
  void loadAllChunksNoMask( vector<string> snapFilenames );
  
  struct cameraParams {
    float zMin;
    float zMax;
    float yMin;
    float yMax;
    float xMin;
    float xMax;
  };
  
  void calcSubCamera( struct cameraParams &cP, int jobNum );
  void convertUthermToKelvin( vector<float> &utherm, vector<float> &nelec );
  
  // data
  string fileBase;
  
  hid_t HDF_Type;
  hid_t HDF_FileID;
  hid_t HDF_GroupID;
  hid_t HDF_DatasetID;
  hid_t HDF_AttributeID;
  
  hid_t HDF_DataspaceID;
  hid_t HDF_MemspaceID;
  
  hsize_t HDF_Dims;
  hsize_t HDF_Dims_2D[2];
};

#endif
