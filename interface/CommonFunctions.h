#ifndef __CommonFunctions__
#define __CommonFunctions__

#include <vector>
#include <stdint.h>
#include <cstddef>

std::vector<int> CrossTalkInv(const std::vector<int>&  Q, const float x1=0.10, const float x2=0.04, bool way=true,float threshold=20,float thresholdSat=25);
bool clusterCleaning(std::vector<int> ampls,  int crosstalkInv=0, uint8_t* exitCode=NULL);


//check function isInside 
//https://github.com/CMS-HSCP/SUSYBSMAnalysis-HSCP/blob/79fa0b3f50969fd92f024fac0ad75a36dbdf389f/test/AnalysisCode/Analysis_CommonFunction.h#L820


#endif
