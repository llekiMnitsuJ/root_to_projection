#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include "TROOT.h"
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

//Use this to load a root file that has the following format:
//   "{prefix}_{radiusCM}_{angleDegrees}.root"
//   e.g.
//   	prefix="proj_";
//   	radiusCM=10;
//   	angleDegrees=90;
//   	"proj_10.0_90.000.root"
//
TFile* loadProjectionFile(const TString prefix, const double radiusCM, const double angleDegrees)
{


	TString filename = prefix +TString::Format("%0.1f", radiusCM)+"_"+TString::Format("%0.3f", angleDegrees)+".root";
	std::cout << "trying to open projection file: " << filename << "\n";
	TFile *f =  (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
	if (!f) {
		f = new TFile(filename);
		std::cout << "opened " << filename << "\n";
	} else {
		std::cout << "unable to open " << filename << "\n";
		f = NULL;
	}

	return f;
}

//this just copies one file to another
void createNewMacFile(TString macFileName, TString newMacFileName, Bool_t overwrite=false)
{
	Int_t result = gSystem->CopyFile(macFileName, newMacFileName, overwrite);
	if(result != 0){
		std::cout << "error when copying file: " << macFileName << " to " << newMacFileName << "\n";
		std::cout << "error code: " << result << "\n";
	}
}

//This will create a lot of mac files!
// It requires that you have generated a list of radii and angles (see printRadiusAndAngle)
//  Example Usage:
//   macro to copy: "myMain.mac"
//   newMacprefix: "myMain_"
//   inputRadiusAngleFile: "myFile.txt"
//   generateAllNewMacFiles("myMain.mac", "myMain_", "myFile.txt");
//
void generateAllNewMacFiles(TString origMacFilename, TString newMacPrefix, TString inputRadiusAngleFile){
	
	//open the inputRadiusAngleFile that contains pairs of radius and angle
	std::ifstream inFile(inputRadiusAngleFile.Data());
	std::string tempString;
	std::vector<std::string> radiusArr;
	std::vector<std::string> angleArr;
	while(std::getline(inFile, tempString))
	{
		std::istringstream iss(tempString);
		std::string radius;
		std::string angle;
		iss >> radius;
		iss >> angle;
		radiusArr.push_back(radius);
		angleArr.push_back(angle);
		cout << "Read values: " << radius << "  " << angle << "\n";
	}
	inFile.close();

	for(Int_t i=0; i < radiusArr.size(); ++i){
		TString myRadius(radiusArr[i].c_str());
		TString myAngle(angleArr[i].c_str());
		TString newFileName = newMacPrefix + myRadius + "_" + myAngle +".mac111";
		std::cout << "copying file # " << i << ":" << newFileName << "\n";
		createNewMacFile(origMacFilename, newFileName);
	}
	
}

void updateRadiusAndAngleInsideMacFile(TString origMacFilename, TString radius, TString angle)
{
	std::ifstream inFile(origMacFilename);
	std::string tempString;
	std::vector<std::string> strArr;
	const std::string replaceRadius("VAR_RADIUS_CM_TO_SET");
	const std::string replaceAngle("VAR_ANGLE_DEGREE_TO_SET");
	while(std::getline(inFile, tempString))
	{
		Int_t i=0;
	}


}


//Use this to generate an input file to generate input macro files with radius and angles to replace
// Example usage to dump output of this function to text file. 
//  	root[0]: printRadiusAndAngle(10, 128); > myFile.txt
//
void printRadiusAndAngle(double radiusCM, const Int_t nProjections, const Bool_t header=false){


	const double  deltaPhi = 360./nProjections;
	double currentPhi =0;
	if (header) std::cout << "radius(cm)	angle(degrees)\n";
	for(Int_t i=0; i < nProjections; ++i){
		std::cout << TString::Format("%0.1f",radiusCM)+" "+TString::Format("%0.3f", currentPhi)+"\n";
		currentPhi+=deltaPhi;
	}
}

//just performs a rotation on the coordinates to get them in the coordinates such that the detector plane
// lies in the 
void transformPointToLocalCoordinatesHack(const double angleDegrees)
{
}



void createProjectionData(TTree* myTree, 
			const double angleDegrees, 
			const double pixSizeMM=4.795, 
			const Int_t nx=128, 
			const Int_t ny=128)
{

}




