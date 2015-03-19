/*
 *
 *
 * Example usage:
 * 	You have a macro for a single SPECThead placement. 
 * 	Replace the radius in the macro with VAR_RADIUS_CM_TO_SET
 * 	Replace the angle of the projection in the macro with VAR_ANGLE_DEGREES_TO_SET
 * 	run
	1)		printRadiusAndAngleToFile(radiusCM, nProjections, TString filename, const Bool_t header=false);
				e.g. printRadiusAndAngleToFile(10.0, 128, "myRadiusAndAngles.txt");
	2)
			generateAllNewMacFiles(TString origMacFilename, TString newMacPrefix, TString inputRadiusAngleFile)
				e.g. generateAllNewMacFiles("orig.mac", "JMnew_", "myRadiusAndAngles.txt");
	3) Run your simulations

	4) Convert root output to projection
		void createProjectionData(TTree* myTree, 
			const double radiusCM, 
			const double angleDegrees,
			const TString filename, 
			const double lowMeV=0.1295,
			const double highMeV=0.1505,
			const double x_startMM=-591./2.,//-crystalWidth/2 in Y direction 591/2 cm
			const double y_startMM=-445./2.,//-crystalLength/2 in Z direction 445/2 cm
			const double x_endMM=591./2., //crystalWidth/2 in Y direction
			const double y_endMM=445./2., //crystalLength/2
			const double pixSizeMM=4.795) 
		Might want to change x_startMM to -537/2, y_startMM to -383/2, x_endMM to 537/2, y_endMM to 383/2.
			e.g. createProjectionData(SinglesSpblurring, 
						10., 
						2.8125, 
						"outputProjection_10_2.8125.txt", 
						0.1295, 
						0.1505, 
						-537./2.,
						-383./2.,
						537./2.,
						383./2.);

	5) Now center and pad the projection data:
		void centerAndPadProjection(const TString projFilename, const Int_t nx=128, const Int_t ny=128)
		e.g. centerAndPadProjection("outputProjection_10_2.8125.txt");

	

			
		

 * 
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include "TROOT.h"
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TMath.h>
#include <TRotation.h>

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

//This opens a *.mac file and replaces the variables VAR_RADIUS_CM_TO_SET and VAR_ANGLE_DEGREE_TO_SET
// with radius and angle respectively.
//  Note: It will overwrite the mac file supplied.
void updateRadiusAndAngleInsideMacFile(TString origMacFilename, TString radiusStr, TString angleStr)
{
	std::ifstream inFile(origMacFilename);
	std::string tempString;
	const std::string radius(radiusStr.Data());
	const std::string angle(angleStr.Data());
	std::vector<std::string> strArr;
	const std::string replaceRadius("VAR_RADIUS_CM_TO_SET");
	const std::string replaceAngle("VAR_ANGLE_DEGREE_TO_SET");

	while(std::getline(inFile, tempString))
	{
	std::cout << tempString << "\n";
		Int_t i=0;
		for(i = tempString.find(replaceRadius, 0); i!=std::string::npos; i=tempString.find(replaceRadius,i))
		{
			tempString.erase(i, replaceRadius.length());
			tempString.insert(i, radius);
			i++;
		}
		for(i = tempString.find(replaceAngle, 0); i!=std::string::npos; i=tempString.find(replaceAngle,i))
		{
			tempString.erase(i,replaceAngle.length());
			tempString.insert(i, angle);
			i++;
		}
		strArr.push_back(tempString);
	}
	inFile.close();

	std::ofstream oFile(origMacFilename);
	for(Int_t i=0; i < strArr.size(); ++i) oFile << strArr[i] << "\n";
	oFile.close();
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
		std::cout << "Read values: " << radius << "  " << angle << "\n";
	}
	inFile.close();

	for(Int_t i=0; i < radiusArr.size(); ++i){
		TString myRadius(radiusArr[i].c_str());
		TString myAngle(angleArr[i].c_str());
		TString newFileName = newMacPrefix + myRadius + "_" + myAngle +".mac";
		std::cout << "copying file # " << i << ":" << newFileName << "\n";
		createNewMacFile(origMacFilename, newFileName);
		updateRadiusAndAngleInsideMacFile(newFileName, myRadius, myAngle);
	}
	
}



//Use this to generate an input file to generate input macro files with radius and angles to replace
// Example usage to dump output of this function to text file. 
//	printRadiusAndAngleToFile(10.0, 128, "myFavorite.txt");
//
void printRadiusAndAngleToFile(double radiusCM, const Int_t nProjections, TString filename, const Bool_t header=false){


	const double  deltaPhi = 360./nProjections;
	double currentPhi =0;
	std::ofstream myFile(filename.Data());

	if(myFile)
	{
		if (header) std::cout << "radius(cm)	angle(degrees)\n";
		for(Int_t i=0; i < nProjections; ++i)
		{
			myFile << TString::Format("%0.1f",radiusCM)+" "+TString::Format("%0.3f", currentPhi)+"\n";
			currentPhi+=deltaPhi;
		}
		myFile.close();
	} else
	{
		std::cout << "unable to open file to write: " << filename << "\n";
	}
}

//just performs a rotation on the coordinates to get them in the coordinates such that the detector plane
// lies in the centered in the YZ plane. 
// This assumes no translation of your detector, only radial extension with rotations around z axis at radius=0.
//
//	outVec: 3 element array representing x,y,z of transformed inVec
//	inVec: 3 element array representing the x,y,z of event in crystal
// 	radiusCM: distance from isocenter to detector face
// 	angleDegrees: the direction of the radius vector from the origin
// 	radiusOffset: distance from the detector face to center of SPECThead (cm)
// 		You can determine this from the Gate simulation aliases. 
// 		This is the = crystalTranslateX+ crystalThick/2 + AluminumCoverThick + CollimatorThick
//
// 	Example usage:
// 	root [21] double outVec[3];
// 	#detector rotated 90 degrees with radius of 10cm (ignoring offset here for example)
// 	root [22] double inVec[3] = {-1.,10.,15.};
// 	root [23] transformPointToLocalCoordinatesHack(myOut,inVec,10., 90.,0,0);
// 	root [24] cout << myOut[0] << " " <<  myOut[1] << " "<< myOut[2] << endl;
// 	-6.12303e-17 1 15
// 	#detector rotated 0 degrees with radius of 10 cm
// 	root [25] double inVec[3] = {10.,0.,0.};
// 	root [26] transformPointToLocalCoordinatesHack(myOut,inVec,10., 0.,0,0);
// 	root [27] cout << myOut[0] << " " <<  myOut[1] << " "<< myOut[2] << endl;
// 	0 0 0
//
//	Alternative implementations...
//		1) perform some direct analysis on the Hits file which has the local coordinates already!
//		2) Modify the Gate source code to also output local coordinates!
//
void transformPointToLocalCoordinatesHack(Float_t* outVec, Float_t* inVec,  
	const double radiusCM, const double angleDegrees, 
	const double SPECTheadTranslate=13.6788, //SPECTheadX/2
	const double crystalTranslate=2.93125, //collimatorThick + aluminumCrystalEncapsulationThick + crystalThick/2
	const double tolerance=1E-6
	)
{
	const double referenceAngle = 0; //corresponds to vector (1,0,0). 90 corresponds to (0,1,0)
	const double fullRadius = radiusCM + SPECTheadTranslate; 
	const double degToRadians = TMath::Pi()/180.;
	const double myAngle = degToRadians*angleDegrees;
	
	//translate so center of SPECThead is at origin
	const double xTrans = inVec[0] - 1.0*fullRadius*cos(myAngle);
	const double yTrans = inVec[1] - 1.0*fullRadius*sin(myAngle);
	const double zTrans = inVec[2] - 0; 
	
	//rotate translated points to YZ plane
	TVector3 myPoint(xTrans, yTrans,zTrans);	
	TRotation m;
	m.RotateZ(-1.0*myAngle);
	myPoint = m*myPoint;
	
	//translate so crystal center is at X=0
	myPoint[0] += (SPECTheadTranslate - crystalTranslate);

	if (abs(myPoint[0]) < tolerance) myPoint[0]=0.;
	if (abs(myPoint[1]) < tolerance) myPoint[1]=0.;
	if (abs(myPoint[2]) < tolerance) myPoint[2]=0.;

	outVec[0] = myPoint[0];
	outVec[1] = myPoint[1];
	outVec[2] = myPoint[2];
}


UInt_t ij_to_index(const Int_t i, const Int_t j, const Int_t nx, const Int_t ny)
{
	if( i >= nx || i < 0 ) std::cout << "warning i=" << i << " is out of range (nx=" << nx << "!\n";
	if( j >= ny || j < 0 ) std::cout << "warning j=" << j << " is out of range (ny=" << ny << "!\n"; 
	return i*ny + j;
}


//nx and ny are the new dimensions. They should be greater than the original projection names.
// Note: this will overwrite the projFilename.
void centerAndPadProjection(const TString projFilename, const Int_t nx=128, const Int_t ny=128)
{
	std::ifstream myFile(projFilename.Data());
	Int_t oldNx, oldNy;
	myFile >> oldNx;
	myFile >> oldNy;
	UInt_t* oldProjArr = new UInt_t[oldNx*oldNy];
	UInt_t temp;
	
	UInt_t totalCounts=0;
	for(Int_t i=0; i < oldNx; ++i)
		for(Int_t j=0; j < oldNy; ++j)
		{
			myFile >> temp;
			oldProjArr[ij_to_index(i,j,oldNx,oldNy)] = temp+1;
			totalCounts+=temp;
		}

	myFile.close();

	std::cout << "total counts read in projection: " << totalCounts << "\n";

	if (nx >= oldNx && ny >= oldNy)
	{
		std::cout << "creating expanded projection data...\n";
		UInt_t* newProjArr = new UInt_t[nx*ny];
		for(UInt_t i=0; i < nx*ny; i++) newProjArr[i]=0;
		Int_t xOffset = round((nx-oldNx)/2.0);
		Int_t yOffset = round((ny-oldNy)/2.0);
		std::cout << "copying old projection data into new projection array\n";
		std::cout << "xOffset, yOffset: " << xOffset << "," << yOffset << "\n";
		for(Int_t i=0; i < oldNx; ++i)
			for(Int_t j=0; j < oldNy; ++j)
			{
				newProjArr[ij_to_index(i+xOffset, j+yOffset, nx,ny)] =
					 oldProjArr[ij_to_index(i,j,oldNx,oldNy)];
			}
		
		//write projection file
		std::ofstream outFile(projFilename.Data());
		outFile << nx << " " << ny << "\n";
		for(Int_t i=0; i< nx; ++i)
		{
			for(Int_t j=0; j< ny; ++j)
			{
				outFile << newProjArr[ij_to_index(i, j, nx, ny)] << " ";
			}
			outFile << "\n";
		}
		outFile.close();
		delete newProjArr;
	} else 
	{
		std::cout << "current projection is larger than propsed expansion!\n";
	}

	delete oldProjArr;
}


//Input a tree that has been energy and spatially blurred by detector...
void createProjectionData(TTree* myTree, 
			const double radiusCM, 
			const double angleDegrees,
			const TString filename, 
			const double lowMeV=0.1295,
			const double highMeV=0.1505,
			const double x_startMM=-537./2.,//-crystalWidth/2 in Y direction 591/2 cm (changed to FOV)
			const double y_startMM=-383./2.,//-crystalLength/2 in Z direction 445/2 cm (changed to FOV)
			const double x_endMM=537./2., //crystalWidth/2 in Y direction (changed to FOV)
			const double y_endMM=383./2., //crystalLength/2 (changed to FOV)
			const double pixSizeMM=4.795) 
{

	Int_t nx,ny, nxMax, nyMax;
	Float_t posx, posy, posz, energy, dx,dy,dz,xmin,xmax,zmin,zmax;

	nxMax = ceil((x_endMM - x_startMM)/pixSizeMM);
	nyMax = ceil((y_endMM - y_startMM)/pixSizeMM);
	std::cout << "nxMax,nyMax: " << nxMax << "," << nyMax << "\n";

	//initialize projection array
	const ULong_t arrSize = nxMax*nyMax;
	UInt_t* projArr = new UInt_t[arrSize];
	for(Long64_t i=0; i< arrSize; ++i) projArr[i] = 0;
	
	Int_t i_index=0, j_index=0;
	Int_t max_i_index=-10000, max_j_index = -10000;
	Int_t min_i_index=10000, min_j_index = 10000;
	
	//now lets loop over the Tree entries
	Long64_t nentries = myTree->GetEntries();
	myTree->SetBranchAddress("globalPosX", &posx);
	myTree->SetBranchAddress("globalPosY", &posy);
	myTree->SetBranchAddress("globalPosZ", &posz);
	myTree->SetBranchAddress("energy", &energy);
	Float_t inVec[3];
	Float_t outVec[3];
	std::cout << "processing " << nentries << " entries in Tree...\n";
	for(Long64_t i=0; i < nentries; ++i)
	{
		myTree->GetEntry(i);
		//check that it falls in our energy window and positions
		if ((energy >= lowMeV) && (energy < highMeV)) 
		{
			//transform so that detector plane is in YZ plane
			inVec[0] = posx;
			inVec[1] = posy;
			inVec[2] = posz;
			transformPointToLocalCoordinatesHack(outVec, inVec, radiusCM, angleDegrees);
			
			//confine counts to within specified region (remember we mapped to y and z plane)
			if ( (outVec[1] >= x_startMM) && (outVec[1] < x_endMM) &&
			     (outVec[2] >= y_startMM) && (outVec[2] < y_endMM))
			{
				//convert the position to pixel position
				i_index = round((outVec[1] - x_startMM)/pixSizeMM);
				j_index = round((outVec[2] - y_startMM)/pixSizeMM);
				//std::cout << i_index << "," << j_index << "\n";
				if (i_index < min_i_index) min_i_index = i_index;
				if (j_index < min_j_index) min_j_index = j_index;
				if (i_index > max_i_index) max_i_index = i_index;
				if (j_index > max_j_index) max_j_index = j_index;
				
				if (i_index < 0) std::cout << "Warning i_index < 0!!!\n";
				if (j_index < 0) std::cout << "Warning j_index < 0!!!\n";
				if(i_index > nxMax) std::cout << "Warning i_index > nxMax!!!\n";
				if(j_index > nyMax) std::cout << "Warning j_index > nyMax!!!\n";
				
				//update counts
				projArr[ij_to_index(i_index, j_index, nxMax,nyMax)]++;
			}
		}
	}

	//std::cout << "min_i_index, max_i_index: " << min_i_index << "," << max_i_index << "\n";
	//std::cout << "min_j_index, max_j_index: " << min_j_index << "," << max_j_index << "\n";
		
	//write projection file
	std::ofstream myFile(filename.Data());
	myFile << nxMax << " " << nyMax << "\n";
	for(Int_t i=0; i< nxMax; ++i)
	{
		for(Int_t j=0; j< nyMax; ++j)
		{
			myFile << projArr[ij_to_index(i, j, nxMax, nyMax)] << " ";
		}
		myFile << "\n";
	}
	myFile.close();
	delete projArr;
}


void generateLocalScript(TString scriptFilename, TString radiusAngleFilename, TString prefix)
{
	
	std::ofstream outFile(scriptFilename.Data());	
	std::ifstream inFile(radiusAngleFilename.Data());
	std::string tempString;
	std::string macFilename;
	while(std::getline(inFile, tempString))
	{
		std::istringstream iss(tempString);
		std::string radius;
		std::string angle;
		iss >> radius;
		iss >> angle;
		macFilename = std::string(prefix.Data()) + radius +"_" + angle +".mac"; 
		outFile << "Gate " << macFilename << " &" << "\n";
	}
	inFile.close();
	outFile.close();
}

//This opens a root file and then calls 
// output file name generated from prefix, radius, angle, and energy window to 4 decimal places.
void generateProjectionFromRootFile(TString inFilename, 
				    TString prefix, Double_t radiusCM, Double_t angleDeg,
				    Double_t lowMeV=0.1295, Double_t highMeV=0.1505)
{

	TFile* f = (TFile*)gROOT->GetListOfFiles()->FindObject(inFilename);
	if (!f) 
	{
		f = new TFile(inFilename);

		TTree* myTree = (TTree*)gDirectory->Get("SinglesSpblurring");
		
		TString radius_str = TString::Format("%0.1f", radiusCM);
		TString angle_str = TString::Format("%0.3f", angleDeg);
		TString lowE_str = TString::Format("%0.4f", lowMeV);
		TString highE_str = TString::Format("%0.4f", highMeV);
		
		TString outFilename = prefix	+ radius_str +"_"
					+ angle_str +"_"
					+ lowE_str + "_" 
					+ highE_str + ".proj";
		
		createProjectionData(myTree, radiusCM, angleDeg, outFilename, lowMeV, highMeV);
		centerAndPadProjection(outFilename);
		delete f;
	} else
	{
		std::cout << "unable to open " << inFilename << "\n";
	}
}

//assumes you have generated root files for each projection
void generateAllProjections(TString prefixRoot, TString inputRadiusAngleFile, 
				TString outPrefix="prim_", Double_t lowMeV=0.1295, Double_t highMeV=0.1505)
{
	//open the inputRadiusAngleFile that contains pairs of radius and angle
	std::ifstream inFile(inputRadiusAngleFile.Data());
	std::string tempString;
	std::vector<Double_t> radiusArr;
	std::vector<Double_t> angleArr;
	while(std::getline(inFile, tempString))
	{
		std::istringstream iss(tempString);
		Double_t radius;
		Double_t angle;
		iss >> radius;
		iss >> angle;
		radiusArr.push_back(radius);
		angleArr.push_back(angle);
		std::cout << "Read values: " << radius << "  " << angle << "\n";
	}
	inFile.close();
	
	//now create projections
	for(Int_t i=0; i < radiusArr.size(); ++i)
	{
		TString rootFilename = prefixRoot 
					+ TString::Format("%0.1f", radiusArr[i]) + "_"
					+ TString::Format("%0.3f", angleArr[i]) + ".root";
		std::cout << "\ngenerating projection for " << rootFilename << "\n";				
		generateProjectionFromRootFile(rootFilename, outPrefix, radiusArr[i], angleArr[i], lowMeV, highMeV);
	}
}

