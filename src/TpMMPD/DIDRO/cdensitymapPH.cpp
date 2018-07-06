#include "cdensitymapPH.h"

// create a map of the entire area with very low resolution that depict the density of tie point (point homologue PH)
// or also Multiplicity of tie point map or Residual map

Im2D_REAL8 aDenoise(3,3,
                    "0 1 0 "
                    "1 2 1 "
                    " 0 1 0"
                    );

cDensityMapPH::cDensityMapPH(int argc,char ** argv)
{

    mOut="TiePoints_DensityMap.tif";
    mDebug=0;
    mExpTxt=0;
    mDir="./";
    mSmoothing=1;
    mMultiplicity=0;
    mResid=0;
    mThreshResid=5;

    ElInitArgMain
            (
                argc,argv,
                LArgMain()   << EAMC(mDir,"Working Directory", eSAM_IsDir)
                << EAMC(mOriPat,"Orientation (xml) list of file, ex 'Ori-Rel/Orientation.*.xml'", eSAM_IsPatFile)
                ,
                LArgMain()  << EAM(mSH,"SH",true, "Set of Homol name")
                << EAM(mExpTxt,"ExpTxt",true, "Are tie points in txt format? default false, means standard binary format is used." )
                << EAM(mFileSH,"FileSH",true, "File of new set of homol format. If provided, argument as 'ExpTxt' and 'SH are not appropriated ",eSAM_IsExistFile )
                << EAM(mOut,"Out",true, "Name of resulting density map, default TiePoints_DensityMap( + SH).tif" )
                << EAM(mGSD,"GSD",true, "Ground Sample Distance of the resulting density map" )
                << EAM(mWidth,"Width",true, "Size [pix] of width resulting density map" )
                << EAM(mSmoothing,"Smooth",true, "Apply Gaussian filter to smooth the result, def true" )
                << EAM(mMultiplicity,"Multi",true, "if true, density map depict not number of tie point but value of max multiplicity. Def false." )
                << EAM(mResid,"Resid",true, "if true, density map depict not number of tie point but value of max reprojection error. Def false." )
                << EAM(mThreshResid,"MaxResid",true, "Tie point are probably not filtered for outliers. Threshold to reject them in case option 'Resid' is true. default 5pixels of reprojection error." )
                << EAM(mDebug,"Debug",true, "Print message in terminal to help debugging." )

                );

    if (!MMVisualMode)
    {
        if (!EAMIsInit(&mOut) && mMultiplicity) mOut="TiePoints_MultiplicityMap.tif";
        if (!EAMIsInit(&mOut) && mResid) mOut="TiePoints_MaxResidualMap.tif";

        mICNM = cInterfChantierNameManipulateur::BasicAlloc(mDir);

        // load orientations
        mOriFL = mICNM->StdGetListOfFile(mOriPat);

        loadPH(); // load PH, convert old to new format is required

        std::string aKey= "NKS-Assoc-Ori2ImGen"  ;
        std::string aTmp1, aNameOri;

        for (auto &aOri : mOriFL){

            // retrieve name of image from name of orientation xml
            SplitDirAndFile(aTmp1,aNameOri,aOri);
            std::string NameIm = mICNM->Assoc1To1(aKey,aNameOri,true);
            mImName.push_back(NameIm);
            // retrieve IdIm
            cCelImTPM * ImTPM=mTPM->CelFromName(NameIm);
            if (ImTPM) {
                // map container of Camera is indexed by the Id of Image (cSetTiePMul)
                mCams[ImTPM->Id()]=CamOrientGenFromFile(aOri,mICNM);
            } else {
                std::cout << "No tie points found for image " << NameIm << ".\n";
            }
        }

        determineMosaicFootprint();
        determineGSD();
        mDM=Im2D_REAL4(mSz.x,mSz.y,0.0);

        populateDensityMap();
        std::string aOut(mDir+mOut);
        if (mSmoothing) ELISE_COPY(mDM.all_pts(),som_masq(mDM.in_proj(),aDenoise)/6,mDM.out());
        Tiff_Im::CreateFromIm(mDM,aOut);
        writeTFW(aOut,Pt2dr(mGSD,-mGSD),Pt2dr(mBoxTerrain.P0().x,mBoxTerrain.P1().y));


        if (mDebug) std::cout << "Result saved in " << aOut << ".\n";

    }
}

void cDensityMapPH::loadPH(){

    if (EAMIsInit(&mFileSH)){
        mTPM = new cSetTiePMul(0);
        //mTPM->SetFilter(mImName);
        mTPM->AddFile(mFileSH);
    }

    // old tie p format; load and convert??
    else {
        std::cout << "Warn, right now only new format of multi tie p is supporte, use FileSH argument.\n";

    }
}


void cDensityMapPH::determineMosaicFootprint(){

    if (mDebug) std::cout << "Determine Density Map footprint.\n";
    double xmin(3.4E+38),ymin(3.4E+38),xmax(0),ymax(0);
    bool first=true;
    for (auto & Cam : mCams)
    {
        double alt=Cam.second->GetAltiSol();
        Box2dr box=Cam.second->BoxTer(alt);
        if (first){
            xmin=box.P0().x;
            ymin=box.P0().y;
            xmax=box.P1().x;
            ymax=box.P1().y;
            first=false;
        } else {
            xmin=ElMin(xmin,box.P0().x);
            ymin=ElMin(ymin,box.P0().y);
            xmax=ElMax(xmax,box.P1().x);
            ymax=ElMax(ymax,box.P1().y);
        }
    }
    mBoxTerrain=Box2dr(Pt2dr(xmin,ymin),Pt2dr(xmax,ymax));
    if (mDebug) std::cout << "Ground box of density map is " << ToString(mBoxTerrain) << ", with means a swath of " << Pt2dr(mBoxTerrain.P1()-mBoxTerrain.P0()) <<" \n";
}

void cDensityMapPH::determineGSD(){

    if (!EAMIsInit(&mGSD) && !EAMIsInit(&mSz)) {
        double aGSDmean(0);
        for (auto& Cam : mCams){
            aGSDmean=aGSDmean+ Cam.second->GlobResol();
        }
        aGSDmean=aGSDmean/mCams.size();
        if (mDebug) std::cout << "mean GSD of images is " << aGSDmean << ".\n";
        // let's give 625 pixels per nadir image in ground geometry
        mGSD=aGSDmean*(((double)mCams[0]->Sz().x)/25.00);

        if (mDebug) std::cout << "Image size is" << ToString(mCams[0]->Sz()) << ".\n";
        if (mDebug) std::cout << "GSD computed for Density map is equal to " << mGSD << ".\n";
    } else if (EAMIsInit(&mWidth)) {
        mGSD= (mBoxTerrain.P1().x-mBoxTerrain.P0().x)/mWidth;
    }
    mSz= Pt2di(round_up((mBoxTerrain.P1().x-mBoxTerrain.P0().x)/mGSD),round_up((mBoxTerrain.P1().y-mBoxTerrain.P0().y)/mGSD));
    if (mDebug) std::cout << "Image size of density map is " << ToString(mSz) << ", which means GSD of " << mGSD <<"\n";
}

void cDensityMapPH::populateDensityMap() {

	// loop on every config of TPM of the set of TPM
	int progBar = ElMax(int(mTPM->VPMul().size() / 10), 1);
	int cnt(0);
	for (auto & aCnf : mTPM->VPMul())
	{
		// retrieve 3D position in model geometry
		if (!mResid) {
			std::vector<Pt3dr> aPts = aCnf->IntersectBundle(mCams);
			// add the points to the density map
			for (auto & Pt : aPts) {
				Pt2dr PosXY(Pt.x, Pt.y);
				Pt2di PosUV = XY2UV(PosXY);
				double aVal = mDM.GetR(PosUV);
				if (!mMultiplicity) { mDM.SetR(PosUV, aVal + 1); }
				else { mDM.SetR(PosUV, ElMax((int)aVal, aCnf->NbIm())); }
			}
		}
		if (mResid) {
			std::vector<double> aResid;
			std::vector<Pt3dr> aPts = aCnf->IntersectBundle(mCams, aResid);
			// add the points to the density map
			int i(0);
			for (auto & Pt : aPts) {
				Pt2dr PosXY(Pt.x, Pt.y);
				Pt2di PosUV = XY2UV(PosXY);
				double aVal = mDM.GetR(PosUV);
				if (aResid.at(i)<mThreshResid) mDM.SetR(PosUV, ElMax(aVal, aResid.at(i)));
				if (mDebug) std::cout << "Reprojection error for this point is equal to " << aResid.at(i) << "\n";
				i++;
			}
		}

		cnt++;
		if (cnt>progBar) {
			std::cout << "-";
			cnt = 0;
		}
	}
	std::cout << "\n";
}

void cDensityMapPH::populateDensityMap4Tests() {
	// Header
	string header = "x y z density multiplicity reprojection_error max_angle// missing RGB";
	std::cout << header << "\n";

	// PAUSE
	//system("PAUSE");

	// Loop on every config of TPM of the set of TPM
	for (auto & aCnf : mTPM->VPMul()) {

		// Initialize residual vector
		std::vector<double> aResid;
		// Retrieve 3D position in model geometry with residual
		std::vector<Pt3dr> aPts = aCnf->IntersectBundle(mCams, aResid);

		// Id of the point
		int i(0);

		// Add the points to the density map
		for (auto & Pt : aPts) {
			Pt2dr PosXY(Pt.x, Pt.y);
			Pt2di PosUV = XY2UV(PosXY);

			double aVal = mDM.GetR(PosUV);

			// Set the value of density in the density map
			mDM.SetR(PosUV, aVal + 1);
		}

		// Try to do things
		for (auto & Pt : aPts) {
			Pt2dr PosXY(Pt.x, Pt.y);
			Pt2di PosUV = XY2UV(PosXY);

			// DEBUG
			//std::cout << "Camera: " << mCams << "\n";
			//std::cout << "Coordinates Pt3dr: " << "x = " << Pt.x << "; y = " << Pt.y << "; z = " << Pt.z << "\n";
			//std::cout << "Coordinates Pt2dr: " << "Pt.x = " << Pt.x << "; Pt.y = " << Pt.y << "\n";
			//std::cout << "Coordinates Pt2di: " << XY2UV(PosXY) << "\n";

			// Density from density map
			double density = mDM.GetR(PosUV);

			// Multiplicity ?
			int multiplicity = aCnf->NbIm();

			// Reprojection error
			double reprojection_error = aResid.at(i);

			// Print features
			//std::cout << "Density = " << aVal << "\n";
			//std::cout << "Multiplicity = " << multiplicity << "\n";
			//std::cout << "Reprojection error " << aResid.at(i) << "\n";

			string line = std::to_string(Pt.x) + " " + std::to_string(Pt.y) + " " + std::to_string(Pt.z) + " " + std::to_string(density) + " " + std::to_string(multiplicity) + " " + std::to_string(reprojection_error) + "\n"; // missing RGB and angles

			std::cout << line;

			i++;
		}
	}
	std::cout << "\n";

	// PAUSE
	// system("PAUSE");
}

Pt2di cDensityMapPH::XY2UV(Pt2dr aVal){
    if (mBoxTerrain.contains(aVal))
    {
        Pt2di aRes((aVal.x-mBoxTerrain.P0().x)/mGSD,(-aVal.y+mBoxTerrain.P1().y)/mGSD);
        return aRes;
    }   else {return Pt2di(0,0);}
}

int main_densityMapPH(int argc,char ** argv)
{
    cDensityMapPH(argc,argv);
    return EXIT_SUCCESS;
}

/***
Les méthodes devraient permettre de pouvoir calculer des angles à partir d'un point 3D et des positions
des sommets de prises de vues.
***/

double sum(vector<double> vect) {
	/***
	* Description: Compute the sum of the elements of a vector
	*
	* Parameters:
	*   - vect : vector
	*
	* Returns: Sum of the elements of the input vector
	*
	***/
	double sum_of_elems;

	std::for_each(vect.begin(), vect.end(), [&](double n) {
		sum_of_elems += n;
	});

	return sum_of_elems;
}

vector<double> unitize(vector<double> vect)
/***
* Description: Compute a unitarized vector from an input vector
*
* Parameters:
*   - vect : vector
*
* Returns: Unitarized input vector
*
***/
{
	double sum_of_elem = sum(vect);

	if (sum_of_elem > 0) {
		for (std::vector<double>::iterator it = vect.begin(); it != vect.end(); ++it) {
			*it = *it / sum_of_elem;
		}
	}

	return vect;
}

double dotProduct(vector<double> u, vector<double> v)
/***
* Description: Compute the scalar product between two vector
*
* Parameters:
*   - u : vector
*   - v : vector
*
* Returns: Scalar product between u and v
*
***/
{
	double dot = 0;
	for (int i = 0; i < v.size(); i++) {
		dot += u[i] * v[i];
	}

	return dot;
}

double length(vector<double> v)
/***
* Description: Calculates the norm of a vector
*
* Parameters:
*   - v : vector
*
* Returns: The norm of the input vector
*
***/
{
	double norm = sqrt(dotProduct(v, v));

	return norm;
}

double angle(vector<double> u, vector<double> v)
/***
* Description: Calculates the angle between two 3D vectors.
*
* Parameters:
*   - u : vector
*   - v : vector
*
* Returns: The computation needed for angle calculation before acos
*          The angle returned is in radians.
*
***/
{
	double elem = dotProduct(u, v) / (length(u) * length(u));

	return elem;
}

double vectorAngle(vector<double> v1, vector<double> v2)
/***
* Description: Calculates the angle between two 3D vectors.
*
* Parameters:
*   - v0 : vector
*   - v1 : vector
*
* Returns: The angle in radians.
*
***/
{
	// Unitize the input vectors
	v1 = unitize(v1);
	v2 = unitize(v2);

	// (v1.v2)/(|v1|.|v2|)
	double elem = angle(v1, v2);

	// Force the dot product of the two input vectors to
	// fall within the domain for inverse cosine, which
	// is -1 <= x <= 1. This will prevent runtime
	// "domain error" math exceptions.
	elem = (elem < -1.0 ? -1.0 : (elem > 1.0 ? 1.0 : elem));

	double angle = acos(elem);

	return angle;
}

double maxInterAngle(int multiplicity, map<int, CamStenope *> mCams, cSetPMul1ConfigTPM *const aCnf, Pt3d<double> &Pt) {
	/***
	* Description:
	*
	* Parameters:
	*   - multiplicity :
	*   - mCams :
	*   - aCnf :
	*   - Pt :
	*
	* Returns:
	*
	***/
	// Initialize vector list containing the directions from a 3D point to a camera camera center
	vector<vector<double>> vectors;

	for (int i = 0; i < multiplicity; i++) {
		// get the optical center of the camera where the point is seen
		Pt3dr camCenter = mCams[aCnf->VIdIm().at(i)]->VraiOpticalCenter();

		double vx = Pt.x - camCenter.x;
		double vy = Pt.y - camCenter.y;
		double vz = Pt.z - camCenter.z;

		vector<double> v;

		v.push_back(vx);
		v.push_back(vy);
		v.push_back(vz);

		vectors.push_back(v);
	}

	vector<double> angles;

	int compt(0);

	for (int indexA = 0; indexA < multiplicity; indexA++) {
		for (int indexB = indexA; indexB < multiplicity; indexB++) {
			if (indexA != indexB) {
				angles.push_back(vectorAngle(vectors[indexA], vectors[indexB]));
			}
		}
	}

	vector<double>::iterator max = std::max_element(begin(angles), end(angles));

	if (*max == 0) {
		std::cout << "INIT DEBUG" << endl;
		std::cout << "VECTORS" << endl;
		for (int j = 0; j < multiplicity; j++) {
			std::cout << vectors[j][0] << " " << vectors[j][1] << " " << vectors[j][2] << endl;
		}
		std::cout << "ANGLES" << endl;
		for (int j = 0; j < multiplicity; j++) {
			std::cout << angles[j] << endl;
		}
		std::cout << "CAMERA POSITION" << endl;
		std::cout << mCams[aCnf->VIdIm().at(compt)]->VraiOpticalCenter() << endl;

		std::cout << "END DEBUG" << endl;

		std::cout << endl;
	}

	return *max;
}

cManipulate_NF_TP::cManipulate_NF_TP(int argc,char ** argv)
{

    mOut="pointsWithFeatures.txt";
    mDebug=0;
    mDir="./";
    mPrintTP_info=0;
    mSavePly=0;
	mWithRadiometry=1;

    ElInitArgMain
            (
                argc,argv,
                LArgMain()   << EAMC(mDir,"Working Directory", eSAM_IsDir)
                << EAMC(mOriPat,"Orientation (xml) list of file, ex 'Ori-Rel/Orientation.*.xml'", eSAM_IsPatFile)
                << EAMC(mFileSH,"File of new set of homol format.", eSAM_IsExistFile )
                ,
                LArgMain()


                << EAM(mOut,"Out",true, "Name of results" )
                << EAM(mDebug,"Debug",true, "Print message in terminal to help debugging." )
                << EAM(mPrintTP_info,"PrintTP",true, "Print tie point info in terminal." )
				<< EAM(mWithRadiometry, "WithRadiometry", true, "Save the radiometric information")
                << EAM(mSavePly, "SavePly", true, "Save the information as ply file.")

                );

    if (!MMVisualMode)
    {
        // this object is used for regular expression manipulation and key of association utilisation
        mICNM = cInterfChantierNameManipulateur::BasicAlloc(mDir);
        // load list of orientation files from regular expression
        mOriFL = mICNM->StdGetListOfFile(mOriPat);
        // load the set of Tie Point
        mTPM = new cSetTiePMul(0);
        mTPM->AddFile(mFileSH);

		// DEBUG
		// Save the tie points to txt format
        mTPM->Save("PMul.txt");

        // Now that we have 1) list of orientation file and 2) Tie point with New format, we recover the ID
        // of each images (used to manipulate tie points) and the name of each image)
        // the map container mCams will contains the index of the image and the camera stenopee object which
        // is the one used for photogrammetric computation like "intersect bundle", which compute 3D position
        // of a tie point from UV coordinates

        // an association key that give the name of an image from the name of an orientation.
        std::string aKey= "NKS-Assoc-Ori2ImGen"  ;
        std::string aTmp1, aNameOri;

        std::cout << "List of orientation files:\n";

        for (auto &aOri : mOriFL){

            std::cout << aOri << "\n";

            // retrieve name of image from name of orientation xml
            SplitDirAndFile(aTmp1,aNameOri,aOri);
            std::string NameIm = mICNM->Assoc1To1(aKey,aNameOri,true);
            mImName.push_back(NameIm);

            // retrieve IdIm
            cCelImTPM * ImTPM=mTPM->CelFromName(NameIm);
            if (ImTPM) {
                // map container of Camera is indexed by the Id of Image (cSetTiePMul)
                mCams[ImTPM->Id()]=CamOrientGenFromFile(aOri,mICNM);
                // map container of image RGB indexed by the Id of image
                mIms[ImTPM->Id()]=new cISR_ColorImg(NameIm);
                std::string tmp("Tmp-MM-Dir/" + NameIm +"_col.tif");
                mIms[ImTPM->Id()]->write(tmp);
            } else {
                std::cout << "No tie points found for image " << NameIm << ".\n";
            }
        }

        // now we are ready to manipulate the set of tie points

        // create a txt output file containing the features on the points
        ofstream aFile;
        aFile.open(mDir + mOut);

        // write file header
        aFile << "Config Point X Y Z R G B mean_reprojection_error multiplicity max_angle\n";

		/*
		// if SavePly option is requested
		if (mSavePly) {
		string *aContentOfPlyFile = "";
		int *countPoints = 0;
		}
		*/

        // loop on every config of TPM of the set of Tie Point Multiple
        int count_Cnf(0); // a counter for the number of tie point configuration in the set of TPM
        for (auto & aCnf : mTPM->VPMul())
        {

            // initialize a vector that will contains all the mean reprojection error for all tie point
            std::vector<double> aResid;
            // do 2 things; compute pseudo intersection of bundle for have 3D position of all tie point of the config and fill the "aResid" vector with mean reprojection error
            std::vector<Pt3dr> aPts=aCnf->IntersectBundle(mCams,aResid);

            // Iterate on each point to have "X Y Z R G B Residual Reprojection_error Max_inter_angle" information
            int i(0);

            for (auto & Pt: aPts){

				// get multiplicity
				int multiplicity = aCnf->NbIm();

				// get max_inter_angle
				double max_inter_angle = maxInterAngle(multiplicity, mCams, aCnf, Pt);

                // r(), g() and b() return u_int1 which are not properly display in terminal, so first a cast in int
                // std::cout << "r = " << (int)image_colors.r() << " ; g = " << (int)image_colors.g() << " ; b = " << (int)image_colors.b() << "\n";

                // Write the coordinates of point Pt in the output file
				aFile << count_Cnf << " " << i << " " << Pt.x << " " << Pt.y << " " << Pt.z;
				
				if (mWithRadiometry) {
					// get the RGB information from UV coordinates
					Pt2di image_coords(aCnf->Pt(i, aCnf->VIdIm().at(0)));
					cISR_Color image_colors = mIms[aCnf->VIdIm().at(0)]->get(image_coords);

					// Write the colorimetric information for point Pt in the output file
					aFile << " " << (int)image_colors.r() << " " << (int)image_colors.g() << " " << (int)image_colors.b();
				}

				// Write the features for the point Pt in the output file
				aFile << " " << aResid.at(i) << " " << multiplicity << " " << max_inter_angle << endl;

				//std::cout << "UV = " << aCnf->Pt(i, aCnf->VIdIm().at(0)) << " ; Image = " << aCnf->VIdIm().at(0) << " ; Name_image = " << mTPM->NameFromId(aCnf->VIdIm().at(0)) << endl;

                if (mPrintTP_info) {
                    std::cout << "Config " << count_Cnf << "Point " << i << "have XYZ position " << Pt << " and mean reprojection error of " << aResid.at(i) << " and multiplicity of " << aCnf->NbIm() << "\n";

                    std::cout << "Radiometry of this point may be extratcted from pixel position UV " << aCnf->Pt(i, aCnf->VIdIm().at(0)) << " of image " << aCnf->VIdIm().at(0) << " which name is " << mTPM->NameFromId(aCnf->VIdIm().at(0)) << "\n";
                }

				/*
				if (mSavePly) {
				// write in the ply file
				aContentOfPlyFile = aContentOfPlyFile + Pt.x + " " + Pt.y + " " + Pt.z;
				if (mWithRadiometry) {
					aContentOfPlyFile = aContentOfPlyFile + (int)image_colors.r() + " " + (int)image_colors.g() + " " + (int)image_colors.b();
				}
				aContentOfPlyFile = aContentOfPlyFile + " " + aResid.at(i) + " " + aCnf->NbIm() << "\n";
				countPoints++;
				}
				*/

                i++;
            }

            if (mDebug) std::cout << "Manipulation of tie point is finished for configuration number " << count_Cnf << "\n\n";

            count_Cnf++;
        }

        // Close the output file
        aFile.close();
        std::cout << mOut << " has been created sucessfully." << "\n";

		/*
		// if SavePly option is requested
		if (mSavePly) {
		// Print success
		int pos = 0;
		int n = (mOut.length() - 4);
		// if SavePly option is requested
		string aContentOfPlyFile = "";
		// Create ply file
		ofstream aPlyFile;
		int pos = 0;
		int n = (mOut.length() - 4);
		aPlyFile.open(mDir + mOut.substr(pos, n) + ".ply");
		// write ply header
		aPlyFile << "ply\nformat ascii 1.0\n";
		aPlyFile << "element vertex %i\n"; // /!\ MISSING THE NUMBER OF POINTS
		aPlyFile << "property float x\n";
		aPlyFile << "property float y\n";
		aPlyFile << "property float z\n";
		aPlyFile << "property uchar red\n";
		aPlyFile << "property uchar green\n";
		aPlyFile << "property uchar blue\n";
		aPlyFile << "property float mean_reprojection_error\n";
		aPlyFile << "property float multiplicity\n";
		aPlyFile << "end_header\n";
		aPlyFile.close();
		std::cout << mOut.substr(pos, n) + ".ply" << " has been created sucessfully." << "\n";
		}
		*/
    }
}

int main_manipulateNF_PH(int argc,char ** argv)
{
    cManipulate_NF_TP(argc,argv);
    return EXIT_SUCCESS;
}

