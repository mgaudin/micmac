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

	// Loop on every config of TPM of the set of TPM
	for (auto & aCnf : mTPM->VPMul()) {

        // Id of the point
        int i(0);

		// Initialize residual vector
		std::vector<double> aResid;
		// Retrieve 3D position in model geometry with residual
		std::vector<Pt3dr> aPts = aCnf->IntersectBundle(mCams, aResid);

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

	// put the angle value in degrees
	double rad2deg = M_PI / 180.0;
	double angle_deg = angle * rad2deg;

	return angle_deg;
}

double maxInterAngle(int multiplicity, map<int, CamStenope *> mCams, cSetPMul1ConfigTPM *const aCnf, Pt3d<double> &Pt) {
	/***
	* Description: Compute the maximum angle of intersection of two bundles for a 3D point.
	*
	* Parameters:
	*   - multiplicity : the numbers of images where a 3D point is seen in
	*   - mCams : map of the cameras of the scene
	*   - aCnf : a configuration of the tie point set
	*   - Pt : a 3D point
	*
	* Returns: maximum angle of intersection of two bundles for a 3D point
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

	// For debug purposes
	int a(0);
	int b(0);

	for (int indexA = 0; indexA < multiplicity; indexA++) {
		for (int indexB = indexA; indexB < multiplicity; indexB++) {
			if (indexA != indexB) {
				angles.push_back(vectorAngle(vectors[indexA], vectors[indexB]));
				// For debug purposes
				int a = indexA;
				int b = indexB;
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
		std::cout << mCams[aCnf->VIdIm().at(a)]->VraiOpticalCenter() << endl;
		std::cout << mCams[aCnf->VIdIm().at(a)]->VraiOpticalCenter() << endl;

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
	mFactElimTieP=5.0;

    ElInitArgMain
            (
                argc,argv,
                LArgMain()
				<< EAMC(mDir,"Working Directory", eSAM_IsDir)
                << EAMC(mOriPat,"Orientation (xml) list of file, ex 'Ori-Rel/Orientation.*.xml'", eSAM_IsPatFile)
                << EAMC(mFileSH,"File of new set of homol format.", eSAM_IsExistFile )
                ,
                LArgMain()
                << EAM(mOut,"Out",true, "Name of results. (default = pointsWithFeatures.txt)" )
                << EAM(mDebug,"Debug",true, "Print message in terminal to help debugging. (default = false)" )
                << EAM(mPrintTP_info,"PrintTP",true, "Print tie point info in terminal. (default = false)" )
				<< EAM(mWithRadiometry, "WithRadiometry", true, "Save the radiometric information. (default = true)")
				<< EAM(mFactElimTieP, "FactElimTieP", true, "Elim Factor for tie point. (default = 5.0)")
                << EAM(mSavePly, "SavePly", true, "Save the information as ply file. (default = false)")

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
                mIms[ImTPM->Id()]= new cISR_ColorImg(NameIm);
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

		// create ply file
		ofstream aPlyFile;
		if (mSavePly) {
			// set n if ply ply writing is needed
			int n = (mOut.length() - 4);
			// open ply file
			aPlyFile.open(mDir + mOut.substr(0, n) + ".ply", ios::out | ios::app);
		}

        // write file header
        aFile << "Config Point X Y Z R G B mean_reprojection_error multiplicity max_angle\n";

		// if SavePly option is requested
		if (mSavePly) {
			int aNbPts = 0;
			for (auto & aCnf : mTPM->VPMul())
			{
				aNbPts += aCnf->NbPts();
			}
			// write ply header
			aPlyFile << "ply\nformat ascii 1.0\n";
			aPlyFile << "element vertex ";
			aPlyFile << aNbPts;
			aPlyFile << "\n";
			aPlyFile << "property float x\n";
			aPlyFile << "property float y\n";
			aPlyFile << "property float z\n";
			aPlyFile << "property uint config\n";
			aPlyFile << "property uint pt_nb_in_config\n";
			if (mWithRadiometry) {
				aPlyFile << "property uchar red\n";
				aPlyFile << "property uchar green\n";
				aPlyFile << "property uchar blue\n";
			}
			aPlyFile << "property float mean_reprojection_error\n";
			aPlyFile << "property float multiplicity\n";
			aPlyFile << "property float max_angle\n";
			aPlyFile << "end_header\n";
		}

        // loop on every config of TPM of the set of Tie Point Multiple
        int count_Cnf(0); // a counter for the number of tie point configuration in the set of TPM
        for (auto & aCnf : mTPM->VPMul())
        {
            // initialize a vector that will contains all the mean reprojection error for all tie point
            std::vector<double> aResid;
            // do 2 things; compute pseudo intersection of bundle to have 3D position of all tie point of the config
			// and fill the "aResid" vector with mean reprojection error
            std::vector<Pt3dr> aPts=aCnf->IntersectBundle(mCams,aResid);

			/*
			// (dev) Try to understand why residuals are so high
			std::cout << "\n";
			for (auto j = aResid.begin(); j != aResid.end(); j++) {
				std::cout << *j << ' ';
			}
			std::cout << "\n";
			std::cout << "\n";
			*/

            // Iterate on each point to have "X Y Z R G B Residual Reprojection_error Max_inter_angle" information
            int i(0);

            for (auto & Pt: aPts){

				// get multiplicity
				int multiplicity = aCnf->NbIm();

				// get max_inter_angle
				double max_inter_angle = maxInterAngle(multiplicity, mCams, aCnf, Pt);

                // Write the coordinates of point Pt in the output file
				aFile << count_Cnf << " " << i << " " << Pt.x << " " << Pt.y << " " << Pt.z;

				// if SavePly option is requested
				if (mSavePly) {
					// write information on point
					aPlyFile << Pt.x << " " << Pt.y << " " << Pt.z << " " << count_Cnf << " " << i;
				}
				
				if (mWithRadiometry) {
					// initialize average color vectors
					vector<int> averageRed, averageGreen, averageBlue;

					// iterate on the images where the point can be seen
					for (int pos = 0; pos < multiplicity; pos++) {
						// get the RGB information from UV coordinates
						Pt2di image_coords(aCnf->Pt(i, aCnf->VIdIm().at(pos)));
						cISR_Color image_colors = mIms[aCnf->VIdIm().at(pos)]->get(image_coords);

						// r(), g() and b() return u_int1 which are not properly display in terminal, so first a cast in int
						// std::cout << "r = " << (int)image_colors.r() << " ; g = " << (int)image_colors.g() << " ; b = " << (int)image_colors.b() << "\n";

						averageRed.push_back((int)image_colors.r());
						averageGreen.push_back((int)image_colors.g());
						averageBlue.push_back((int)image_colors.b());
					}

					// compute the average value for each channel
					int aRed = std::accumulate(averageRed.begin(), averageRed.end(), 0) / averageRed.size();
					int aGreen = std::accumulate(std::begin(averageGreen), std::end(averageGreen), 0) / averageGreen.size();
					int aBlue = std::accumulate(std::begin(averageBlue), std::end(averageBlue), 0) / averageBlue.size();

					// Write the colorimetric information for point Pt in the output file
					aFile << " " << aRed << " " << aGreen << " " << aBlue;

					// if SavePly option is requested
					if (mSavePly) {
						// write information on point
						aPlyFile << " " << aRed << " " << aGreen << " " << aBlue;
					}
				}

				// Write the features for the point Pt in the output file
				aFile << " " << aResid.at(i) << " " << multiplicity << " " << max_inter_angle << endl;

				// if SavePly option is requested
				if (mSavePly) {
					// write information on point
					aPlyFile << " " << aResid.at(i) << " " << multiplicity << " " << max_inter_angle << endl;
				}

				//std::cout << "UV = " << aCnf->Pt(i, aCnf->VIdIm().at(0)) << " ; Image = " << aCnf->VIdIm().at(0) << " ; Name_image = " << mTPM->NameFromId(aCnf->VIdIm().at(0)) << endl;

                if (mPrintTP_info) {
                    std::cout << "Config " << count_Cnf << "Point " << i << "have XYZ position " << Pt << " and mean reprojection error of " << aResid.at(i) << " and multiplicity of " << aCnf->NbIm() << "\n";

                    std::cout << "Radiometry of this point may be extratcted from pixel position UV " << aCnf->Pt(i, aCnf->VIdIm().at(0)) << " of image " << aCnf->VIdIm().at(0) << " which name is " << mTPM->NameFromId(aCnf->VIdIm().at(0)) << "\n";
                }

                i++;
            }

            if (mDebug) std::cout << "Manipulation of tie point is finished for configuration number " << count_Cnf << "\n\n";

            count_Cnf++;
        }

        // Close the output file
        aFile.close();
		// print sucess
        std::cout << mOut << " has been created sucessfully." << "\n";
		
		// if SavePly option is requested
		if (mSavePly) {
			// set n if ply ply writing is needed
			int n = (mOut.length() - 4);
			// print sucess
			std::cout << mOut.substr(0, n) + ".ply" << " has been created sucessfully." << "\n";
		}
    }
}

int main_manipulateNF_PH(int argc,char ** argv)
{
    cManipulate_NF_TP(argc,argv);
    return EXIT_SUCCESS;
}

