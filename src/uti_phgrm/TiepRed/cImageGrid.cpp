/*Header-MicMac-eLiSe-25/06/2007

    MicMac : Multi Image Correspondances par Methodes Automatiques de Correlation
    eLiSe  : ELements of an Image Software Environnement

    www.micmac.ign.fr


    Copyright : Institut Geographique National
    Author : Marc Pierrot Deseilligny
    Contributors : Gregoire Maillet, Didier Boldo.

[1] M. Pierrot-Deseilligny, N. Paparoditis.
    "A multiresolution and optimization-based image matching approach:
    An application to surface reconstruction from SPOT5-HRS stereo imagery."
    In IAPRS vol XXXVI-1/W41 in ISPRS Workshop On Topographic Mapping From Space
    (With Special Emphasis on Small Satellites), Ankara, Turquie, 02-2006.

[2] M. Pierrot-Deseilligny, "MicMac, un lociel de mise en correspondance
    d'images, adapte au contexte geograhique" to appears in
    Bulletin d'information de l'Institut Geographique National, 2007.

Francais :

   MicMac est un logiciel de mise en correspondance d'image adapte
   au contexte de recherche en information geographique. Il s'appuie sur
   la bibliotheque de manipulation d'image eLiSe. Il est distibue sous la
   licences Cecill-B.  Voir en bas de fichier et  http://www.cecill.info.


English :

    MicMac is an open source software specialized in image matching
    for research in geographic information. MicMac is built on the
    eLiSe image library. MicMac is governed by the  "Cecill-B licence".
    See below and http://www.cecill.info.

Header-MicMac-eLiSe-25/06/2007*/

#include "TiepRed.h"

#if (!BUG_PUSH_XML_TIEP)


struct mPComparerClass {
  bool operator() (cPMulTiepRed * i,cPMulTiepRed * j) { return (i->Gain() > j->Gain());}
} mPComparer;

/**********************************************************************/
/*                                                                    */
/*                         cImageGrid                                 */
/*                                                                    */
/**********************************************************************/
cImageGrid::cImageGrid(cImageTiepRed * aImage, Box2dr & aBox , int aNumCellsX, int aNumCellsY, int aNumImages) :
		mImage (aImage),
		mBox (aBox),
		mNumCellsX (aNumCellsX),
		mNumCellsY (aNumCellsY)
{
	std::size_t numCells = aNumCellsX*aNumCellsY;
	mCellsPoints.resize(numCells); // We reserve space for the cells
	if (mImage->ImageId() == 0){ // This is a grid for th emaster image, in this case we also need counters
		mCellsImageCounters.resize(numCells); // we reserve space for the cells
		for (std::size_t i = 0; i < numCells ; i++){
			mCellsImageCounters[i].resize(aNumImages-1, 0); // for each cell we initialize the counters to 0, there is a counter for all the images except the current one (master)
		}
	}
}

cImageTiepRed &  cImageGrid::Image() {return *mImage;}

void cImageGrid::Add(cPMulTiepRed * aMultiPoint)
{
	const Pt2df & aPt = aMultiPoint->MergedHomolPoints()->GetVal(mImage->ImageId()); // Get the merge structure of the multi-point and from it get the 2D position of the tiePoint in the current image
	int cellIndex = CellIndex(aPt); // Get the index of the grid cell in which the tiePoint in the current image
	mCellsPoints[cellIndex].push_back(aMultiPoint);

	if (mImage->ImageId() == 0){
		const std::vector<U_INT2> & vecInd = aMultiPoint->MergedHomolPoints()->VecInd(); // We check in which images this point is present. It should be visible in 0 (the master image) and at least anther image
		for (std::size_t k = 0; k < vecInd.size(); k++){ // For all the images where the point is present
			if (vecInd[k] != 0) mCellsImageCounters.at(cellIndex).at(vecInd[k]-1)++; // We update the counters
		}
	}
}

void cImageGrid::SortCells()
{
	// Sort the the vector of multi-points in each cell. It uses the mPComparer
	for (std::vector<std::vector<cPMulTiepRed *> >::iterator itP=mCellsPoints.begin(); itP!=mCellsPoints.end();  itP++){
		std::vector<cPMulTiepRed *> & mPS = (*itP);
		std::sort (mPS.begin(), mPS.end(), mPComparer);
	}
}

int cImageGrid::CellIndex(const Pt2df & aPoint)
{
	int xpos = int((aPoint.x - mBox._p0.x) * mNumCellsX / (mBox._p1.x - mBox._p0.x));
	int ypos = int((aPoint.y - mBox._p0.y) * mNumCellsY / (mBox._p1.y - mBox._p0.y));

	if (xpos == mNumCellsX) xpos -= 1; // If it is in the edge of the box (in the maximum side) we need to put in the last tile
	if (ypos == mNumCellsY) ypos -= 1;

	return (ypos * mNumCellsX) + xpos;
}

int cImageGrid::NumPointsCell(int aCellIndex)
{
	// Gets the number of points in cell, but only points which are not selected fro removal
	int n = 0;
	std::vector<cPMulTiepRed *> & mPS = mCellsPoints[aCellIndex];
	for (std::vector<cPMulTiepRed *>::iterator itP=mPS.begin(); itP!=mPS.end();  itP++){
		if ((*itP)->Removed() == false){
			n++;
		}
	}
	return n;
}

int cImageGrid::NumPointsCellImage(int aCellIndex, int aImageId){
	return mCellsImageCounters[aCellIndex][aImageId-1];
}

std::vector<cPMulTiepRed *> & cImageGrid::CellPoints(int aCellIndex){
	return mCellsPoints[aCellIndex];
}

void cImageGrid::Remove(cPMulTiepRed * aMultiPoint, int aCellIndex){
	aMultiPoint->Remove();
	if (mImage->ImageId() == 0){
		const std::vector<U_INT2> & vecInd = aMultiPoint->MergedHomolPoints()->VecInd(); // We check in which images this point is present. It should be visible in 0 (the master image) and at least anther image
		for (std::size_t k = 0; k < vecInd.size(); k++){ // For all the images where the point is present
			if (vecInd[k] != 0) mCellsImageCounters[aCellIndex][vecInd[k]-1]--;
		}
	}
}
#endif


/*Footer-MicMac-eLiSe-25/06/2007

Ce logiciel est un programme informatique servant à la mise en
correspondances d'images pour la reconstruction du relief.

Ce logiciel est régi par la licence CeCILL-B soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL-B telle que diffusée par le CEA, le CNRS et l'INRIA
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargement,  à l'utilisation,  à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant
donné sa spécificité de logiciel libre, qui peut le rendre complexe à
manipuler et qui le réserve donc à des développeurs et des professionnels
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
logiciel à leurs besoins dans des conditions permettant d'assurer la
sécurité de leurs systèmes et ou de leurs données et, plus généralement,
à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
pris connaissance de la licence CeCILL-B, et que vous en avez accepté les
termes.
Footer-MicMac-eLiSe-25/06/2007*/
