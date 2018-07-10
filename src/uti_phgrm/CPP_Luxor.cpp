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


#include "StdAfx.h"


/*typedef enum
{
	eBBA,      //bundle block adjustment  ie Tapas
	eSBBA,     //structure-less BBA 	  ie Martini
	eSBBAFus   //SBBA and BBA intertwined ie Martini et Tapas

}eTypeOriVid;
*/


/* to do :

fusion
- faire coherent les sorties de BBA et SBBA pourque ; peut-etre faudrait re-ecrire un peu pour que InOri soit mis plus globalement?
- 


*/

void BanniereLuxor()
{
    std::cout <<  "\n";
    std::cout <<  " *********************************************\n";
    std::cout <<  " *     L-inear                               *\n";
    std::cout <<  " *     U-ltra                                *\n";
    std::cout <<  " *     X (Q)uick (S)liding                   *\n";
    std::cout <<  " *     OR-ientation                          *\n";
    std::cout <<  " *********************************************\n\n";

}



class cAppliLuxor
{
	public :
		cAppliLuxor(int argc, char** argv);

			
	private :
	 
		int                        mWSz;        //sliding window size
		int				           mNbWinAll;	//number of all windows
		int				           mWPas;		//sliding window step (how many frames added at a time)
		int				           mBBAPas;		//bba step
		std::list< std::string >   mWName;      //list of all frames

		std::string        mInOri;
		std::string        mInCal;
		std::string  	   mDir;
		std::string        mIms;
		std::string        mSH;
		std::string        mOut;

        bool DoBASCULE;
        bool DoLBA; 

		cInterfChantierNameManipulateur 			* mICNM; 
	    const cInterfChantierNameManipulateur::tSet * mSetIm;	
		int 										  mNbIm;

		bool        			mModeHelp;	
		std::string 			mStrType;	
		eTypeOriVid 			mType;

		void        ReadType    (const std::string & aType);	
		std::string MakeFenName (const int aNum);
		void        CalculFen   ();


		void DoBBAGlob     ();
		void DoBBA         ();
		void DoBBA         (const std::string &,bool Init=false);
		void DoSBBA        ();
		void DoSBBA        (const std::string & aNameCur,
                            const std::string & aNamePrev="");
		void BasculeMorito (const std::string & aNameCur,
                            const std::string & aNamePrev);
		void LBACampari    (const std::list<std::string> );	
		 
		
};

std::string cAppliLuxor::MakeFenName (const int aNum)
{
	return "Fen-" + ToString(aNum) + ".xml";
}

void cAppliLuxor::CalculFen()
{

	//calculate the windows
	std::vector< std::list< std::string >> aWVTmp;

	std::cout << "fenetres sous traitement: \n";
	for (int aW=0; aW<mNbWinAll; aW++)
	{
		cListOfName              aXml;
		std::list< std::string > aImInWL;

		for (int aIm=0; aIm<mWSz; aIm++)
		{	
			int aIGlob ;   
			aIGlob = aW * mWPas + aIm;

		
			if (aIGlob<mNbIm)
				aImInWL.push_back((*mSetIm)[aIGlob]);
			
			std::cout << " Num=" << aIGlob << " " << (*mSetIm)[aIGlob] << " \n" ;

		}
		std::cout << ".." << aW+1 << "/" << mNbWinAll << "\n" ;
		
		aWVTmp.push_back(aImInWL);

		aXml.Name() = aImInWL;
		MakeFileXML(aXml,MakeFenName(aW));
		mWName.push_back(MakeFenName(aW));
	}

getchar();
	//add a "tail" to each window    --- move the xml save outside
	/*if (mNbWinVois)	
	{
		for (int aW=mNbWinAll-1; aW>=0; aW--)
		{
			aWinsV->at(aW) 	
		}
	}
	*/
	
	

	
}

void cAppliLuxor::LBACampari(const std::list<std::string> )
{

}

void cAppliLuxor::BasculeMorito (const std::string & aNameCur,
                                 const std::string & aNamePrev)
{
	std::string aCom = MMBinFile("mm3d Morito ") +
                        "Ori-Martini" + StdPrefix(aNamePrev) + "/Ori.*xml " +
                        "Ori-Martini" + StdPrefix(aNameCur) + "/Ori.*xml " +
                        "Martini" + StdPrefix(aNameCur) ; 
    
    std::cout << "CMD=" << aCom << "\n";                        
	TopSystem(aCom.c_str());
     

}


void cAppliLuxor::DoBBA(const std::string & aName,bool Init)
{


	std::string aCom = MMBinFile("mm3d Tapas Figee ") 
                       + "NKS-Set-OfFile@" + aName 
                       + " RefineAll=0 Out=" + mOut;

	if (EAMIsInit(&mInOri) && Init)	
		aCom += " InOri=" + mInOri;
	else if (!Init)
	{ 
		mInOri = mOut;
		aCom += " InOri=" + mInOri;

	}

	if (EAMIsInit(&mInCal))	
		aCom += " InCal=" + mInCal;

	//io figee
	//aCom += " LibFoc=0 DegRadMax=0 LibPP=0" ;  

	if (EAMIsInit(&mSH))	
		aCom += " SH=" + mSH;


	std::cout << "aCom=" << aCom << "\n";

	TopSystem(aCom.c_str());
}

void cAppliLuxor::DoBBA()
{
	int i=0;

	for (auto aW : mWName)
	{
		if (i==0)
			DoBBA(aW,true);
		else
			DoBBA(aW);
		i++;
	}	
}

void cAppliLuxor::DoSBBA(const std::string & aName,const std::string & aNamePrev)
{

	std::string aCom = MMBinFile("mm3d Martini ") 
                       + "NKS-Set-OfFile@" + aName;


    if (aNamePrev=="")
        aCom += " OriOut=Martini" + StdPrefix(aName);
    else
    {
        mInOri = "Martini" + StdPrefix(aNamePrev);
       
        aCom += " InOri=" + mInOri +
                " OriOut=Martini" + StdPrefix(aName);

    }


	if (EAMIsInit(&mInCal))	
		aCom += " OriCalib=" + mInCal;

	if (EAMIsInit(&mSH))	
		aCom += " SH=" + mSH;

	std::cout << "aCom=" << aCom << "\n";

	TopSystem(aCom.c_str());
	
}

/* Two options : 
   a- DoSBBA DoBASCULE=false - Martini with    InOri from prev sliding window to keep the repere
   b- DoSBBA DoBASCULE=true  - Martini without InOri but with Morito    
    then, LBA with Campari if desired */
void cAppliLuxor::DoSBBA()
{
	int i      =0;

    std::list<std::string>::iterator aW = mWName.begin();
    for ( ; aW!=mWName.end(); aW++)
    {
        if (i==0)
            DoSBBA(*aW);
        else if(DoBASCULE)
        {
            DoSBBA(*aW);
            BasculeMorito (*aW,*std::prev(aW,1));
        }
        else
			DoSBBA(*aW,*std::prev(aW,i));

        if (DoLBA)
        {
            if ( (i+1)%mBBAPas == 0)
            {
                std::list<std::string> aLBAWAll;
                for (int aLBAWTmp=mBBAPas; aLBAWTmp>0; aLBAWTmp--)
                {
                    std::cout << "ewwww="<<  i << " " << *aW << " " << "\n";
                    std::cout << "ewwwwwwww="<< aLBAWTmp << " " << *std::prev(aW,aLBAWTmp) << "\n";
                    aLBAWAll.push_back( *std::prev(aW,aLBAWTmp) );
                }
        getchar();
                LBACampari(aLBAWAll);
                
            }
        }
                
        i++;
    }
}





cAppliLuxor::cAppliLuxor(int argc, char** argv) :
	mWSz(5),
	mWPas(2),
	mBBAPas(2),
	mInOri(""),
	mInCal(""),
	mSH(""),
	mOut("Luxor"),
    DoBASCULE(false),
    DoLBA(false)
{
	std::string aPattern;

	ElInitArgMain
	(
		argc, argv,
		LArgMain() << EAMC(mStrType,"Orientation mode (enum values)")
		           << EAMC(aPattern,"Pattern of images")
				   << EAMC(mWSz,"Processing window size"),
		LArgMain() << EAM (mWPas,"M",true,"Motion of the processing window in frames; Def=2")  
				   << EAM (mBBAPas,"BbaPas",true,"BBa step; Def=2 every other processing window")  
				   << EAM (mInOri,"InOri",true,"Input external orientation")  
				   << EAM (mInCal,"InCal",true,"Input internal orientation")  
				   << EAM (DoBASCULE,"Basc",true,"Do Bascule with Morito, Def=false")  
				   << EAM (DoLBA,"LBA",true,"Do LBA with Campari, Def=false")  
				   << EAM (mSH,"SH",true,"Homol prefix") 
				   << EAM (mOut,"Out",true,"Output orientation") 
	);

	
#if (ELISE_windows)
      replace( aPattern.begin(), aPattern.end(), '\\', '/' );
#endif
    SplitDirAndFile(mDir,mIms,aPattern);
    StdCorrecNameOrient(mInOri,mDir);

	ReadType(mStrType);



	mICNM = cInterfChantierNameManipulateur::BasicAlloc(mDir);
    mSetIm = mICNM->Get(mIms);
    mNbIm = (int)mSetIm->size();


	mNbWinAll = std::floor(double((mNbIm-1) - (mWSz-1))/mWPas) +1;

	CalculFen();

	//InitDefValFromType
	switch (mType)
	{
		case eBBA :
			DoBBA();
			break;

		case eSBBA :
			DoSBBA();
			break;

		case eSBBAFus :
			break;

		case eUndefVal :
			break;
	}
}

void cAppliLuxor::ReadType(const std::string & aType)
{
    mStrType = aType;
    StdReadEnum(mModeHelp,mType,mStrType,eUndefVal);
/*	eTypeMalt   xType;
    StdReadEnum(mModeHelp,xType,mStrType,eNbTypesMNE);*/
}


int Luxor_main(int argc,char** argv)
{
	cAppliLuxor anAppli(argc,argv);

    BanniereLuxor();

	return EXIT_SUCCESS;
}

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
