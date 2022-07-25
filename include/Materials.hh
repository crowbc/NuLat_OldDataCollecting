//MATERIALS HEADER
//Contains information on all materials to be used in simulations
// phase this out and define materials using NIST manager in DetectorConstruction.cc

#ifndef Materials_H

#define Materials_H 1

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

class Materials
{
public:
	G4Element* H;
	G4Element* Be;
	G4Element* B;
	G4Element* C;
	G4Element* O;
	G4Element* N;
	G4Element* F;
	G4Element* Si;
	G4Element* Ca;
	G4Element* Fe;
	G4Element* Ni;
	G4Element* Cu;
	G4Element* Mo;
	G4Element* In;
	G4Element* Pb;
	G4Element* enriched6Li;
	G4Element* Al;
	
	G4Isotope* isoLi6;
	
	G4Material* pseudocumene;
	G4Material* teflonFEP;
	G4Material* acrylic;
	G4Material* LENSscint;
	G4Material* LENSbufferScint;
	G4Material* air;
	G4Material* lead;
	G4Material* limestone;
	G4Material* borosilicateGlass;
	G4Material* BoronPVT_1percent;
	G4Material* Li6PVT_1TenthsOfAPercent;
	G4Material* Li6PVT_2TenthsOfAPercent;
	G4Material* Li6PVT_5TenthsOfAPercent;
	G4Material* Li6PVT_10TenthsOfAPercent;
	G4Material* Li6PVT_15TenthsOfAPercent;
	G4Material* muMetal;
	G4Material* vacuum;
	G4Material* BeCuPhotoCathode;
	G4Material* matAl;
	
public: 
	Materials()
	{
		G4NistManager *nist = G4NistManager::Instance();
		
		// Useful Constants
		const G4double hcnm=1.239841939*keV;
        	//Photon energy range for energy dependent material responses
        	const size_t numEntries = 182;
        	G4double photonEnergy[numEntries] = {
        	  2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV, 2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV, 2.341*eV, 2.386*eV, //10
        	  2.433*eV, 2.481*eV, 2.487*eV, 2.496*eV, 2.506*eV, 2.516*eV, 2.524*eV, 2.531*eV, 2.539*eV, 2.547*eV, //20
        	  2.554*eV, 2.561*eV, 2.569*eV, 2.577*eV, 2.586*eV, 2.595*eV, 2.605*eV, 2.614*eV, 2.622*eV, 2.630*eV, //30
        	  2.638*eV, 2.646*eV, 2.653*eV, 2.660*eV, 2.669*eV, 2.676*eV, 2.681*eV, 2.688*eV, 2.693*eV, 2.698*eV, //40
        	  2.703*eV, 2.706*eV, 2.711*eV, 2.718*eV, 2.723*eV, 2.731*eV, 2.742*eV, 2.755*eV, 2.768*eV, 2.782*eV, //50
        	  2.793*eV, 2.803*eV, 2.811*eV, 2.819*eV, 2.829*eV, 2.837*eV, 2.845*eV, 2.853*eV, 2.860*eV, 2.867*eV, //60
        	  2.875*eV, 2.882*eV, 2.888*eV, 2.894*eV, 2.900*eV, 2.907*eV, 2.913*eV, 2.919*eV, 2.924*eV, 2.930*eV, //70
        	  2.937*eV, 2.942*eV, 2.948*eV, 2.954*eV, 2.960*eV, 2.968*eV, 2.976*eV, 2.983*eV, 2.991*eV, 3.001*eV, //80
        	  3.008*eV, 3.017*eV, 3.028*eV, 3.038*eV, 3.048*eV, 3.055*eV, 3.070*eV, 3.087*eV, 3.103*eV, 3.121*eV, //90
        	  3.138*eV, 3.155*eV, 3.173*eV, 3.191*eV, 3.220*eV, 3.250*eV, 3.281*eV, 3.313*eV, 3.344*eV, 3.375*eV, //100
        	  3.403*eV, 3.439*eV, 3.479*eV, 3.522*eV, 3.566*eV, 3.611*eV, 3.644*eV, 3.684*eV, 3.731*eV, 3.780*eV, //110
        	  3.831*eV, 3.868*eV, 3.892*eV, 3.910*eV, 3.921*eV, 3.934*eV, 3.946*eV, 3.957*eV, 3.970*eV, 3.994*eV, //120
        	  4.044*eV, 4.102*eV, 4.160*eV, 4.202*eV, 4.236*eV, 4.267*eV, 4.298*eV, 4.328*eV, 4.357*eV, 4.387*eV, //130
        	  4.422*eV, 4.455*eV, 4.494*eV, 4.563*eV, 4.607*eV, 4.616*eV, 4.624*eV, 4.627*eV, 4.628*eV, 4.633*eV, //140
        	  4.640*eV, 4.642*eV, 4.649*eV, 4.656*eV, 4.661*eV, 4.666*eV, 4.678*eV, 4.685*eV, 4.692*eV, 4.699*eV, //150
        	  4.706*eV, 4.713*eV, 4.720*eV, 4.727*eV, 4.740*eV, 4.751*eV, 4.763*eV, 4.775*eV, 4.788*eV, 4.798*eV, //160
        	  4.813*eV, 4.828*eV, 4.840*eV, 4.853*eV, 4.869*eV, 4.886*eV, 4.905*eV, 4.928*eV, 4.953*eV, 5.015*eV, //170
        	  5.099*eV, 5.143*eV, 5.174*eV, 5.202*eV, 5.235*eV, 5.265*eV, 5.294*eV, 5.330*eV, 5.413*eV, 5.493*eV, //180
        	  5.556*eV, 5.611*eV}; //182
		
	//***************************
	//*        Materials        *
	//***************************
	
	//&&&&&&&&&&&&&&&&&&&&&&&&&
	//&  Element Data Base    &
	//&&&&&&&&&&&&&&&&&&&&&&&&&
	
		G4double aH = 1.00797 * g/mole;
		G4double aisoLi6 = 6.015   * g/mole;
		G4double aBe = 9.012   * g/mole;
		G4double aB = 10.81    * g/mole;
		G4double aC = 12.01115 * g/mole;
		G4double aO = 15.9994  * g/mole;
		G4double aN = 14.0067  * g/mole;
		G4double aF = 18.9984  * g/mole;
		G4double aSi = 28.085   * g/mole;
		G4double aCa = 40.078   * g/mole;
		G4double aFe = 55.845   * g/mole;
		G4double aNi = 58.6934  * g/mole;
		G4double aCu = 63.546   * g/mole;
		G4double aMo = 95.95    * g/mole;
		G4double aIn = 114.82    * g/mole;
		G4double aPb = 207.19    * g/mole;
		G4double aAl = 26.98     * g/mole;
		
		G4int z;  // atomic number
		
		//z = 1;
		H  = nist->FindOrBuildElement("H");//new G4Element("Hydrogen", "H" , z, aH);
		//z = 4;
		Be  = nist->FindOrBuildElement("Be");//new G4Element("Beryllium", "Be" , z, aBe);
		//z = 5;
		B  = nist->FindOrBuildElement("B");//new G4Element("Boron", "B", z, aB);
		//z = 6;
		C  = nist->FindOrBuildElement("C");//new G4Element("Carbon", "C" , z, aC);
		//z = 7;
		N  = nist->FindOrBuildElement("N");//new G4Element("Nitrogen", "N" , z, aN);
		//z = 8;
		O  = nist->FindOrBuildElement("O");//new G4Element("Oxygen", "O" , z, aO);
		//z = 9;
		F  = nist->FindOrBuildElement("F");//new G4Element("Fluorine", "F" , z, aF);
		//z = 13;
		Al = nist->FindOrBuildElement("Al");//new G4Element("Aluminium", "Al", z, aAl);
		//z = 14;
		Si  = nist->FindOrBuildElement("Si");//new G4Element("Silicon", "Si" , z, aSi);
		//z = 20;
		Ca = nist->FindOrBuildElement("Ca");//new G4Element( "Calcium", "Ca", z, aCa);
		//z = 26;
		Fe = nist->FindOrBuildElement("Fe");//new G4Element( "Iron", "Fe", z, aFe);
		//z = 28;
		Ni = nist->FindOrBuildElement("Ni");//new G4Element( "Nickel", "Ni", z, aNi);
		//z = 29;
		Cu = nist->FindOrBuildElement("Cu");//new G4Element( "Copper", "Cu", z, aCu);
		//z = 42;
		Mo = nist->FindOrBuildElement("Mo");//new G4Element( "Molybdenum", "Mo", z, aMo);
		//z = 49;
		In  = nist->FindOrBuildElement("In");//new G4Element("Indium", "In" , z, aIn);
		//z = 82;
		Pb = nist->FindOrBuildElement("Pb");//new G4Element( "Lead", "Pb", z, aPb);
		
		z = 3;
		isoLi6  = new G4Isotope("Li6", z, 6, aisoLi6);
		enriched6Li = new G4Element("enriched Li6", "Li6", 1);
		enriched6Li->AddIsotope(isoLi6, 100*perCent);
		
		//&&&&&&&&&&&&&&&&&&&&&&&&&
		//&  Molecule Data Base   &
		//&&&&&&&&&&&&&&&&&&&&&&&&&
		
		G4double density;
		
		//Pseudocumene
		density = 0.875*g/cm3;
		pseudocumene = new G4Material("Pseudocumene",density,2);
		pseudocumene->AddElement(C, 9);
		pseudocumene->AddElement(H, 12);
		
		//TeflonFEP
		density = 2.150*g/cm3;
		teflonFEP = new G4Material("Teflon FEP",density,2);
		teflonFEP->AddElement(F, 10);
		teflonFEP->AddElement(C, 5);
		
		// Acrylic Material properties
		density = 1.180*g/cm3;
		acrylic = new G4Material("Acrylic",density,3);
		acrylic->AddElement(O, 2);
		acrylic->AddElement(H, 8);
		acrylic->AddElement(C, 5);
		
		// Acrylic Optical properties
		// see 2012 J. Phys.: Conf. Ser. 356 012049
		// Assumes that it is thin film acrylic (idealistic)
		G4double refractiveIndexAcrylic[numEntries];
		G4double absorptionAcrylic[numEntries];
		G4cout << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl;
		for(size_t i=0;i<numEntries;i++)
		{
			G4double wavelength = ((hcnm/photonEnergy[i])*nm);//*static_cast<G4double>(pow(10.,-3.));//why divide by 1000?
			refractiveIndexAcrylic[i] = 1.492;
	         
			//sqrt( 1.866120+
			 	//0.2085454*pow(wavelength,2.)+0.4806770*pow(wavelength,-2.)-
				//0.1840693*pow(wavelength,-4.)+0.03424849*pow(wavelength,-6.)-
				//0.002340796*pow(wavelength,-8.));
			absorptionAcrylic[i] = 10.0*m;  //absorptions will need to be tuned to data
			//if(refractiveIndexAcrylic[i]<1.49244 || refractiveIndexAcrylic[i]!=refractiveIndexAcrylic[i])
			// this is removing the nan from Cauchy equation... look for emperical measurments
			//{
				//refractiveIndexAcrylic[i]=1.49244;
			//}
			G4cout << wavelength << "nm: n = " << refractiveIndexAcrylic[i] << G4endl;
		}
		G4cout << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl << G4endl;
		
		G4MaterialPropertiesTable *acrylicMat = new G4MaterialPropertiesTable();
		acrylicMat->AddProperty("RINDEX",photonEnergy,refractiveIndexAcrylic,numEntries);
		acrylicMat->AddProperty("ABSLENGTH", photonEnergy, absorptionAcrylic, numEntries);
		acrylic->SetMaterialPropertiesTable(acrylicMat);
		
		// Lead
		density = 11.35*g/cm3;
		lead = new G4Material( "Lead", density, 1);
		lead->AddElement(Pb, 1);
		
		// Limestone
		density = 2360*kg/m3;
		limestone = new G4Material("Limestone",density,3);
		limestone->AddElement(Ca, 2);
		limestone->AddElement(C, 8);
		limestone->AddElement(O, 5);
		 
		// Vacuum
		vacuum = new G4Material( "vacuum", 1., 1.008*g/mole, 1.e-25*g/cm3,kStateGas, 
					273*kelvin, 3.8e-18*pascal );
		
	 	// Artificial Photocathode material density set such that 20nm ->1mm
		BeCuPhotoCathode = new G4Material( "BeCuPhotoCathode", 5.0/50000.0 * g/cm3, 2 );
		BeCuPhotoCathode->AddElement( Be, 1 );
		BeCuPhotoCathode->AddElement( Cu, 1 );
		// BeCuPhotoCathode Optical Properties
		G4double refractiveIndexBeCuPhotoCathode[numEntries];
		G4double absorptionBeCuPhotoCathode[numEntries];
		for(size_t i=0; i<numEntries; i++){
			refractiveIndexBeCuPhotoCathode[i] = 2.7;
			absorptionBeCuPhotoCathode[i] = 0.001*mm;
		}
		G4MaterialPropertiesTable *BeCuPhotoCathodeMat = new G4MaterialPropertiesTable();
		BeCuPhotoCathodeMat->AddProperty("RINDEX", photonEnergy, refractiveIndexBeCuPhotoCathode, numEntries);
		BeCuPhotoCathodeMat->AddProperty("ABSLENGTH", photonEnergy, absorptionBeCuPhotoCathode, numEntries);
		BeCuPhotoCathode->SetMaterialPropertiesTable( BeCuPhotoCathodeMat );
		
		// BorosilicateGlass Material Properties
		borosilicateGlass = new G4Material( "BorosilicateGlass", 2.65*g/cm3, 2 );
		borosilicateGlass->AddElement( Si, 1 );
		borosilicateGlass->AddElement( O, 2 );
		// BorosilicateGlass Optical Properties
		G4double refractiveIndexBorosilicateGlass[numEntries];
		G4double absorptionBorosilicateGlass[numEntries];
		for(size_t i=0; i<numEntries; i++){
			refractiveIndexBorosilicateGlass[i] = 1.55;
			absorptionBorosilicateGlass[i] = 10*m;
		}
		G4MaterialPropertiesTable *borosilicateGlassMat = new G4MaterialPropertiesTable();
		borosilicateGlassMat->AddProperty("RINDEX", photonEnergy, refractiveIndexBorosilicateGlass, numEntries);
		borosilicateGlassMat->AddProperty("ABSLENGTH", photonEnergy, absorptionBorosilicateGlass, numEntries);
		borosilicateGlass->SetMaterialPropertiesTable( borosilicateGlassMat );
		
		// Aluminium Material Properties
		matAl = new G4Material("Aluminium", 2.7*g/cm3, 1);
		matAl->AddElement( Al, 1);
		
		G4double reflectivityAlMetal[numEntries];
		G4double refractiveIndexAl[numEntries];
		G4double absorptionAlMetal[numEntries];
		G4double imgrefractiveIndexAl[numEntries];          
		for(size_t i=0; i<numEntries; i++){
			reflectivityAlMetal[i] = 0.95;
			absorptionAlMetal[i] = 0.1*m;
			refractiveIndexAl[i] = 1.5;
			imgrefractiveIndexAl[i] = 4.836;          
			
		}
		
		G4MaterialPropertiesTable *Al_SurfaceMPT = new G4MaterialPropertiesTable();
		//Al_SurfaceMPT->AddProperty("REFLECTIVITY", photonEnergy, reflectivityAlMetal, numEntries);   //overall
		Al_SurfaceMPT->AddProperty("ABSLENGTH", photonEnergy, absorptionAlMetal, numEntries );
		//Al_SurfaceMPT->AddProperty("RINDEX", photonEnergy, refractiveIndexAl, numEntries );
		//Al_SurfaceMPT->AddProperty("IMAGINARYRINDEX", photonEnergy, imgrefractiveIndexAl, numEntries );          
		matAl->SetMaterialPropertiesTable(Al_SurfaceMPT);
		
		//&&&&&&&&&&&&&&&&&&&&&&&&&
		//&  Mixture Data Base    &
		//&&&&&&&&&&&&&&&&&&&&&&&&&
		
		density = 0.950*g/cm3;
		LENSscint = new G4Material("LENSscint",density,2);
		LENSscint->AddElement(In, 8*perCent);
		LENSscint->AddMaterial(pseudocumene, 92*perCent);
		
		density = 0.875*g/cm3;
		LENSbufferScint = new G4Material("bufferScint",density,1);
		LENSbufferScint->AddMaterial(pseudocumene, 100*perCent);
		
		// Air Material Properies
		density = 1.290*mg/cm3;
		air = nist->FindOrBuildMaterial("G4_AIR");//new G4Material("air",density,2);
		//air->AddElement(N, 70*perCent);
		//air->AddElement(O, 30*perCent);
		// Air Optics
		G4double refractiveIndexAir[numEntries];
		G4double absorptionAir[numEntries];
		for(size_t i=0; i<numEntries; i++){
			refractiveIndexAir[i] = 1.003;
			absorptionAir[i] = 100*m;
		}
		G4MaterialPropertiesTable *airMat = new G4MaterialPropertiesTable();
		airMat->AddProperty("RINDEX", photonEnergy, refractiveIndexAir, numEntries );
		airMat->AddProperty("ABSLENGTH", photonEnergy, absorptionAir, numEntries );
		air->SetMaterialPropertiesTable(airMat);
		
		// 1% by wt Boron loaded polyvinyltoluene
		G4double CAtomsPerCC=4.62, HAtomsPerCC=5.16, BAtomsPerCC=1.14/100/.199;
		G4double CMassPerCC=CAtomsPerCC*aC, HMassPerCC=HAtomsPerCC*aH, BMassPerCC=BAtomsPerCC*aB;
		
		G4double C_fractionmass = CMassPerCC/(CMassPerCC + HMassPerCC + BMassPerCC);
		G4double H_fractionmass = HMassPerCC/(CMassPerCC + HMassPerCC + BMassPerCC);
		G4double B_fractionmass = BMassPerCC/(CMassPerCC + HMassPerCC + BMassPerCC);
		
		density = 1.021*g/cm3;
		BoronPVT_1percent = new G4Material("BoronPVT_1percent",density, 3);
		BoronPVT_1percent->AddElement(H, H_fractionmass);
		BoronPVT_1percent->AddElement(C, C_fractionmass);
		BoronPVT_1percent->AddElement(B, B_fractionmass);
		
		G4cout << BoronPVT_1percent << G4endl;
		
		G4double C_StandardFractionMass = C_fractionmass * 1.01;
		G4double H_StandardFractionMass = H_fractionmass * 1.01;
		
		// Li 6 Generic optical properties
		// scintillator
		// see Journal of Colloid and Interface Science Volume 118, Issue 2, August 1987, Pages 314–325
		// gives Sellmeier Coefficients
		
		G4double refractiveIndexPlastic[numEntries];
		G4double absorptionPlastic[numEntries];
		G4double rayleigh_scattering[numEntries];
		
		for(size_t i=0; i<numEntries; i++)  {
		// refractive index assumes PVT, absorption is something to be tuned
			G4double wavelength = (hcnm/photonEnergy[i])*nm;
			refractiveIndexPlastic[i] = 1.58;
			absorptionPlastic[i]=90.*cm;
			rayleigh_scattering[i] =90.0*cm;
		}
		
		// Fast scintillation properties, no reference and were here with the original file
		G4double fastComponentPlastic[numEntries] = {
			0.000,  0.000,  0.000,  0.000,  0.000,  0.010,  0.020,  0.035,  0.050,  0.060, //10
			0.070,  0.085,  0.090,  0.095,  0.098,  0.100,  0.110,  0.120,  0.130,  0.140, //20
			0.150,  0.160,  0.170,  0.180,  0.200,  0.220,  0.240,  0.250,  0.270,  0.290, //30  
			0.300,  0.320,  0.340,  0.350,  0.360,  0.390,  0.400,  0.420,  0.430,  0.440, //40  
			0.445,  0.450,  0.460,  0.470,  0.480,  0.500,  0.550,  0.600,  0.630,  0.700, //50  
			0.730,  0.750,  0.800,  0.830,  0.850,  0.870,  0.900,  0.920,  0.940,  0.950, //60  
			0.960,  0.970,  0.980,  0.985,  0.990,  0.995,  1.000,  1.000,  1.000,  0.995, //70  
			0.990,  0.985,  0.980,  0.970,  0.960,  0.930,  0.900,  0.870,  0.850,  0.800, //80 
			0.700,  0.600,  0.500,  0.400,  0.300,  0.220,  0.130,  0.070,  0.010,  0.000, //90  
			0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //100
			0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //110  
			0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //120  
			0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //130  
			0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //140  
			0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //150  
			0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //160  
			0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //170  
			0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //180  
			0.000,  0.000}; //182
		
		G4MaterialPropertiesTable *PVT_MPT = new G4MaterialPropertiesTable();      
		PVT_MPT->AddProperty("RINDEX", photonEnergy, refractiveIndexPlastic, numEntries);
		PVT_MPT->AddProperty("ABSLENGTH", photonEnergy, absorptionPlastic, numEntries); 
		PVT_MPT->AddProperty("FASTCOMPONENT",photonEnergy, fastComponentPlastic, numEntries, true);
		// add a slow compontent at some point
		PVT_MPT->AddConstProperty("SCINTILLATIONYIELD",10000./MeV);//10000./MeV);
		PVT_MPT->AddConstProperty("RESOLUTIONSCALE", 1.0, true);    
		PVT_MPT->AddConstProperty("FASTTIMECONSTANT", 2.1*ns, true);    
		PVT_MPT->AddConstProperty("YIELDRATIO", 1., true);
		//PVT_MPT->AddProperty("RAYLEIGH",photonEnergy,rayleigh_scattering,numEntries); 
		
		// 0.1% by wt 6Li loaded polyvinyltoluene  
		G4double Li6_fractionmass = 0.001;
		H_fractionmass = C_StandardFractionMass/(1.0+Li6_fractionmass);
		C_fractionmass = H_StandardFractionMass/(1.0+Li6_fractionmass);
		
		density = 1.021*g/cm3;
		Li6PVT_1TenthsOfAPercent = new G4Material("Li6PVT_1TenthsOfAPercent",density, 3);
		Li6PVT_1TenthsOfAPercent->AddElement(H, H_fractionmass);
		Li6PVT_1TenthsOfAPercent->AddElement(C, C_fractionmass);
		Li6PVT_1TenthsOfAPercent->AddElement(enriched6Li, Li6_fractionmass);
		Li6PVT_1TenthsOfAPercent->SetMaterialPropertiesTable(PVT_MPT);
		Li6PVT_1TenthsOfAPercent->GetIonisation()->SetBirksConstant(0.01045*cm/MeV);
		
		G4cout << Li6PVT_1TenthsOfAPercent << G4endl << G4endl << G4endl;
		
		// 0.2% by wt 6Li loaded polyvinyltoluene  
		Li6_fractionmass = 0.002;
		H_fractionmass = C_StandardFractionMass/(1.0+Li6_fractionmass);
		C_fractionmass = H_StandardFractionMass/(1.0+Li6_fractionmass);
		
		density = 1.021*g/cm3;
		Li6PVT_2TenthsOfAPercent = new G4Material("Li6PVT_1TenthsOfAPercent",density, 3);
		Li6PVT_2TenthsOfAPercent->AddElement(H, H_fractionmass);
		Li6PVT_2TenthsOfAPercent->AddElement(C, C_fractionmass);
		Li6PVT_2TenthsOfAPercent->AddElement(enriched6Li, Li6_fractionmass);
		Li6PVT_2TenthsOfAPercent->SetMaterialPropertiesTable(PVT_MPT);
		Li6PVT_2TenthsOfAPercent->GetIonisation()->SetBirksConstant(0.01045*cm/MeV);
		
		G4cout << Li6PVT_2TenthsOfAPercent << G4endl << G4endl << G4endl;
		
		// 0.5% by wt 6Li loaded polyvinyltoluene  
		Li6_fractionmass = 0.005;
		H_fractionmass = C_StandardFractionMass/(1.0+Li6_fractionmass);
		C_fractionmass = H_StandardFractionMass/(1.0+Li6_fractionmass);
		
		density = 1.021*g/cm3;
		Li6PVT_5TenthsOfAPercent = new G4Material("Li6PVT_5TenthsOfAPercent",density, 3);
		Li6PVT_5TenthsOfAPercent->AddElement(H, H_fractionmass);
		Li6PVT_5TenthsOfAPercent->AddElement(C, C_fractionmass);
		Li6PVT_5TenthsOfAPercent->AddElement(enriched6Li, Li6_fractionmass);
		Li6PVT_5TenthsOfAPercent->SetMaterialPropertiesTable(PVT_MPT);
		Li6PVT_5TenthsOfAPercent->GetIonisation()->SetBirksConstant(0.01045*cm/MeV);
		
		G4cout << Li6PVT_5TenthsOfAPercent << G4endl << G4endl << G4endl;
		
		// 1.0% by wt 6Li loaded polyvinyltoluene  
		Li6_fractionmass = 0.010;
		H_fractionmass = C_StandardFractionMass/(1.0+Li6_fractionmass);
		C_fractionmass = H_StandardFractionMass/(1.0+Li6_fractionmass);
		
		density = 1.021*g/cm3;
		Li6PVT_10TenthsOfAPercent = new G4Material("Li6PVT_10TenthsOfAPercent",density, 3);
		Li6PVT_10TenthsOfAPercent->AddElement(H, H_fractionmass);
		Li6PVT_10TenthsOfAPercent->AddElement(C, C_fractionmass);
		Li6PVT_10TenthsOfAPercent->AddElement(enriched6Li, Li6_fractionmass);
		Li6PVT_10TenthsOfAPercent->SetMaterialPropertiesTable(PVT_MPT);
		Li6PVT_10TenthsOfAPercent->GetIonisation()->SetBirksConstant(0.01045*cm/MeV);
		
		G4cout << Li6PVT_10TenthsOfAPercent << G4endl << G4endl << G4endl;
		
		// 1.5% by wt 6Li loaded polyvinyltoluene  
		Li6_fractionmass = 0.015;
		H_fractionmass = C_StandardFractionMass/(1.0+Li6_fractionmass);
		C_fractionmass = H_StandardFractionMass/(1.0+Li6_fractionmass);
		
		density = 1.021*g/cm3;
		Li6PVT_15TenthsOfAPercent = new G4Material("Li6PVT_15TenthsOfAPercent",density, 3);
		Li6PVT_15TenthsOfAPercent->AddElement(H, H_fractionmass);
		Li6PVT_15TenthsOfAPercent->AddElement(C, C_fractionmass);
		Li6PVT_15TenthsOfAPercent->AddElement(enriched6Li, Li6_fractionmass);
		Li6PVT_15TenthsOfAPercent->SetMaterialPropertiesTable(PVT_MPT);
		Li6PVT_15TenthsOfAPercent->GetIonisation()->SetBirksConstant(0.01045*cm/MeV);
		
		G4cout << Li6PVT_15TenthsOfAPercent << G4endl << G4endl << G4endl;
		// materials
		// Water, Steel, Quartz, Stainless steel (304, 316, etc), tungsten, teflon (all others), polyester, polyethyliene..
		
		// MuMetal  
		muMetal = new G4Material( "muMetal",8.77*g/cm3,4);
		muMetal->AddElement(Ni,77.*perCent);
		muMetal->AddElement(Fe,16.*perCent);
		muMetal->AddElement(Cu,5.*perCent);
		muMetal->AddElement(Mo,2.*perCent);   
		// Mu-metal Surface Properties
		G4double reflectivityMuMetal[numEntries];
		G4double specularLopeMuMetal[numEntries];
		G4double specularSpikeMuMetal[numEntries];
		G4double backscatterMuMetal[numEntries];
		G4double efficiencyMuMetal[numEntries];
		for(size_t i=0; i<numEntries; i++){
			specularLopeMuMetal[i] = 0.0;
			specularSpikeMuMetal[i] = 0.;
			backscatterMuMetal[i] = 0.0;
			reflectivityMuMetal[i] = 0.0;
			efficiencyMuMetal[i] = 0.;
		}
		G4MaterialPropertiesTable *muMetal_SurfaceMPT = new G4MaterialPropertiesTable();
		muMetal_SurfaceMPT->AddProperty("REFLECTIVITY", photonEnergy, 
						reflectivityMuMetal, numEntries);// overall
		muMetal_SurfaceMPT->AddProperty("SPECULARLOBECONSTANT", photonEnergy, 
						specularLopeMuMetal, numEntries);// normal of micro-facet
		muMetal_SurfaceMPT->AddProperty("SPECULARSPIKECONSTANT", photonEnergy, 
						specularSpikeMuMetal, numEntries);// mirror
		muMetal_SurfaceMPT->AddProperty("BACKSCATTERCONSTANT", photonEnergy,
						backscatterMuMetal, numEntries);// backscatter
		muMetal_SurfaceMPT->AddProperty("EFFICIENCY",photonEnergy,efficiencyMuMetal,
						numEntries);
		muMetal->SetMaterialPropertiesTable(muMetal_SurfaceMPT);
		
	}
};

#endif
