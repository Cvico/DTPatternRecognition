//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file DTGeometry.cc
/// \brief Implementation of the DTGeometry class

#include "DTGeometry.hh"
#include "SuperLayerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"
#include "G4SDChargedFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4ThreadLocal G4UniformMagField* DTGeometry::fMagneticField = nullptr;
G4ThreadLocal G4FieldManager* DTGeometry::fFieldMgr = nullptr;

using namespace DTSim;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DTGeometry::DTGeometry()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true)
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DTGeometry::~DTGeometry()
{

}

G4VPhysicalVolume* DTGeometry::Construct()
{
  // Construct materials
  DefineMaterials();
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DTGeometry::ConstructSDandField()
{
  // sensitive detectors -----------------------------------------------------
  auto sdManager = G4SDManager::GetSDMpointer();
  sdManager->SetVerboseLevel(1);
  // 
  // Scorers
  //

  // declare cells as a MultiFunctionalDetector scorer
  //  
  G4String SDname;

  auto superlayer1 = new SuperLayerSD(SDname="/sl1");
  sdManager->AddNewDetector(superlayer1);
  fSuperLayer1LV->SetSensitiveDetector(superlayer1);

  auto superlayer2 = new SuperLayerSD(SDname="/sl2");
  sdManager->AddNewDetector(superlayer2);
  fSuperLayer2LV->SetSensitiveDetector(superlayer2);

  auto superlayer3 = new SuperLayerSD(SDname="/sl3");
  sdManager->AddNewDetector(superlayer2);
  fSuperLayer3LV->SetSensitiveDetector(superlayer3);

  //
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.

  // magnetic field ----------------------------------------------------------
  fMagneticField = new G4UniformMagField(G4ThreeVector(0., 2*tesla, 0.));
  fFieldMgr = new G4FieldManager();
  fFieldMgr->SetDetectorField(fMagneticField);
  fFieldMgr->CreateChordFinder(fMagneticField);
  fYokeLV->SetFieldManager(fFieldMgr, true);

  // Register the field and its manager for deleting
  G4AutoDelete::Register(fMagneticField);
  G4AutoDelete::Register(fFieldMgr);  
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DTGeometry::DefineMaterials()
{
  auto nistManager = G4NistManager::Instance();

  // Air 
  auto air = nistManager->FindOrBuildMaterial("G4_AIR");
  G4double air_density = air->GetDensity();

  // Gas mixture
  auto gas_mixture = new G4Material("GasMixture", air_density*g/cm3, 2);
  gas_mixture->AddMaterial(nistManager->FindOrBuildMaterial("G4_Ar"), 85*perCent);
  gas_mixture->AddMaterial(nistManager->FindOrBuildMaterial("G4_CARBON_DIOXIDE"), 15*perCent);  
  // Aluminium - Honeycomb (GAP)
  nistManager->FindOrBuildMaterial("G4_Al");
  
  // iron Yoke
  nistManager->FindOrBuildMaterial("G4_Fe");
  
  // Vacuum "Galactic"
  nistManager->FindOrBuildMaterial("G4_Galactic");

  // Vacuum "Air with low density"
  // auto air = G4Material::GetMaterial("G4_AIR");
  // G4double density = 1.0e-5*air->GetDensity();
  // nistManager
  //   ->BuildMaterialWithNewDensity("Air_lowDensity", "G4_AIR", density);

  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DTGeometry::DefineVolumes()
{
  // Geometry parameters
  G4int  nofLayers = 4;
  G4int  nofCells = 50;
  G4int  nofSuperLayers = 3; 
  G4double  cellThickness = 13.*mm;
  G4double  cellWidth     = 42.*mm;
  G4double  gapThickness  = 23.5*cm;
  G4double  yokeThickness = 50*cm;  // The other flavour is 63cm

  auto  layerThickness = cellThickness;
  auto  superlayerThickness = nofLayers * layerThickness;
  auto  dtThickness = superlayerThickness * nofSuperLayers + gapThickness + yokeThickness;
  auto  dtWidth = cellWidth * nofCells;

  auto  worldSizeXY = 1.2 * dtWidth;
  auto  worldSizeZ  = 1.2 * dtThickness; 
  
  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("G4_Galactic");
  auto yokeMaterial = G4Material::GetMaterial("G4_Fe");
  auto cellMaterial = G4Material::GetMaterial("GasMixture");
  auto gapMaterial = G4Material::GetMaterial("G4_Al");
  
  // POSITIONS OF DIFFERENT SUBDETECTORS...
  auto zpos_yoke = -dtThickness/2+yokeThickness/2;
  auto zpos_superlayer1 = -dtThickness/2+yokeThickness+superlayerThickness/2;
  auto zpos_honeycomb = -dtThickness/2+yokeThickness+superlayerThickness+gapThickness/2;
  auto zpos_superlayer2 = -dtThickness/2+yokeThickness+superlayerThickness+gapThickness+superlayerThickness/2;
  auto zpos_superlayer3 = -dtThickness/2+yokeThickness+superlayerThickness+gapThickness+superlayerThickness+superlayerThickness/2;

  if ( ! defaultMaterial || ! yokeMaterial || ! cellMaterial || !gapMaterial) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("DTGeometry::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
   
  //     
  // World
  //
  auto worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // DT chamber
  //  
  auto dtchamberS
    = new G4Box("DTChamber",     // its name
                 dtWidth/2, dtWidth/2, dtThickness/2); // its size
                         
  auto dtchamberLV
    = new G4LogicalVolume(
                 dtchamberS,    // its solid
                 defaultMaterial, // its material
                 "DTChamber");  // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 dtchamberLV,          // its logical volume                         
                 "DTChamber",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //
  // Iron Yoke
  //
  auto yokeS 
    = new G4Box("Yoke",            // its name
                 dtWidth/2, dtWidth/2, yokeThickness/2); // its size
                         
  fYokeLV
    = new G4LogicalVolume(
                 yokeS,        // its solid
                 yokeMaterial, // its material
                 "Yoke");          // its name
                                   
   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., zpos_yoke), //  its position
                 fYokeLV,       // its logical volume                         
                 "Yoke",           // its name
                 dtchamberLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                                 
  // Super Layer
  //
  auto superLayer1S 
    = new G4Box("SuperLayer1",           // its name
                 dtWidth/2, dtWidth/2, superlayerThickness/2); // its size
  
  fSuperLayer1LV
    = new G4LogicalVolume(
                 superLayer1S,   // its solid
                 cellMaterial,  // its material
                 "SuperLayer1"); // its name
  
  auto superLayer2S 
    = new G4Box("SuperLayer2",           // its name
                 dtWidth/2, dtWidth/2, superlayerThickness/2); // its size
  
  fSuperLayer2LV
    = new G4LogicalVolume(
                 superLayer2S,   // its solid
                 cellMaterial,  // its material
                 "SuperLayer2"); // its name

  auto superLayer3S 
    = new G4Box("SuperLayer3",           // its name
                 dtWidth/2, dtWidth/2, superlayerThickness/2); // its size
  
  fSuperLayer3LV
    = new G4LogicalVolume(
                 superLayer3S,   // its solid
                 cellMaterial,  // its material
                 "SuperLayer3"); // its name

  buildSuperLayer(fSuperLayer1LV, cellMaterial, nofLayers, nofCells, layerThickness, cellWidth, dtWidth);
  buildSuperLayer(fSuperLayer2LV, cellMaterial, nofLayers, nofCells, layerThickness, cellWidth, dtWidth);
  buildSuperLayer(fSuperLayer3LV, cellMaterial, nofLayers, nofCells, layerThickness, cellWidth, dtWidth);

  new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0.,  zpos_superlayer1), //  its position
                  fSuperLayer1LV,          // its logical volume
                  "SuperLayer",    // its name
                  dtchamberLV,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps

  //                               
  // Honeycomb
  //
  auto gapS 
    = new G4Box("Honeycomb",            // its name
                 dtWidth/2, dtWidth/2, gapThickness/2); // its size
                         
  auto gapLV
    = new G4LogicalVolume(
                 gapS,        // its solid
                 gapMaterial, // its material
                 "Honeycomb");          // its name
                                   
   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0.,  zpos_honeycomb), //  its position
                 gapLV,       // its logical volume                         
                 "Honeycomb", // its name
                 dtchamberLV, // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  // 
  // Now SL-2 (90º rotation) 
  // 

  // Rotate logical volumes by 90 degrees in X-axis
  G4RotationMatrix* sl2_rot = new G4RotationMatrix;
  sl2_rot->rotateZ(90.*deg);

  new G4PVPlacement(sl2_rot,G4ThreeVector(0.,0.,zpos_superlayer2),fSuperLayer2LV,
                      "SuperLayer2",dtchamberLV,
                      false,1,fCheckOverlaps);
  
  // Now SL-3 
  new G4PVPlacement(0,G4ThreeVector(0.,0.,zpos_superlayer3),fSuperLayer3LV,
                      "SuperLayer3",dtchamberLV,
                      false,2,fCheckOverlaps);


  //
  // print parameters
  //
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The DT chamber is made of " << G4endl
    << yokeThickness/mm << "mm of " << yokeMaterial->GetName() << G4endl
    << nofSuperLayers << " superlayers of " << nofLayers << " layers with "  
    << nofCells << " cells each.  The cells are " << cellThickness/mm << "mm of " << cellMaterial->GetName() << G4endl
    << "The gap between superlayers 1 and 2 is " << gapThickness/mm << "mm of " << gapMaterial->GetName() << G4endl
    << "------------------------------------------------------------" << G4endl;
  
  //                                        
  // Visualization attributes
  //
  G4VisAttributes invisible(G4VisAttributes::GetInvisible());
  G4VisAttributes invisibleBlue(false, G4Colour::Blue());
  G4VisAttributes invisibleGreen(false, G4Colour::Green());
  G4VisAttributes invisibleYellow(false, G4Colour::Yellow());
  G4VisAttributes blue(G4Colour::Blue());
  G4VisAttributes cgray(G4Colour::Gray());
  G4VisAttributes green(G4Colour::Green());
  G4VisAttributes red(G4Colour::Red());
  G4VisAttributes yellow(G4Colour::Yellow());

  worldLV->SetVisAttributes (invisible);

  fYokeLV->SetVisAttributes(red);
  fSuperLayer1LV->SetVisAttributes(blue);
  fSuperLayer2LV->SetVisAttributes(blue);
  fSuperLayer3LV->SetVisAttributes(blue);
  gapLV->SetVisAttributes(green);
  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DTGeometry::buildSuperLayer(G4LogicalVolume* superLayerLV, G4Material* cellMaterial, G4int nofLayers, G4int nofCells, G4double layerThickness, G4double cellWidth, G4double dtWidth)
{

  //                               
  // Layer
  //
  auto layerS 
    = new G4Box("Layer",           // its name
                 dtWidth/2, dtWidth/2, layerThickness/2); // its size                      
  auto layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 cellMaterial,  // its material
                 "Layer");         // its name
  
  // 
  // Cell 
  // 
  auto cellS 
    = new G4Box("Cell",            // its name
                 cellWidth/2, dtWidth/2, layerThickness/2); // its size
  
  auto cellLV
    = new G4LogicalVolume(
                 cellS,        // its solid
                 cellMaterial, // its material
                 "Cell");          // its name

  // Fill each layer with cells
  for (auto l=0;l<nofLayers;l++) {
    G4double zlayer = (l-nofLayers/2)*layerThickness;
    for (auto c=0;c<nofCells;c++) {
      G4double shift = (l%2)*cellWidth/2;
      G4double xcell = (c-nofCells/2)*cellWidth+shift;
      new G4PVPlacement(0,G4ThreeVector(xcell,0.,0.),cellLV,
                        "Cell",layerLV,
                        false,c,fCheckOverlaps);
    }
    new G4PVPlacement(0,G4ThreeVector(0.,0.,zlayer),layerLV,
                        "Layer",superLayerLV,
                        false,l,fCheckOverlaps);

  }

  G4VisAttributes blue(G4Colour::Blue());
  G4VisAttributes cgray(G4Colour::Gray());

  cellLV->SetVisAttributes(cgray);
  layerLV->SetVisAttributes(blue);

}