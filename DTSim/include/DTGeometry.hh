//
//

#ifndef DTGeometry_h
#define DTGeometry_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"

#include <vector>


class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;
class G4UniformMagField;
class G4Material;

/// Detector construction

class DTGeometry : public G4VUserDetectorConstruction
{
  public:
    DTGeometry();
    virtual ~DTGeometry();
    
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;
    
    void ConstructMaterials();
    
    void buildSuperLayer(G4LogicalVolume* superLayerLV, G4Material* cellMaterial, G4int nofLayers, G4int nofCells, G4double layerThickness, G4double cellWidth, G4double dtWidth);

  private:
    
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();


    // data members
    //
    G4LogicalVolume* fYokeLV = nullptr;
    G4LogicalVolume* fSuperLayer1LV = nullptr;
    G4LogicalVolume* fSuperLayer2LV = nullptr;
    G4LogicalVolume* fSuperLayer3LV = nullptr;

    static G4ThreadLocal G4UniformMagField* fMagneticField;
    static G4ThreadLocal G4FieldManager* fFieldMgr;

    G4bool fCheckOverlaps = true; // option to activate checking of volumes overlaps

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
