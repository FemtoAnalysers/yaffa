#include "TGeoManager.h"
#include "TGeoMatrix.h"

void chat1() {
    // Create the geometry manager
    TGeoManager *geom = new TGeoManager("world", "The world");

    // Create a top volume
    // TGeoVolume *top = geom->MakeBox("Top", nullptr, 10, 10, 10);
    // geom->SetTopVolume(top);

    // Apply a rotation to change the z-axis to horizontal
    // For example, rotate 90 degrees around the x-axis
    TGeoRotation *rot = new TGeoRotation("rot", 90, 0, 0);
    TGeoCombiTrans *trans = new TGeoCombiTrans();
    trans->SetRotation(rot);


        TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98, 13, 2.7);
    //--- define some media
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum", 1, matVacuum);
    TGeoMedium *Al = new TGeoMedium("Aluminium", 2, matAl);

    // Add a box to the top volume with the new transformation
    TGeoVolume *box = geom->MakeBox("Box", nullptr, 2, 2, 2);
    geom->SetTopVolume(box);
    // top->AddNode(box, 1, trans);
    TGeoTube *tube2 = new TGeoTube(50, 50 + 0.1, 258);
    TGeoVolume *tubeV2 = new TGeoVolume("volTube", tube2, Al);
    box->AddNode(tubeV2, 1234);

    // Close the geometry
    geom->CloseGeometry();



    // Draw the geometry
    box->Draw("ogl");
}