#include "TGeoManager.h"
#include "TGeoTube.h"
#include "TGeoBBox.h"
#include "TMath.h"

// Elements
TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98, 13, 2.7);

void MakeLayer(TGeoVolume *top, double radius, double thickness, double length) {
    // Define the media
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum", 1, matVacuum);
    TGeoMedium *Al = new TGeoMedium("Aluminium", 2, matAl);

    // Create the layer
    double width = 4.54 * 2;

    int nStaves = (int) (radius * 2 * 3.14 / width + 1.5);
    double overlap = width - 2 * radius * std::sin(TMath::Pi() / nStaves);

    double theta = 2 * TMath::Pi() / nStaves;
    double theta1 = std::atan(width / 2 / radius);
    double st = std::sin(theta);
    double ct = std::cos(theta);
    double theta2 = std::atan((radius * st - width / 2 * ct)/(radius * ct + width / 2 * st));
    double overlap2 = (theta1 - theta2) * radius;

    printf("--> Creating layer ar r=%.2f cm with %d staves, each %.2f cm wide. Overlap = %.3f  ov2 = %.3f\n", radius, nStaves, width, overlap, overlap2);

    int counter = 0;
    for (int iStave = 0; iStave < nStaves; iStave++) {
        TGeoBBox* stave = new TGeoBBox(Form("stave_X_%d", iStave), width / 2, thickness/2, length/2);
        TGeoVolume* staveVol = new TGeoVolume("volTube", stave, Al);

        // Make transformation as a composition of a rotation, translation and another rotation
        TGeoCombiTrans *trans = new TGeoCombiTrans();

        double theta = 360. * iStave / nStaves;
        TGeoRotation *rot2 = new TGeoRotation("rot2", theta + 90 + 2, 0, 0);
        trans->SetRotation(rot2);

        // translation
        trans->SetTranslation(radius * std::cos(2. * TMath::Pi() * iStave / nStaves), radius * std::sin(2 * TMath::Pi() * iStave / nStaves), 0);

        top->AddNode(staveVol, counter++, trans);
    }
}

void test_geom() {
    TGeoManager *geom = new TGeoManager("Assemblies", "Geometry using assemblies");

    // Rotate the axes
    TGeoRotation *rot = new TGeoRotation("rot", 90, 0, 0);
    TGeoCombiTrans *trans = new TGeoCombiTrans();
    trans->SetRotation(rot);

    // Make top volume
    TGeoVolume *top = geom->MakeBox("top", nullptr, 1000, 1000, 1000);
    geom->SetTopVolume(top);

    MakeLayer(top, 30, 0.03, 264);
    MakeLayer(top, 45.5, 0.03, 264);
    MakeLayer(top, 60, 0.03, 264);
    MakeLayer(top, 80, 0.03, 264);

    //--- close the geometry
    geom->CloseGeometry();


    geom->SetVisLevel(1);
    geom->SetVisOption(2);
    geom->SetClipping(true);
    // top->Raytrace();
    // gGeoManager->GetGeomPainter()->SetDrawOption("x");

    // geom->SetDrawOption("x");
    top->Draw("ogl");

    // TVirtualViewer3D *view = (TVirtualViewer3D*) gPad->GetViewer3D();
    // view->ShowAxis();
}
