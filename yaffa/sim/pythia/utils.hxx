/*
Collection of different functions used in pythia simulations.
*/

#include <set>

#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "Pythia8/Pythia.h"

/*
Return true if the particle satisfies all the following criteria:
- is a pion, kaon or proton
- is a final-state particle (status code > 0)
- is within the TPC acceptance (|eta| < 0.8)
- has pT > 0.3 GeV/c
*/
// inline bool IsDetectableInTPC(Pythia8::Particle p, std::set<int> pdgs = {211, 321, 2212}) {
//     return pdgs.find(std::abs(p.id())) != pdgs.end() &&
//            p->GetStatusCode() > 0 &&
//            p->Pt() > 0.3 &&
//            std::abs(p->Eta()) < 0.8;
// }


bool IsPrimary(const Pythia8::Event &event, const int &iPart) {
    Pythia8::Particle part = event[iPart];

    int m1 = part.mother1();
    int m2 = part.mother2();

    // Based on the definition in Pythia, see https://pythia.org/latest-manual/ParticleProperties.html
    return m1 > 0 && m2 > 0 && m1 < m2 && 81 <= std::abs(part.status()) && std::abs(part.status() <= 86);
}



bool IsInAcc(const Pythia8::Event &event, const int &iPart) {
    Pythia8::Particle p = event[iPart];

    int abspdg = std::abs(p.id());
    if (abspdg == 211 || abspdg == 2212 || abspdg == 321) { // pion, proton and kaon
        return std::abs(p.eta()) < 0.8;
    } else if (abspdg == 3122) { // Lambda
        Pythia8::Particle d1 = p.daughter1();
        Pythia8::Particle d2 = p.daughter2();

        return IsInAcc(event, p.daughter1()) && IsInAcc(event, p.daughter1());
    } else {
        printf("Error: Acceptance check for pdg=%d is not implemented. Exit!\n", p.id());
        // exit(1);
    }
    return false;
}


float RelativePairMomentum(const TLorentzVector &PartOne, const TLorentzVector &PartTwo) {
  TLorentzVector trackSum = PartOne + PartTwo;

  float beta = trackSum.Beta();
  float betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
  float betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
  float betaz = beta * cos(trackSum.Theta());

  TLorentzVector PartOneCMS = PartOne;
  TLorentzVector PartTwoCMS = PartTwo;

  PartOneCMS.Boost(-betax, -betay, -betaz);
  PartTwoCMS.Boost(-betax, -betay, -betaz);

  TLorentzVector trackRelK = PartOneCMS - PartTwoCMS;

  return 0.5 * trackRelK.P();
}

float RelativePairMomentum(const Pythia8::Particle &part1, const Pythia8::Particle &part2) {
    TLorentzVector p1 = TLorentzVector(part1.px(), part1.py(), part1.pz(), part1.e());
    TLorentzVector p2 = TLorentzVector(part2.px(), part2.py(), part2.pz(), part2.e());

    return RelativePairMomentum(p1, p2);
}


std::set<int> SelectInDaughterTree(const Pythia8::Event &event, const int &start, const std::vector<int> &pdgs) {
    std::set<int> result = {};

    if (start < 3) return result;

    Pythia8::Particle part = event[start];
    if (std::find(pdgs.begin(), pdgs.end(), std::abs(part.id())) != pdgs.end()) {
        result.insert(start);
    } else {
        for (int iDau = part.daughter1(); iDau <= part.daughter2(); iDau++) {
            auto r = SelectInDaughterTree(event, iDau, pdgs);
            result.insert(r.begin(), r.end());
        }
    }
    return result;
}
