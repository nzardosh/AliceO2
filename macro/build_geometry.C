// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "TGeoManager.h"
#include "TString.h"
#include "TSystem.h"

#include "DetectorsPassive/Cave.h"
#include "DetectorsPassive/Magnet.h"
#include "DetectorsPassive/Dipole.h"
#include "DetectorsPassive/Compensator.h"
#include "DetectorsPassive/Absorber.h"
#include "DetectorsPassive/Shil.h"
#include "DetectorsPassive/Hall.h"
#include "DetectorsPassive/Pipe.h"
#include <Field/MagneticField.h>
#include <TPCSimulation/Detector.h>
#include <ITSSimulation/Detector.h>
#include <MFTSimulation/Detector.h>
#include <MCHSimulation/Detector.h>
#include <MIDSimulation/Detector.h>
#include <EMCALSimulation/Detector.h>
#include <TOFSimulation/Detector.h>
#include <TRDSimulation/Detector.h>
#include <FT0Simulation/Detector.h>
#include <FV0Simulation/Detector.h>
#include <FDDSimulation/Detector.h>
#include <HMPIDSimulation/Detector.h>
#include <PHOSSimulation/Detector.h>
#include <CPVSimulation/Detector.h>
#include <ZDCSimulation/Detector.h>
#include <DetectorsPassive/Cave.h>
#include <DetectorsPassive/FrameStructure.h>
#include <SimConfig/SimConfig.h>
#include "FairRunSim.h"
#include <FairLogger.h>
#include <algorithm>
#endif

#ifdef ENABLE_UPGRADES
#include <ITS3Simulation/Detector.h>
#endif

void finalize_geometry(FairRunSim* run);

bool isActivated(std::string s)
{
  // access user configuration for list of wanted modules
  auto& modulelist = o2::conf::SimConfig::Instance().getActiveDetectors();
  auto active = std::find(modulelist.begin(), modulelist.end(), s) != modulelist.end();
  if (active) {
    LOG(INFO) << "Activating " << s << " module";
  }
  return active;
}

// a "factory" like macro to instantiate the O2 geometry
void build_geometry(FairRunSim* run = nullptr)
{
  bool geomonly = (run == nullptr);

  // minimal macro to test setup of the geometry
  auto& confref = o2::conf::SimConfig::Instance();

  TString dir = getenv("VMCWORKDIR");
  TString geom_dir = dir + "/Detectors/Geometry/";
  gSystem->Setenv("GEOMPATH", geom_dir.Data());

  TString tut_configdir = dir + "/Detectors/gconfig";
  gSystem->Setenv("CONFIG_DIR", tut_configdir.Data());

  // Create simulation run if it does not exist
  if (run == nullptr) {
    run = new FairRunSim();
    run->SetOutputFile("foo.root"); // Output file
    run->SetName("TGeant3");        // Transport engine
  }
  // Create media
  run->SetMaterials("media.geo"); // Materials

  // we need a field to properly init the media
  int fld = confref.getConfigData().mField, fldAbs = std::abs(fld);
  float fldCoeff;
  o2::field::MagFieldParam::BMap_t fldType;
  switch (fldAbs) {
    case 5:
      fldType = o2::field::MagFieldParam::k5kG;
      fldCoeff = fld > 0 ? 1. : -1;
      break;
    case 0:
      fldType = o2::field::MagFieldParam::k5kG;
      fldCoeff = 0;
      break;
    case 2:
      fldType = o2::field::MagFieldParam::k2kG;
      fldCoeff = fld > 0 ? 1. : -1;
      break;
    default:
      LOG(FATAL) << "Field option " << fld << " is not supported, use +-2, +-5 or 0";
  };

  auto field = new o2::field::MagneticField("Maps", "Maps", fldCoeff, fldCoeff, fldType);
  run->SetField(field);

  // Create geometry
  // we always need the cave
  o2::passive::Cave* cave = new o2::passive::Cave("CAVE");
  // adjust size depending on content
  cave->includeZDC(isActivated("ZDC"));
  // the experiment hall (cave)
  cave->SetGeometryFileName("cave.geo");
  run->AddModule(cave);

  // the experimental hall
  if (isActivated("HALL")) {
    auto hall = new o2::passive::Hall("HALL", "Experimental Hall");
    run->AddModule(hall);
  }

  // the magnet
  if (isActivated("MAG")) {
    // the frame structure to support other detectors
    auto magnet = new o2::passive::Magnet("MAG", "L3 Magnet");
    run->AddModule(magnet);
  }

  // the dipole
  if (isActivated("DIPO")) {
    auto dipole = new o2::passive::Dipole("DIPO", "Alice Dipole");
    run->AddModule(dipole);
  }

  // the compensator dipole
  if (isActivated("COMP")) {
    run->AddModule(new o2::passive::Compensator("COMP", "Alice Compensator Dipole"));
  }

  // beam pipe
  if (isActivated("PIPE")) {
#ifdef ENABLE_UPGRADES
    run->AddModule(new o2::passive::Pipe("PIPE", "Beam pipe", 1.6f, 0.05));
#else
    run->AddModule(new o2::passive::Pipe("PIPE", "Beam pipe"));
#endif
  }

  // the absorber
  if (isActivated("ABSO")) {
    // the frame structure to support other detectors
    auto abso = new o2::passive::Absorber("ABSO", "Absorber");
    run->AddModule(abso);
  }

  // the shil
  if (isActivated("SHIL")) {
    auto shil = new o2::passive::Shil("SHIL", "Small angle beam shield");
    run->AddModule(shil);
  }

  if (isActivated("TOF") || isActivated("TRD") || isActivated("FRAME")) {
    // the frame structure to support other detectors
    auto frame = new o2::passive::FrameStructure("FRAME", "Frame");
    run->AddModule(frame);
  }

  if (isActivated("TOF")) {
    // TOF
    auto tof = new o2::tof::Detector(true);
    run->AddModule(tof);
  }

  if (isActivated("TRD")) {
    // TRD
    auto trd = new o2::trd::Detector(true);
    run->AddModule(trd);
  }

  if (isActivated("TPC")) {
    // tpc
    auto tpc = new o2::tpc::Detector(true);
    run->AddModule(tpc);
  }
#ifdef ENABLE_UPGRADES
  if (isActivated("IT3")) {
    // ITS3
    auto its3 = new o2::its3::Detector(kTRUE);
    run->AddModule(its3);
  }
#endif

  if (isActivated("ITS")) {
    // its
    auto its = new o2::its::Detector(true);
    run->AddModule(its);
  }

  if (isActivated("MFT")) {
    // mft
    auto mft = new o2::mft::Detector();
    run->AddModule(mft);
  }

  if (isActivated("MCH")) {
    // mch
    run->AddModule(new o2::mch::Detector(true));
  }

  if (isActivated("MID")) {
    // mid
    run->AddModule(new o2::mid::Detector(true));
  }

  if (isActivated("EMC")) {
    // emcal
    run->AddModule(new o2::emcal::Detector(true));
  }

  if (isActivated("PHS")) {
    // phos
    run->AddModule(new o2::phos::Detector(true));
  }

  if (isActivated("CPV")) {
    // cpv
    run->AddModule(new o2::cpv::Detector(true));
  }

  if (isActivated("FT0")) {
    // FIT-T0
    run->AddModule(new o2::ft0::Detector(true));
  }

  if (isActivated("FV0")) {
    // FIT-V0
    run->AddModule(new o2::fv0::Detector(true));
  }

  if (isActivated("FDD")) {
    // FIT-FDD
    run->AddModule(new o2::fdd::Detector(true));
  }

  if (isActivated("HMP")) {
    // HMP
    run->AddModule(new o2::hmpid::Detector(true));
  }

  if (isActivated("ZDC")) {
    // ZDC
    run->AddModule(new o2::zdc::Detector(true));
  }

  if (geomonly) {
    run->Init();
  }
}
