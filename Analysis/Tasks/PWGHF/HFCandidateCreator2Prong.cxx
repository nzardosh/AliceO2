// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file hfcandidatecreator2prong.cxx
/// \brief Reconstruction of heavy-flavour 2-prong decay candidates
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "Analysis/SecondaryVertexHF.h"
#include "Analysis/SecondaryVertexHFML.h"
#include "Analysis/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Analysis/RecoDecay.h"
#include "PID/PIDResponse.h"
#include <cmath>
#include <array>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using std::array;

/// Reconstruction of heavy-flavour 2-prong decay candidates
struct HFCandidateCreator2Prong {
  Produces<aod::HfCandBase> rowCandidateBase;
  Produces<aod::VertexHFML> secondaryVertexHFML;
  //Produces<aod::HfCandProng2Base> rowCandidateProng2Base; // TODO split table
  Configurable<double> magneticField{"d_bz", 5.0, "magnetic field"};
  Configurable<bool> b_propdca{"b_propdca", true, "create tracks version propagated to PCA"};
  Configurable<double> d_maxr{"d_maxr", 200., "reject PCA's above this radius"};
  Configurable<double> d_maxdzini{"d_maxdzini", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> d_minparamchange{"d_minparamchange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> d_minrelchi2change{"d_minrelchi2change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<bool> b_dovalplots{"b_dovalplots", true, "do validation plots"};
  Configurable<bool> b_doSecondaryVertexHFML{"b_doSecondaryVertexHFML", false, "make secondary vertex hf ml table"};
  OutputObj<TH1F> hmass2{TH1F("hmass2", "2-track inv mass", 500, 0., 5.0)};

  double massPi = RecoDecay::getMassPDG(kPiPlus);
  double massK = RecoDecay::getMassPDG(kKPlus);
  double massPiK{0};
  double massKPi{0};

  void process(aod::Collision const& collision,
               aod::HfTrackIndexProng2 const& rowsTrackIndexProng2,
               aod::BigTracks const& tracks)
  {
    // 2-prong vertex fitter
    o2::vertexing::DCAFitterN<2> df;
    df.setBz(magneticField);
    df.setPropagateToPCA(b_propdca);
    df.setMaxR(d_maxr);
    df.setMaxDZIni(d_maxdzini);
    df.setMinParamChange(d_minparamchange);
    df.setMinRelChi2Change(d_minrelchi2change);
    df.setUseAbsDCA(true);

    // loop over pairs of track indeces
    for (const auto& rowTrackIndexProng2 : rowsTrackIndexProng2) {
      auto trackParVarPos1 = getTrackParCov(rowTrackIndexProng2.index0());
      auto trackParVarNeg1 = getTrackParCov(rowTrackIndexProng2.index1());

      // reconstruct the 2-prong secondary vertex
      if (df.process(trackParVarPos1, trackParVarNeg1) == 0)
        continue;
      const auto& secondaryVertex = df.getPCACandidate();
      auto chi2PCA = df.getChi2AtPCACandidate();
      auto covMatrixPCA = df.calcPCACovMatrix().Array();
      auto trackParVar0 = df.getTrack(0);
      auto trackParVar1 = df.getTrack(1);

      // get track momenta
      array<float, 3> pvec0;
      array<float, 3> pvec1;
      trackParVar0.getPxPyPzGlo(pvec0);
      trackParVar1.getPxPyPzGlo(pvec1);

      // calculate invariant masses
      auto arrayMomenta = array{pvec0, pvec1};
      massPiK = RecoDecay::M(arrayMomenta, array{massPi, massK});
      massKPi = RecoDecay::M(arrayMomenta, array{massK, massPi});

      // get track impact parameters
      // This modifies track momenta!
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      o2::dataformats::DCA impactParameter0;
      o2::dataformats::DCA impactParameter1;
      trackParVar0.propagateToDCA(primaryVertex, magneticField, &impactParameter0);
      trackParVar1.propagateToDCA(primaryVertex, magneticField, &impactParameter1);

      // get uncertainty of the decay length
      double phi, theta;
      getPointDirection(array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertex, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

      // fill candidate table rows
      rowCandidateBase(collision.posX(), collision.posY(), collision.posZ(),
                       secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                       errorDecayLength, errorDecayLengthXY,
                       chi2PCA, //);
                                //rowCandidateProng2Base( // TODO split table
                       pvec0[0], pvec0[1], pvec0[2],
                       pvec1[0], pvec1[1], pvec1[2],
                       impactParameter0.getY(), impactParameter1.getY(),
                       std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()),
                       rowTrackIndexProng2.index0Id(), rowTrackIndexProng2.index1Id());
      if (b_doSecondaryVertexHFML) {
        auto track1 = rowTrackIndexProng2.index0();
        auto track2 = rowTrackIndexProng2.index1();
        //auto trackParVarCov1 = getTrackParCov(trackParVar0);
        //auto trackParVarCov2 = getTrackParCov(trackParVar1);
        secondaryVertexHFML(track1.collisionId(),
                            track1.x(), track1.alpha(),
                            track1.y(), track1.z(), track1.snp(), track1.tgl(), track1.signed1Pt(),
                            track1.phi(), track1.eta(),
                            std::sqrt(trackParVarPos1.getSigmaY2()), std::sqrt(trackParVarPos1.getSigmaZ2()), std::sqrt(trackParVarPos1.getSigmaSnp2()),
                            std::sqrt(trackParVarPos1.getSigmaTgl2()), std::sqrt(trackParVarPos1.getSigma1Pt2()),
                            trackParVarPos1.getSigmaZY(),
                            trackParVarPos1.getSigmaSnpY(), trackParVarPos1.getSigmaSnpZ(),
                            trackParVarPos1.getSigmaTglY(), trackParVarPos1.getSigmaTglZ(), trackParVarPos1.getSigmaTglSnp(),
                            trackParVarPos1.getSigma1PtY(), trackParVarPos1.getSigma1PtZ(), trackParVarPos1.getSigma1PtSnp(), trackParVarPos1.getSigma1PtTgl(),
                            track2.x(), track2.alpha(),
                            track2.y(), track2.z(), track2.snp(), track2.tgl(), track2.signed1Pt(),
                            track2.phi(), track2.eta(),
                            std::sqrt(trackParVarNeg1.getSigmaY2()), std::sqrt(trackParVarNeg1.getSigmaZ2()), std::sqrt(trackParVarNeg1.getSigmaSnp2()),
                            std::sqrt(trackParVarNeg1.getSigmaTgl2()), std::sqrt(trackParVarNeg1.getSigma1Pt2()),
                            trackParVarNeg1.getSigmaZY(),
                            trackParVarNeg1.getSigmaSnpY(), trackParVarNeg1.getSigmaSnpZ(),
                            trackParVarNeg1.getSigmaTglY(), trackParVarNeg1.getSigmaTglZ(), trackParVarNeg1.getSigmaTglSnp(),
                            trackParVarNeg1.getSigma1PtY(), trackParVarNeg1.getSigma1PtZ(), trackParVarNeg1.getSigma1PtSnp(), trackParVarNeg1.getSigma1PtTgl(),
                            secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                            errorDecayLength, chi2PCA,
                            impactParameter0.getY(), impactParameter0.getZ(),
                            std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter0.getSigmaZ2()),
                            impactParameter0.getSigmaYZ(),
                            impactParameter1.getY(), impactParameter1.getZ(),
                            std::sqrt(impactParameter1.getSigmaY2()), std::sqrt(impactParameter1.getSigmaZ2()),
                            impactParameter1.getSigmaYZ(),
                            trackParVar0.getX(), trackParVar0.getAlpha(),
                            trackParVar0.getY(), trackParVar0.getZ(), trackParVar0.getSnp(), trackParVar0.getTgl(), trackParVar0.getQ2Pt(), //is this correct?
                            trackParVar0.getPhi(), trackParVar0.getEta(),
                            //std::sqrt(trackParVarCov1.getSigmaY2()),std::sqrt(trackParVarCov1.getSigmaZ2()),std::sqrt(trackParVarCov1.getSigmaSnp2()),
                            //std::sqrt(trackParVarCov1.getSigmaTgl2()),std::sqrt(trackParVarCov1.getSigma1Pt2()),
                            //trackParVarCov1.getSigmaZY(),
                            //trackParVarCov1.getSigmaSnpY(),trackParVarCov1.getSigmaSnpZ(),
                            //trackParVarCov1.getSigmaTglY(),trackParVarCov1.getSigmaTglZ(),trackParVarCov1.getSigmaTglSnp(),
                            //trackParVarCov1.getSigma1PtY(), trackParVarCov1.getSigma1PtZ(), trackParVarCov1.getSigma1PtSnp(),trackParVarCov1.getSigma1PtTgl(),
                            trackParVar1.getX(), trackParVar1.getAlpha(),
                            trackParVar1.getY(), trackParVar1.getZ(), trackParVar1.getSnp(), trackParVar1.getTgl(), trackParVar1.getQ2Pt(), //is this correct?
                            trackParVar1.getPhi(), trackParVar1.getEta()                                                                    //,
                            //std::sqrt(trackParVarCov2.getSigmaY2()),std::sqrt(trackParVarCov2.getSigmaZ2()),std::sqrt(trackParVarCov2.getSigmaSnp2()),
                            //std::sqrt(trackParVarCov2.getSigmaTgl2()),std::sqrt(trackParVarCov2.getSigma1Pt2()),
                            //trackParVarCov2.getSigmaZY(),
                            //trackParVarCov2.getSigmaSnpY(),trackParVarCov2.getSigmaSnpZ(),
                            //trackParVarCov2.getSigmaTglY(),trackParVarCov2.getSigmaTglZ(),trackParVarCov2.getSigmaTglSnp(),
                            //trackParVarCov2.getSigma1PtY(), trackParVarCov2.getSigma1PtZ(), trackParVarCov2.getSigma1PtSnp(),trackParVarCov2.getSigma1PtTgl()
        );
      }
      // fill histograms
      if (b_dovalplots) {
        hmass2->Fill(massPiK);
        hmass2->Fill(massKPi);
      }
    }
  }
};

/// Extends the base table with expression columns.
struct HFCandidateCreator2ProngExpressions {
  Spawns<aod::HfCandProng2Ext> rowCandidateProng2;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    adaptAnalysisTask<HFCandidateCreator2Prong>("hf-cand-creator-2prong"),
    adaptAnalysisTask<HFCandidateCreator2ProngExpressions>("hf-cand-creator-2prong-expressions")};
}
