// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file SecondaryVertexHFML.h
/// \brief Definitions of tables for Secondary Vertex ML training
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN

#ifndef O2_ANALYSIS_SECONDARYVERTEXHFML_H_
#define O2_ANALYSIS_SECONDARYVERTEXHFML_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace vertexhfml
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
// track params before secondary verrtex propagation
DECLARE_SOA_COLUMN(XTrack1, xTrack1, float);
DECLARE_SOA_COLUMN(AlphaTrack1, alphaTrack1, float);
DECLARE_SOA_COLUMN(YTrack1, yTrack1, float);
DECLARE_SOA_COLUMN(ZTrack1, zTrack1, float);
DECLARE_SOA_COLUMN(SnpTrack1, snpTrack1, float);
DECLARE_SOA_COLUMN(TglTrack1, tglTrack1, float);
DECLARE_SOA_COLUMN(Signed1PtTrack1, signed1PtTrack1, float);
DECLARE_SOA_COLUMN(PhiTrack1, phiTrack1, float);
DECLARE_SOA_COLUMN(EtaTrack1, etaTrack1, float);
DECLARE_SOA_COLUMN(SigmaYTrack1, sigmaYTrack1, float);
DECLARE_SOA_COLUMN(SigmaZTrack1, sigmaZTrack1, float);
DECLARE_SOA_COLUMN(SigmaSnpTrack1, sigmaSnpTrack1, float);
DECLARE_SOA_COLUMN(SigmaTglTrack1, sigmaTglTrack1, float);
DECLARE_SOA_COLUMN(SigmaSigned1PtTrack1, sigmaSigned1PtTrack1, float);
DECLARE_SOA_COLUMN(RhoZYTrack1, rhoZYTrack1, int8_t);
DECLARE_SOA_COLUMN(RhoSnpYTrack1, rhoSnpYTrack1, int8_t);
DECLARE_SOA_COLUMN(RhoSnpZTrack1, rhoSnpZTrack1, int8_t);
DECLARE_SOA_COLUMN(RhoTglYTrack1, rhoTglYTrack1, int8_t);
DECLARE_SOA_COLUMN(RhoTglZTrack1, rhoTglZTrack1, int8_t);
DECLARE_SOA_COLUMN(RhoTglSnpTrack1, rhoTglSnpTrack1, int8_t);
DECLARE_SOA_COLUMN(Rho1PtYTrack1, rho1PtYTrack1, int8_t);
DECLARE_SOA_COLUMN(Rho1PtZTrack1, rho1PtZTrack1, int8_t);
DECLARE_SOA_COLUMN(Rho1PtSnpTrack1, rho1PtSnpTrack1, int8_t);
DECLARE_SOA_COLUMN(Rho1PtTglTrack1, rho1PtTglTrack1, int8_t);
DECLARE_SOA_COLUMN(XTrack2, xTrack2, float);
DECLARE_SOA_COLUMN(AlphaTrack2, alphaTrack2, float);
DECLARE_SOA_COLUMN(YTrack2, yTrack2, float);
DECLARE_SOA_COLUMN(ZTrack2, zTrack2, float);
DECLARE_SOA_COLUMN(SnpTrack2, snpTrack2, float);
DECLARE_SOA_COLUMN(TglTrack2, tglTrack2, float);
DECLARE_SOA_COLUMN(Signed1PtTrack2, signed1PtTrack2, float);
DECLARE_SOA_COLUMN(PhiTrack2, phiTrack2, float);
DECLARE_SOA_COLUMN(EtaTrack2, etaTrack2, float);
DECLARE_SOA_COLUMN(SigmaYTrack2, sigmaYTrack2, float);
DECLARE_SOA_COLUMN(SigmaZTrack2, sigmaZTrack2, float);
DECLARE_SOA_COLUMN(SigmaSnpTrack2, sigmaSnpTrack2, float);
DECLARE_SOA_COLUMN(SigmaTglTrack2, sigmaTglTrack2, float);
DECLARE_SOA_COLUMN(SigmaSigned1PtTrack2, sigmaSigned1PtTrack2, float);
DECLARE_SOA_COLUMN(RhoZYTrack2, rhoZYTrack2, int8_t);
DECLARE_SOA_COLUMN(RhoSnpYTrack2, rhoSnpYTrack2, int8_t);
DECLARE_SOA_COLUMN(RhoSnpZTrack2, rhoSnpZTrack2, int8_t);
DECLARE_SOA_COLUMN(RhoTglYTrack2, rhoTglYTrack2, int8_t);
DECLARE_SOA_COLUMN(RhoTglZTrack2, rhoTglZTrack2, int8_t);
DECLARE_SOA_COLUMN(RhoTglSnpTrack2, rhoTglSnpTrack2, int8_t);
DECLARE_SOA_COLUMN(Rho1PtYTrack2, rho1PtYTrack2, int8_t);
DECLARE_SOA_COLUMN(Rho1PtZTrack2, rho1PtZTrack2, int8_t);
DECLARE_SOA_COLUMN(Rho1PtSnpTrack2, rho1PtSnpTrack2, int8_t);
DECLARE_SOA_COLUMN(Rho1PtTglTrack2, rho1PtTglTrack2, int8_t);

// secondary vertex
DECLARE_SOA_COLUMN(XSecondaryVertex, xSecondaryVertex, float);
DECLARE_SOA_COLUMN(YSecondaryVertex, ySecondaryVertex, float);
DECLARE_SOA_COLUMN(ZSecondaryVertex, zSecondaryVertex, float);
DECLARE_SOA_COLUMN(ErrorDecayLength, errorDecayLength, float);
DECLARE_SOA_COLUMN(Chi2PCA, chi2PCA, float);

// tracks after propagation
DECLARE_SOA_COLUMN(YDCA1, yDCA1, float);
DECLARE_SOA_COLUMN(ZDCA1, zDCA1, float);
DECLARE_SOA_COLUMN(SigmaYDCA1, sigmaYDCA1, float);
DECLARE_SOA_COLUMN(SigmaZDCA1, sigmaZDCA1, float);
DECLARE_SOA_COLUMN(RHOYZDCA1, rhoYZDCA1, float);
DECLARE_SOA_COLUMN(YDCA2, yDCA2, float);
DECLARE_SOA_COLUMN(ZDCA2, zDCA2, float);
DECLARE_SOA_COLUMN(SigmaYDCA2, sigmaYDCA2, float);
DECLARE_SOA_COLUMN(SigmaZDCA2, sigmaZDCA2, float);
DECLARE_SOA_COLUMN(RHOYZDCA2, rhoYZDCA2, float);

// tracks after propagation to secondary vertex
DECLARE_SOA_COLUMN(XTrackPropag1, xTrackPropag1, float);
DECLARE_SOA_COLUMN(AlphaTrackPropag1, alphaTrackPropag1, float);
DECLARE_SOA_COLUMN(YTrackPropag1, yTrackPropag1, float);
DECLARE_SOA_COLUMN(ZTrackPropag1, zTrackPropag1, float);
DECLARE_SOA_COLUMN(SnpTrackPropag1, snpTrackPropag1, float);
DECLARE_SOA_COLUMN(TglTrackPropag1, tglTrackPropag1, float);
DECLARE_SOA_COLUMN(Signed1PtTrackPropag1, signed1PtTrackPropag1, float);
DECLARE_SOA_COLUMN(PhiTrackPropag1, phiTrackPropag1, float);
DECLARE_SOA_COLUMN(EtaTrackPropag1, etaTrackPropag1, float);
DECLARE_SOA_COLUMN(SigmaYTrackPropag1, sigmaYTrackPropag1, float);
DECLARE_SOA_COLUMN(SigmaZTrackPropag1, sigmaZTrackPropag1, float);
DECLARE_SOA_COLUMN(SigmaSnpTrackPropag1, sigmaSnpTrackPropag1, float);
DECLARE_SOA_COLUMN(SigmaTglTrackPropag1, sigmaTglTrackPropag1, float);
DECLARE_SOA_COLUMN(SigmaSigned1PtTrackPropag1, sigmaSigned1PtTrackPropag1, float);
DECLARE_SOA_COLUMN(RhoZYTrackPropag1, rhoZYTrackPropag1, int8_t);
DECLARE_SOA_COLUMN(RhoSnpYTrackPropag1, rhoSnpYTrackPropag1, int8_t);
DECLARE_SOA_COLUMN(RhoSnpZTrackPropag1, rhoSnpZTrackPropag1, int8_t);
DECLARE_SOA_COLUMN(RhoTglYTrackPropag1, rhoTglYTrackPropag1, int8_t);
DECLARE_SOA_COLUMN(RhoTglZTrackPropag1, rhoTglZTrackPropag1, int8_t);
DECLARE_SOA_COLUMN(RhoTglSnpTrackPropag1, rhoTglSnpTrackPropag1, int8_t);
DECLARE_SOA_COLUMN(Rho1PtYTrackPropag1, rho1PtYTrackPropag1, int8_t);
DECLARE_SOA_COLUMN(Rho1PtZTrackPropag1, rho1PtZTrackPropag1, int8_t);
DECLARE_SOA_COLUMN(Rho1PtSnpTrackPropag1, rho1PtSnpTrackPropag1, int8_t);
DECLARE_SOA_COLUMN(Rho1PtTglTrackPropag1, rho1PtTglTrackPropag1, int8_t);
DECLARE_SOA_COLUMN(XTrackPropag2, xTrackPropag2, float);
DECLARE_SOA_COLUMN(AlphaTrackPropag2, alphaTrackPropag2, float);
DECLARE_SOA_COLUMN(YTrackPropag2, yTrackPropag2, float);
DECLARE_SOA_COLUMN(ZTrackPropag2, zTrackPropag2, float);
DECLARE_SOA_COLUMN(SnpTrackPropag2, snpTrackPropag2, float);
DECLARE_SOA_COLUMN(TglTrackPropag2, tglTrackPropag2, float);
DECLARE_SOA_COLUMN(Signed1PtTrackPropag2, signed1PtTrackPropag2, float);
DECLARE_SOA_COLUMN(PhiTrackPropag2, phiTrackPropag2, float);
DECLARE_SOA_COLUMN(EtaTrackPropag2, etaTrackPropag2, float);
DECLARE_SOA_COLUMN(SigmaYTrackPropag2, sigmaYTrackPropag2, float);
DECLARE_SOA_COLUMN(SigmaZTrackPropag2, sigmaZTrackPropag2, float);
DECLARE_SOA_COLUMN(SigmaSnpTrackPropag2, sigmaSnpTrackPropag2, float);
DECLARE_SOA_COLUMN(SigmaTglTrackPropag2, sigmaTglTrackPropag2, float);
DECLARE_SOA_COLUMN(SigmaSigned1PtTrackPropag2, sigmaSigned1PtTrackPropag2, float);
DECLARE_SOA_COLUMN(RhoZYTrackPropag2, rhoZYTrackPropag2, int8_t);
DECLARE_SOA_COLUMN(RhoSnpYTrackPropag2, rhoSnpYTrackPropag2, int8_t);
DECLARE_SOA_COLUMN(RhoSnpZTrackPropag2, rhoSnpZTrackPropag2, int8_t);
DECLARE_SOA_COLUMN(RhoTglYTrackPropag2, rhoTglYTrackPropag2, int8_t);
DECLARE_SOA_COLUMN(RhoTglZTrackPropag2, rhoTglZTrackPropag2, int8_t);
DECLARE_SOA_COLUMN(RhoTglSnpTrackPropag2, rhoTglSnpTrackPropag2, int8_t);
DECLARE_SOA_COLUMN(Rho1PtYTrackPropag2, rho1PtYTrackPropag2, int8_t);
DECLARE_SOA_COLUMN(Rho1PtZTrackPropag2, rho1PtZTrackPropag2, int8_t);
DECLARE_SOA_COLUMN(Rho1PtSnpTrackPropag2, rho1PtSnpTrackPropag2, int8_t);
DECLARE_SOA_COLUMN(Rho1PtTglTrackPropag2, rho1PtTglTrackPropag2, int8_t);

} // namespace vertexhfml
DECLARE_SOA_TABLE(VertexHFML, "AOD", "VERTEXHFML", vertexhfml::CollisionId,
                  vertexhfml::XTrack1, vertexhfml::AlphaTrack1,
                  vertexhfml::YTrack1, vertexhfml::ZTrack1, vertexhfml::SnpTrack1, vertexhfml::TglTrack1, vertexhfml::Signed1PtTrack1,
                  vertexhfml::PhiTrack1, vertexhfml::EtaTrack1,
                  vertexhfml::SigmaYTrack1, vertexhfml::SigmaZTrack1, vertexhfml::SigmaSnpTrack1, vertexhfml::SigmaTglTrack1, vertexhfml::SigmaSigned1PtTrack1,
                  vertexhfml::RhoZYTrack1,
                  vertexhfml::RhoSnpYTrack1, vertexhfml::RhoSnpZTrack1,
                  vertexhfml::RhoTglYTrack1, vertexhfml::RhoTglZTrack1, vertexhfml::RhoTglSnpTrack1,
                  vertexhfml::Rho1PtYTrack1, vertexhfml::Rho1PtZTrack1, vertexhfml::Rho1PtSnpTrack1, vertexhfml::Rho1PtTglTrack1,
                  vertexhfml::XTrack2, vertexhfml::AlphaTrack2,
                  vertexhfml::YTrack2, vertexhfml::ZTrack2, vertexhfml::SnpTrack2, vertexhfml::TglTrack2, vertexhfml::Signed1PtTrack2,
                  vertexhfml::PhiTrack2, vertexhfml::EtaTrack2,
                  vertexhfml::SigmaYTrack2, vertexhfml::SigmaZTrack2, vertexhfml::SigmaSnpTrack2, vertexhfml::SigmaTglTrack2, vertexhfml::SigmaSigned1PtTrack2,
                  vertexhfml::RhoZYTrack2,
                  vertexhfml::RhoSnpYTrack2, vertexhfml::RhoSnpZTrack2,
                  vertexhfml::RhoTglYTrack2, vertexhfml::RhoTglZTrack2, vertexhfml::RhoTglSnpTrack2,
                  vertexhfml::Rho1PtYTrack2, vertexhfml::Rho1PtZTrack2, vertexhfml::Rho1PtSnpTrack2, vertexhfml::Rho1PtTglTrack2,
                  vertexhfml::XSecondaryVertex, vertexhfml::YSecondaryVertex, vertexhfml::ZSecondaryVertex,
                  vertexhfml::ErrorDecayLength, vertexhfml::Chi2PCA,
                  vertexhfml::YDCA1, vertexhfml::ZDCA1,
                  vertexhfml::SigmaYDCA1, vertexhfml::SigmaZDCA1,
                  vertexhfml::RHOYZDCA1,
                  vertexhfml::YDCA2, vertexhfml::ZDCA2,
                  vertexhfml::SigmaYDCA2,
                  vertexhfml::SigmaZDCA2, vertexhfml::RHOYZDCA2,
                  vertexhfml::XTrackPropag1, vertexhfml::AlphaTrackPropag1,
                  vertexhfml::YTrackPropag1, vertexhfml::ZTrackPropag1, vertexhfml::SnpTrackPropag1, vertexhfml::TglTrackPropag1, vertexhfml::Signed1PtTrackPropag1,
                  vertexhfml::PhiTrackPropag1, vertexhfml::EtaTrackPropag1,
                  //vertexhfml::SigmaYTrackPropag1,vertexhfml::SigmaZTrackPropag1,vertexhfml::SigmaSnpTrackPropag1,vertexhfml::SigmaTglTrackPropag1,vertexhfml::SigmaSigned1PtTrackPropag1,
                  //vertexhfml::RhoZYTrackPropag1,
                  //vertexhfml::RhoSnpYTrackPropag1,vertexhfml::RhoSnpZTrackPropag1,
                  //vertexhfml::RhoTglYTrackPropag1,vertexhfml::RhoTglZTrackPropag1,vertexhfml::RhoTglSnpTrackPropag1,
                  //vertexhfml::Rho1PtYTrackPropag1,vertexhfml::Rho1PtZTrackPropag1,vertexhfml::Rho1PtSnpTrackPropag1,vertexhfml::Rho1PtTglTrackPropag1,
                  vertexhfml::XTrackPropag2, vertexhfml::AlphaTrackPropag2,
                  vertexhfml::YTrackPropag2, vertexhfml::ZTrackPropag2, vertexhfml::SnpTrackPropag2, vertexhfml::TglTrackPropag2, vertexhfml::Signed1PtTrackPropag2,
                  vertexhfml::PhiTrackPropag2, vertexhfml::EtaTrackPropag2 //,
                  //vertexhfml::SigmaYTrackPropag2,vertexhfml::SigmaZTrackPropag2,vertexhfml::SigmaSnpTrackPropag2,vertexhfml::SigmaTglTrackPropag2,vertexhfml::SigmaSigned1PtTrackPropag2,
                  //vertexhfml::RhoZYTrackPropag2,
                  //vertexhfml::RhoSnpYTrackPropag2,vertexhfml::RhoSnpZTrackPropag2,
                  //vertexhfml::RhoTglYTrackPropag2,vertexhfml::RhoTglZTrackPropag2,vertexhfml::RhoTglSnpTrackPropag2,
                  //vertexhfml::Rho1PtYTrackPropag2,vertexhfml::Rho1PtZTrackPropag2,vertexhfml::Rho1PtSnpTrackPropag2,vertexhfml::Rho1PtTglTrackPropag2
);
} // namespace o2::aod

#endif // O2_ANALYSIS_SECONDARYVERTEXHFML_H_
