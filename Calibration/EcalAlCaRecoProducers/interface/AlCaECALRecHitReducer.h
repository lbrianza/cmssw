#ifndef _ALCAECALRECHITREDUCER_H
#define _ALCAECALRECHITREDUCER_H

// -*- C++ -*-
//
// Package:    AlCaECALRecHitReducer
// Class:      AlCaECALRecHitReducer
// 
/**\class AlCaECALRecHitReducer AlCaECALRecHitReducer.cc Calibration/EcalAlCaRecoProducers/src/AlCaECALRecHitReducer.cc

 Description: Example of a producer of AlCa electrons

 Implementation:
     <Notes on implementation>

*/
//
// Original Author:  Lorenzo AGOSTINO
//         Created:  Mon Jul 17 18:07:01 CEST 2006
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//!
//! class declaration
//!

class AlCaECALRecHitReducer : public edm::EDProducer {
 public:
  explicit AlCaECALRecHitReducer(const edm::ParameterSet&);
  ~AlCaECALRecHitReducer();
  
  //! producer
  virtual void produce(edm::Event &, const edm::EventSetup&);

   private:
      // ----------member data ---------------------------

  
  edm::InputTag ebRecHitsLabel_;
  edm::InputTag eeRecHitsLabel_;
  edm::InputTag esRecHitsLabel_;
  edm::InputTag electronLabel_;
  std::string alcaBarrelHitsCollection_;
  std::string alcaEndcapHitsCollection_;
  std::string alcaPreshowerHitsCollection_;
  int etaSize_;
  int phiSize_;
  float weight_;
  int esNstrips_;
  int esNcolumns_;

  bool selectByEleNum_;
  int minEleNumber_;
  double minElePt_;

};

#endif
