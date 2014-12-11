#ifndef ElectronRecalibSuperClusterAssociator_h
#define ElectronRecalibSuperClusterAssociator_h
  
//
// Package:         RecoEgamma/EgammaElectronProducers
// Class:           ElectronRecalibSuperClusterAssociator
// 
// Description:   
  
  
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>

class PixelMatchElectronAlgo;

class ElectronRecalibSuperClusterAssociator : public edm::EDProducer
{
 public:

  explicit ElectronRecalibSuperClusterAssociator(const edm::ParameterSet& conf);

  virtual ~ElectronRecalibSuperClusterAssociator();

  virtual void produce(edm::Event& e, const edm::EventSetup& c);

 private:

  edm::InputTag superClusterCollectionEB_;
  edm::InputTag superClusterCollectionEE_;

  std::string outputLabel_;
   
  edm::InputTag electronSrc_;

};
#endif
