// -*- C++ -*-
//
//

// user include files
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "Calibration/EcalCalibAlgos/interface/ElectronRecalibSuperClusterAssociator.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCoreFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include <iostream>

#include "DataFormats/Math/interface/deltaR.h"
#define DEBUG 

using namespace reco;
using namespace edm;
 
ElectronRecalibSuperClusterAssociator::ElectronRecalibSuperClusterAssociator(const edm::ParameterSet& iConfig) 
{
#ifdef DEBUG
  std::cout<< "ElectronRecalibSuperClusterAssociator::ElectronRecalibSuperClusterAssociator" << std::endl;
#endif

  //register your products
  produces<GsfElectronCollection>();
  produces<GsfElectronCoreCollection>() ;
  //  produces<SuperClusterCollection>();
  
  superClusterCollectionEB_ = iConfig.getParameter<edm::InputTag > ("superClusterCollectionEB");
  superClusterCollectionEE_ = iConfig.getParameter<edm::InputTag > ("superClusterCollectionEE");

  outputLabel_ = iConfig.getParameter<std::string>("outputLabel");
  electronSrc_ = iConfig.getParameter<edm::InputTag > ("electronSrc");

#ifdef DEBUG
  std::cout<< "ElectronRecalibSuperClusterAssociator::ElectronRecalibSuperClusterAssociator::end" << std::endl;
#endif
}

ElectronRecalibSuperClusterAssociator::~ElectronRecalibSuperClusterAssociator()
{
}

// ------------ method called to produce the data  ------------
void ElectronRecalibSuperClusterAssociator::produce(edm::Event& e, const edm::EventSetup& iSetup) 
{
  // Create the output collections   
  std::auto_ptr<GsfElectronCollection> pOutEle(new GsfElectronCollection);
  std::auto_ptr<GsfElectronCoreCollection> pOutEleCore(new GsfElectronCoreCollection);
  //  std::auto_ptr<SuperClusterCollection> pOutNewEndcapSC(new SuperClusterCollection);

  //Get Hybrid SuperClusters
  Handle<reco::SuperClusterCollection> superClusterEBHandle;
  e.getByLabel(superClusterCollectionEB_, superClusterEBHandle);
  if (!superClusterEBHandle.isValid()) {
    std::cerr << "Error! can't get the product SuperClusterCollection "<< std::endl;
  }
  const reco::SuperClusterCollection* scCollection = superClusterEBHandle.product();
  
#ifdef DEBUG
  std::cout<<"scCollection->size()"<<scCollection->size()<<std::endl;
#endif
  
  Handle<reco::SuperClusterCollection> superClusterEEHandle;
  e.getByLabel(superClusterCollectionEE_, superClusterEEHandle);
  if (!superClusterEEHandle.isValid()) {
    std::cerr << "Error! can't get the product IslandSuperClusterCollection "<< std::endl;
  }
  const reco::SuperClusterCollection* scIslandCollection = superClusterEEHandle.product();
  
#ifdef DEBUG
  std::cout<<"scEECollection->size()"<<scIslandCollection->size()<<std::endl;
#endif

  // Get Electrons
  Handle<edm::View<reco::GsfElectron> > pElectrons;
  e.getByLabel(electronSrc_, pElectrons);
  if (!pElectrons.isValid()) {
    std::cerr << "Error! can't get the product ElectronCollection "<< std::endl;
  }
  const edm::View<reco::GsfElectron>* electronCollection = pElectrons.product();

  GsfElectronCoreRefProd rEleCore=e.getRefBeforePut<GsfElectronCoreCollection>();
  edm::Ref<GsfElectronCoreCollection>::key_type idxEleCore = 0;
  
  for(edm::View<reco::GsfElectron>::const_iterator eleIt = electronCollection->begin(); eleIt != electronCollection->end(); eleIt++)
    {
      float DeltaRMineleSCbarrel(0.15); //initial minDeltaR
      float DeltaRMineleSCendcap(0.15); 
      const reco::SuperCluster* nearestSCbarrel=0;
      const reco::SuperCluster* nearestSCendcap=0;
      int iscRef=-1, iscRefendcap=-1;
      int iSC=0;


#ifdef DEBUG
      std::cout << "+++++++++++" << std::endl;
      std::cout << &(*eleIt->gsfTrack()) << std::endl;
      std::cout << eleIt->core()->ecalDrivenSeed() << std::endl;
      std::cout << "+++++++++++" << std::endl;
#endif
      
      // first loop is on EB superClusters
      iSC=0;
      for(reco::SuperClusterCollection::const_iterator scIt = scCollection->begin();
	  scIt != scCollection->end(); scIt++, iSC++){

	double DeltaReleSC = sqrt(reco::deltaR2(eleIt->eta(), eleIt->phi(),
						scIt->eta(), scIt->phi()));

	if(DeltaReleSC<DeltaRMineleSCbarrel)
	  {
	    DeltaRMineleSCbarrel = DeltaReleSC;
	    nearestSCbarrel = &*scIt;
	    iscRef = iSC;
	  }
#ifdef DEBUG	
	std::cout << "EB: " << scIt - scCollection->begin() << " " << iSC << " " << iscRef 
		  << "\t" << std::setprecision(4) << scIt->energy() 
		  << " " << scIt->eta() << " " << scIt->phi() 
		  << " " << eleIt->eta() << " " << eleIt->phi() 
		  << "\t" << DeltaRMineleSCbarrel
		  << std::endl;
#endif
      }
      
      // second loop is on EE superClusters
      iSC=0;
      for(reco::SuperClusterCollection::const_iterator scIt = scIslandCollection->begin();
	  scIt != scIslandCollection->end(); scIt++, iSC++){
#ifdef DEBUG	
	std::cout << "EE: " << scIt - scCollection->begin() << " " << iSC << " " << iscRef 
		  << "\t" << std::setprecision(4) << scIt->energy() 
		  << " " << scIt->eta() << " " << scIt->phi() 
		  << " " << eleIt->eta() << " " << eleIt->phi() 
		  << "\t" << DeltaRMineleSCbarrel 
		  << std::endl;
#endif
	
	double DeltaReleSC = sqrt(reco::deltaR2(eleIt->eta(), eleIt->phi(),
						scIt->eta(), scIt->phi()));

	if(DeltaReleSC<DeltaRMineleSCendcap)
	  {
	    DeltaRMineleSCendcap = DeltaReleSC;
	    nearestSCendcap = &*scIt;
	    iscRefendcap = iSC;
	  }
      }
      ////////////////////////      
      //      if(eleIt->isEB()) assert(DeltaRMineleSCbarrel < DeltaRMineleSCendcap);
      //else assert(DeltaRMineleSCbarrel > DeltaRMineleSCendcap);
      if(eleIt->isEB() && DeltaRMineleSCbarrel > DeltaRMineleSCendcap) continue;
      
      if(eleIt->isEB() && nearestSCbarrel){
#ifdef DEBUG
	std::cout << "Starting Association is with EB superCluster "<< std::endl;
#endif  
 	reco::GsfElectronCore newEleCore(*(eleIt->core()));
#ifdef DEBUG
	std::cout << "newEleCore created "<< std::endl;
#endif  
	
	newEleCore.setGsfTrack(eleIt->gsfTrack());
	reco::SuperClusterRef scRef(reco::SuperClusterRef(superClusterEBHandle, iscRef));
	newEleCore.setSuperCluster(scRef);
#ifdef DEBUG
	std::cout << "SC ref associated: " << scRef->energy() << std::endl;
	std::cout << "Creating newEleCoreRef with idxEleCore: " << idxEleCore << std::endl;
#endif
	reco::GsfElectronCoreRef newEleCoreRef(reco::GsfElectronCoreRef(rEleCore, idxEleCore++));
#ifdef DEBUG
	std::cout << "Created " << std::endl;
#endif
	pOutEleCore->push_back(newEleCore);
#ifdef DEBUG
	
	std::cout << "pushed core: " 
	  //<< newEleCoreRef->ecalDrivenSeed() 
		  << std::endl;
#endif
	//reco::GsfElectron newEle(*eleIt,newEleCoreRef);
	reco::GsfElectron newEle(*eleIt,eleIt->core());
	//	newEle.setCorrectedEcalEnergy(eleIt->p4().energy()*(nearestSCbarrel->energy()/eleIt->ecalEnergy()));
	//newEle.setCorrectedEcalEnergyError(eleIt->ecalEnergyError()*(nearestSCbarrel->energy()/eleIt->ecalEnergy()));
#ifdef DEBUG
	std::cout << "new Electron created" << std::endl;
#endif
	pOutEle->push_back(newEle);
#ifdef DEBUG
	std::cout << "Association is with EB superCluster "<< std::endl;
#endif  
      }  

      if(!(eleIt->isEB()) && nearestSCendcap)
	{
#ifdef DEBUG
	std::cout << "Starting Association is with EE superCluster "<< std::endl;
#endif  

// 	float preshowerEnergy=eleIt->superCluster()->preshowerEnergy(); 
// #ifdef DEBUG
// 	  std::cout << "preshowerEnergy"<< preshowerEnergy << std::endl;
// #endif
// 	  /// fixme : should have a vector of ptr of ref, to avoid copying
// 	  CaloClusterPtrVector newBCRef;
// 	  for (CaloCluster_iterator bcRefIt=nearestSCendcap->clustersBegin();bcRefIt!=nearestSCendcap->clustersEnd();++bcRefIt){
// 	    CaloClusterPtr cPtr(*bcRefIt);
// 	    newBCRef.push_back(cPtr);
// 	  }
	 

// 	  //	  reco::SuperCluster newSC(nearestSCendcap->energy() + preshowerEnergy, nearestSCendcap->position() , nearestSCendcap->seed(),newBCRef , preshowerEnergy );
// 	  reco::SuperCluster newSC(nearestSCendcap->energy(), nearestSCendcap->position() , nearestSCendcap->seed(),newBCRef , preshowerEnergy );
// 	  pOutNewEndcapSC->push_back(newSC);
// 	  reco::SuperClusterRef scRef(reco::SuperClusterRef(rSC, idxSC ++));
	  
	  reco::GsfElectronCore newEleCore(*(eleIt->core()));
	  newEleCore.setGsfTrack(eleIt->gsfTrack());
	  reco::SuperClusterRef scRef(reco::SuperClusterRef(superClusterEBHandle, iscRefendcap));
	  newEleCore.setSuperCluster(scRef);
	  reco::GsfElectronCoreRef newEleCoreRef(reco::GsfElectronCoreRef(rEleCore, idxEleCore ++));
	  pOutEleCore->push_back(newEleCore);
	  reco::GsfElectron newEle(*eleIt,newEleCoreRef);
	  //,CaloClusterPtr(),
	  //				  TrackRef(),GsfTrackRefVector());
	  //TrackRef(),TrackBaseRef(), GsfTrackRefVector());
	  newEle.setCorrectedEcalEnergy(eleIt->p4().energy()*(nearestSCendcap->energy()/eleIt->ecalEnergy()));
	  newEle.setCorrectedEcalEnergyError(eleIt->ecalEnergyError()*(nearestSCendcap->energy()/eleIt->ecalEnergy()));
	//	std::cout << "FROM REF " << newEle.superCluster().key() << std::endl;
	  pOutEle->push_back(newEle);
// 	  reco::GsfElectronCore newEleCore(*(eleIt->core()));
// 	  newEleCore.setGsfTrack(eleIt->gsfTrack());
// 	  newEleCore.setSuperCluster(scRef);
// 	  reco::GsfElectronCoreRef newEleCoreRef(reco::GsfElectronCoreRef(rEleCore, idxEleCore ++));
// 	  pOutEleCore->push_back(newEleCore);
// 	  reco::GsfElectron newEle(*eleIt,newEleCoreRef,CaloClusterPtr(),
// //				  TrackRef(),GsfTrackRefVector());
// 				  TrackRef(),TrackBaseRef(), GsfTrackRefVector());
// 	  newEle.setCorrectedEcalEnergy(eleIt->p4().energy()*(newSC.energy()/eleIt->ecalEnergy()),eleIt->ecalEnergyError()*(newSC.energy()/eleIt->ecalEnergy())); 
// 	  pOutEle->push_back(newEle);

#ifdef DEBUG
	std::cout << "Association is with EE superCluster "<< std::endl;
#endif  
      }  
    
    }
  
  
  
#ifdef DEBUG
  std::cout << "Filled new electrons  " << pOutEle->size() << std::endl;
  std::cout << "Filled new electronsCore  " << pOutEleCore->size() << std::endl;
  //  std::cout << "Filled new endcapSC  " << pOutNewEndcapSC->size() << std::endl;
#endif  
  
  // put result into the Event

  e.put(pOutEle);
  e.put(pOutEleCore);
  //  e.put(pOutNewEndcapSC);
  
}

DEFINE_FWK_MODULE(ElectronRecalibSuperClusterAssociator);
