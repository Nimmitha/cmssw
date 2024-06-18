// #ifndef _miniAODtrigRate_h
// #define _miniAODtrigRate_h

// // system include files
// #include <memory>

// // user include files
// #include "trigFireRateAnalyzer/miniAODtrigRate/plugins/miniAODtrigRate.h"
// #include "FWCore/Framework/interface/Frameworkfwd.h"
// #include "FWCore/Framework/interface/one/EDAnalyzer.h"

// #include "FWCore/Framework/interface/Event.h"
// #include "FWCore/Framework/interface/MakerMacros.h"

// #include "FWCore/ParameterSet/interface/ParameterSet.h"
// #include "FWCore/Utilities/interface/InputTag.h"

// #include "DataFormats/PatCandidates/interface/Muon.h"

// // class declaration
// //
// // If the analyzer does not use TFileService, please remove
// // the template argument to the base class so the class inherits
// // from  edm::one::EDAnalyzer<>
// // This will improve performance in multithreaded jobs.

// class miniAODtrigRate : public edm::one::EDAnalyzer<edm::one::SharedResources> {
// public:
//   explicit miniAODtrigRate(const edm::ParameterSet&);
//   ~miniAODtrigRate() override;

//   static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

// private:
//   void beginJob() override;
//   void analyze(const edm::Event&, const edm::EventSetup&) override;
//   void endJob() override;

//   // ----------member data ---------------------------
//   edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
//   edm::EDGetTokenT<edm::TriggerResults> triggerBits_;


// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//   edm::ESGetToken<SetupData, SetupRecord> setupToken_;
// #endif
// };

// #endif