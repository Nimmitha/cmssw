//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 15 17:42:38 2020 by ROOT version 6.08/06
// from TTree ntuple/ntuple
// found on file: Rootuple_MC_2017-MiniAOD_ZY1Y1_All-5K.root
//////////////////////////////////////////////////////////

#ifndef AnalysisWithIso_h
#define AnalysisWithIso_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class AnalysisWithIso {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          nB;
   vector<float>   *Run;
   vector<float>   *LumiBlock;
   vector<float>   *Event;
   vector<float>   *FourL_mass;
   vector<float>   *FourL_px;
   vector<float>   *FourL_py;
   vector<float>   *FourL_pz;
   vector<float>   *FourL_pt;
   vector<float>   *FourL_eta;
   vector<float>   *FourL_phi;
   vector<float>   *FourL_VtxProb;
   vector<float>   *FourL_PVx;
   vector<float>   *FourL_PVy;
   vector<float>   *FourL_PVz;
   vector<float>   *FourL_PVxError;
   vector<float>   *FourL_PVyError;
   vector<float>   *FourL_PVzError;
   vector<bool>    *Y_TriggerPath;
   vector<float>   *Y_TriggerPt1;
   vector<float>   *Y_TriggerEta1;
   vector<float>   *Y_TriggerPhi1;
   vector<float>   *Y_TriggerPt2;
   vector<float>   *Y_TriggerEta2;
   vector<float>   *Y_TriggerPhi2;
   vector<float>   *Y_TriggerPt3;
   vector<float>   *Y_TriggerEta3;
   vector<float>   *Y_TriggerPhi3;
   vector<float>   *Y_TriggerPt4;
   vector<float>   *Y_TriggerEta4;
   vector<float>   *Y_TriggerPhi4;
   vector<float>   *Y_TriggerPt5;
   vector<float>   *Y_TriggerEta5;
   vector<float>   *Y_TriggerPhi5;
   vector<float>   *Y_mass;
   vector<float>   *Y_px;
   vector<float>   *Y_py;
   vector<float>   *Y_pz;
   vector<float>   *Y_pt;
   vector<float>   *Y_eta;
   vector<float>   *Y_phi;
   vector<float>   *Y_rapidity;

   vector<float>   *Y_eta1;
   vector<float>   *Y_pt1;
   vector<float>   *Y_eta2;
   vector<float>   *Y_pt2;   

   vector<float>   *Z_mass;
   vector<float>   *Z_px;
   vector<float>   *Z_py;
   vector<float>   *Z_pz;
   vector<float>   *Z_pt;
   vector<float>   *Z_eta;
   vector<float>   *Z_phi;
   vector<float>   *Z_rapidity;

   vector<float>   *Z_pt1;
   vector<float>   *Z_eta1;
   vector<float>   *Z_pt2;
   vector<float>   *Z_eta2;

   vector<float>   *B_J3_mass;
   vector<float>   *B_J3_px;
   vector<float>   *B_J3_py;
   vector<float>   *B_J3_pz;
   vector<float>   *B_J3_pt;
   vector<float>   *B_J3_eta;
   vector<float>   *B_J3_phi;
   vector<float>   *B_J3_rapidity;
   vector<float>   *B_J4_mass;
   vector<float>   *B_J4_px;
   vector<float>   *B_J4_py;
   vector<float>   *B_J4_pz;
   vector<float>   *B_J4_pt;
   vector<float>   *B_J4_eta;
   vector<float>   *B_J4_phi;
   
   vector<float>   *B_J4_rapidity;
   vector<float>   *B_J1_VtxPx;
   vector<float>   *B_J1_VtxPy;
   vector<float>   *B_J1_VtxPz;
   vector<float>   *B_J1_VtxPt;
   vector<float>   *B_J1_VtxEta;
   vector<float>   *B_J1_VtxPhi;
   vector<float>   *B_J1_VtxRapidity;
   vector<float>   *B_J1_VtxMass;
   vector<float>   *B_J1_PVx;
   vector<float>   *B_J1_PVy;
   vector<float>   *B_J1_PVz;
   vector<float>   *B_J1_PVxError;
   vector<float>   *B_J1_PVyError;
   vector<float>   *B_J1_PVzError;
   vector<float>   *B_J2_VtxPx;
   vector<float>   *B_J2_VtxPy;
   vector<float>   *B_J2_VtxPz;
   vector<float>   *B_J2_VtxPt;
   vector<float>   *B_J2_VtxEta;
   vector<float>   *B_J2_VtxPhi;
   vector<float>   *B_J2_VtxRapidity;
   vector<float>   *B_J2_VtxMass;
   vector<float>   *B_J2_PVx;
   vector<float>   *B_J2_PVy;
   vector<float>   *B_J2_PVz;
   vector<float>   *B_J2_PVxError;
   vector<float>   *B_J2_PVyError;
   vector<float>   *B_J2_PVzError;
   vector<float>   *B_J3_VtxPx;
   vector<float>   *B_J3_VtxPy;
   vector<float>   *B_J3_VtxPz;
   vector<float>   *B_J3_VtxPt;
   vector<float>   *B_J3_VtxEta;
   vector<float>   *B_J3_VtxPhi;
   vector<float>   *B_J3_VtxRapidity;
   vector<float>   *B_J3_VtxMass;
   vector<float>   *B_J3_PVx;
   vector<float>   *B_J3_PVy;
   vector<float>   *B_J3_PVz;
   vector<float>   *B_J3_PVxError;
   vector<float>   *B_J3_PVyError;
   vector<float>   *B_J3_PVzError;
   vector<float>   *B_J4_VtxPx;
   vector<float>   *B_J4_VtxPy;
   vector<float>   *B_J4_VtxPz;
   vector<float>   *B_J4_VtxPt;
   vector<float>   *B_J4_VtxEta;
   vector<float>   *B_J4_VtxPhi;
   vector<float>   *B_J4_VtxRapidity;
   vector<float>   *B_J4_VtxMass;
   vector<float>   *B_J4_PVx;
   vector<float>   *B_J4_PVy;
   vector<float>   *B_J4_PVz;
   vector<float>   *B_J4_PVxError;
   vector<float>   *B_J4_PVyError;
   vector<float>   *B_J4_PVzError;
   vector<float>   *B_Mu1_px;
   vector<float>   *B_Mu1_py;
   vector<float>   *B_Mu1_pz;
   vector<float>   *B_Mu1_pt;
   vector<float>   *B_Mu1_eta;
   vector<float>   *B_Mu1_phi;
   vector<int>     *B_Mu1_charge;
   vector<bool>    *B_Mu1_soft;
   vector<bool>    *B_Mu1_tight;
   vector<bool>    *B_Mu1_loose;
   vector<float>   *B_Mu1_IsoTrack;

//   vector<float>   *Y_trackiso1;

   vector<float>   *B_Mu1_IsoHcal;
   vector<float>   *B_Mu1_IsoEcal;
   vector<float>   *B_Mu1_IsoCalo;
   vector<float>   *B_Mu2_px;
   vector<float>   *B_Mu2_py;
   vector<float>   *B_Mu2_pz;
   vector<float>   *B_Mu2_pt;
   vector<float>   *B_Mu2_eta;
   vector<float>   *B_Mu2_phi;
   vector<int>     *B_Mu2_charge;
   vector<bool>    *B_Mu2_soft;
   vector<bool>    *B_Mu2_tight;
   vector<bool>    *B_Mu2_loose;
   vector<float>   *B_Mu2_IsoTrack;
   vector<float>   *B_Mu2_IsoHcal;
   vector<float>   *B_Mu2_IsoEcal;
   vector<float>   *B_Mu2_IsoCalo;
   vector<float>   *B_Mu3_px;
   vector<float>   *B_Mu3_py;
   vector<float>   *B_Mu3_pz;
   vector<float>   *B_Mu3_pt;
   vector<float>   *B_Mu3_eta;
   vector<float>   *B_Mu3_phi;
   vector<bool>    *B_Mu3_soft;
   vector<bool>    *B_Mu3_tight;
   vector<bool>    *B_Mu3_loose;
   vector<int>     *B_Mu3_charge;
   vector<float>   *B_Mu3_IsoTrack;
   vector<float>   *B_Mu3_IsoHcal;
   vector<float>   *B_Mu3_IsoEcal;
   vector<float>   *B_Mu3_IsoCalo;
   vector<float>   *B_Mu4_px;
   vector<float>   *B_Mu4_py;
   vector<float>   *B_Mu4_pz;
   vector<float>   *B_Mu4_pt;
   vector<float>   *B_Mu4_eta;
   vector<float>   *B_Mu4_phi;
   vector<bool>    *B_Mu4_soft;
   vector<bool>    *B_Mu4_tight;
   vector<bool>    *B_Mu4_loose;
   vector<int>     *B_Mu4_charge;
   vector<float>   *B_Mu4_IsoTrack;
   vector<float>   *B_Mu4_IsoHcal;
   vector<float>   *B_Mu4_IsoEcal;
   vector<float>   *B_Mu4_IsoCalo;

   vector<float>   *Y_VtxProb;
   vector<float>   *Z_VtxProb;

   vector<float>   *B_J3_VtxProb;
   vector<float>   *B_J4_VtxProb;
   vector<float>   *Z_xyP1;
   vector<float>   *Z_xyM1;
   vector<float>   *Z_zP1;
   vector<float>   *Z_zM1;
   vector<float>   *Z_xyP2;
   vector<float>   *Z_xyM2;
   vector<float>   *Z_zP2;
   vector<float>   *Z_zM2;
   vector<float>   *mu1mC2;
   vector<int>     *mu1mNHits;
   vector<int>     *mu1mNPHits;
   vector<float>   *mu1pC2;
   vector<int>     *mu1pNHits;
   vector<int>     *mu1pNPHits;
   vector<float>   *mu2mC2;
   vector<int>     *mu2mNHits;
   vector<int>     *mu2mNPHits;
   vector<float>   *mu2pC2;
   vector<int>     *mu2pNHits;
   vector<int>     *mu2pNPHits;
   vector<float>   *B_M1_pt;
   vector<float>   *B_M1_eta;
   vector<float>   *B_M1_phi;
   vector<float>   *B_M1_px;
   vector<float>   *B_M1_py;
   vector<float>   *B_M1_pz;
   vector<float>   *B_M2_pt;
   vector<float>   *B_M2_eta;
   vector<float>   *B_M2_phi;
   vector<float>   *B_M2_px;
   vector<float>   *B_M2_py;
   vector<float>   *B_M2_pz;
   vector<float>   *B_M3_pt;
   vector<float>   *B_M3_eta;
   vector<float>   *B_M3_phi;
   vector<float>   *B_M3_px;
   vector<float>   *B_M3_py;
   vector<float>   *B_M3_pz;
   vector<float>   *B_M4_pt;
   vector<float>   *B_M4_eta;
   vector<float>   *B_M4_phi;
   vector<float>   *B_M4_px;
   vector<float>   *B_M4_py;
   vector<float>   *B_M4_pz;
   vector<float>   *Z_GenMuonPt;
   vector<float>   *Z_GenMuonEta;
   vector<float>   *Z_GenMuonPhi;
   vector<float>   *Y_GenMuonPt;
   vector<float>   *Y_GenMuonEta;
   vector<float>   *Y_GenMuonPhi;

   // vector<float>   *Y_looseHEEPwp1;
   vector<float>   *Z_loose1;
   // vector<bool>   *Y_looseHZZwp1;
   // vector<bool>   *Y_looseHZZwp2;
   vector<bool>   *Y_mvaIsoWP90_1;
   vector<bool>   *Y_mvaIsoWP90_2;

   vector<bool>   *Y_looseCutBase1;
   vector<bool>   *Y_looseCutBase2;
   vector<bool>   *Z_soft1;
   vector<bool>   *Z_soft2;

   // List of branches

   // TBranch        *b_Y_looseHEEPwp1;
   TBranch        *b_Z_loose1;
   TBranch        *b_Y_mvaIsoWP90_1;
   TBranch        *b_Y_mvaIsoWP90_2;
 //  TBranch        *Y_looseHZZwp1;

   TBranch        *b_nB;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_FourL_mass;   //!
   TBranch        *b_FourL_px;   //!
   TBranch        *b_FourL_py;   //!
   TBranch        *b_FourL_pz;   //!
   TBranch        *b_FourL_pt;   //!
   TBranch        *b_FourL_eta;   //!
   TBranch        *b_FourL_phi;   //!
   TBranch        *b_FourL_VtxProb;   //!
   TBranch        *b_FourL_PVx;   //!
   TBranch        *b_FourL_PVy;   //!
   TBranch        *b_FourL_PVz;   //!
   TBranch        *b_FourL_PVxError;   //!
   TBranch        *b_FourL_PVyError;   //!
   TBranch        *b_FourL_PVzError;   //!
   TBranch        *b_Y_TriggerPath;   //!
   TBranch        *b_Y_TriggerPt1;   //!
   TBranch        *b_Y_TriggerEta1;   //!
   TBranch        *b_Y_TriggerPhi1;   //!
   TBranch        *b_Y_TriggerPt2;   //!
   TBranch        *b_Y_TriggerEta2;   //!
   TBranch        *b_Y_TriggerPhi2;   //!
   TBranch        *b_Y_TriggerPt3;   //!
   TBranch        *b_Y_TriggerEta3;   //!
   TBranch        *b_Y_TriggerPhi3;   //!
   TBranch        *b_Y_TriggerPt4;   //!
   TBranch        *b_Y_TriggerEta4;   //!
   TBranch        *b_Y_TriggerPhi4;   //!
   TBranch        *b_Y_TriggerPt5;   //!
   TBranch        *b_Y_TriggerEta5;   //!
   TBranch        *b_Y_TriggerPhi5;   //!
   TBranch        *b_Y_mass;   //!
   TBranch        *b_Y_px;   //!
   TBranch        *b_Y_py;   //!
   TBranch        *b_Y_pz;   //!
   TBranch        *b_Y_pt;   //!
   TBranch        *b_Y_eta;   //!
   TBranch        *b_Y_phi;   //!
   TBranch        *b_Y_rapidity;   //!
   TBranch        *b_Z_mass;   //!
   TBranch        *b_Z_px;   //!
   TBranch        *b_Z_py;   //!
   TBranch        *b_Z_pz;   //!
   TBranch        *b_Z_pt;   //!
   TBranch        *b_Z_eta;   //!
   TBranch        *b_Z_phi;   //!
   TBranch        *b_Z_rapidity;   //!

   TBranch        *Z_eta1_b;   //!
   TBranch        *Z_pt1_b;
   TBranch        *Z_eta2_b;   //!
   TBranch        *Z_pt2_b;
   TBranch        *Y_eta1_b;   //!
   TBranch        *Y_pt1_b;
   TBranch        *Y_eta2_b;   //!
   TBranch        *Y_pt2_b;

   // TBranch        *b_B_J3_mass;   //!
   // TBranch        *b_B_J3_px;   //!
   // TBranch        *b_B_J3_py;   //!
   // TBranch        *b_B_J3_pz;   //!
   // TBranch        *b_B_J3_pt;   //!
   // TBranch        *b_B_J3_eta;   //!
   // TBranch        *b_B_J3_phi;   //!
   // TBranch        *b_B_J3_rapidity;   //!
   // TBranch        *b_B_J4_mass;   //!
   // TBranch        *b_B_J4_px;   //!
   // TBranch        *b_B_J4_py;   //!
   // TBranch        *b_B_J4_pz;   //!
   // TBranch        *b_B_J4_pt;   //!
   // TBranch        *b_B_J4_eta;   //!
   // TBranch        *b_B_J4_phi;   //!
   // TBranch        *b_B_J4_rapidity;   //!
   // TBranch        *b_B_J1_VtxPx;   //!
   // TBranch        *b_B_J1_VtxPy;   //!
   // TBranch        *b_B_J1_VtxPz;   //!
   // TBranch        *b_B_J1_VtxPt;   //!
   // TBranch        *b_B_J1_VtxEta;   //!
   // TBranch        *b_B_J1_VtxPhi;   //!
   // TBranch        *b_B_J1_VtxRapidity;   //!
   // TBranch        *b_B_J1_VtxMass;   //!
   // TBranch        *b_B_J1_PVx;   //!
   // TBranch        *b_B_J1_PVy;   //!
   // TBranch        *b_B_J1_PVz;   //!
   // TBranch        *b_B_J1_PVxError;   //!
   // TBranch        *b_B_J1_PVyError;   //!
   // TBranch        *b_B_J1_PVzError;   //!
   // TBranch        *b_B_J2_VtxPx;   //!
   // TBranch        *b_B_J2_VtxPy;   //!
   // TBranch        *b_B_J2_VtxPz;   //!
   // TBranch        *b_B_J2_VtxPt;   //!
   // TBranch        *b_B_J2_VtxEta;   //!
   // TBranch        *b_B_J2_VtxPhi;   //!
   // TBranch        *b_B_J2_VtxRapidity;   //!
   // TBranch        *b_B_J2_VtxMass;   //!
   // TBranch        *b_B_J2_PVx;   //!
   // TBranch        *b_B_J2_PVy;   //!
   // TBranch        *b_B_J2_PVz;   //!
   // TBranch        *b_B_J2_PVxError;   //!
   // TBranch        *b_B_J2_PVyError;   //!
   // TBranch        *b_B_J2_PVzError;   //!
   // TBranch        *b_B_J3_VtxPx;   //!
   // TBranch        *b_B_J3_VtxPy;   //!
   // TBranch        *b_B_J3_VtxPz;   //!
   // TBranch        *b_B_J3_VtxPt;   //!
   // TBranch        *b_B_J3_VtxEta;   //!
   // TBranch        *b_B_J3_VtxPhi;   //!
   // TBranch        *b_B_J3_VtxRapidity;   //!
   // TBranch        *b_B_J3_VtxMass;   //!
   // TBranch        *b_B_J3_PVx;   //!
   // TBranch        *b_B_J3_PVy;   //!
   // TBranch        *b_B_J3_PVz;   //!
   // TBranch        *b_B_J3_PVxError;   //!
   // TBranch        *b_B_J3_PVyError;   //!
   // TBranch        *b_B_J3_PVzError;   //!
   // TBranch        *b_B_J4_VtxPx;   //!
   // TBranch        *b_B_J4_VtxPy;   //!
   // TBranch        *b_B_J4_VtxPz;   //!
   // TBranch        *b_B_J4_VtxPt;   //!
   // TBranch        *b_B_J4_VtxEta;   //!
   // TBranch        *b_B_J4_VtxPhi;   //!
   // TBranch        *b_B_J4_VtxRapidity;   //!
   // TBranch        *b_B_J4_VtxMass;   //!
   // TBranch        *b_B_J4_PVx;   //!
   // TBranch        *b_B_J4_PVy;   //!
   // TBranch        *b_B_J4_PVz;   //!
   // TBranch        *b_B_J4_PVxError;   //!
   // TBranch        *b_B_J4_PVyError;   //!
   // TBranch        *b_B_J4_PVzError;   //!
   // TBranch        *b_B_Mu1_px;   //!
   // TBranch        *b_B_Mu1_py;   //!
   // TBranch        *b_B_Mu1_pz;   //!
   // TBranch        *b_B_Mu1_pt;   //!
   // TBranch        *b_B_Mu1_eta;   //!
   // TBranch        *b_B_Mu1_phi;   //!
   // TBranch        *b_B_Mu1_charge;   //!
   // TBranch        *b_B_Mu1_soft;   //!
   // TBranch        *b_B_Mu1_tight;   //!
   // TBranch        *b_B_Mu1_loose;   //!
   // TBranch        *b_B_Mu1_IsoTrack;   //!

   TBranch        *b_Y_trackiso1;
   // TBranch        *b_Y_looseHZZwp1;
   // TBranch        *b_Y_looseHZZwp2;

   TBranch        *b_Y_looseCutBase1;
   TBranch        *b_Y_looseCutBase2;
   TBranch        *b_Z_soft1;
   TBranch        *b_Z_soft2;

   // TBranch        *b_B_Mu1_IsoHcal;   //!
   // TBranch        *b_B_Mu1_IsoEcal;   //!
   // TBranch        *b_B_Mu1_IsoCalo;   //!
   // TBranch        *b_B_Mu2_px;   //!
   // TBranch        *b_B_Mu2_py;   //!
   // TBranch        *b_B_Mu2_pz;   //!
   // TBranch        *b_B_Mu2_pt;   //!
   // TBranch        *b_B_Mu2_eta;   //!
   // TBranch        *b_B_Mu2_phi;   //!
   // TBranch        *b_B_Mu2_charge;   //!
   // TBranch        *b_B_Mu2_soft;   //!
   // TBranch        *b_B_Mu2_tight;   //!
   // TBranch        *b_B_Mu2_loose;   //!
   // TBranch        *b_B_Mu2_IsoTrack;   //!
   // TBranch        *b_B_Mu2_IsoHcal;   //!
   // TBranch        *b_B_Mu2_IsoEcal;   //!
   // TBranch        *b_B_Mu2_IsoCalo;   //!
   // TBranch        *b_B_Mu3_px;   //!
   // TBranch        *b_B_Mu3_py;   //!
   // TBranch        *b_B_Mu3_pz;   //!
   // TBranch        *b_B_Mu3_pt;   //!
   // TBranch        *b_B_Mu3_eta;   //!
   // TBranch        *b_B_Mu3_phi;   //!
   // TBranch        *b_B_Mu3_soft;   //!
   // TBranch        *b_B_Mu3_tight;   //!
   // TBranch        *b_B_Mu3_loose;   //!
   // TBranch        *b_B_Mu3_charge;   //!
   // TBranch        *b_B_Mu3_IsoTrack;   //!
   // TBranch        *b_B_Mu3_IsoHcal;   //!
   // TBranch        *b_B_Mu3_IsoEcal;   //!
   // TBranch        *b_B_Mu3_IsoCalo;   //!
   // TBranch        *b_B_Mu4_px;   //!
   // TBranch        *b_B_Mu4_py;   //!
   // TBranch        *b_B_Mu4_pz;   //!
   // TBranch        *b_B_Mu4_pt;   //!
   // TBranch        *b_B_Mu4_eta;   //!
   // TBranch        *b_B_Mu4_phi;   //!
   // TBranch        *b_B_Mu4_soft;   //!
   // TBranch        *b_B_Mu4_tight;   //!
   // TBranch        *b_B_Mu4_loose;   //!
   // TBranch        *b_B_Mu4_charge;   //!
   // TBranch        *b_B_Mu4_IsoTrack;   //!
   // TBranch        *b_B_Mu4_IsoHcal;   //!
   // TBranch        *b_B_Mu4_IsoEcal;   //!
   // TBranch        *b_B_Mu4_IsoCalo;   //!

   TBranch        *b_Y_VtxProb;   //!
   TBranch        *b_Z_VtxProb;   //!

   // TBranch        *b_B_J3_VtxProb;   //!
   // TBranch        *b_B_J4_VtxProb;   //!
   // TBranch        *b_Z_xyP1;   //!
   // TBranch        *b_Z_xyM1;   //!
   // TBranch        *b_Z_zP1;   //!
   // TBranch        *b_Z_zM1;   //!
   // TBranch        *b_Z_xyP2;   //!
   // TBranch        *b_Z_xyM2;   //!
   // TBranch        *b_Z_zP2;   //!
   // TBranch        *b_Z_zM2;   //!
   // TBranch        *b_mu1mC2;   //!
   // TBranch        *b_mu1mNHits;   //!
   // TBranch        *b_mu1mNPHits;   //!
   // TBranch        *b_mu1pC2;   //!
   // TBranch        *b_mu1pNHits;   //!
   // TBranch        *b_mu1pNPHits;   //!
   // TBranch        *b_mu2mC2;   //!
   // TBranch        *b_mu2mNHits;   //!
   // TBranch        *b_mu2mNPHits;   //!
   // TBranch        *b_mu2pC2;   //!
   // TBranch        *b_mu2pNHits;   //!
   // TBranch        *b_mu2pNPHits;   //!
   // TBranch        *b_B_M1_pt;   //!
   // TBranch        *b_B_M1_eta;   //!
   // TBranch        *b_B_M1_phi;   //!
   // TBranch        *b_B_M1_px;   //!
   // TBranch        *b_B_M1_py;   //!
   // TBranch        *b_B_M1_pz;   //!
   // TBranch        *b_B_M2_pt;   //!
   // TBranch        *b_B_M2_eta;   //!
   // TBranch        *b_B_M2_phi;   //!
   // TBranch        *b_B_M2_px;   //!
   // TBranch        *b_B_M2_py;   //!
   // TBranch        *b_B_M2_pz;   //!
   // TBranch        *b_B_M3_pt;   //!
   // TBranch        *b_B_M3_eta;   //!
   // TBranch        *b_B_M3_phi;   //!
   // TBranch        *b_B_M3_px;   //!
   // TBranch        *b_B_M3_py;   //!
   // TBranch        *b_B_M3_pz;   //!
   // TBranch        *b_B_M4_pt;   //!
   // TBranch        *b_B_M4_eta;   //!
   // TBranch        *b_B_M4_phi;   //!
   // TBranch        *b_B_M4_px;   //!
   // TBranch        *b_B_M4_py;   //!
   // TBranch        *b_B_M4_pz;   //!
   // TBranch        *b_Z_GenMuonPt;   //!
   // TBranch        *b_Z_GenMuonEta;   //!
   // TBranch        *b_Z_GenMuonPhi;   //!
   // TBranch        *b_Y_GenMuonPt;   //!
   // TBranch        *b_Y_GenMuonEta;   //!
   // TBranch        *b_Y_GenMuonPhi;   //!

   AnalysisWithIso(TTree *tree=0);
   virtual ~AnalysisWithIso();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnalysisWithIso_cxx
AnalysisWithIso::AnalysisWithIso(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/fourmuonBothJY/SingleElectronEgamma_Run2_v2.root");
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/ZeeYuu/ZeeYuu_10k_Apr23Preselection.root");//UL MC 10k
     
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/Isolationredo/Y2toY1_Zee_2p5k.root"); 
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/Isolationredo/Y3toY1_Zee_3k.root"); 


     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/ZeeYuu/ZeeYuuUL2018_MCtestUL_1000_Apr23_addedHZZID.root");
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/singleElectronMar23/SingleElectron_Run2UL_2point8to12.root");//UL SingleElectron run2
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/singleElectronMar23/SingleElectron_Run2UL_2point8to12_Triggers.root");

     TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject("ZmmMCOut.root");

     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/singleElectronMar23/SingleElectron-UL2017_Data_2point8to12_all.root");
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/singleElectronMar23/SingleElectron-UL2017_Data_2point8to12_TriggerL1.root");

     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/singleElectronMar23/SingleElectron-UL2016_Data_2point8to12_most.root");
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/singleElectronMar23/Egamma_2018UL_march2023Apr_Data_2point8to12_Triggers.root");

     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/singleElectronMar23/SingleElectron_Run2UL_2point8to12_TriggerOffical.root");
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/singleElectronMar23/SingleElectron_Run2UL_2point8to12_L1.root");

     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/ZeeYuu/SimulationHZeeYuu-2018_s12.root");// MC 1k
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/ZeeYuu/ZeeYuuUL2018_MCtestUL_1000_Apr23_tirgger27.root");// UL MC 1k trigger 27
     // TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("SingleElectron_Run2017B-UL2017_take2_Data_UltraLegacy.root"); // MC
      if (!f || !f->IsOpen()) {
         //f = new TFile("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/fourmuonBothJY/SingleElectronEgamma_Run2_v2.root");
       //f = new TFile("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/ZeeYuu/ZeeYuuUL2018_MCtestUL_1000_Apr23.root"); //UL MC 1k
       //f = new TFile("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/ZeeYuu/ZeeYuu_10k_Apr23Preselection.root"); //UL MC 10k
       
       //f = new TFile("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/Isolationredo/Y2toY1_Zee_2p5k.root");
       //f = new TFile("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/Isolationredo/Y3toY1_Zee_3k.root");


       //f = new TFile("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/ZeeYuu/ZeeYuuUL2018_MCtestUL_1000_Apr23_addedHZZID.root");
       //f = new TFile("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/singleElectronMar23/SingleElectron_Run2UL_2point8to12.root"); //UL SingleElectron run2
       //f = new TFile("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/singleElectronMar23/SingleElectron_Run2UL_2point8to12_Triggers.root");

         f = new TFile("ZmmMCOut.root");

         //f = new TFile("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/singleElectronMar23/SingleElectron-UL2017_Data_2point8to12_all.root");
         //f = new TFile("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/singleElectronMar23/SingleElectron-UL2017_Data_2point8to12_TriggerL1.root");

         //f = new TFile("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/singleElectronMar23/SingleElectron-UL2016_Data_2point8to12_most.root");
         //f = new TFile("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/singleElectronMar23/Egamma_2018UL_march2023Apr_Data_2point8to12_Triggers.root");

         //f = new TFile("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/singleElectronMar23/Egamma_2018UL_march2023Apr_Data_2point8to12_all.root");
         //f = new TFile("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/singleElectronMar23/SingleElectron_Run2UL_2point8to12_L1.root");

         //f = new TFile("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/ZeeYuu/SimulationHZeeYuu-2018_s12.root"); // MC 1k
         //f = new TFile("~/../../../mnt/d/Documents/PhD/Research/workingDirectory/ZeeYuu/ZeeYuuUL2018_MCtestUL_1000_Apr23_tirgger27.root"); // UL MC 1k trigger 27
         // f = new TFile("~/looseCuts/singleElectron/SimulationHZeeYuu-2018-10k.root"); // MC
          
      }
      f->GetObject("ntuple",tree);

   }
   Init(tree);
}

AnalysisWithIso::~AnalysisWithIso()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalysisWithIso::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalysisWithIso::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AnalysisWithIso::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Run = 0;
   LumiBlock = 0;
   Event = 0;
   FourL_mass = 0;
   FourL_px = 0;
   FourL_py = 0;
   FourL_pz = 0;
   FourL_pt = 0;
   FourL_eta = 0;
   FourL_phi = 0;
   FourL_VtxProb = 0;
   FourL_PVx = 0;
   FourL_PVy = 0;
   FourL_PVz = 0;
   FourL_PVxError = 0;
   FourL_PVyError = 0;
   FourL_PVzError = 0;
   Y_TriggerPath = 0;
   Y_TriggerPt1 = 0;
   Y_TriggerEta1 = 0;
   Y_TriggerPhi1 = 0;
   Y_TriggerPt2 = 0;
   Y_TriggerEta2 = 0;
   Y_TriggerPhi2 = 0;
   Y_TriggerPt3 = 0;
   Y_TriggerEta3 = 0;
   Y_TriggerPhi3 = 0;
   Y_TriggerPt4 = 0;
   Y_TriggerEta4 = 0;
   Y_TriggerPhi4 = 0;
   Y_TriggerPt5 = 0;
   Y_TriggerEta5 = 0;
   Y_TriggerPhi5 = 0;
   Y_mass = 0;
   Y_px = 0;
   Y_py = 0;
   Y_pz = 0;
   Y_pt = 0;
   Y_eta = 0;
   Y_phi = 0;
   Y_rapidity = 0;
   Z_mass = 0;
   Z_px = 0;
   Z_py = 0;
   Z_pz = 0;
   Z_pt = 0;
   Z_eta = 0;
   Z_phi = 0;
   Z_rapidity = 0;

   Z_eta1 = 0;
   Z_pt1 = 0;
   Z_eta2 = 0;
   Z_pt2 = 0;
   Y_pt1 = 0;
   Y_eta1 = 0;
   Y_pt2 = 0;
   Y_eta2 = 0;

   Z_mass = 0;
   B_J3_px = 0;
   B_J3_py = 0;
   B_J3_pz = 0;
   B_J3_pt = 0;
   B_J3_eta = 0;
   B_J3_phi = 0;
   B_J3_rapidity = 0;
   B_J4_mass = 0;
   B_J4_px = 0;
   B_J4_py = 0;
   B_J4_pz = 0;
   B_J4_pt = 0;
   B_J4_eta = 0;
   B_J4_phi = 0;
   B_J4_rapidity = 0;
   B_J1_VtxPx = 0;
   B_J1_VtxPy = 0;
   B_J1_VtxPz = 0;
   B_J1_VtxPt = 0;
   B_J1_VtxEta = 0;
   B_J1_VtxPhi = 0;
   B_J1_VtxRapidity = 0;
   B_J1_VtxMass = 0;
   B_J1_PVx = 0;
   B_J1_PVy = 0;
   B_J1_PVz = 0;
   B_J1_PVxError = 0;
   B_J1_PVyError = 0;
   B_J1_PVzError = 0;
   B_J2_VtxPx = 0;
   B_J2_VtxPy = 0;
   B_J2_VtxPz = 0;
   B_J2_VtxPt = 0;
   B_J2_VtxEta = 0;
   B_J2_VtxPhi = 0;
   B_J2_VtxRapidity = 0;
   B_J2_VtxMass = 0;
   B_J2_PVx = 0;
   B_J2_PVy = 0;
   B_J2_PVz = 0;
   B_J2_PVxError = 0;
   B_J2_PVyError = 0;
   B_J2_PVzError = 0;
   B_J3_VtxPx = 0;
   B_J3_VtxPy = 0;
   B_J3_VtxPz = 0;
   B_J3_VtxPt = 0;
   B_J3_VtxEta = 0;
   B_J3_VtxPhi = 0;
   B_J3_VtxRapidity = 0;
   B_J3_VtxMass = 0;
   B_J3_PVx = 0;
   B_J3_PVy = 0;
   B_J3_PVz = 0;
   B_J3_PVxError = 0;
   B_J3_PVyError = 0;
   B_J3_PVzError = 0;
   B_J4_VtxPx = 0;
   B_J4_VtxPy = 0;
   B_J4_VtxPz = 0;
   B_J4_VtxPt = 0;
   B_J4_VtxEta = 0;
   B_J4_VtxPhi = 0;
   B_J4_VtxRapidity = 0;
   B_J4_VtxMass = 0;
   B_J4_PVx = 0;
   B_J4_PVy = 0;
   B_J4_PVz = 0;
   B_J4_PVxError = 0;
   B_J4_PVyError = 0;
   B_J4_PVzError = 0;
   B_Mu1_px = 0;
   B_Mu1_py = 0;
   B_Mu1_pz = 0;
   B_Mu1_pt = 0;
   B_Mu1_eta = 0;
   B_Mu1_phi = 0;
   B_Mu1_charge = 0;
   B_Mu1_soft = 0;
   B_Mu1_tight = 0;
   B_Mu1_loose = 0;
   B_Mu1_IsoTrack = 0;

//   Y_trackiso1 = 0;
   // Y_looseHEEPwp1 = 0;
   Z_loose1 = 0;

   // Y_looseHZZwp1 = 0;
   // Y_looseHZZwp2 = 0;
   
   Y_mvaIsoWP90_1 = 0;
   Y_mvaIsoWP90_2 = 0;

   Y_looseCutBase1 = 0;
   Y_looseCutBase2 = 0;

   B_Mu1_IsoHcal = 0;
   B_Mu1_IsoEcal = 0;
   B_Mu1_IsoCalo = 0;
   B_Mu2_px = 0;
   B_Mu2_py = 0;
   B_Mu2_pz = 0;
   B_Mu2_pt = 0;
   B_Mu2_eta = 0;
   B_Mu2_phi = 0;
   B_Mu2_charge = 0;
   B_Mu2_soft = 0;
   B_Mu2_tight = 0;
   B_Mu2_loose = 0;
   B_Mu2_IsoTrack = 0;
   B_Mu2_IsoHcal = 0;
   B_Mu2_IsoEcal = 0;
   B_Mu2_IsoCalo = 0;
   B_Mu3_px = 0;
   B_Mu3_py = 0;
   B_Mu3_pz = 0;
   B_Mu3_pt = 0;
   B_Mu3_eta = 0;
   B_Mu3_phi = 0;
   B_Mu3_soft = 0;
   B_Mu3_tight = 0;
   B_Mu3_loose = 0;
   B_Mu3_charge = 0;
   B_Mu3_IsoTrack = 0;
   B_Mu3_IsoHcal = 0;
   B_Mu3_IsoEcal = 0;
   B_Mu3_IsoCalo = 0;
   B_Mu4_px = 0;
   B_Mu4_py = 0;
   B_Mu4_pz = 0;
   B_Mu4_pt = 0;
   B_Mu4_eta = 0;
   B_Mu4_phi = 0;
   B_Mu4_soft = 0;
   B_Mu4_tight = 0;
   B_Mu4_loose = 0;
   B_Mu4_charge = 0;
   B_Mu4_IsoTrack = 0;
   B_Mu4_IsoHcal = 0;
   B_Mu4_IsoEcal = 0;
   B_Mu4_IsoCalo = 0;

   Y_VtxProb = 0;
   Z_VtxProb = 0;

   Z_soft1 = 0;
   Z_soft2 = 0;

   B_J3_VtxProb = 0;
   B_J4_VtxProb = 0;
   Z_xyP1 = 0;
   Z_xyM1 = 0;
   Z_zP1 = 0;
   Z_zM1 = 0;
   Z_xyP2 = 0;
   Z_xyM2 = 0;
   Z_zP2 = 0;
   Z_zM2 = 0;
   mu1mC2 = 0;
   mu1mNHits = 0;
   mu1mNPHits = 0;
   mu1pC2 = 0;
   mu1pNHits = 0;
   mu1pNPHits = 0;
   mu2mC2 = 0;
   mu2mNHits = 0;
   mu2mNPHits = 0;
   mu2pC2 = 0;
   mu2pNHits = 0;
   mu2pNPHits = 0;
   B_M1_pt = 0;
   B_M1_eta = 0;
   B_M1_phi = 0;
   B_M1_px = 0;
   B_M1_py = 0;
   B_M1_pz = 0;
   B_M2_pt = 0;
   B_M2_eta = 0;
   B_M2_phi = 0;
   B_M2_px = 0;
   B_M2_py = 0;
   B_M2_pz = 0;
   B_M3_pt = 0;
   B_M3_eta = 0;
   B_M3_phi = 0;
   B_M3_px = 0;
   B_M3_py = 0;
   B_M3_pz = 0;
   B_M4_pt = 0;
   B_M4_eta = 0;
   B_M4_phi = 0;
   B_M4_px = 0;
   B_M4_py = 0;
   B_M4_pz = 0;
   Z_GenMuonPt = 0;
   Z_GenMuonEta = 0;
   Z_GenMuonPhi = 0;
   Y_GenMuonPt = 0;
   Y_GenMuonEta = 0;
   Y_GenMuonPhi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nB", &nB, &b_nB);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("FourL_mass", &FourL_mass, &b_FourL_mass);
   fChain->SetBranchAddress("FourL_px", &FourL_px, &b_FourL_px);
   fChain->SetBranchAddress("FourL_py", &FourL_py, &b_FourL_py);
   fChain->SetBranchAddress("FourL_pz", &FourL_pz, &b_FourL_pz);
   fChain->SetBranchAddress("FourL_pt", &FourL_pt, &b_FourL_pt);
   fChain->SetBranchAddress("FourL_eta", &FourL_eta, &b_FourL_eta);
   fChain->SetBranchAddress("FourL_phi", &FourL_phi, &b_FourL_phi);
   fChain->SetBranchAddress("FourL_VtxProb", &FourL_VtxProb, &b_FourL_VtxProb);
   fChain->SetBranchAddress("FourL_PVx", &FourL_PVx, &b_FourL_PVx);
   fChain->SetBranchAddress("FourL_PVy", &FourL_PVy, &b_FourL_PVy);
   fChain->SetBranchAddress("FourL_PVz", &FourL_PVz, &b_FourL_PVz);
   fChain->SetBranchAddress("FourL_PVxError", &FourL_PVxError, &b_FourL_PVxError);
   fChain->SetBranchAddress("FourL_PVyError", &FourL_PVyError, &b_FourL_PVyError);
   fChain->SetBranchAddress("FourL_PVzError", &FourL_PVzError, &b_FourL_PVzError);
   fChain->SetBranchAddress("Y_TriggerPath", &Y_TriggerPath, &b_Y_TriggerPath);
   fChain->SetBranchAddress("Y_TriggerPt1", &Y_TriggerPt1, &b_Y_TriggerPt1);
   fChain->SetBranchAddress("Y_TriggerEta1", &Y_TriggerEta1, &b_Y_TriggerEta1);
   fChain->SetBranchAddress("Y_TriggerPhi1", &Y_TriggerPhi1, &b_Y_TriggerPhi1);
   fChain->SetBranchAddress("Y_TriggerPt2", &Y_TriggerPt2, &b_Y_TriggerPt2);
   fChain->SetBranchAddress("Y_TriggerEta2", &Y_TriggerEta2, &b_Y_TriggerEta2);
   fChain->SetBranchAddress("Y_TriggerPhi2", &Y_TriggerPhi2, &b_Y_TriggerPhi2);
   fChain->SetBranchAddress("Y_TriggerPt3", &Y_TriggerPt3, &b_Y_TriggerPt3);
   fChain->SetBranchAddress("Y_TriggerEta3", &Y_TriggerEta3, &b_Y_TriggerEta3);
   fChain->SetBranchAddress("Y_TriggerPhi3", &Y_TriggerPhi3, &b_Y_TriggerPhi3);
   fChain->SetBranchAddress("Y_TriggerPt4", &Y_TriggerPt4, &b_Y_TriggerPt4);
   fChain->SetBranchAddress("Y_TriggerEta4", &Y_TriggerEta4, &b_Y_TriggerEta4);
   fChain->SetBranchAddress("Y_TriggerPhi4", &Y_TriggerPhi4, &b_Y_TriggerPhi4);
   fChain->SetBranchAddress("Y_TriggerPt5", &Y_TriggerPt5, &b_Y_TriggerPt5);
   fChain->SetBranchAddress("Y_TriggerEta5", &Y_TriggerEta5, &b_Y_TriggerEta5);
   fChain->SetBranchAddress("Y_TriggerPhi5", &Y_TriggerPhi5, &b_Y_TriggerPhi5);
   fChain->SetBranchAddress("Y_mass", &Y_mass, &b_Y_mass);
   fChain->SetBranchAddress("Y_px", &Y_px, &b_Y_px);
   fChain->SetBranchAddress("Y_py", &Y_py, &b_Y_py);
   fChain->SetBranchAddress("Y_pz", &Y_pz, &b_Y_pz);
   fChain->SetBranchAddress("Y_pt", &Y_pt, &b_Y_pt);
   fChain->SetBranchAddress("Y_eta", &Y_eta, &b_Y_eta);
   fChain->SetBranchAddress("Y_phi", &Y_phi, &b_Y_phi);
   fChain->SetBranchAddress("Y_rapidity", &Y_rapidity, &b_Y_rapidity);
   fChain->SetBranchAddress("Z_mass", &Z_mass, &b_Z_mass);
   fChain->SetBranchAddress("Z_px", &Z_px, &b_Z_px);
   fChain->SetBranchAddress("Z_py", &Z_py, &b_Z_py);
   fChain->SetBranchAddress("Z_pz", &Z_pz, &b_Z_pz);
   fChain->SetBranchAddress("Z_pt", &Z_pt, &b_Z_pt);
   fChain->SetBranchAddress("Z_eta", &Z_eta, &b_Z_eta);
   fChain->SetBranchAddress("Z_phi", &Z_phi, &b_Z_phi);
   fChain->SetBranchAddress("Z_rapidity", &Z_rapidity, &b_Z_rapidity);

   fChain->SetBranchAddress("Y_eta1", &Y_eta1, &Y_eta1_b);
   fChain->SetBranchAddress("Y_pt1", &Y_pt1, &Y_pt1_b);
   fChain->SetBranchAddress("Y_eta2", &Y_eta2, &Y_eta2_b);
   fChain->SetBranchAddress("Y_pt2", &Y_pt2, &Y_pt2_b);
   fChain->SetBranchAddress("Z_pt1", &Z_pt1, &Z_pt1_b);
   fChain->SetBranchAddress("Z_eta1", &Z_eta1, &Z_eta1_b);  
   fChain->SetBranchAddress("Z_pt2", &Z_pt2, &Z_pt2_b);
   fChain->SetBranchAddress("Z_eta2", &Z_eta2, &Z_eta2_b);

   // fChain->SetBranchAddress("Y_looseHEEPwp1", &Y_looseHEEPwp1, &b_Y_looseHEEPwp1);
   fChain->SetBranchAddress("Z_loose1", &Z_loose1, &b_Z_loose1);
   // fChain->SetBranchAddress("Y_looseHZZwp1", &Y_looseHZZwp1, &b_Y_looseHZZwp1);
   // fChain->SetBranchAddress("Y_looseHZZwp2", &Y_looseHZZwp2, &b_Y_looseHZZwp2);
   fChain->SetBranchAddress("Y_mvaIsoWP90_1", &Y_mvaIsoWP90_1, &b_Y_mvaIsoWP90_1);
   fChain->SetBranchAddress("Y_mvaIsoWP90_2", &Y_mvaIsoWP90_2, &b_Y_mvaIsoWP90_2);


   fChain->SetBranchAddress("Y_looseCutBase1", &Y_looseCutBase1, &b_Y_looseCutBase1);
   fChain->SetBranchAddress("Y_looseCutBase2", &Y_looseCutBase2, &b_Y_looseCutBase2);
 
   fChain->SetBranchAddress("Z_soft1", &Z_soft1, &b_Z_soft1);
   fChain->SetBranchAddress("Z_soft2", &Z_soft2, &b_Z_soft2);

   // fChain->SetBranchAddress("B_J3_mass", &B_J3_mass, &b_B_J3_mass);
   // fChain->SetBranchAddress("B_J3_px", &B_J3_px, &b_B_J3_px);
   // fChain->SetBranchAddress("B_J3_py", &B_J3_py, &b_B_J3_py);
   // fChain->SetBranchAddress("B_J3_pz", &B_J3_pz, &b_B_J3_pz);
   // fChain->SetBranchAddress("B_J3_pt", &B_J3_pt, &b_B_J3_pt);
   // fChain->SetBranchAddress("B_J3_eta", &B_J3_eta, &b_B_J3_eta);
   // fChain->SetBranchAddress("B_J3_phi", &B_J3_phi, &b_B_J3_phi);
   // fChain->SetBranchAddress("B_J3_rapidity", &B_J3_rapidity, &b_B_J3_rapidity);
   // fChain->SetBranchAddress("B_J4_mass", &B_J4_mass, &b_B_J4_mass);
   // fChain->SetBranchAddress("B_J4_px", &B_J4_px, &b_B_J4_px);
   // fChain->SetBranchAddress("B_J4_py", &B_J4_py, &b_B_J4_py);
   // fChain->SetBranchAddress("B_J4_pz", &B_J4_pz, &b_B_J4_pz);
   // fChain->SetBranchAddress("B_J4_pt", &B_J4_pt, &b_B_J4_pt);
   // fChain->SetBranchAddress("B_J4_eta", &B_J4_eta, &b_B_J4_eta);
   // fChain->SetBranchAddress("B_J4_phi", &B_J4_phi, &b_B_J4_phi);
   // fChain->SetBranchAddress("B_J4_rapidity", &B_J4_rapidity, &b_B_J4_rapidity);
   // fChain->SetBranchAddress("B_J1_VtxPx", &B_J1_VtxPx, &b_B_J1_VtxPx);
   // fChain->SetBranchAddress("B_J1_VtxPy", &B_J1_VtxPy, &b_B_J1_VtxPy);
   // fChain->SetBranchAddress("B_J1_VtxPz", &B_J1_VtxPz, &b_B_J1_VtxPz);
   // fChain->SetBranchAddress("B_J1_VtxPt", &B_J1_VtxPt, &b_B_J1_VtxPt);
   // fChain->SetBranchAddress("B_J1_VtxEta", &B_J1_VtxEta, &b_B_J1_VtxEta);
   // fChain->SetBranchAddress("B_J1_VtxPhi", &B_J1_VtxPhi, &b_B_J1_VtxPhi);
   // fChain->SetBranchAddress("B_J1_VtxRapidity", &B_J1_VtxRapidity, &b_B_J1_VtxRapidity);
   // fChain->SetBranchAddress("B_J1_VtxMass", &B_J1_VtxMass, &b_B_J1_VtxMass);
   // fChain->SetBranchAddress("B_J1_PVx", &B_J1_PVx, &b_B_J1_PVx);
   // fChain->SetBranchAddress("B_J1_PVy", &B_J1_PVy, &b_B_J1_PVy);
   // fChain->SetBranchAddress("B_J1_PVz", &B_J1_PVz, &b_B_J1_PVz);
   // fChain->SetBranchAddress("B_J1_PVxError", &B_J1_PVxError, &b_B_J1_PVxError);
   // fChain->SetBranchAddress("B_J1_PVyError", &B_J1_PVyError, &b_B_J1_PVyError);
   // fChain->SetBranchAddress("B_J1_PVzError", &B_J1_PVzError, &b_B_J1_PVzError);
   // fChain->SetBranchAddress("B_J2_VtxPx", &B_J2_VtxPx, &b_B_J2_VtxPx);
   // fChain->SetBranchAddress("B_J2_VtxPy", &B_J2_VtxPy, &b_B_J2_VtxPy);
   // fChain->SetBranchAddress("B_J2_VtxPz", &B_J2_VtxPz, &b_B_J2_VtxPz);
   // fChain->SetBranchAddress("B_J2_VtxPt", &B_J2_VtxPt, &b_B_J2_VtxPt);
   // fChain->SetBranchAddress("B_J2_VtxEta", &B_J2_VtxEta, &b_B_J2_VtxEta);
   // fChain->SetBranchAddress("B_J2_VtxPhi", &B_J2_VtxPhi, &b_B_J2_VtxPhi);
   // fChain->SetBranchAddress("B_J2_VtxRapidity", &B_J2_VtxRapidity, &b_B_J2_VtxRapidity);
   // fChain->SetBranchAddress("B_J2_VtxMass", &B_J2_VtxMass, &b_B_J2_VtxMass);
   // fChain->SetBranchAddress("B_J2_PVx", &B_J2_PVx, &b_B_J2_PVx);
   // fChain->SetBranchAddress("B_J2_PVy", &B_J2_PVy, &b_B_J2_PVy);
   // fChain->SetBranchAddress("B_J2_PVz", &B_J2_PVz, &b_B_J2_PVz);
   // fChain->SetBranchAddress("B_J2_PVxError", &B_J2_PVxError, &b_B_J2_PVxError);
   // fChain->SetBranchAddress("B_J2_PVyError", &B_J2_PVyError, &b_B_J2_PVyError);
   // fChain->SetBranchAddress("B_J2_PVzError", &B_J2_PVzError, &b_B_J2_PVzError);
   // fChain->SetBranchAddress("B_J3_VtxPx", &B_J3_VtxPx, &b_B_J3_VtxPx);
   // fChain->SetBranchAddress("B_J3_VtxPy", &B_J3_VtxPy, &b_B_J3_VtxPy);
   // fChain->SetBranchAddress("B_J3_VtxPz", &B_J3_VtxPz, &b_B_J3_VtxPz);
   // fChain->SetBranchAddress("B_J3_VtxPt", &B_J3_VtxPt, &b_B_J3_VtxPt);
   // fChain->SetBranchAddress("B_J3_VtxEta", &B_J3_VtxEta, &b_B_J3_VtxEta);
   // fChain->SetBranchAddress("B_J3_VtxPhi", &B_J3_VtxPhi, &b_B_J3_VtxPhi);
   // fChain->SetBranchAddress("B_J3_VtxRapidity", &B_J3_VtxRapidity, &b_B_J3_VtxRapidity);
   // fChain->SetBranchAddress("B_J3_VtxMass", &B_J3_VtxMass, &b_B_J3_VtxMass);
   // fChain->SetBranchAddress("B_J3_PVx", &B_J3_PVx, &b_B_J3_PVx);
   // fChain->SetBranchAddress("B_J3_PVy", &B_J3_PVy, &b_B_J3_PVy);
   // fChain->SetBranchAddress("B_J3_PVz", &B_J3_PVz, &b_B_J3_PVz);
   // fChain->SetBranchAddress("B_J3_PVxError", &B_J3_PVxError, &b_B_J3_PVxError);
   // fChain->SetBranchAddress("B_J3_PVyError", &B_J3_PVyError, &b_B_J3_PVyError);
   // fChain->SetBranchAddress("B_J3_PVzError", &B_J3_PVzError, &b_B_J3_PVzError);
   // fChain->SetBranchAddress("B_J4_VtxPx", &B_J4_VtxPx, &b_B_J4_VtxPx);
   // fChain->SetBranchAddress("B_J4_VtxPy", &B_J4_VtxPy, &b_B_J4_VtxPy);
   // fChain->SetBranchAddress("B_J4_VtxPz", &B_J4_VtxPz, &b_B_J4_VtxPz);
   // fChain->SetBranchAddress("B_J4_VtxPt", &B_J4_VtxPt, &b_B_J4_VtxPt);
   // fChain->SetBranchAddress("B_J4_VtxEta", &B_J4_VtxEta, &b_B_J4_VtxEta);
   // fChain->SetBranchAddress("B_J4_VtxPhi", &B_J4_VtxPhi, &b_B_J4_VtxPhi);
   // fChain->SetBranchAddress("B_J4_VtxRapidity", &B_J4_VtxRapidity, &b_B_J4_VtxRapidity);
   // fChain->SetBranchAddress("B_J4_VtxMass", &B_J4_VtxMass, &b_B_J4_VtxMass);
   // fChain->SetBranchAddress("B_J4_PVx", &B_J4_PVx, &b_B_J4_PVx);
   // fChain->SetBranchAddress("B_J4_PVy", &B_J4_PVy, &b_B_J4_PVy);
   // fChain->SetBranchAddress("B_J4_PVz", &B_J4_PVz, &b_B_J4_PVz);
   // fChain->SetBranchAddress("B_J4_PVxError", &B_J4_PVxError, &b_B_J4_PVxError);
   // fChain->SetBranchAddress("B_J4_PVyError", &B_J4_PVyError, &b_B_J4_PVyError);
   // fChain->SetBranchAddress("B_J4_PVzError", &B_J4_PVzError, &b_B_J4_PVzError);
   // fChain->SetBranchAddress("B_Mu1_px", &B_Mu1_px, &b_B_Mu1_px);
   // fChain->SetBranchAddress("B_Mu1_py", &B_Mu1_py, &b_B_Mu1_py);
   // fChain->SetBranchAddress("B_Mu1_pz", &B_Mu1_pz, &b_B_Mu1_pz);
   // fChain->SetBranchAddress("B_Mu1_pt", &B_Mu1_pt, &b_B_Mu1_pt);
   // fChain->SetBranchAddress("B_Mu1_eta", &B_Mu1_eta, &b_B_Mu1_eta);
   // fChain->SetBranchAddress("B_Mu1_phi", &B_Mu1_phi, &b_B_Mu1_phi);
   // fChain->SetBranchAddress("B_Mu1_charge", &B_Mu1_charge, &b_B_Mu1_charge);
   // fChain->SetBranchAddress("B_Mu1_soft", &B_Mu1_soft, &b_B_Mu1_soft);
   // fChain->SetBranchAddress("B_Mu1_tight", &B_Mu1_tight, &b_B_Mu1_tight);
   // fChain->SetBranchAddress("B_Mu1_loose", &B_Mu1_loose, &b_B_Mu1_loose);
   // fChain->SetBranchAddress("B_Mu1_IsoTrack", &B_Mu1_IsoTrack, &b_B_Mu1_IsoTrack);

//    fChain->SetBranchAddress("Y_trackiso1", &Y_trackiso1, &b_Y_trackiso1);

   // fChain->SetBranchAddress("B_Mu1_IsoHcal", &B_Mu1_IsoHcal, &b_B_Mu1_IsoHcal);
   // fChain->SetBranchAddress("B_Mu1_IsoEcal", &B_Mu1_IsoEcal, &b_B_Mu1_IsoEcal);
   // fChain->SetBranchAddress("B_Mu1_IsoCalo", &B_Mu1_IsoCalo, &b_B_Mu1_IsoCalo);
   // fChain->SetBranchAddress("B_Mu2_px", &B_Mu2_px, &b_B_Mu2_px);
   // fChain->SetBranchAddress("B_Mu2_py", &B_Mu2_py, &b_B_Mu2_py);
   // fChain->SetBranchAddress("B_Mu2_pz", &B_Mu2_pz, &b_B_Mu2_pz);
   // fChain->SetBranchAddress("B_Mu2_pt", &B_Mu2_pt, &b_B_Mu2_pt);
   // fChain->SetBranchAddress("B_Mu2_eta", &B_Mu2_eta, &b_B_Mu2_eta);
   // fChain->SetBranchAddress("B_Mu2_phi", &B_Mu2_phi, &b_B_Mu2_phi);
   // fChain->SetBranchAddress("B_Mu2_charge", &B_Mu2_charge, &b_B_Mu2_charge);
   // fChain->SetBranchAddress("B_Mu2_soft", &B_Mu2_soft, &b_B_Mu2_soft);
   // fChain->SetBranchAddress("B_Mu2_tight", &B_Mu2_tight, &b_B_Mu2_tight);
   // fChain->SetBranchAddress("B_Mu2_loose", &B_Mu2_loose, &b_B_Mu2_loose);
   // fChain->SetBranchAddress("B_Mu2_IsoTrack", &B_Mu2_IsoTrack, &b_B_Mu2_IsoTrack);
   // fChain->SetBranchAddress("B_Mu2_IsoHcal", &B_Mu2_IsoHcal, &b_B_Mu2_IsoHcal);
   // fChain->SetBranchAddress("B_Mu2_IsoEcal", &B_Mu2_IsoEcal, &b_B_Mu2_IsoEcal);
   // fChain->SetBranchAddress("B_Mu2_IsoCalo", &B_Mu2_IsoCalo, &b_B_Mu2_IsoCalo);
   // fChain->SetBranchAddress("B_Mu3_px", &B_Mu3_px, &b_B_Mu3_px);
   // fChain->SetBranchAddress("B_Mu3_py", &B_Mu3_py, &b_B_Mu3_py);
   // fChain->SetBranchAddress("B_Mu3_pz", &B_Mu3_pz, &b_B_Mu3_pz);
   // fChain->SetBranchAddress("B_Mu3_pt", &B_Mu3_pt, &b_B_Mu3_pt);
   // fChain->SetBranchAddress("B_Mu3_eta", &B_Mu3_eta, &b_B_Mu3_eta);
   // fChain->SetBranchAddress("B_Mu3_phi", &B_Mu3_phi, &b_B_Mu3_phi);
   // fChain->SetBranchAddress("B_Mu3_soft", &B_Mu3_soft, &b_B_Mu3_soft);
   // fChain->SetBranchAddress("B_Mu3_tight", &B_Mu3_tight, &b_B_Mu3_tight);
   // fChain->SetBranchAddress("B_Mu3_loose", &B_Mu3_loose, &b_B_Mu3_loose);
   // fChain->SetBranchAddress("B_Mu3_charge", &B_Mu3_charge, &b_B_Mu3_charge);
   // fChain->SetBranchAddress("B_Mu3_IsoTrack", &B_Mu3_IsoTrack, &b_B_Mu3_IsoTrack);
   // fChain->SetBranchAddress("B_Mu3_IsoHcal", &B_Mu3_IsoHcal, &b_B_Mu3_IsoHcal);
   // fChain->SetBranchAddress("B_Mu3_IsoEcal", &B_Mu3_IsoEcal, &b_B_Mu3_IsoEcal);
   // fChain->SetBranchAddress("B_Mu3_IsoCalo", &B_Mu3_IsoCalo, &b_B_Mu3_IsoCalo);
   // fChain->SetBranchAddress("B_Mu4_px", &B_Mu4_px, &b_B_Mu4_px);
   // fChain->SetBranchAddress("B_Mu4_py", &B_Mu4_py, &b_B_Mu4_py);
   // fChain->SetBranchAddress("B_Mu4_pz", &B_Mu4_pz, &b_B_Mu4_pz);
   // fChain->SetBranchAddress("B_Mu4_pt", &B_Mu4_pt, &b_B_Mu4_pt);
   // fChain->SetBranchAddress("B_Mu4_eta", &B_Mu4_eta, &b_B_Mu4_eta);
   // fChain->SetBranchAddress("B_Mu4_phi", &B_Mu4_phi, &b_B_Mu4_phi);
   // fChain->SetBranchAddress("B_Mu4_soft", &B_Mu4_soft, &b_B_Mu4_soft);
   // fChain->SetBranchAddress("B_Mu4_tight", &B_Mu4_tight, &b_B_Mu4_tight);
   // fChain->SetBranchAddress("B_Mu4_loose", &B_Mu4_loose, &b_B_Mu4_loose);
   // fChain->SetBranchAddress("B_Mu4_charge", &B_Mu4_charge, &b_B_Mu4_charge);
   // fChain->SetBranchAddress("B_Mu4_IsoTrack", &B_Mu4_IsoTrack, &b_B_Mu4_IsoTrack);
   // fChain->SetBranchAddress("B_Mu4_IsoHcal", &B_Mu4_IsoHcal, &b_B_Mu4_IsoHcal);
   // fChain->SetBranchAddress("B_Mu4_IsoEcal", &B_Mu4_IsoEcal, &b_B_Mu4_IsoEcal);
   // fChain->SetBranchAddress("B_Mu4_IsoCalo", &B_Mu4_IsoCalo, &b_B_Mu4_IsoCalo);
   fChain->SetBranchAddress("Y_VtxProb", &Y_VtxProb, &b_Y_VtxProb);
   fChain->SetBranchAddress("Z_VtxProb", &Z_VtxProb, &b_Z_VtxProb);
   // fChain->SetBranchAddress("B_J3_VtxProb", &B_J3_VtxProb, &b_B_J3_VtxProb);
   // fChain->SetBranchAddress("B_J4_VtxProb", &B_J4_VtxProb, &b_B_J4_VtxProb);
   // fChain->SetBranchAddress("Z_xyP1", &Z_xyP1, &b_Z_xyP1);
   // fChain->SetBranchAddress("Z_xyM1", &Z_xyM1, &b_Z_xyM1);
   // fChain->SetBranchAddress("Z_zP1", &Z_zP1, &b_Z_zP1);
   // fChain->SetBranchAddress("Z_zM1", &Z_zM1, &b_Z_zM1);
   // fChain->SetBranchAddress("Z_xyP2", &Z_xyP2, &b_Z_xyP2);
   // fChain->SetBranchAddress("Z_xyM2", &Z_xyM2, &b_Z_xyM2);
   // fChain->SetBranchAddress("Z_zP2", &Z_zP2, &b_Z_zP2);
   // fChain->SetBranchAddress("Z_zM2", &Z_zM2, &b_Z_zM2);
   // fChain->SetBranchAddress("mu1mC2", &mu1mC2, &b_mu1mC2);
   // fChain->SetBranchAddress("mu1mNHits", &mu1mNHits, &b_mu1mNHits);
   // fChain->SetBranchAddress("mu1mNPHits", &mu1mNPHits, &b_mu1mNPHits);
   // fChain->SetBranchAddress("mu1pC2", &mu1pC2, &b_mu1pC2);
   // fChain->SetBranchAddress("mu1pNHits", &mu1pNHits, &b_mu1pNHits);
   // fChain->SetBranchAddress("mu1pNPHits", &mu1pNPHits, &b_mu1pNPHits);
   // fChain->SetBranchAddress("mu2mC2", &mu2mC2, &b_mu2mC2);
   // fChain->SetBranchAddress("mu2mNHits", &mu2mNHits, &b_mu2mNHits);
   // fChain->SetBranchAddress("mu2mNPHits", &mu2mNPHits, &b_mu2mNPHits);
   // fChain->SetBranchAddress("mu2pC2", &mu2pC2, &b_mu2pC2);
   // fChain->SetBranchAddress("mu2pNHits", &mu2pNHits, &b_mu2pNHits);
   // fChain->SetBranchAddress("mu2pNPHits", &mu2pNPHits, &b_mu2pNPHits);
   // fChain->SetBranchAddress("B_M1_pt", &B_M1_pt, &b_B_M1_pt);
   // fChain->SetBranchAddress("B_M1_eta", &B_M1_eta, &b_B_M1_eta);
   // fChain->SetBranchAddress("B_M1_phi", &B_M1_phi, &b_B_M1_phi);
   // fChain->SetBranchAddress("B_M1_px", &B_M1_px, &b_B_M1_px);
   // fChain->SetBranchAddress("B_M1_py", &B_M1_py, &b_B_M1_py);
   // fChain->SetBranchAddress("B_M1_pz", &B_M1_pz, &b_B_M1_pz);
   // fChain->SetBranchAddress("B_M2_pt", &B_M2_pt, &b_B_M2_pt);
   // fChain->SetBranchAddress("B_M2_eta", &B_M2_eta, &b_B_M2_eta);
   // fChain->SetBranchAddress("B_M2_phi", &B_M2_phi, &b_B_M2_phi);
   // fChain->SetBranchAddress("B_M2_px", &B_M2_px, &b_B_M2_px);
   // fChain->SetBranchAddress("B_M2_py", &B_M2_py, &b_B_M2_py);
   // fChain->SetBranchAddress("B_M2_pz", &B_M2_pz, &b_B_M2_pz);
   // fChain->SetBranchAddress("B_M3_pt", &B_M3_pt, &b_B_M3_pt);
   // fChain->SetBranchAddress("B_M3_eta", &B_M3_eta, &b_B_M3_eta);
   // fChain->SetBranchAddress("B_M3_phi", &B_M3_phi, &b_B_M3_phi);
   // fChain->SetBranchAddress("B_M3_px", &B_M3_px, &b_B_M3_px);
   // fChain->SetBranchAddress("B_M3_py", &B_M3_py, &b_B_M3_py);
   // fChain->SetBranchAddress("B_M3_pz", &B_M3_pz, &b_B_M3_pz);
   // fChain->SetBranchAddress("B_M4_pt", &B_M4_pt, &b_B_M4_pt);
   // fChain->SetBranchAddress("B_M4_eta", &B_M4_eta, &b_B_M4_eta);
   // fChain->SetBranchAddress("B_M4_phi", &B_M4_phi, &b_B_M4_phi);
   // fChain->SetBranchAddress("B_M4_px", &B_M4_px, &b_B_M4_px);
   // fChain->SetBranchAddress("B_M4_py", &B_M4_py, &b_B_M4_py);
   // fChain->SetBranchAddress("B_M4_pz", &B_M4_pz, &b_B_M4_pz);
   // fChain->SetBranchAddress("Z_GenMuonPt", &Z_GenMuonPt, &b_Z_GenMuonPt);
   // fChain->SetBranchAddress("Z_GenMuonEta", &Z_GenMuonEta, &b_Z_GenMuonEta);
   // fChain->SetBranchAddress("Z_GenMuonPhi", &Z_GenMuonPhi, &b_Z_GenMuonPhi);
   // fChain->SetBranchAddress("Y_GenMuonPt", &Y_GenMuonPt, &b_Y_GenMuonPt);
   // fChain->SetBranchAddress("Y_GenMuonEta", &Y_GenMuonEta, &b_Y_GenMuonEta);
   // fChain->SetBranchAddress("Y_GenMuonPhi", &Y_GenMuonPhi, &b_Y_GenMuonPhi);
   Notify();
}

Bool_t AnalysisWithIso::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnalysisWithIso::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalysisWithIso::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnalysisWithIso_cxx
