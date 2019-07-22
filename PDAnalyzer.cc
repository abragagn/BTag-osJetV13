#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>

#include "PDAnalyzer.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TFile.h"

// additional features
#include "PDSecondNtupleWriter.h"   // second ntuple
#include "PDMuonVar.cc"
#include "PDSoftMuonMvaEstimator.cc"
#include "AlbertoUtil.cc"
#include "OSMuonMvaTag.cc"

using namespace std;

/*
pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2018Lists/BsToJpsiPhi_2018_DCAP.list hist.root -v outputFile ntu.root -v histoMode RECREATE -v use_gen t -n 10000
*/

PDAnalyzer::PDAnalyzer() {

    std::cout << "new PDAnalyzer" << std::endl;

    // user parameters are set as names associated to a string, 
    // default values can be set in the analyzer class contructor

    setUserParameter( "verbose", "f" );

    setUserParameter( "process", "BsJPsiPhi" );
    setUserParameter( "useHLT", "true" );

    //Jet Parameters
    setUserParameter( "minPtJet", "10" );
    setUserParameter( "cutBTag", "0.15" ); 
    // DeepJet loose = 0.0494, medium = 0.2770, tight = 0.7264
    setUserParameter( "jetDrCut", "0.5" ); // min dR wrt signal B

    setUserParameter( "jetDzCut", "1.0" ); // max trk dZ wrt PV
    
    setUserParameter( "ptCut", "40.0" ); //needed for paolo's code for unknow reasons

}


PDAnalyzer::~PDAnalyzer() {
}



void PDAnalyzer::beginJob() {

    PDAnalyzerUtil::beginJob();

    getUserParameter( "verbose", verbose );

    getUserParameter( "useHLT", useHLT );
    getUserParameter( "process", process );
 
    getUserParameter( "minPtJet", minPtJet );
    getUserParameter( "cutBTag", cutBTag );
    getUserParameter( "jetDrCut", jetDrCut );

    getUserParameter( "jetDzCut", jetDzCut );
    
    getUserParameter( "ptCut", ptCut ); //needed for paolo's code for unknow reasons

// additional features

    tWriter = new PDSecondNtupleWriter; // second ntuple
    tWriter->open( getUserParameter("outputFile"), "RECREATE" ); // second ntuple

    inizializeMuonMvaReader();

    if(process=="BsJPsiPhi") SetBsMassRange(5.20, 5.65);
    if(process=="BuJPsiK") SetBuMassRange(5.1, 5.65);

    nTot = 0;
    nPass = 0;

    return;

}


void PDAnalyzer::book() {

    // putting "autoSavedObject" in front of the histo creation 
    // it's automatically marked for saving on file; the option 
    // is uneffective when not using the full utility

    return;

}


void PDAnalyzer::reset() {
// automatic reset
    autoReset();
    return;
}


bool PDAnalyzer::analyze( int entry, int event_file, int event_tot ) {

    if ( verbose ) {
        cout << " +++++++++++++++++++++++++++ " << endl;
        cout << "entry: "
             << entry << " " << event_file << " " << event_tot << endl;
        cout << "run: " << runNumber << " , "
             << "evt: " << eventNumber << endl;
    }
    else {

        if ( (!(event_tot%10) && event_tot<100 ) || 
     (!(event_tot %100) && event_tot<1000 ) || 
     (!(event_tot %1000) && event_tot<10000 ) || 
     (!(event_tot %10000) && event_tot<100000 ) || 
     (!(event_tot %100000) && event_tot<1000000 ) || 
     (!(event_tot %1000000) && event_tot<10000000 ) )
            cout << " == at event " << event_file << " " << event_tot << endl;
    }
// additional features
    computeMuonVar();
    inizializeOsMuonTagVars();
    tWriter->Reset();
    convSpheCart(jetPt, jetEta, jetPhi, jetPx, jetPy, jetPz);
    convSpheCart(muoPt, muoEta, muoPhi, muoPx, muoPy, muoPz);
    convSpheCart(trkPt, trkEta, trkPhi, trkPx, trkPy, trkPz);
    convSpheCart(pfcPt, pfcEta, pfcPhi, pfcPx, pfcPy, pfcPz);

    if( !((process=="BsJPsiPhi")||(process=="BuJPsiK")) ) {
        cout<<"!$!#$@$% PROCESS NAME WRONG"<<endl;
        return false;
    }

//------------------------------------------------HLT---------------------------------------

    bool jpsimu = false;
    bool jpsitktk = false;
    bool jpsitk = false;

    if(hlt(PDEnumString::HLT_Dimuon0_Jpsi3p5_Muon2_v)||hlt(PDEnumString::HLT_Dimuon0_Jpsi_Muon_v)) jpsimu = true;
    if(hlt(PDEnumString::HLT_DoubleMu4_JpsiTrkTrk_Displaced_v)) jpsitktk =  true;
    if(hlt(PDEnumString::HLT_DoubleMu4_JpsiTrk_Displaced_v)) jpsitk = true;

    if(useHLT && process=="BsJPsiPhi" && !jpsitktk) return false;
    if(useHLT && process=="BuJPsiK" && !jpsitk) return false;
    if( jpsimu ) return false;

    SetJpsiTrkTrkCut();

//------------------------------------------------SEARCH FOR SS---------------------------------------

    int ssbSVT = GetCandidate(process);
    if(ssbSVT<0) return false;

    bool isTight = false;
    int ssbSVTtight = GetTightCandidate(process);
    if(ssbSVTtight>=0){
        isTight = true;
        ssbSVT = ssbSVTtight;
    }

    int iJPsi = (subVtxFromSV(ssbSVT)).at(0);
    vector <int> tkJpsi = tracksFromSV(iJPsi);
    vector <int> tkSsB = tracksFromSV(ssbSVT);

    TLorentzVector tB = GetTLorentzVecFromJpsiX(ssbSVT);

    //generation information

    vector <int> ListLongLivedB;
    vector <int> ListB;
    int genBindex = -1;
    int ssbLund = 0;
    int tagMix = -1;
    float evtWeight = 1;

    if(use_gen){
        for( uint i=0 ; i<genId->size() ; ++i ){
            if(TagMixStatus( i ) == 2) continue;
            if( IsB(i) ) ListB.push_back(i);
            uint Code = abs(genId->at(i));
            if( Code == 511 || Code == 521 || Code == 531 || Code == 541 || Code == 5122 ) ListLongLivedB.push_back(i);
        }

        genBindex = GetClosestGenLongLivedB( tB.Eta(), tB.Phi(), tB.Pt(), &ListLongLivedB);
        if(genBindex<0) return false;

        ssbLund = genId->at(genBindex);
        if((process=="BsJPsiPhi") && (abs(ssbLund)!=531)) return false;
        if((process=="BuJPsiK") && (abs(ssbLund)!=521)) return false;

        tagMix = TagMixStatus( genBindex );
        if(tagMix == 2) return false;
        if(tagMix == 1) ssbLund*=-1;

        for(auto it:ListLongLivedB){
            if(it == genBindex) continue;
            if(abs(genId->at(it)) == abs(ssbLund)) evtWeight = 2;
        }


    }else{
        if(process=="BsJPsiPhi") ssbLund = ((double)rand() / (RAND_MAX)) < 0.5 ? +531 : -531; //this should not be used
        if(process=="BuJPsiK"){
            for( auto it:tkSsB ){
                if( it == tkJpsi[0] || it == tkJpsi[1] ) continue;
                ssbLund = trkCharge->at(it) > 0 ? +521 : -521;
            }
        }
    }

    int ssbPVT = GetBestPV(ssbSVT, tB);
    if(ssbPVT < 0) return false;

    setVtxOsMuonTag(ssbSVT, ssbPVT);
    nTot += evtWeight;

    bool osmuon = false;
    if(selectOsMuon()>=0) osmuon = true;

    //Mario's Electron Selection
    bool osele = false;
    for(int iEle = 0; iEle<nElectrons; ++iEle){
        if(elePt->at(iEle) < 2.5) continue;
        if(fabs(eleEta->at(iEle)) > 2.4) continue;
        if(fabs(dZele(iEle, ssbPVT)) > 0.5) continue;
        if( deltaR(eleEta->at(iEle),elePhi->at(iEle),tB.Eta(),tB.Phi()) < 0.4) continue;
        // if(std::find(tkSsB.begin(), tkSsB.end(), eletk) != tkSsB.end()) continue;
        float idvalue = -2;
        for (int iUserInfo = 0; iUserInfo < nUserInfo; ++iUserInfo){
            if(useObjType->at(iUserInfo) == PDEnumString::recElectron && useObjIndex->at(iUserInfo) == iEle){
                if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Fall17NoIsoV2Values){
                    idvalue = useInfoValue->at(iUserInfo);
                }
            }
        }
        if(idvalue < -0.999) continue;
        osele = true;
    }


    //FILLING SS
    (tWriter->ssbMass) = svtMass->at(ssbSVT);
    (tWriter->ssbIsTight) = isTight;
    (tWriter->ssbLund) = ssbLund;

    (tWriter->ssbPt) = tB.Pt();
    (tWriter->ssbEta) = tB.Eta();
    (tWriter->ssbPhi) = tB.Phi();

    (tWriter->ssbPVTx) = pvtX->at(ssbPVT);
    (tWriter->ssbPVTy) = pvtY->at(ssbPVT);
    (tWriter->ssbPVTz) = pvtZ->at(ssbPVT);

    (tWriter->hltJpsiMu) = jpsimu;
    (tWriter->hltJpsiTrkTrk) = jpsitktk;
    (tWriter->hltJpsiTrk) = jpsitk;

    (tWriter->evtWeight) = evtWeight;
    (tWriter->evtNumber) = event_tot;
    (tWriter->isOSMuon) = osmuon;
    (tWriter->isOSEle) = osele;

//-----------------------------------------JET-----------------------------------------

    int bestJet = -1;
    float bestJetTag = cutBTag;

    //SELECTION
    for (int iJet = 0; iJet<nJets; ++iJet){

        if(goodJet(iJet)!=true) continue;
        if(abs(jetEta->at(iJet))>2.5) continue;
        if(jetPt->at(iJet)<minPtJet) continue;

        float bTag = GetJetProbb(iJet);
        if(bTag < cutBTag) continue;

        if( deltaR(jetEta->at(iJet),jetPhi->at(iJet),tB.Eta(),tB.Phi()) < jetDrCut ) continue;

        vector <int> jet_tks = tracksFromJet( iJet );

        bool skip = false; 
        for(auto it:jet_tks){
            if(std::find(tkSsB.begin(),tkSsB.end(),it) != tkSsB.end()){
                skip = true;
                break;
            }
        }
        if(skip) continue; // skip jet if contains signal tracks

        int nTrkNearPV = 0;
        for(auto it:jet_tks){
            if(fabs(dZ(it, ssbPVT))>=jetDzCut) continue;
            if( !isTrkHighPurity(it) ) continue;
            nTrkNearPV++;
        }

        if(nTrkNearPV < 2) continue; // skip jet if less than 2 tracks are near the PV

        if(bTag>bestJetTag){ // take the jet with the highest btag prob
            bestJetTag = bTag;
            bestJet = iJet;
        }
    }

    //TAG
    if(bestJet < 0){
        (tWriter->isOSJet) = false;
        tWriter->fill();
        return true;
    }

    (tWriter->isOSJet) = true;

    nPass += evtWeight;

    //indices
    int iJet = bestJet;
    vector <int> jet_tks = tracksFromJet( iJet );

    vector<int> selectedJetTracks;
    for(auto it:jet_tks){
        if(fabs(dZ(it, ssbPVT))>=jetDzCut) continue;
        if( !isTrkHighPurity(it) ) continue;
        selectedJetTracks.push_back(it);
    }

    //GENINFO
    int jetAncestor = GetJetAncestor( iJet, &ListB );

    //TAGGING VARIABLES
    //Separation
    float jetDzB=0;
    for(auto it:selectedJetTracks)
        jetDzB += fabs(dZ(it, ssbPVT));

    jetDzB/=selectedJetTracks.size();

    //Jet Charge
    float jet_charge = GetListCharge(&selectedJetTracks, 1.0);

    (tWriter->jetPt) = jetPt->at(iJet);
    (tWriter->jetEta) = jetEta->at(iJet);
    (tWriter->jetPhi) = jetPhi->at(iJet);
    (tWriter->jetCharge) = jet_charge;
    (tWriter->jetProbb) = GetJetProbb(iJet);
    (tWriter->jetDrB) = deltaR(jetEta->at(iJet),jetPhi->at(iJet),tB.Eta(),tB.Phi());
    (tWriter->jetDzB) = jetDzB;

    //------------------------------------------------TAG------------------------------------------------

    (tWriter->jetHasAncestor) = jetAncestor;

    //------------------------------------------------TRACKS------------------------------------------------

    for(auto it:jet_tks){
        if( !isTrkHighPurity(it)) continue;
        float dxy = dSign(it, jetPx->at(iJet), jetPy->at(iJet))*abs(trkDxy->at(it));
        float dz = dZ(it, ssbPVT);
        (tWriter->trkPt)->push_back(trkPt->at(it));
        (tWriter->trkEta)->push_back(trkEta->at(it));
        (tWriter->trkPhi)->push_back(trkPhi->at(it));
        (tWriter->trkCharge)->push_back(trkCharge->at(it));
        (tWriter->trkDxySigned)->push_back( dxy );
        (tWriter->trkDxy)->push_back( trkDxy->at(it) );
        (tWriter->trkDz)->push_back( dz );
        (tWriter->trkDxyz)->push_back( sqrt(pow(dxy,2)+pow(dz,2)) );
        (tWriter->trkExy)->push_back( trkExy->at(it) );
        (tWriter->trkEz)->push_back( trkEz->at(it) );
        (tWriter->trkDrJet)->push_back(deltaR(jetEta->at(iJet),jetPhi->at(iJet),trkEta->at(it),trkPhi->at(it)));        
        (tWriter->trkIsInJet)->push_back( 1 );
    }

    // all tracks collinear with the jet 
    for (int it = 0; it<nTracks; ++it){

        if( deltaR(jetEta->at(iJet), jetPhi->at(iJet), trkEta->at(it), trkPhi->at(it)) > 0.5 ) continue;
        if( fabs(dZ(it, ssbPVT)) > jetDzCut ) continue;
        if(std::find(jet_tks.begin(), jet_tks.end(), it) != jet_tks.end()) continue;
        if( !isTrkHighPurity(it) ) continue;

        float dxy = dSign(it, jetPx->at(iJet), jetPy->at(iJet))*abs(trkDxy->at(it));
        float dz = dZ(it, ssbPVT);
        (tWriter->trkPt)->push_back(trkPt->at(it));
        (tWriter->trkEta)->push_back(trkEta->at(it));
        (tWriter->trkPhi)->push_back(trkPhi->at(it));
        (tWriter->trkCharge)->push_back(trkCharge->at(it));
        (tWriter->trkDxySigned)->push_back( dxy );
        (tWriter->trkDxy)->push_back( trkDxy->at(it) );
        (tWriter->trkDz)->push_back( dz );
        (tWriter->trkDxyz)->push_back( sqrt(pow(dxy,2)+pow(dz,2)) );
        (tWriter->trkExy)->push_back( trkExy->at(it) );
        (tWriter->trkEz)->push_back( trkEz->at(it) );
        (tWriter->trkDrJet)->push_back(deltaR(jetEta->at(iJet),jetPhi->at(iJet),trkEta->at(it),trkPhi->at(it)));        
        (tWriter->trkIsInJet)->push_back( 0 );
    }

    tWriter->fill();

    return true;

}


void PDAnalyzer::endJob() {


// additional features

    tWriter->close();   // second ntuple

    cout<<nTot<<"   "<<nPass<<endl;
    cout<<100.*(float)nPass/nTot<<endl;

    return;
}


void PDAnalyzer::save() {
#   if UTIL_USE == FULL
    // explicit saving not necessary for "autoSavedObjects"
    autoSave();
#elif UTIL_USE == BARE
    // explicit save histos when not using the full utility

#endif

    return;
}


// ======MY FUNCTIONS===============================================================================
// =====================================================================================
int PDAnalyzer::GetJetAncestor( unsigned int iJet, vector<int> *GenList )
{
    double drb = 0.4;
    double dpb = minPtJet; 
    int best = -1;

    for(auto it:*GenList){

        float dr = deltaR(jetEta->at(iJet), jetPhi->at(iJet), genEta->at(it), genPhi->at(it));
        float pt = genPt->at(it);

        if( dr > drb ) continue;
        if( pt < dpb) continue;

        best = it;
        dpb = pt;

    }

    return best;
}
// =====================================================================================
int PDAnalyzer::GetClosesJet( float pt, float eta, float phi )
{
    double drb = 0.5;
    int best = -1;
    for( int i=0; i<nJets; ++i ){
       float dr = deltaR(eta, phi, jetEta->at(i), jetPhi->at(i));
       if( dr > drb ) continue;
       best = (int) i;
       drb = dr;
    }
    return best;
}


double PDAnalyzer::dZele(const int iEle, const int iVtx)
{
    number px;
    number py;
    number pz;

    if(eleGsfPx->size())
    {
        px = eleGsfPx->at(iEle);
        py = eleGsfPy->at(iEle);
        pz = eleGsfPz->at(iEle);
    }
    else if(eleGsfPt->size())
    {
        TVector3 pEle;
        pEle.SetPtEtaPhi(eleGsfPt->at(iEle), eleGsfEta->at(iEle), eleGsfPhi->at(iEle));
        px = pEle.Px();
        py = pEle.Py();
        pz = pEle.Pz();
    }
    else  if(eleGsfPxAtVtx->size())
    {
        px = eleGsfPxAtVtx->at(iEle);
        py = eleGsfPyAtVtx->at(iEle);
        pz = eleGsfPzAtVtx->at(iEle);
    }
    else if(eleGsfPtAtVtx->size())
    {
        TVector3 pEle;
        pEle.SetPtEtaPhi(eleGsfPtAtVtx->at(iEle), eleGsfEtaAtVtx->at(iEle), eleGsfPhiAtVtx->at(iEle));
        px = pEle.Px();
        py = pEle.Py();
        pz = pEle.Pz();
    }
    else if(elePx->size())
    {
        px = elePx->at(iEle);
        py = elePy->at(iEle);
        pz = elePz->at(iEle);
    }
    else
    {
        TVector3 pEle;
        pEle.SetPtEtaPhi(elePt->at(iEle), eleEta->at(iEle), elePhi->at(iEle));
        px = pEle.Px();
        py = pEle.Py();
        pz = pEle.Pz();
    }

    number pq = modSqua(px, py, 0);
    return (((px * (pvtX->at(iVtx) - bsX))) + (py * (pvtY->at(iVtx) - bsY)) * pz/pq) +
            eleGsfDz->at(iEle) + bsZ - pvtZ->at(iVtx);
}
