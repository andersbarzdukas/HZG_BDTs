#!/usr/bin/env python3
'''
Function that generates n-tuples for MVA training
'''
import ROOT
ROOT.ROOT.EnableImplicitMT(0)
def write_ntuples(filenames, cuts, out_name, defines=[], tree_name='tree', branches=()):
  '''Generate ROOT n-tuple from existing n-tuple
  
  Parameters:
  filenames - list of filenames of signal ROOT n-tuples
  cuts - list of cuts expressed as strings in order they should be applied
  out_name - output filename of n-tuple
  defines - list of 2-tuples describing new branches to define in the format (name, expr)
            note that these must be defineable before cuts
  tree_name - name of tree in ROOT file
  branches - tuple of branches to save; if empty all branches are saved
  '''
  filenames_vec = ROOT.std.vector('string')()
  print(filenames)
  for filename in filenames:
    filenames_vec.push_back(filename)

  df = ROOT.RDataFrame('tree',filenames_vec)
  print(df)

  for define in defines:
    print(define)
    df = df.Define(define[0],define[1])

  for cut in cuts:
    df = df.Filter(cut)

  if (branches == ()):
    df.Snapshot(tree_name,'ntuples/'+out_name)

  else:
    df.Snapshot(tree_name,'ntuples/'+out_name,branches)
  print('Wrote ntuples/'+out_name)

ROOT.gInterpreter.Declare("""
template <class C>
using RVec = ROOT::VecOps::RVec<C>;
#include "TROOT.h"
#include "TMVA/Reader.h"

float getMVA(float jj_deta,float jj_dphi,float Zyjj_dr,float pt_bal,float zep_var,float mjj,float drmin_yj,float j1_pt,float j2_pt,int event){

    TMVA::Reader *rdr =  new TMVA::Reader("Silent");
    rdr->TMVA::Reader::AddVariable("jj_deta", &jj_deta);
    rdr->TMVA::Reader::AddVariable("jj_dphi", &jj_dphi);
    rdr->TMVA::Reader::AddVariable("Zyjj_dr", &Zyjj_dr);
    rdr->TMVA::Reader::AddVariable("pt_bal",  &pt_bal);
    rdr->TMVA::Reader::AddVariable("zep_var", &zep_var);
    rdr->TMVA::Reader::AddVariable("mjj", &mjj);
    rdr->TMVA::Reader::AddVariable("drmin_yj", &drmin_yj);
    rdr->TMVA::Reader::AddVariable("j1_pt", &j1_pt);
    rdr->TMVA::Reader::AddVariable("j2_pt", &j2_pt);
    rdr -> TMVA::Reader::BookMVA("BDT","/net/cms27/cms27r0/abarzdukas/ANTools/HZG_Start/HZG_VBF_BDT/VBF_BDT/dataset/weights/TMVAClassification_VBF_test_noypt_BDTG.weights.xml"); 

    float score = rdr -> EvaluateMVA("BDT");
    delete rdr;
    return score;
}

RVec<int> truth_matched(RVec<float> mc_pt, RVec<float> mc_eta, RVec<float> mc_phi,RVec<float> mc_m, RVec<float> mc_id, RVec<float> mc_momidx, 
                        RVec<float> jet_pt, RVec<float> jet_eta, RVec<float> jet_phi, RVec<float> jet_m,RVec<bool> jet_isgood){
   RVec<int> jet_idx = {-1,-1}; 
   TLorentzVector jet,mc;

   for(int idx = 0; idx < jet_pt.size(); idx++){
        if(!jet_isgood[idx]){continue;}

        //IF YOU DONT WANT TRUTH-MATCHED JETS
        if(jet_idx[0]==-1){jet_idx[0]=idx;continue;}
          jet_idx[1]=idx;break;
          
       
        jet.SetPtEtaPhiM(jet_pt[idx],jet_eta[idx],jet_phi[idx],jet_m[idx]);

        for(int mc_idx = 0; mc_idx < mc_pt.size(); mc_idx++){
          if( (fabs(mc_id[mc_idx]) >10) ||(mc_momidx[mc_idx]==-1)||( mc_id[ mc_momidx[mc_idx] ] == 21) ){ continue; }
          mc.SetPtEtaPhiM(mc_pt[mc_idx],mc_eta[mc_idx],mc_phi[mc_idx],mc_m[mc_idx]);

          if(mc.DeltaR(jet) > 0.1 || fabs(mc.Pt() - jet.Pt()) > 0.1*jet.Pt()){ continue; }

          if(jet_idx[0]==-1){
            jet_idx[0] = idx;
            break;
          }
          jet_idx[1] = idx;
          break;

        } 

     if(jet_idx[1] > -1){break;}
    
   }
    
   return jet_idx;
}

bool truth_matched_bool(RVec<float> mc_pt, RVec<float> mc_eta, RVec<float> mc_phi,RVec<float> mc_m, RVec<float> mc_id, RVec<float> mc_momidx, 
                        RVec<float> jet_pt, RVec<float> jet_eta, RVec<float> jet_phi, RVec<float> jet_m,RVec<bool> jet_isgood){
    RVec<int> jet_idx = truth_matched(mc_pt,mc_eta,mc_phi,mc_m,mc_id,mc_momidx,jet_pt,jet_eta,jet_phi,jet_m,jet_isgood);
    if(jet_idx[0] == -1 || jet_idx[1] == -1){ return false; }
    return true;
}

float get_dr(float eta1, float phi1, float eta2, float phi2) {
  const double PI = 3.1415;
  double dphi = fmod(fabs(phi2-phi1), 2.*PI);
  dphi = dphi>PI ? 2.*PI-dphi : dphi;
  double deta = fabs(eta1-eta2);
  return sqrt(deta*deta+dphi*dphi);
}

float get_max_dr(RVec<float> photon_eta, RVec<float> photon_phi, 
    RVec<float> el_eta, RVec<float> el_phi, RVec<float> mu_eta,
    RVec<float> mu_phi, RVec<int> ll_lepid, RVec<int> ll_i1,
    RVec<int> ll_i2) {
  float dr1, dr2;
  if (ll_lepid[0]==11) {
    dr1 = get_dr(photon_eta[0],photon_phi[0],el_eta[ll_i1[0]],el_phi[ll_i1[0]]);
    dr2 = get_dr(photon_eta[0],photon_phi[0],el_eta[ll_i2[0]],el_phi[ll_i2[0]]);
    return dr1 > dr2 ? dr1 : dr2;
  }
  dr1 = get_dr(photon_eta[0],photon_phi[0],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]]);
  dr2 = get_dr(photon_eta[0],photon_phi[0],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]]);
  return dr1 > dr2 ? dr1 : dr2;
}

float get_l1_rapidity(RVec<float> el_pt, RVec<float> el_eta, 
    RVec<float> mu_pt, RVec<float> mu_eta, RVec<int> ll_lepid, 
    RVec<int> ll_i1, RVec<int> ll_i2) {
  if (ll_lepid[0]==11) {
    return (el_pt[ll_i1[0]] > el_pt[ll_i2[0]]) ? el_eta[ll_i1[0]] : el_eta[ll_i2[0]];
  }
  return (mu_pt[ll_i1[0]] > mu_pt[ll_i2[0]]) ? mu_eta[ll_i1[0]] : mu_eta[ll_i2[0]];
}

float get_l2_rapidity(RVec<float> el_pt, RVec<float> el_eta, 
    RVec<float> mu_pt, RVec<float> mu_eta, RVec<int> ll_lepid, 
    RVec<int> ll_i1, RVec<int> ll_i2) {
  if (ll_lepid[0]==11) {
    return (el_pt[ll_i1[0]] > el_pt[ll_i2[0]]) ? el_eta[ll_i2[0]] : el_eta[ll_i1[0]];
  }
  return (mu_pt[ll_i1[0]] > mu_pt[ll_i2[0]]) ? mu_eta[ll_i2[0]] : mu_eta[ll_i1[0]];
}

float get_dphi(float phi1, float phi2){
  const double PI = 3.1415;
  double dphi = fmod(fabs(phi2-phi1), 2.*PI);
    dphi = dphi>PI ? 2.*PI-dphi : dphi;

  return dphi;
}

float jj_dphi(RVec<float> jet_phi){ return get_dphi(jet_phi[0],jet_phi[1]); }



float jj_dphi(RVec<float> mc_pt, RVec<float> mc_eta, RVec<float> mc_phi,RVec<float> mc_m, RVec<float> mc_id, RVec<float> mc_momidx, 
              RVec<float> jet_pt, RVec<float> jet_eta, RVec<float> jet_phi, RVec<float> jet_m, RVec<bool> jet_isgood){ 
    RVec<int> jet_idx = truth_matched(mc_pt,mc_eta,mc_phi,mc_m,mc_id,mc_momidx,jet_pt,jet_eta,jet_phi,jet_m,jet_isgood);
    return get_dphi(jet_phi[jet_idx[0]],jet_phi[jet_idx[1]]);
}



float get_Zy_dphi(RVec<float> ll_phi, RVec<float> photon_phi){ return get_dphi(ll_phi[0],photon_phi[0]);  }

float get_Zyjj_dphi(RVec<float> ll_pt,RVec<float> ll_eta,RVec<float> ll_phi,RVec<float> ll_m, 
                RVec<float> photon_pt, RVec<float> photon_eta, RVec<float> photon_phi,
                float dijet_pt, float dijet_eta, float dijet_phi, float dijet_m){

  TLorentzVector y,Z,jj;
  y.SetPtEtaPhiM(photon_pt[0],photon_eta[0],photon_phi[0],0);
  Z.SetPtEtaPhiM(ll_pt[0],ll_eta[0],ll_phi[0],ll_m[0]);
  jj.SetPtEtaPhiM(dijet_pt,dijet_eta,dijet_phi,dijet_m);
  return get_dphi((Z+y).Phi(),jj.Phi());  
               
}


float get_Zyjj_dphi(RVec<float> mc_pt, RVec<float> mc_eta, RVec<float> mc_phi,RVec<float> mc_m, RVec<float> mc_id, RVec<float> mc_momidx,
                    RVec<float> ll_pt,RVec<float> ll_eta,RVec<float> ll_phi,RVec<float> ll_m, 
                    RVec<float> photon_pt, RVec<float> photon_eta, RVec<float> photon_phi,
                    RVec<float> jet_pt, RVec<float> jet_eta, RVec<float> jet_phi, RVec<float> jet_m, RVec<bool> jet_isgood){

  RVec<int> jet_idx = truth_matched(mc_pt,mc_eta,mc_phi,mc_m,mc_id,mc_momidx,jet_pt,jet_eta,jet_phi,jet_m,jet_isgood);
  TLorentzVector y,Z,jj;
  jj.SetPtEtaPhiM(jet_pt[jet_idx[0]],jet_eta[jet_idx[0]],jet_phi[jet_idx[0]],jet_m[jet_idx[0]]);
  Z.SetPtEtaPhiM(jet_pt[jet_idx[1]],jet_eta[jet_idx[1]],jet_phi[jet_idx[1]],jet_m[jet_idx[1]]);
  jj = Z+jj;

  y.SetPtEtaPhiM(photon_pt[0],photon_eta[0],photon_phi[0],0);
  Z.SetPtEtaPhiM(ll_pt[0],ll_eta[0],ll_phi[0],ll_m[0]);
  return get_dphi((Z+y).Phi(),jj.Phi());  
 
 }


float get_Zyjj_dr(RVec<float> ll_pt,RVec<float> ll_eta,RVec<float> ll_phi,RVec<float> ll_m, 
                RVec<float> photon_pt, RVec<float> photon_eta, RVec<float> photon_phi,
                float dijet_pt, float dijet_eta, float dijet_phi, float dijet_m){

  TLorentzVector y,Z,jj;
  y.SetPtEtaPhiM(photon_pt[0],photon_eta[0],photon_phi[0],0);
  Z.SetPtEtaPhiM(ll_pt[0],ll_eta[0],ll_phi[0],ll_m[0]);
  jj.SetPtEtaPhiM(dijet_pt,dijet_eta,dijet_phi,dijet_m);
  return (Z+y).DeltaR(jj);  
               
}

float get_Zyjj_dr(RVec<float> mc_pt, RVec<float> mc_eta, RVec<float> mc_phi,RVec<float> mc_m, RVec<float> mc_id, RVec<float> mc_momidx,
      RVec<float> ll_pt,RVec<float> ll_eta,RVec<float> ll_phi,RVec<float> ll_m, 
      RVec<float> photon_pt, RVec<float> photon_eta, RVec<float> photon_phi,
      RVec<float> jet_pt, RVec<float> jet_eta, RVec<float> jet_phi, RVec<float> jet_m, RVec<bool> jet_isgood){


      RVec<int> jet_idx = truth_matched(mc_pt,mc_eta,mc_phi,mc_m,mc_id,mc_momidx,jet_pt,jet_eta,jet_phi,jet_m,jet_isgood);
      TLorentzVector y,Z,jj;
      jj.SetPtEtaPhiM(jet_pt[jet_idx[0]],jet_eta[jet_idx[0]],jet_phi[jet_idx[0]],jet_m[jet_idx[0]]);
      Z.SetPtEtaPhiM(jet_pt[jet_idx[1]],jet_eta[jet_idx[1]],jet_phi[jet_idx[1]],jet_m[jet_idx[1]]);
      jj = Z+jj;

      y.SetPtEtaPhiM(photon_pt[0],photon_eta[0],photon_phi[0],0);
      Z.SetPtEtaPhiM(ll_pt[0],ll_eta[0],ll_phi[0],ll_m[0]);
      return (Z+y).DeltaR(jj);  

}

float get_deta(float eta1, float eta2){
  return fabs(eta1-eta2);
}

float jj_deta(RVec<float> jet_eta){ return get_deta(jet_eta[0],jet_eta[1]); }


float jj_deta(RVec<float> mc_pt, RVec<float> mc_eta, RVec<float> mc_phi,RVec<float> mc_m, RVec<float> mc_id, RVec<float> mc_momidx, 
              RVec<float> jet_pt, RVec<float> jet_eta, RVec<float> jet_phi, RVec<float> jet_m, RVec<bool> jet_isgood){ 

    RVec<int> jet_idx = truth_matched(mc_pt,mc_eta,mc_phi,mc_m,mc_id,mc_momidx,jet_pt,jet_eta,jet_phi,jet_m,jet_isgood);
    return get_deta(jet_eta[jet_idx[0]],jet_eta[jet_idx[1]]);
}



float zep_var(RVec<float> ll_pt,RVec<float> ll_eta,RVec<float> ll_phi, RVec<float> ll_m,
              RVec<float> jet_pt,RVec<float> jet_eta,RVec<float> jet_phi, RVec<float> jet_m,
              RVec<float> photon_pt,RVec<float> photon_eta,RVec<float> photon_phi){
 
  TLorentzVector y,Z,j1,j2;
  y.SetPtEtaPhiM(photon_pt[0],photon_eta[0],photon_phi[0],0);
  Z.SetPtEtaPhiM(ll_pt[0],ll_eta[0],ll_phi[0],ll_m[0]);
  j1.SetPtEtaPhiM(jet_pt[0],jet_eta[0],jet_phi[0],jet_m[0]);
  j2.SetPtEtaPhiM(jet_pt[1],jet_eta[1],jet_phi[1],jet_m[1]);
  //return get_deta(photon_eta[0],0.5*(get_deta(jet_eta[0],-1.0*jet_eta[1])) );
  return fabs(((y+Z).Rapidity() - 0.5*(j1.Rapidity() + j2.Rapidity()))/fabs(j1.Rapidity() - j2.Rapidity() ));

}

float zep_var(RVec<float> mc_pt, RVec<float> mc_eta, RVec<float> mc_phi,RVec<float> mc_m, RVec<float> mc_id, RVec<float> mc_momidx, 
              RVec<float> ll_pt,RVec<float> ll_eta,RVec<float> ll_phi, RVec<float> ll_m,
              RVec<float> jet_pt,RVec<float> jet_eta,RVec<float> jet_phi, RVec<float> jet_m, RVec<bool> jet_isgood,
              RVec<float> photon_pt,RVec<float> photon_eta,RVec<float> photon_phi){

  RVec<int> jet_idx = truth_matched(mc_pt,mc_eta,mc_phi,mc_m,mc_id,mc_momidx,jet_pt,jet_eta,jet_phi,jet_m,jet_isgood);
  TLorentzVector y,Z,j1,j2;
  y.SetPtEtaPhiM(photon_pt[0],photon_eta[0],photon_phi[0],0);
  Z.SetPtEtaPhiM(ll_pt[0],ll_eta[0],ll_phi[0],ll_m[0]);
  j1.SetPtEtaPhiM(jet_pt[jet_idx[0]],jet_eta[jet_idx[0]],jet_phi[jet_idx[0]],jet_m[jet_idx[0]]);
  j2.SetPtEtaPhiM(jet_pt[jet_idx[1]],jet_eta[jet_idx[1]],jet_phi[jet_idx[1]],jet_m[jet_idx[1]]);
  //return get_deta(photon_eta[0],0.5*(get_deta(jet_eta[0],-1.0*jet_eta[1])) );
  return fabs(((y+Z).Rapidity() - 0.5*(j1.Rapidity() + j2.Rapidity()))/fabs(j1.Rapidity() - j2.Rapidity() ));

}


float drmax_yj(RVec<float> photon_eta, RVec<float> photon_phi,RVec<float> jet_eta,RVec<float> jet_phi){
    float dr1,dr2;
    dr1 = get_dr(photon_eta[0],photon_phi[0],jet_eta[0],jet_phi[0]);
    dr2 = get_dr(photon_eta[0],photon_phi[0],jet_eta[1],jet_phi[0]);
    return dr1 > dr2 ? dr1 : dr2;
}


float drmax_yj(RVec<float> mc_pt, RVec<float> mc_eta, RVec<float> mc_phi,RVec<float> mc_m, RVec<float> mc_id, RVec<float> mc_momidx, 
              RVec<float> jet_pt, RVec<float> jet_eta, RVec<float> jet_phi, RVec<float> jet_m, RVec<bool> jet_isgood,RVec<float> photon_eta, RVec<float> photon_phi){ 
  RVec<int> jet_idx = truth_matched(mc_pt,mc_eta,mc_phi,mc_m,mc_id,mc_momidx,jet_pt,jet_eta,jet_phi,jet_m,jet_isgood);
  float dr1,dr2;
  dr1 = get_dr(photon_eta[0],photon_phi[0],jet_eta[jet_idx[0]],jet_phi[jet_idx[0]]);
  dr2 = get_dr(photon_eta[0],photon_phi[0],jet_eta[jet_idx[1]],jet_phi[jet_idx[1]]);
  return dr1 > dr2 ? dr1 : dr2;
}

float drmin_yj(RVec<float> photon_eta, RVec<float> photon_phi,RVec<float> jet_eta, RVec<float> jet_phi){
    float dr1,dr2;
    dr1 = get_dr(photon_eta[0],photon_phi[0],jet_eta[0],jet_phi[0]);
    dr2 = get_dr(photon_eta[0],photon_phi[0],jet_eta[1],jet_phi[1]);
    return dr1 < dr2 ? dr1 : dr2;
}

float drmin_yj(RVec<float> mc_pt, RVec<float> mc_eta, RVec<float> mc_phi,RVec<float> mc_m, RVec<float> mc_id, RVec<float> mc_momidx, 
              RVec<float> jet_pt, RVec<float> jet_eta, RVec<float> jet_phi, RVec<float> jet_m, RVec<bool> jet_isgood,RVec<float> photon_eta, RVec<float> photon_phi){ 
  RVec<int> jet_idx = truth_matched(mc_pt,mc_eta,mc_phi,mc_m,mc_id,mc_momidx,jet_pt,jet_eta,jet_phi,jet_m,jet_isgood);
  float dr1,dr2;
  dr1 = get_dr(photon_eta[0],photon_phi[0],jet_eta[jet_idx[0]],jet_phi[jet_idx[0]]);
  dr2 = get_dr(photon_eta[0],photon_phi[0],jet_eta[jet_idx[1]],jet_phi[jet_idx[1]]);
  return dr1 < dr2 ? dr1 : dr2;
}

float drmin_l1j(RVec<int> ll_lepid, RVec<float> el_eta, RVec<float> el_phi,RVec<float> mu_eta, RVec<float> mu_phi, RVec<float> jet_eta, RVec<float> jet_phi){
    float dr1,dr2;
    if(ll_lepid[0]==11){ 
      dr1 = get_dr(el_eta[0],el_phi[0],jet_eta[0],jet_phi[0]); 
      dr2 = get_dr(el_eta[0],el_phi[0],jet_eta[1],jet_phi[1]);
    } else {
      dr1 = get_dr(mu_eta[0],mu_phi[0],jet_eta[0],jet_phi[0]); 
      dr2 = get_dr(mu_eta[0],mu_phi[0],jet_eta[1],jet_phi[1]);
    }
    return dr1 < dr2 ? dr1 : dr2;
}

float drmin_l2j(RVec<int> ll_lepid, RVec<float> el_eta, RVec<float> el_phi,RVec<float> mu_eta, RVec<float> mu_phi, RVec<float> jet_eta, RVec<float> jet_phi){
    float dr1,dr2;
    if(ll_lepid[0]==11){ 
      dr1 = get_dr(el_eta[1],el_phi[1],jet_eta[0],jet_phi[0]); 
      dr2 = get_dr(el_eta[1],el_phi[1],jet_eta[1],jet_phi[1]);
    } else {
      dr1 = get_dr(mu_eta[1],mu_phi[1],jet_eta[0],jet_phi[0]); 
      dr2 = get_dr(mu_eta[1],mu_phi[1],jet_eta[1],jet_phi[1]);
    }
    
    return dr1 < dr2 ? dr1 : dr2;
}

float dr_jj(RVec<float> jet_eta, RVec<float> jet_phi){
    float dr1;
    return get_dr(jet_eta[0],jet_phi[0],jet_eta[1],jet_phi[1]);
}

float dr_jj(RVec<float> mc_pt, RVec<float> mc_eta, RVec<float> mc_phi,RVec<float> mc_m, RVec<float> mc_id, RVec<float> mc_momidx, 
            RVec<float> jet_pt, RVec<float> jet_eta, RVec<float> jet_phi, RVec<float> jet_m, RVec<bool> jet_isgood){ 
    RVec<int> jet_idx = truth_matched(mc_pt,mc_eta,mc_phi,mc_m,mc_id,mc_momidx,jet_pt,jet_eta,jet_phi,jet_m,jet_isgood);
    return get_dr(jet_eta[jet_idx[0]],jet_phi[jet_idx[0]],jet_eta[jet_idx[1]],jet_phi[jet_idx[1]]);
}


float pt_bal_func(RVec<float> ll_pt, RVec<float> ll_eta, RVec<float> ll_phi,RVec<float> ll_m, 
               RVec<float> jet_pt,RVec<float> jet_eta, RVec<float> jet_phi, RVec<float> jet_m, 
               RVec<float> photon_pt, RVec<float> photon_eta, RVec<float> photon_phi){

  TLorentzVector Z,y,j1,j2;
  Z.SetPtEtaPhiM(ll_pt[0],ll_eta[0],ll_phi[0],ll_m[0]);
  y.SetPtEtaPhiM(photon_pt[0],photon_eta[0],photon_phi[0],0.0);
  j1.SetPtEtaPhiM(jet_pt[0],jet_eta[0],jet_phi[0],jet_m[0]);
  j2.SetPtEtaPhiM(jet_pt[1],jet_eta[1],jet_phi[1],jet_m[1]);
  return ( (Z+y+j1+j2).Pt()/(Z.Pt() + y.Pt() + j1.Pt() + j2.Pt()) ); 
  }

float pt_bal_func(RVec<float> mc_pt, RVec<float> mc_eta, RVec<float> mc_phi,RVec<float> mc_m, RVec<float> mc_id, RVec<float> mc_momidx, 
                  RVec<float> ll_pt, RVec<float> ll_eta, RVec<float> ll_phi,RVec<float> ll_m, 
                  RVec<float> jet_pt,RVec<float> jet_eta, RVec<float> jet_phi, RVec<float> jet_m, RVec<bool> jet_isgood, 
                  RVec<float> photon_pt, RVec<float> photon_eta, RVec<float> photon_phi){
    RVec<int> jet_idx = truth_matched(mc_pt,mc_eta,mc_phi,mc_m,mc_id,mc_momidx,jet_pt,jet_eta,jet_phi,jet_m,jet_isgood);
    TLorentzVector Z,y,j1,j2;
    Z.SetPtEtaPhiM(ll_pt[0],ll_eta[0],ll_phi[0],ll_m[0]);
    y.SetPtEtaPhiM(photon_pt[0],photon_eta[0],photon_phi[0],0.0);
    j1.SetPtEtaPhiM(jet_pt[jet_idx[0]],jet_eta[jet_idx[0]],jet_phi[jet_idx[0]],jet_m[jet_idx[0]]);
    j2.SetPtEtaPhiM(jet_pt[jet_idx[1]],jet_eta[jet_idx[1]],jet_phi[jet_idx[1]],jet_m[jet_idx[1]]);
    return ( (Z+y+j1+j2).Pt()/(Z.Pt() + y.Pt() + j1.Pt() + j2.Pt()) ); 
}

float pTt_funcVec(RVec<float> ll_pt, RVec<float> ll_eta, RVec<float> ll_phi,RVec<float> ll_m, RVec<float> photon_pt, RVec<float> photon_eta, RVec<float> photon_phi){
    TLorentzVector Z,y;
    Z.SetPtEtaPhiM(ll_pt[0],ll_eta[0],ll_phi[0],ll_m[0]);
    y.SetPtEtaPhiM(photon_pt[0],photon_eta[0],photon_phi[0],0.0);
    return fabs( (Z+y).Vect().Cross( (Z-y).Vect() ).Pz()/( (Z+y).Pt() ) );
}

float get_weight(float w_lumi ,float w_year){//RVec<float> llphoton_l1_masserr,RVec<float> llphoton_l2_masserr,RVec<float> llphoton_ph_masserr, float weight, bool isNotSig){
  float dm = 1.0;
  /*if(isNotSig){dm=1;} else {
  float dml1,dml2,dmph;
  dml1 = llphoton_l1_masserr[0];
  dml2 = llphoton_l2_masserr[0];
  dmph = llphoton_ph_masserr[0];
  dm = sqrt(dml1 * dml1 + dml2 * dml2 + dmph * dmph);
  }*/

  //return weight;
  //if(SampleType == 2016) { return w_lumi*36.32264/dm;}
  //if(SampleType == 2017) { return w_lumi*41.52756/dm;}
  //if(SampleType == 2018) { return w_lumi*59.67377/dm;}
  //return weight/dm;
  return w_year*w_lumi;
  }



bool signal_lead_electron_pt(RVec<float> el_pt, RVec<float> el_sig){
  for (unsigned iPart = 0; iPart<el_pt.size(); iPart++) {
    if (el_sig.at(iPart)) {
      return (el_pt.at(iPart) > 25);
    }
  }
  return false;
}

bool signal_sublead_electron_pt(RVec<float> el_pt, RVec<float> el_sig){
  bool sublead = false;
  for (unsigned iPart = 0; iPart<el_pt.size(); iPart++) {
    if (el_sig.at(iPart)) {
      if (sublead == false) sublead = true;
      else return (el_pt.at(iPart) > 15);
    }
  }
  return false;
}


bool signal_lead_muon_pt(RVec<float> mu_pt, RVec<float> mu_sig){
  for (unsigned iPart = 0; iPart<mu_pt.size(); iPart++) {
    if (mu_sig.at(iPart)) {
      return (mu_pt.at(iPart) > 20);
    }
  }
  return false;
}

bool signal_sublead_muon_pt(RVec<float> mu_pt, RVec<float> mu_sig){
  bool sublead = false;
  for (unsigned iPart = 0; iPart<mu_pt.size(); iPart++) {
    if (mu_sig.at(iPart)) {
      if (sublead == false) sublead = true;
        else return (mu_pt.at(iPart) > 10);
    }
  }
  return false;
}

double j_pt(RVec<float>jet_pt,RVec<float> jet_isgood,int which=0){
  int count=0;
  for(unsigned iPart = 0; iPart<jet_pt.size(); iPart++){
    if(!jet_isgood.at(iPart)){continue;}
    if(count<which){count++;continue;}
      return jet_pt.at(iPart);
  }
  return -1.0;
}

""")

if __name__=='__main__':
  ROOT.EnableImplicitMT()
  #parameters
  #TODO: move generic function to /lib/

  w_year = [19.5,16.8,41.52756,59.67377]
  year = ["2016","2016APV","2017","2018"]
                
  filters_year = ['1','1','!((jet_pt[0]<50&&fabs(jet_eta[0])<3.139&&fabs(jet_eta[0])>2.65)||(jet_pt[1]<50&&fabs(jet_eta[0])<3.139&&fabs(jet_eta[0])>2.65))','1']

  for sel in range(len(year)):
    defines = [    
              #KINEMATIC BDT VARIABLES   
             ('photon_mva','photon_idmva[0]'),
             ('min_dR','photon_drmin[0]'),
             ('max_dR','get_max_dr(photon_eta,photon_phi,el_eta,el_phi,mu_eta,mu_phi,ll_lepid,ll_i1,ll_i2)'),
             ('pt_mass','llphoton_pt[0]/llphoton_m[0]'),
             ('cosTheta','llphoton_cosTheta[0]'),
             ('costheta','llphoton_costheta[0]'),
             ('phi','llphoton_psi[0]'),
             ('photon_res','photon_energyErr[0]/photon_pt[0]'),
             ('photon_prap','photon_eta[0]'),
             ('y_pt','photon_pt[0]'),      
             ('l1_rapidity','get_l1_rapidity(el_pt,el_eta,mu_pt,mu_eta,ll_lepid,ll_i1,ll_i2)'),
             ('l2_rapidity','get_l2_rapidity(el_pt,el_eta,mu_pt,mu_eta,ll_lepid,ll_i1,ll_i2)'),

             #DIJET BDT VARIABLEs
             ('j1_pt','j_pt(jet_pt,jet_isgood,0)'),
             ('j2_pt','j_pt(jet_pt,jet_isgood,1)'),
             ('jj_dphi','dijet_dphi'),
             ('jj_deta','dijet_eta'),
             ('Zy_dphi','llphoton_dphi[0]'),
             ('Zyjj_dphi','llphoton_dijet_dphi[0]'),
             ('Zyjj_dr','get_Zyjj_dr(ll_pt, ll_eta, ll_phi, ll_m, photon_pt, photon_eta, photon_phi, dijet_pt, dijet_eta, dijet_phi, dijet_m)'),
             ('zep_var','photon_zeppenfeld[0]'),
             ('pTt','llphoton_pTt2[0]'),
             ('pt_bal','llphoton_dijet_balance[0]'),
             ('mjj','dijet_m'),

             #MISC. VARIABLES
             ('mlly','llphoton_m[0]'),
             ('drmin_yj','photon_jet_mindr[0]'),
             ('drmax_yj','drmax_yj(mc_pt,mc_eta,mc_phi,mc_mass,mc_id,mc_momidx,jet_pt,jet_eta,jet_phi,jet_m,jet_isgood,photon_eta,photon_phi)'),
             ('dr_jj','dr_jj(mc_pt,mc_eta,mc_phi,mc_mass,mc_id,mc_momidx,jet_pt,jet_eta,jet_phi,jet_m,jet_isgood)'),
             ('drmin_l1j','drmin_l1j(ll_lepid,el_eta,el_phi,mu_eta,mu_phi,jet_eta,jet_phi)'),
             ('drmin_l2j','drmin_l2j(ll_lepid,el_eta,el_phi,mu_eta,mu_phi,jet_eta,jet_phi)'),
             ('el_lead_pt'   ,'signal_lead_electron_pt(el_pt,el_sig)'),
             ('el_sublead_pt','signal_sublead_electron_pt(el_pt,el_sig)'),
             ('mu_lead_pt'   ,'signal_lead_muon_pt(mu_pt,mu_sig)'),
             ('mu_sublead_pt','signal_sublead_muon_pt(mu_pt,mu_sig)'),
             ('binning_bdt','getMVA(jj_deta,jj_dphi,Zyjj_dr,pt_bal,zep_var,mjj,drmin_yj,j1_pt,j2_pt,event)'),
             ('tm_jets','truth_matched_bool(mc_pt,mc_eta,mc_phi,mc_mass,mc_id,mc_momidx,jet_pt,jet_eta,jet_phi,jet_m,jet_isgood)'),
             ('event_number','event'),
             ('wgt','get_weight(w_lumi,' + str(w_year[sel]) + ')')]
   
    branches = ('photon_mva','min_dR','max_dR','pt_mass','cosTheta','costheta', 'phi','photon_res','photon_prap','y_pt','l1_rapidity','l2_rapidity',
               'j1_pt','j2_pt','jj_deta','jj_dphi','Zyjj_dphi','Zyjj_dr',"dr_jj",'pTt','Zy_dphi','drmin_yj','drmax_yj',"drmin_l1j","drmin_l2j",'pt_bal','zep_var',
               'mlly','mjj','binning_bdt','event_number','weight')

#             ('kinMVA','getMVA(photon_mva,min_dR,max_dR,pt_mass,cosTheta,costheta,phi,photon_res,photon_prap,l1_rapidity,l2_rapidity)'),

    #make n-tuples
    cuts = ['trig_double_el || trig_double_mu', 'llphoton_m.size()>0 && photon_pt.size()>0','(mu_lead_pt && mu_sublead_pt)||(el_lead_pt && el_sublead_pt)',
            'use_event', '(photon_drmin[0]>0.4)&&(photon_id80[0])&&(photon_pt[0]/mlly>=15.0/110.0)',
            '(ll_m[0]>80 && ll_m[0]<100)&&(photon_elveto[0])',
            'mlly>100 && mlly<180', 'njet>1&&mjj>0.0&&nlep==2&&nbdfm==0']
          


    names = 'photon_bdt_ntuples' # notm'
    base_dir  = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'
    pico_type = '/mc/merged_zgmc_llg/'
   

    sig_samples_VBF = ['*VBF*HToZG*125_*.root']                              
    sig_samples_ggF = ['*GluGlu*HToZG*125_*.root']
    sig_samples_VH = ['*WplusH_HToZG_WToAll_M-125*.root','*WminusH_HToZG_WToAll_M-125*.root','*ZH_HToZG_ZToAll*M-125*.root']                              

    bkg_samples_tt = ['*_TGJets*.root','*_TTGJets*.root','*_TTTo2L2Nu*.root']                  
    bkg_samples_DY = ['*DYJetsToLL*amcatnloFXFX*.root','*ZGToLLG*.root']                 
    bkg_samples_EWK= ['*ZGamma2JToGamma2L2J_EWK*.root']                  

    #Make signal ntuples ggF
    write_ntuples([base_dir + year[sel] + pico_type + sig for sig in sig_samples_ggF],  cuts,  
                   names + '_sig_ggF_' + year[sel] + '.root',  defines,  'tree',  branches)

    #Make signal ntuples VBF
    write_ntuples([base_dir + year[sel] + pico_type + sig for sig in sig_samples_VBF],  cuts,  
                   names + '_sig_VBF_' + year[sel] + '.root',  defines,  'tree',  branches)

    #Make signal ntuples VH
    write_ntuples([base_dir + year[sel] + pico_type + sig for sig in sig_samples_VH],  cuts,  
                   names + '_sig_VH_' + year[sel] + '.root',  defines,  'tree',  branches)

    #Make bkg_tt ntuples               
    write_ntuples([base_dir + year[sel] + pico_type  + bkg  for bkg in bkg_samples_tt],  cuts,
                   names + '_bkg_tt' + year[sel] + '.root',  defines,  'tree',  branches)

    #Make bkg_DY ntuples               
    write_ntuples([base_dir + year[sel] + pico_type  + bkg  for bkg in bkg_samples_DY],  cuts,
                   names + '_bkg_DY' + year[sel] + '.root',  defines,  'tree',  branches)

    #Make bkg_DY ntuples               
    write_ntuples([base_dir + year[sel] + pico_type  + bkg  for bkg in bkg_samples_EWK],  cuts,
                   names + '_bkg_EWKZG' + year[sel] + '.root',  defines,  'tree',  branches)


#List of potential samples to use
#DYJetsToLL EWKZ2Jets GluGluHToGG GluGluHToTauTau GluGluHToWW GluGluHToZZ TGJets_TuneCP5 TTGJets TTTo2L2Nu W ZG ZZ 

