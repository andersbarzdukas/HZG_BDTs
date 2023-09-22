#!/usr/bin/env python3
'''
Function that generates n-tuples for MVA training
'''
import ROOT

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
  for filename in filenames:
    filenames_vec.push_back(filename)
  df = ROOT.RDataFrame('tree',filenames_vec)
  for define in defines:
    df = df.Define(define[0],define[1])
  for cut in cuts:
    df = df.Filter(cut)
  if (branches == ()):
    df.Snapshot(tree_name,'./ntuples/'+out_name)
  else:
    df.Snapshot(tree_name,'./ntuples/'+out_name,branches)
  print('Wrote test_ntuples/'+out_name)

ROOT.gInterpreter.Declare("""
template <class C>
using RVec = ROOT::VecOps::RVec<C>;


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

float H_t(RVec<float> photon_pt,RVec<float> photon_eta,RVec<float> photon_phi,
          RVec<float> el_pt,RVec<float> el_eta,RVec<float> el_phi,
          RVec<float> mu_pt,RVec<float> mu_eta,RVec<float> mu_phi,RVec<int> ll_lepid,RVec<int> ll_i1,RVec<int> ll_i2,
          RVec<float> jet_pt,RVec<float> jet_eta,RVec<float> jet_phi,RVec<float> jet_m){

    TLorentzVector tot,temp;
    tot.SetPtEtaPhiM(photon_pt[0],photon_eta[0],photon_phi[0],0.0);

    for(unsigned int idx = 0; idx < jet_pt.size(); idx++){
      temp.SetPtEtaPhiM(jet_pt[idx],jet_eta[idx],jet_phi[idx],jet_m[idx]);
      tot = tot + temp;
    }

    for(unsigned int idx = 0; idx < el_pt.size(); idx++){
      temp.SetPtEtaPhiM(el_pt[idx],el_eta[idx],el_phi[idx],0.000511);
      tot = tot + temp;
    }

    for(unsigned int idx = 0; idx < mu_pt.size(); idx++){
      temp.SetPtEtaPhiM(mu_pt[idx],mu_eta[idx],mu_phi[idx],0.1057);
      tot = tot + temp;
    }

    return tot.Et();

}



float get_weight(float w_lumi ,float w_year, RVec<float> llphoton_l1_masserr,RVec<float> llphoton_l2_masserr,RVec<float> llphoton_ph_masserr, bool isNotSig){
    float dm = 1.0;
    if(isNotSig){dm=1;} else {
      float dml1,dml2,dmph;
      dml1 = llphoton_l1_masserr[0];
      dml2 = llphoton_l2_masserr[0];
      dmph = llphoton_ph_masserr[0];
      dm = sqrt(dml1 * dml1 + dml2 * dml2 + dmph * dmph);
    }

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

double m_lly(RVec<float> llphoton_m, RVec<float> ll_pt,RVec<float> ll_eta,RVec<float> ll_phi,RVec<float> ll_m,RVec<float> photon_pt,RVec<float> photon_eta,RVec<float> photon_phi){
  if(llphoton_m.size()>0){return llphoton_m[0];}
  TLorentzVector ll,y;
  ll.SetPtEtaPhiM(ll_pt[0],ll_eta[0],ll_phi[0],ll_m[0]);
  y.SetPtEtaPhiM(photon_pt[0],photon_eta[0],photon_phi[0],0);
  return (ll+y).M();
}



""")

if __name__=='__main__':
  ROOT.EnableImplicitMT()
  #parameters
  #TODO: move generic function to /lib/
 
  w_year = [19.5,16.8,41.52756,59.67377]
  year = ["2016","2016APV","2017","2018"]

  triggers = ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"]

  #make n-tuples
  pt_cuts=['1']#['1','photon_pt[0]<35.0','photon_pt[0]>35.0']
  mH_cuts=['1','mlly > 120 && mlly < 130'] 
  pt_names = ['']#['_noptcut','_l35','_g35'] #['_l35','_g35']
  mH_names = ['','_mH'] 
  #eta_names = ['_noetacut'] #['_barrel','_endcap']
  #idmva_names = ['_noidmvacut'] #['_g06','_l06','_g04','_l04']
  #idmva_names = [[["_g06","_l06"],["_g04","_l04"]],[['_g04','_l04'],["_g02","_l02"]]]

  #branches = ('photon_mva','min_dR','max_dR','pt_mass','cosTheta','costheta','met_calo','met_phi','Ht','phi','photon_res','photon_eta','y_pt','y_pt_deco','l1_rapidity','l2_rapidity','wgt')

  for i in range(len(pt_names)):
    for j in range(len(mH_names)):
      for sel in range(len(year)): #[i][j])):

        defines = [('photon_mva','photon_idmva[0]'),
             ('min_dR','photon_drmin[0]'),
             ('max_dR','get_max_dr(photon_eta,photon_phi,el_eta,el_phi,mu_eta,mu_phi,ll_lepid,ll_i1,ll_i2)'),
             ('Ht','H_t(photon_pt,photon_eta,photon_phi,el_pt,el_eta,el_phi,mu_pt,mu_eta,mu_phi,ll_lepid,ll_i1,ll_i2,jet_pt,jet_eta,jet_phi,jet_m)'),
             ('pt_mass','llphoton_pt[0]/llphoton_m[0]'),
             ('cosTheta','llphoton_cosTheta[0]'),
             ('costheta','llphoton_costheta[0]'),
             ('phi','llphoton_psi[0]'),
             ('photon_res','photon_pterr[0]/photon_pt[0]'),
             ('photon_prap','photon_eta[0]'),
             ('y_pt','photon_pt[0]'),     
             ('y_pt_deco','photon_pt[0]/llphoton_m[0]'),      
             ('l1_rapidity','get_l1_rapidity(el_pt,el_eta,mu_pt,mu_eta,ll_lepid,ll_i1,ll_i2)'),
             ('l2_rapidity','get_l2_rapidity(el_pt,el_eta,mu_pt,mu_eta,ll_lepid,ll_i1,ll_i2)'),
             ('el_lead_pt'   ,'signal_lead_electron_pt(el_pt,el_sig)'),
             ('el_sublead_pt','signal_sublead_electron_pt(el_pt,el_sig)'),
             ('mu_lead_pt'   ,'signal_lead_muon_pt(mu_pt,mu_sig)'),
             ('mu_sublead_pt','signal_sublead_muon_pt(mu_pt,mu_sig)'),
             ('mlly','llphoton_m[0]'),
             ('wgt','get_weight(w_lumi,' + str(w_year[sel]) + ',llphoton_l1_masserr,llphoton_l2_masserr,llphoton_ph_masserr,true)')]
 
        branches = ('photon_mva','min_dR','max_dR','pt_mass','cosTheta','costheta','Ht','phi','photon_res','photon_prap','y_pt','y_pt_deco','l1_rapidity','l2_rapidity','wgt')
        #'met_calo','met_phi',
        #make n-tuples
        cuts = [triggers[sel], 'llphoton_m.size()>0 && photon_pt.size()>0','(mu_lead_pt && mu_sublead_pt)||(el_lead_pt && el_sublead_pt)',
            'stitch_dy||(type/1000!=6)||(type != 1600)', '(photon_drmin[0]>0.4)&&(photon_id80[0])&&(photon_pt[0]/mlly>=15.0/110.0)',
            '(ll_m[0]>80 && ll_m[0]<100)&&(photon_elveto[0])',
            'mlly>100 && mlly<180', 'njet<2 && nlep<3',mH_cuts[j],pt_cuts[i]]
          

        names = 'ggF_ntuples' + pt_names[i] + mH_names[j]
        base_dir  = '/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'
        pico_type = '/mc/merged_zgmc_llg/'
        sig_samples = ['*HToZG*125_*.root']                              
        bkg_samples = ['*DYJetsToLL*amcatnloFXFX*.root','*ZGToLLG*.root','*_TGJets*.root','*_TTGJets*.root','*_TTTo2L2Nu*.root']                  
        print([base_dir + year[sel] + pico_type + sig for sig in sig_samples])
        print([base_dir + year[sel] + pico_type  + bkg  for bkg in bkg_samples])

        #write_ntuples(filenames, cuts, out_name, defines=[], tree_name='tree', branches=())
        #Make signal ntuples
        write_ntuples([base_dir + year[sel] + pico_type + sig for sig in sig_samples],  cuts,  
                 (names + '_sig_' + year[sel] + '.root'),  defines,  'tree',  branches)

        #Make bkg ntuples               
        write_ntuples([base_dir + year[sel] + pico_type  + bkg  for bkg in bkg_samples],  cuts,
                 (names + '_bkg_' + year[sel] + '.root'),  defines,  'tree',  branches)



