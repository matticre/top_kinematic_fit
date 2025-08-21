// chiudo sull'angolo non sull'impulso
   
#include <iostream>
#include <vector>
#include <string>
#include "TRandom3.h"
#include <TGenPhaseSpace.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TMath.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMinuit.h" 
#include "TLatex.h"

TRandom3 rnd;

using namespace std;

//Plotting
//Risoluzioni da ricavare
TH1D hNuResPhi("hPhiRes","Risoluzione Phi neutrino prefit",50,-4,4);        // istogramma pt del neutrino
TH1D hNuResPt("hNuResPt","Risoluzione Pt neutrino prefit",50,-50,50);       // istogramma phi del neutrino 

//Massa top
TH1D hMtopl("m_{t}","Massa del top postfit",30,80,260);        // istogramma massa top leptonica
TH1D hMtop_pre("m_{t}^{prefit}","Massa del top prefit",30,60,260);  // istogramma massa top adronica

//Massa W
TH1D hMW("m_{W}","Massa del W postfit",30,70,90);                    // istogramma massa del w dopo il fit
TH1D hMWpre("m_{W}^{prefit}","Massa del W prefit",30,-10,160);            // istogramma massa del w prefit

TH1D hThenu("hThenu","Distribuzione #theta_{#nu} postfit",30,-0.5,M_PI+0.5);       // istogramma theta del neutrino

//impulso trasverso
TH1D hNuPtSim("hNuPtSim","PtSim neutrino prefit",50,0,100);
TH1D hNuPtReal("hNuPtRec","PtRec neutrino prefit",50,0,100);

//Pull delle variabili pre e post fit - sono piantate
TH1D hPullep("hPullep","#Delta(p) elettrone",60,-6,6);         // pull impulso elettrone
TH1D hPullephi("hPullephi","#Delta(#phi) elettrone",60,-6,6);   // pull phi elettrone
TH1D hPullethe("hPullethe","#Delta(#theta) elettrone",60,-6,6); // pull theta elettrone
TH1D hPulljp("hPulljp","#Delta(p) b-jet",60,-6,6);         // pull impulso jet
TH1D hPulljphi("hPulljphi","#Delta(#phi) b-jet ",60,-6,6);  // pull phi jet
TH1D hPulljthe("hPulljthe","#Delta(#theta) b-jet",60,-6,6); // pull the jet

TH1D hThere("hThere","Pull theta postfit",60,-6,6);     // pull theta (vero-fit) del neutrino

//Componenti x e y dell'impulso del neutrino
TH1D hPnux("hPnux","p_{x}^{#nu} - p_{x}^{#nu}_{sim}",60,-30,30);  
TH1D hPnuy("hPnuy","p_{y}^{#nu} - p_{y}^{#nu}_{sim}",60,-30,30); 

//Massa W
TH1D hMWRe("MW sim","Massa del W simulata",80,50,110);  // istogramma massa del w dopo il fit
TH1D hPulltop("Pull top","Pull top",40,-6,6);  // istogramma massa del w dopo il fit


//Particelle
TDatabasePDG pdg;
TParticlePDG *top = pdg.GetParticle(6);
TParticlePDG *bot = pdg.GetParticle(5);
TParticlePDG *em  = pdg.GetParticle(11);
TParticlePDG *Wp  = pdg.GetParticle(24);
TParticlePDG *nu  = pdg.GetParticle(12);

// Parametri simulazione
const int N_EV = 1e4;   // Numero di eventi

// Array di storage
namespace data{
    vector<double> theta_nu(3);     //salvo i dati per la distribuzione di theta
}

struct reso{
    double ep   = 0.001;
    double ephi = 0.02;
    double ethe = 0.02;
    double jp   = 0.15;
    double jphi = 0.10;
    double jthe = 0.10;
    double np   = 0;    
    double nphi = 0; 
    double efit;
    double jfit;

    vector<double> data, edata;     //salvo i dati iniziali per il fit
    vector<double> topsim, topfit;
};

reso  res;

class Particle: public TLorentzVector {
    
    public:
        using TLorentzVector::TLorentzVector;
        
        Particle(TParticlePDG p, TLorentzVector t):
            TLorentzVector(t), m_code(p.PdgCode()), m_charge(p.Charge())
            {};
        
        Particle(TParticlePDG *p, TLorentzVector* t): 
            TLorentzVector(*t), m_code(p->PdgCode()), m_charge(p->Charge())
            {};

        int Charge(){
            return m_charge;
        }

        int PdgCode(){
            return m_code;
        }
        
    private:
        int m_code;
        float m_charge;
};

TLorentzVector returnV (double p, double phi, double theta, double m){
  
    double px = p*sin(theta)*cos(phi);
    double py = p*sin(theta)*sin(phi);
    double pz = p*cos(theta);
    double p0 = sqrt(px*px+py*py+pz*pz+m*m);

    TLorentzVector P(px,py,pz,p0);

    return P;
}

//funzione chi quadro
void fcn (int &npar, double *gin, double &f, double *par, int flag){
    f = 0.0;

    for (int i=0; i< res.data.size()-1; i++){
        f += pow((res.data[i]-par[i])/res.edata[i],2);
    }

    TLorentzVector P_e   = returnV(par[3], par[4], par[5], em->Mass());
    TLorentzVector P_nu  = returnV(par[6], par[7], par[8], 0);

    double M   = (P_e+P_nu).Mag();
    double eMW = Wp->Width();
    
    //vincolo cinematico
    //f += pow((Wp->Mass()-M)/eMW,2);
    f += -2*log(TMath::BreitWignerRelativistic(M,Wp->Mass(),eMW));
}

//applico la risoluzione
void ApplyResolution (Particle &track){

    double m   = sqrt(track.T()*track.T()-track.Rho()*track.Rho());
    double P   = track.Rho();
    double phi = track.Phi();
    double the = track.Theta();

    double eP, ethe, ephi;

    if(abs(track.PdgCode())!=12){

        if(abs(track.PdgCode())==11){   //elettrone
            res.efit = res.ep*P*P;

            P   = rnd.Gaus(P, res.ep*P*P);
            phi = rnd.Gaus(phi, res.ephi);
            the = rnd.Gaus(the, res.ethe);

        } else if (abs(track.PdgCode())==5) {  //quark b
            res.jfit = res.jp*P;

            P   = rnd.Gaus(P,   res.jp*P);
            phi = rnd.Gaus(phi, res.jphi);
            the = rnd.Gaus(the, res.jthe);
        }
    }

    track.SetRho(P);
    track.SetPhi(phi);
    track.SetTheta(the);
    track.SetE(sqrt(P*P+m*m));
}

void Reso(vector<Particle>& recoparticles){
    for (int i=0; i<recoparticles.size(); i++){
        ApplyResolution(recoparticles[i]);
    }
}

//simulazione evento
void Event(TGenPhaseSpace &tgps, vector<Particle> &decayparticle, double &wlep){
    
    double t_mass  = top->Mass();
    double t_width = top->Width();
    double w_mass  = Wp->Mass();
    double w_width = Wp->Width();

    double t_sim   = sqrt(t_mass*t_width*tan(M_PI*(rnd.Rndm()-0.5))+t_mass*t_mass);
    double w_sim   = sqrt(w_mass*w_width*tan(M_PI*(rnd.Rndm()-0.5))+w_mass*w_mass);
    
    double bw_masses[2]  = {w_sim, bot->Mass()}; //array di masse del primo decadimento
    double enu_masses[2] = {em->Mass(), 0};      //array di masse del secondo decadimento

    TLorentzVector t1_in(0,0,0,t_sim);
    
    // ttbar ->  W+ b W- bbar -> e+ nu b + qqbar- bbar
    tgps.SetDecay(t1_in, 2, bw_masses);
    wlep = tgps.Generate();        
    Particle W1(Wp, tgps.GetDecay(0));
    Particle bot1(bot, tgps.GetDecay(1));

    // faccio decadere il W positivo in positrone e nu
    tgps.SetDecay(W1, 2, enu_masses);
    wlep *= tgps.Generate();
    Particle pos(em, tgps.GetDecay(0));
    Particle nu_r(nu,tgps.GetDecay(1));

    decayparticle.push_back(bot1);
    decayparticle.push_back(pos);
    decayparticle.push_back(nu_r);        
} 
    
//programma fit cinematico
vector<TLorentzVector> Kinematic_Fit(TMinuit *minuit, vector<Particle> &recoparticles, vector<Particle> &decayparticles){

    //jet
    res.data.push_back(recoparticles[0].P());
    res.data.push_back(recoparticles[0].Phi());
    res.data.push_back(recoparticles[0].Theta());
    //elettrone
    res.data.push_back(recoparticles[1].P());
    res.data.push_back(recoparticles[1].Phi());
    res.data.push_back(recoparticles[1].Theta());
    //neutrino
    res.data.push_back(recoparticles[2].P());
    res.data.push_back(recoparticles[2].Phi());
    res.data.push_back(recoparticles[2].Theta());

    res.edata.push_back(res.jfit);
    res.edata.push_back(res.jphi);
    res.edata.push_back(res.jthe);

    res.edata.push_back(res.efit);
    res.edata.push_back(res.ephi);
    res.edata.push_back(res.ethe);
    
    res.edata.push_back(res.np);
    res.edata.push_back(res.nphi);

    minuit->DefineParameter(0,"bP"  ,res.data[0],0.01,0.,0.);
    minuit->DefineParameter(1,"bPhi",res.data[1],0.01,0.,0.);
    minuit->DefineParameter(2,"bThe",res.data[2],0.01,0.,0.);
    
    minuit->DefineParameter(3,"eP"  ,res.data[3],0.01,0.,0.);
    minuit->DefineParameter(4,"ePhi",res.data[4],0.01,0.,0.);
    minuit->DefineParameter(5,"eThe",res.data[5],0.01,0.,0.);
    
    minuit->DefineParameter(6,"nuP"  ,res.data[6],0.01,0.,0.);
    minuit->DefineParameter(7,"nuPhi",res.data[7],0.01,0.,0.);
    minuit->DefineParameter(8,"nuThe",res.data[8],0.01,0.,0.);

    minuit->Command("MIGRAD");

    double pb,  epb,  phib,  ephib,  theb,  etheb;  //variabili b       
    double pe,  epe,  phie,  ephie,  thee,  ethee;  //variabili e    
    double pnu, epnu, phinu, ephinu, thenu, ethenu; //variabili nu

    minuit->GetParameter(0,pb,epb);
    minuit->GetParameter(1,phib,ephib);
    minuit->GetParameter(2,theb,etheb);

    minuit->GetParameter(3,pe,epe);
    minuit->GetParameter(4,phie,ephie);
    minuit->GetParameter(5,thee,ethee);

    minuit->GetParameter(6,pnu,epnu);
    minuit->GetParameter(7,phinu,ephinu);
    minuit->GetParameter(8,thenu,ethenu);

    TLorentzVector P_b  = returnV(pb,phib,theb,bot->Mass()); //b-jet
    TLorentzVector P_e  = returnV(pe,phie,thee,em->Mass());  //elettrone/positrone
    TLorentzVector P_nu = returnV(pnu,phinu,thenu,0);        //neutrino    

    //Plotting dei pull
    hPulljp.Fill((recoparticles[0].P()-pb)/epb);            // pull impulso jet
    hPulljphi.Fill((recoparticles[0].Phi()-phib)/ephib);    // pull phi jet
    hPulljthe.Fill((recoparticles[0].Theta()-theb)/etheb);  // pull the jet

    hPullep.Fill((recoparticles[1].P()-pe)/epe);            // pull impulso elettrone
    hPullephi.Fill((recoparticles[1].Phi()-phie)/ephie);    // pull phi elettrone
    hPullethe.Fill((recoparticles[1].Theta()-thee)/ethee);  // pull theta elettrone

    hThere.Fill((recoparticles[2].Theta()-thenu)/ethenu);        
    
    return {P_b, P_e, P_nu};
}

int main(){
    TApplication app("app",NULL,NULL);

    //Setting generators
    TGenPhaseSpace tgps;    
    rnd.SetSeed(12345);
    
    //Inizializzo minuit
    TMinuit *minuit = new TMinuit(9);
    minuit->SetPrintLevel(-1); 
    minuit->SetFCN(fcn);

    vector<Particle> decayparticles, recoparticles;

    //calcolo l'ampiezza le distribuzioni di phi_nu e p_nu
    for (int i=0; i<N_EV; i++){   
        if (i%10000==0){
            cout << i << endl;
        }
        double wlep;
        Event(tgps, decayparticles, wlep);
        recoparticles = decayparticles;
        Reso(recoparticles);    
        TLorentzVector nu_r = -recoparticles[0]-recoparticles[1];
    
        hNuResPhi.Fill(nu_r.Phi()-decayparticles[2].Phi());
        hNuResPt.Fill(nu_r.Pt()-decayparticles[2].Pt());
        
        hPnux.Fill(nu_r.X()-decayparticles[2].X());
        hPnuy.Fill(nu_r.Y()-decayparticles[2].Y());
        
        hNuPtSim.Fill(decayparticles[2].Pt());
        hNuPtReal.Fill(nu_r.Pt());

        decayparticles.clear();
        recoparticles.clear();
    }

    res.np   = hNuResPt.GetStdDev();
    res.nphi = hNuResPhi.GetStdDev();

    for (int i=0; i<N_EV; i++){   

        if(i%10000==0){
            cout << i << endl;
        }

        double wlep; //peso relativo al decadimento
        Event(tgps, decayparticles, wlep);
        recoparticles = decayparticles;
        Reso(recoparticles);
        hMWRe.Fill((recoparticles[1]+recoparticles[2]).Mag());

        res.topsim.push_back((decayparticles[0]+decayparticles[1]+decayparticles[2]).Mag());

        //calcolo il quadrivettore del neutrino con le informazioni che ho
        TLorentzVector p_nu  = -recoparticles[0]-recoparticles[1];
        double pt_nu  = p_nu.Pt();
        double phi_nu = p_nu.Phi();
        TLorentzVector P_nu_pre = returnV(pt_nu, phi_nu, decayparticles[2].Theta(), 0.0);
        
        double mtop_pre = (P_nu_pre+recoparticles[0]+recoparticles[1]).Mag();
        double mw_pre   = (P_nu_pre+recoparticles[1]).Mag();

        //eseguo il fit cinematico
        vector<TLorentzVector> vettori = Kinematic_Fit(minuit, recoparticles, decayparticles);
        double mtoplep = (vettori[0]+vettori[1]+vettori[2]).Mag();
        double mw      = (vettori[1]+vettori[2]).Mag();  

        res.topfit.push_back(mtoplep);

        hMtop_pre.Fill(mtop_pre,wlep);
        hMWpre.Fill(mw_pre,wlep);
        hMW.Fill(mw,wlep);
        hMtopl.Fill(mtoplep,wlep);
        hThenu.Fill(vettori[2].Theta(),wlep);

        decayparticles.clear();
        recoparticles.clear(); 
        res.data.clear(); res.edata.clear(); 
    }
    
    cout << "La massa del top è: " << top->Mass() << endl;
    cout << "La massa del top misurata è: " << hMtopl.GetMean() << " +- " << hMtopl.GetMeanError() << endl;




    //Plotting

    TCanvas c2("c2","",5,5,1500,900);
    c2.Divide(2,2);
    c2.cd(1); hNuResPt.Draw();
    c2.cd(2); hNuResPhi.Draw();
    c2.cd(3); hMWpre.Draw();
    c2.cd(4); hThenu.Draw();

    TCanvas c1("c1","",5,5,1500,900);
    c1.Divide(2,2);
    c1.cd(1); hMtopl.Draw();
    c1.cd(2); hMW.Draw();
    c1.cd(3); hMtop_pre.Draw();
    c1.cd(4); hThere.Draw();

    TCanvas cnu("cnu","",5,5,1500,900);
    cnu.Divide(2,1);
    cnu.cd(1);  hPnux.Draw();   hPnux.SetXTitle("#Delta(p_{x}) (GeV)");
    cnu.cd(2);  hPnuy.Draw();   hPnuy.SetXTitle("#Delta(p_{y}) (GeV)");

    TCanvas theta("theta","",5,5,1500,9000);
    hThere.Draw();  

/*     TCanvas massa_w("massa_w","Massa W",5,5,1500,900);
    massa_w.Divide(2,1);
    massa_w.cd(1); hMWpre.Draw();   hMWpre.SetXTitle("m_{W} (GeV)");
    massa_w.cd(2); hMW.Draw();      hMW.SetXTitle("m_{W} (GeV)");

    TCanvas massa_top("massa_top","Massa top",5,5,1500,900);
    massa_top.Divide(2,1);
    massa_top.cd(1); hMtop_pre.Draw(); hMtop_pre.SetXTitle("m_{t} (GeV)");
    massa_top.cd(2); hMtopl.Draw();    hMtopl.SetXTitle("m_{t} (GeV)"); */

/*     TCanvas pull("Pull","Pull",5,5,1500,900);
    pull.Divide(3,2);
    
    pull.cd(1); hPullep.Draw();     hPullep.SetXTitle("(|p_{e}^{fit}|-|p_{e}^{prefit}|)/#sigma_{fit} (GeV)");
    pull.cd(2); hPullephi.Draw();   hPullephi.SetXTitle("(#phi_{e}^{fit}-#phi_{e}^{prefit})/#sigma_{fit}");
    pull.cd(3); hPullethe.Draw();   hPullethe.SetXTitle("(#theta_{e}^{fit}-#theta_{e}^{prefit})/#sigma_{fit}");
    pull.cd(4); hPulljp.Draw();     hPulljp.SetXTitle("(|p_{b}^{fit}|-|p_{b}^{prefit}|)/#sigma_{fit} (GeV)");
    pull.cd(5); hPulljphi.Draw();   hPulljphi.SetXTitle("(#phi_{b}^{fit}-#phi_{b}^{prefit})/#sigma_{fit}");
    pull.cd(6); hPulljthe.Draw();   hPulljthe.SetXTitle("(#theta_{e}^{fit}-#theta_{e}^{prefit})/#sigma_{fit}"); */

/*     TCanvas reso("Reso","Risoluzioni",5,5,1500,900);
    reso.Divide(2,1);
    reso.cd(1); hNuResPhi.Draw(); hNuResPhi.SetXTitle("#phi");
    reso.cd(2); hNuResPt.Draw();  hNuResPt.SetXTitle("P_{T} (GeV)"); */
/* 
    TCanvas massa_w_chi("massa_w_chi","Massa W",5,5,1500,900);
    hMW.Draw();      hMW.SetXTitle("m_{W} (GeV)");

 */

    for (int i=0; i<res.topsim.size(); i++){
        hPulltop.Fill((res.topfit[i]-res.topsim[i])/hMtopl.GetStdDev());
    }

    TCanvas c6;
    hPulltop.Draw(); hPulltop.SetXTitle("(m_{t}^{fit}-m_{t}^{sim})/#sigma");

    app.Run(true);

    return 0;
}
