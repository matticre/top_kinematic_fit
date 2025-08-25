/**
 * @file main.cpp
 * @brief This program simulates top quark decay and performs a kinematic fit
 * to determine the top quark mass.
 *
 * It uses the ROOT framework for event generation, simulation, data analysis,
 * and visualization. The core of the project is a chi-square minimization
 * to reconstruct the kinematics of the decay products, including the undetected
 * neutrino.
 */

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

/// Global random number generator.
TRandom3 rnd;

using namespace std;

//-----------------------------------------------------------------------------
// HISTOGRAMS AND PLOTTING
//-----------------------------------------------------------------------------

/// Histograms for plotting resolutions and kinematic quantities.

/// Resolution of the neutrino's phi angle before the kinematic fit.
TH1D hNuResPhi("hPhiRes","Neutrino Phi Resolution Prefit",50,-4,4);
/// Resolution of the neutrino's transverse momentum before the kinematic fit.
TH1D hNuResPt("hNuResPt","Neutrino Pt Resolution Prefit",50,-50,50);

/// Top quark mass after the kinematic fit.
TH1D hMtopl("m_{t}","Top Mass Postfit",30,80,260);
/// Top quark mass before the kinematic fit.
TH1D hMtop_pre("m_{t}^{prefit}","Top Mass Prefit",30,60,260);

/// W boson mass after the kinematic fit.
TH1D hMW("m_{W}","W Mass Postfit",30,70,90);
/// W boson mass before the kinematic fit.
TH1D hMWpre("m_{W}^{prefit}","W Mass Prefit",30,-10,160);

/// Distribution of the neutrino's polar angle after the fit.
TH1D hThenu("hThenu","Neutrino #theta Distribution Postfit",30,-0.5,M_PI+0.5);

/// Simulated neutrino transverse momentum.
TH1D hNuPtSim("hNuPtSim","Simulated Neutrino Pt Prefit",50,0,100);
/// Reconstructed neutrino transverse momentum before the fit.
TH1D hNuPtReal("hNuPtRec","Reconstructed Neutrino Pt Prefit",50,0,100);

/// Pull distributions for pre- and post-fit variables.
TH1D hPullep("hPullep","Electron #Delta(p) Pull",60,-6,6);
TH1D hPullephi("hPullephi","Electron #Delta(#phi) Pull",60,-6,6);
TH1D hPullethe("hPullethe","Electron #Delta(#theta) Pull",60,-6,6);
TH1D hPulljp("hPulljp","B-jet #Delta(p) Pull",60,-6,6);
TH1D hPulljphi("hPulljphi","B-jet #Delta(#phi) Pull ",60,-6,6);
TH1D hPulljthe("hPulljthe","B-jet #Delta(#theta) Pull",60,-6,6);

/// Pull of the neutrino's theta (true-fit).
TH1D hThere("hThere","Neutrino #theta Pull Postfit",60,-6,6);

/// Difference between reconstructed and simulated neutrino px.
TH1D hPnux("hPnux","p_{x}^{#nu} - p_{x}^{#nu}_{sim}",60,-30,30);
/// Difference between reconstructed and simulated neutrino py.
TH1D hPnuy("hPnuy","p_{y}^{#nu} - p_{y}^{#nu}_{sim}",60,-30,30);

/// Simulated W boson mass.
TH1D hMWRe("MW sim","Simulated W Mass",80,50,110);
/// Pull distribution for the top quark mass.
TH1D hPulltop("Pull top","Top Pull",40,-6,6);

//-----------------------------------------------------------------------------
// PARTICLE DEFINITIONS AND GLOBAL PARAMETERS
//-----------------------------------------------------------------------------

/// Particle database from ROOT.
TDatabasePDG pdg;
TParticlePDG *top = pdg.GetParticle(6);
TParticlePDG *bot = pdg.GetParticle(5);
TParticlePDG *em  = pdg.GetParticle(11);
TParticlePDG *Wp  = pdg.GetParticle(24);
TParticlePDG *nu  = pdg.GetParticle(12);

/// Number of events to simulate.
const int N_EV = 1e4;

/// Namespace to store data.
namespace data{
    vector<double> theta_nu(3);
}

/// Structure to hold resolution parameters and data for the fit.
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

    vector<double> data, edata;
    vector<double> topsim, topfit;
};

/// Global instance of the resolution parameters.
reso  res;

//-----------------------------------------------------------------------------
// CLASSES AND FUNCTIONS
//-----------------------------------------------------------------------------

/**
 * @class Particle
 * @brief Custom class for particles, inheriting from TLorentzVector.
 *
 * This class adds particle properties like PDG code and charge.
 */
class Particle: public TLorentzVector {

    public:
        using TLorentzVector::TLorentzVector;

        Particle(TParticlePDG p, TLorentzVector t):
            TLorentzVector(t), m_code(p.PdgCode()), m_charge(p.Charge())
            {};

        Particle(TParticlePDG *p, TLorentzVector* t):
            TLorentzVector(*t), m_code(p->PdgCode()), m_charge(p->Charge())
            {};

        /**
         * @brief Returns the particle's charge.
         */
        int Charge(){
            return m_charge;
        }

        /**
         * @brief Returns the particle's PDG code.
         */
        int PdgCode(){
            return m_code;
        }

    private:
        int m_code;
        float m_charge;
};

/**
 * @brief Creates a TLorentzVector from momentum, angles, and mass.
 * @param p The momentum magnitude.
 * @param phi The azimuthal angle.
 * @param theta The polar angle.
 * @param m The mass of the particle.
 * @return A TLorentzVector with the specified parameters.
 */
TLorentzVector returnV (double p, double phi, double theta, double m){

    double px = p * sin(theta)*cos(phi);
    double py = p * sin(theta)*sin(phi);
    double pz = p * cos(theta);
    double p0 = sqrt(px*px+py*py+pz*pz+m*m);

    TLorentzVector P(px,py,pz,p0);

    return P;
}

/**
 * @brief The chi-square minimization function for TMinuit.
 * @param npar Number of parameters.
 * @param gin Unused.
 * @param f The calculated chi-square value.
 * @param par Array of parameters to be fitted.
 * @param flag Unused.
 *
 * This function calculates the generalized chi-square value by summing the
 * squared differences between measured and fitted parameters, and adds a
 * Breit-Wigner term for the W boson mass constraint.
 */
void fcn (int &npar, double *gin, double &f, double *par, int flag){
    f = 0.0;

    for (int i=0; i< res.data.size()-1; i++){
        f += pow((res.data[i]-par[i])/res.edata[i],2);
    }

    TLorentzVector P_e   = returnV(par[3], par[4], par[5], em->Mass());
    TLorentzVector P_nu  = returnV(par[6], par[7], par[8], 0);

    double M   = (P_e+P_nu).Mag();
    double eMW = Wp->Width();

    // The Breit-Wigner constraint for the W boson mass.
    f += -2*log(TMath::BreitWignerRelativistic(M,Wp->Mass(),eMW));
}

/**
 * @brief Applies detector resolution (Gaussian smearing) to a particle's kinematics.
 * @param track The particle to which the resolution is applied.
 *
 * The function modifies the momentum, phi, and theta of the particle based
 * on predefined resolution values for electrons and b-jets.
 */
void ApplyResolution (Particle &track){

    double m   = sqrt(track.T()*track.T()-track.Rho()*track.Rho());
    double P   = track.Rho();
    double phi = track.Phi();
    double the = track.Theta();

    double eP, ethe, ephi;

    if(abs(track.PdgCode())!=12){

        if(abs(track.PdgCode())==11){
            res.efit = res.ep*P*P;

            P   = rnd.Gaus(P, res.ep*P*P);
            phi = rnd.Gaus(phi, res.ephi);
            the = rnd.Gaus(the, res.ethe);

        } else if (abs(track.PdgCode())==5) {
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

/**
 * @brief Applies detector resolution to a vector of particles.
 * @param recoparticles The vector of particles to be smeared.
 */
void Reso(vector<Particle>& recoparticles){
    for (int i=0; i<recoparticles.size(); i++){
        ApplyResolution(recoparticles[i]);
    }
}

/**
 * @brief Simulates a single top quark decay event.
 * @param tgps The TGenPhaseSpace object for decay generation.
 * @param decayparticle A vector to store the generated particles.
 * @param wlep The weight of the decay.
 *
 * This function generates the kinematics of the final state particles from
 * the top quark decay using a phase space generator.
 */
void Event(TGenPhaseSpace &tgps, vector<Particle> &decayparticle, double &wlep){

    double t_mass  = top->Mass();
    double t_width = top->Width();
    double w_mass  = Wp->Mass();
    double w_width = Wp->Width();

    double t_sim   = sqrt(t_mass*t_width*tan(M_PI*(rnd.Rndm()-0.5))+t_mass*t_mass);
    double w_sim   = sqrt(w_mass*w_width*tan(M_PI*(rnd.Rndm()-0.5))+w_mass*w_mass);

    double bw_masses[2]  = {w_sim, bot->Mass()};
    double enu_masses[2] = {em->Mass(), 0};

    TLorentzVector t1_in(0,0,0,t_sim);

    tgps.SetDecay(t1_in, 2, bw_masses);
    wlep = tgps.Generate();
    Particle W1(Wp, tgps.GetDecay(0));
    Particle bot1(bot, tgps.GetDecay(1));

    tgps.SetDecay(W1, 2, enu_masses);
    wlep *= tgps.Generate();
    Particle pos(em, tgps.GetDecay(0));
    Particle nu_r(nu,tgps.GetDecay(1));

    decayparticle.push_back(bot1);
    decayparticle.push_back(pos);
    decayparticle.push_back(nu_r);
}

/**
 * @brief Performs the kinematic fit using TMinuit.
 * @param minuit The TMinuit object for minimization.
 * @param recoparticles The measured particles after detector simulation.
 * @param decayparticles The true simulated particles.
 * @return A vector of TLorentzVectors for the fitted particles.
 *
 * This function sets up the parameters for the fit, runs the minimization,
 * and returns the fitted kinematic quantities.
 */
vector<TLorentzVector> Kinematic_Fit(TMinuit *minuit, const vector<Particle> &recoparticles, const vector<Particle> &decayparticles){

    // Setting up initial data for the fit.
    //b-jet
    res.data.push_back(recoparticles[0].P());
    res.data.push_back(recoparticles[0].Phi());
    res.data.push_back(recoparticles[0].Theta());
    //electron
    res.data.push_back(recoparticles[1].P());
    res.data.push_back(recoparticles[1].Phi());
    res.data.push_back(recoparticles[1].Theta());
    //neutrino
    res.data.push_back(recoparticles[2].P());
    res.data.push_back(recoparticles[2].Phi());
    res.data.push_back(recoparticles[2].Theta());

    // Setting up errors for the fit parameters.
    res.edata.push_back(res.jfit);
    res.edata.push_back(res.jphi);
    res.edata.push_back(res.jthe);

    res.edata.push_back(res.efit);
    res.edata.push_back(res.ephi);
    res.edata.push_back(res.ethe);

    res.edata.push_back(res.np);
    res.edata.push_back(res.nphi);

    // Defining parameters for TMinuit.
    minuit->DefineParameter(0,"bP"  ,res.data[0],0.01,0.,0.);
    minuit->DefineParameter(1,"bPhi",res.data[1],0.01,0.,0.);
    minuit->DefineParameter(2,"bThe",res.data[2],0.01,0.,0.);

    minuit->DefineParameter(3,"eP"  ,res.data[3],0.01,0.,0.);
    minuit->DefineParameter(4,"ePhi",res.data[4],0.01,0.,0.);
    minuit->DefineParameter(5,"eThe",res.data[5],0.01,0.,0.);

    minuit->DefineParameter(6,"nuP"  ,res.data[6],0.01,0.,0.);
    minuit->DefineParameter(7,"nuPhi",res.data[7],0.01,0.,0.);
    minuit->DefineParameter(8,"nuThe",res.data[8],0.01,0.,0.);

    // Execute the minimization.
    minuit->Command("MIGRAD");

    // Retrieving fitted parameters.
    double pb,  epb,  phib,  ephib,  theb,  etheb;
    double pe,  epe,  phie,  ephie,  thee,  ethee;
    double pnu, epnu, phinu, ephinu, thenu, ethenu;

    minuit->GetParameter(0,pb,epb);
    minuit->GetParameter(1,phib,ephib);
    minuit->GetParameter(2,theb,etheb);

    minuit->GetParameter(3,pe,epe);
    minuit->GetParameter(4,phie,ephie);
    minuit->GetParameter(5,thee,ethee);

    minuit->GetParameter(6,pnu,epnu);
    minuit->GetParameter(7,phinu,ephinu);
    minuit->GetParameter(8,thenu,ethenu);

    // Creating TLorentzVectors from the fitted parameters.
    TLorentzVector P_b  = returnV(pb,phib,theb,bot->Mass());
    TLorentzVector P_e  = returnV(pe,phie,thee,em->Mass());
    TLorentzVector P_nu = returnV(pnu,phinu,thenu,0);

    // Filling pull histograms.
    hPulljp.Fill((recoparticles[0].P()-pb)/epb);
    hPulljphi.Fill((recoparticles[0].Phi()-phib)/ephib);
    hPulljthe.Fill((recoparticles[0].Theta()-theb)/etheb);

    hPullep.Fill((recoparticles[1].P()-pe)/epe);
    hPullephi.Fill((recoparticles[1].Phi()-phie)/ephie);
    hPullethe.Fill((recoparticles[1].Theta()-thee)/ethee);

    hThere.Fill((recoparticles[2].Theta()-thenu)/ethenu);

    return {P_b, P_e, P_nu};
}

/**
 * @brief Calculates the neutrino resolutions by simulating events.
 *
 * This function runs the first loop to simulate events and fill the
 * histograms used to determine the neutrino's resolution for the kinematic fit.
 */
void CalculateResolutions(TGenPhaseSpace &tgps) {
    vector<Particle> decayparticles, recoparticles;
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
}

/**
 * @brief Runs the main simulation and kinematic fit loop.
 *
 * This function performs the core logic of the program, generating events,
 * applying resolutions, performing the kinematic fit, and filling the
 * final result histograms.
 */
void RunSimulationAndFit(TGenPhaseSpace &tgps, TMinuit *minuit) {
    vector<Particle> decayparticles, recoparticles;
    for (int i=0; i<N_EV; i++){

        if(i%10000==0){
            cout << i << endl;
        }

        double wlep;
        Event(tgps, decayparticles, wlep);
        recoparticles = decayparticles;
        Reso(recoparticles);
        hMWRe.Fill((recoparticles[1]+recoparticles[2]).Mag());

        res.topsim.push_back((decayparticles[0]+decayparticles[1]+decayparticles[2]).Mag());

        // Calculate pre-fit neutrino four-vector and masses.
        TLorentzVector p_nu  = -recoparticles[0]-recoparticles[1];
        double pt_nu  = p_nu.Pt();
        double phi_nu = p_nu.Phi();
        TLorentzVector P_nu_pre = returnV(pt_nu, phi_nu, decayparticles[2].Theta(), 0.0);

        double mtop_pre = (P_nu_pre+recoparticles[0]+recoparticles[1]).Mag();
        double mw_pre   = (P_nu_pre+recoparticles[1]).Mag();

        // Perform the kinematic fit.
        vector<TLorentzVector> vettori = Kinematic_Fit(minuit, recoparticles, decayparticles);
        double mtoplep = (vettori[0]+vettori[1]+vettori[2]).Mag();
        double mw      = (vettori[1]+vettori[2]).Mag();

        res.topfit.push_back(mtoplep);

        // Fill histograms.
        hMtop_pre.Fill(mtop_pre,wlep);
        hMWpre.Fill(mw_pre,wlep);
        hMW.Fill(mw,wlep);
        hMtopl.Fill(mtoplep,wlep);
        hThenu.Fill(vettori[2].Theta(),wlep);

        decayparticles.clear();
        recoparticles.clear();
        res.data.clear(); res.edata.clear();
    }
}


/**
 * @brief The main function that orchestrates the simulation and analysis.
 */
int main(){
    TApplication app("app",NULL,NULL);

    // Setting up generators and TMinuit.
    TGenPhaseSpace tgps;
    rnd.SetSeed(12345);

    TMinuit *minuit = new TMinuit(9);
    minuit->SetPrintLevel(-1);
    minuit->SetFCN(fcn);

    vector<Particle> decayparticles, recoparticles;

    CalculateResolutions(tgps);

    // Set the neutrino resolution for the fit based on the first loop.
    res.np   = hNuResPt.GetStdDev();
    res.nphi = hNuResPhi.GetStdDev();

    RunSimulationAndFit(tgps, minuit);

    // Print final results.
    cout << "La massa del top è: " << top->Mass() << endl;
    cout << "La massa del top misurata è: " << hMtopl.GetMean() << " +- " << hMtopl.GetMeanError() << endl;


    // Drawing canvases with histograms.
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


    // Fill the pull histogram for the top quark mass.
    for (int i=0; i<res.topsim.size(); i++){
        hPulltop.Fill((res.topfit[i]-res.topsim[i])/hMtopl.GetStdDev());
    }

    // Draw the top quark mass pull.
    TCanvas c6;
    hPulltop.Draw(); hPulltop.SetXTitle("(m_{t}^{fit}-m_{t}^{sim})/#sigma");

    // Run the ROOT application.
    app.Run(true);

    delete minuit;
    return 0;
}