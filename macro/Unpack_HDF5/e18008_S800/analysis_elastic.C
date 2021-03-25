//#include "spline.h"  //cubic interpolation
//#include "TInverseMap.hh"
#include <unistd.h>

TGraph *graphtable;

//change MTDC gate S800_timeMTDCXf.at(k)>215 && etc..
//change the ObjCorr1C1 ObjCorr1C2
//change the Ekin_proj to the correct value,
//check the energy loss function
//check calibration/correction coef on the dE, and ToF S800.
//If plots make no sens first check that the angles are well calculated (e.g. for real data GetThetaPhi(..,..,1)

static Double_t proton_mass = 1.0078250322 * 931.494;
static Double_t proj_mass = 14.008596359 * 931.494;
static Double_t target_mass = 2.01410177812 * 931.494;
static Double_t recoil_mass = 14.00307400443 * 931.494;
static Double_t he2_mass = 2.0 * proton_mass;
static Double_t Ekin_proj = 105.16 * 14.008596359;//115

static Double_t corrGainE1up = 1.;
static Double_t corrGainE1down = 1.;
static double x0_corr_tof = 0.0;
static double afp_corr_tof = 0.0;
static double rf_offset = 19320.0;
static double afp_corr_dE = 0.0;
static double x0_corr_dE = 0.0;
static double bta_corr = 1.0;
static double coeff_hodo[32] = {100000000, 1164, 887, 100000000, 824, 1728, 826,
	100000000, 941, 962, 875, 915, 875, 835, 1079, 842, 855, 870, 891, 1786, 960,
	940, 825, 625, 990, 2195, 805, 870, 575, 1004, 1164, 748 };

std::vector <double> get_invmap_vars(TInverseMap *inv_map, double x0, double y0, double afp, double bfp)
	{
		int order;
		double sinb_sina;
		vector <double> outPuts;

		order = 5;

		outPuts.push_back(inv_map->Ata(order, x0, afp, y0, bfp));//ata
		outPuts.push_back(inv_map->Bta(order, x0, afp, y0, bfp) * bta_corr);//bta
		outPuts.push_back(inv_map->Yta(order, x0, afp, y0, bfp) * 1000.);//yta
		outPuts.push_back(inv_map->Dta(order, x0, afp, y0, bfp));//dta
		outPuts.push_back(atan(sqrt(pow(tan(outPuts.at(0)), 2) + pow(tan(outPuts.at(1)), 2))));//theta_lab
		outPuts.push_back(atan(tan(outPuts.at(1))/tan(outPuts.at(0))));//phi
		if (outPuts.at(0) < 0)
		outPuts.at(5) = 3.141592653 + outPuts.at(5);
		else if (outPuts.at(1) < 0)
		outPuts.at(5) = 2*3.141592653 + outPuts.at(5);

		return outPuts;
	}

std::vector <double> get_invmap_vars(TInverseMap *inv_map, double x0, double y0, double afp, double bfp, double z)
	{
		int order;
		double sinb_sina;
		vector <double> outPuts;

		order = 5;

		outPuts.push_back(inv_map->Ata(order, x0, afp, y0, bfp, z));//ata
		outPuts.push_back(inv_map->Bta(order, x0, afp, y0, bfp, z) * bta_corr);//bta
		outPuts.push_back(inv_map->Yta(order, x0, afp, y0, bfp, z) * 1000.);//yta
		outPuts.push_back(inv_map->Dta(order, x0, afp, y0, bfp, z));//dta
		outPuts.push_back(atan(sqrt(pow(tan(outPuts.at(0)), 2) + pow(tan(outPuts.at(1)), 2))));//theta_lab
		outPuts.push_back(atan(tan(outPuts.at(1))/tan(outPuts.at(0))));//phi
		if (outPuts.at(0) < 0)
		outPuts.at(5) = 3.141592653 + outPuts.at(5);
		else if (outPuts.at(1) < 0)
		outPuts.at(5) = 2*3.141592653 + outPuts.at(5);

		return outPuts;
	}


Double_t FindAngleBetweenTracks(const TVector3 &vec1,const TVector3 &vec2)
    {
        Double_t ang = vec1.Angle(vec2);
      	return ang;
    }

double eloss_approx(double rangep){//rangep in mm, elossp in MeV; get this approximated function by fitting the GEANT4 E vs R obtained with the simulation
		double elossp;
		if(rangep<131.) elossp =  0.0833284+ 0.0118757*rangep -7.3698e-05*pow(rangep,2)+  2.42758e-07*pow(rangep,3);
		if(rangep>=131.) elossp =  0.36519+ 0.004761*rangep -4.28371e-06*pow(rangep,2)+ 2.35727e-09*pow(rangep,3);
		return elossp;
	}


void SetERtable(){//fit of the GEANT4 E vs R obtained from the simulation with the function model given by LISE++
		//ifstream fER("p_in_d_LISE_1atm.dat");//from LISE++, with this one comment the convertion um to mm
	//	ifstream fER("p_in_d_350torr.txt");//from LISE++,
		ifstream fER("eLossTables/p_in_d_530torr.txt");//from LISE++,
		//ifstream fER("p_in_d_730torr.txt");//from LISE++,
		//ifstream fER("p_in_D2_LISE_530torr.dat");//from LISE++
		Double_t l1=0, l2=0, l3=0, l4=0, l5=0;
		Int_t model=1;
		vector <vector<Double_t>> Energy_Range;

		for (string line; getline(fER, line);) {
			stringstream parse_die(line);
			vector<Double_t> iRE;
			parse_die >> l1 >> l2 >> l3 >> l4 >> l5;
			iRE.push_back(l1);//E in MeV
			iRE.push_back(l2/1000.);//range in mm, model 1
			iRE.push_back(l3/1000.);//range in mm, model 2
			iRE.push_back(l4/1000.);//range in mm, model 3
			iRE.push_back(l5/1000.);//range in mm, model 4
			Energy_Range.push_back(iRE);
		}
		fER.close();
		Int_t v_size = Energy_Range.size();
		Double_t X[v_size];
		Double_t Y[v_size];
		for(Int_t i=0; i<v_size; i++){
			X[i]=Energy_Range.at(i).at(0);
			Y[i]=Energy_Range.at(i).at(model);//*0.738 for LISE++ Eloss to match with GEANT4
			//cout<<X[i]<<" "<<Y[i]<<endl;
		}
		graphtable = new TGraph(v_size,Y,X);
	}


TGraph *aGraph1;
Double_t aGraph1toFunc(double *xx, double *)
{
   return aGraph1->Eval(xx[0]);
}
TGraph *aGraph2;
Double_t aGraph2toFunc(double *xx, double *)
{
   return aGraph2->Eval(xx[0]);
}


void analysis_elastic(int runNumberS800=2063, int runNumberATTPC=63)
	{

		SetERtable();

		FairRunAna* run = new FairRunAna(); //Forcing a dummy run
		ATd2HeAnalysis *d2heana = new ATd2HeAnalysis ();

		// std::string digiFileName = "run_2021_0030.root";//merged.root
		TString digiFileName = TString::Format("/home/juan/NSCL/elastics/run_%04d_%04d.root", runNumberS800, runNumberATTPC);
		//TString digiFileName = TString::Format("run_%04d_%04d.root", runNumberS800, runNumberATTPC);

		TFile* file = new TFile(digiFileName,"READ");

		TTree* tree = (TTree*) file -> Get("cbmsim");
		Int_t nEvents = tree -> GetEntries();

		S800Calc *s800cal = new S800Calc();
		TBranch *bS800cal = tree->GetBranch("s800cal");
	  bS800cal->SetAddress(&s800cal);


		TTreeReader reader("cbmsim", file);
		// TTreeReaderValue<S800Calc> s800Calc(reader, "s800cal");
		TTreeReaderValue<TClonesArray> ransacArray(reader, "ATRansac");

		TFile* outfile;
		// TString  outFileNameHead = "attpcana_merg_newRansac.root";
		//TString outFileNameHead = TString::Format("/mnt/analysis/e18008/rootAna/runAnalyzed_%04d_%04d.root", runNumberS800, runNumberATTPC);
		TString outFileNameHead = TString::Format("runAnalyzed_%04d_%04d.root", runNumberS800, runNumberATTPC);
		//TString digiFileName = TString::Format("/mnt/analysis/e18008/rootMerg/jcz_test/run_%04d_%04d.root", runNumberS800, runNumberATTPC);
		//TString outFileNameHead = TString::Format("runAnalyzed_%04d_%04d.root", runNumberS800, runNumberATTPC);

		outfile   = TFile::Open(outFileNameHead.Data(),"recreate");
		outfile->cd();//New Bragg
	  outfile->mkdir("Bragg");//New Bragg
		outfile->mkdir("nogate");
		outfile->mkdir("gated");

		TH1F* scatteringAngle = new TH1F("scatteringAngle","scatteringAngle",1000,0,200);
		TH1F* energy = new TH1F("enegy","energy",100,0,40);
		TH2F* ang_vs_energy = new TH2F("ang_vs_energy","ang_vs_energy,",100,0,200,100,0,40);
		//TH2D *tracks_z_r = new TH2D ("tracks_z_r", "ZvsR", 500, -100, 1000, 500, 0, 300);
		//TH2D *tracks_x_y = new TH2D ("tracks_x_y", "XvsY", 500, -300, 300, 500, -300, 300);
		TH1D *theta_r_he2_reco = new TH1D ("theta_r_he2_reco", "theta 2He", 1800, 0, 180);
		TH1D *kin_r_he2_reco = new TH1D ("kin_r_he2_reco", "Energy 2He", 100, 0, 5);
		TH1D *phi_r_he2_reco = new TH1D ("phi_r_he2_reco", "phi 2He", 3600, -180, 180);
		TH2D *theta_kin_he2_reco = new TH2D ("theta_kin_he2_reco", "Kin vs Theta 2He", 1800, 0, 180, 100, 0, 5);
		TH1D *thetacm_he2_reco = new TH1D ("thetacm_he2_reco", "thetacm_he2", 200, 0, 20);
		TH1D *Ex_reco = new TH1D ("Ex_reco", "Ex_reco", 350, -5, 30);
		TH2D *thetacm_Ex_he2_reco = new TH2D ("thetacm_Ex_he2_reco", "thetacm_Ex_he2", 200, 0, 20, 350, -5, 30);
		TH1D *ex_he2_reco = new TH1D ("ex_he2_reco", "ex_he2", 100, 0, 10);
		TH1D *epsilon_pp_reco = new TH1D ("epsilon_pp_reco", "#epsilon_{pp}", 100, 0, 10);
		TH2D *thetacm_epsilon_pp_reco = new TH2D ("thetacm_epsilon_pp_reco", "#theta_{cm} #epsilon_{pp} #^{2}He", 200, 0, 20, 100, 0, 10);
		//----- S800

		// TH2D *tof_dE = new TH2D ("tof_dE", "tof_dE", 250, 1000, 1500, 250, 0, 500);//PID
		TH2D *XfpObj_tof = new TH2D ("XfpObj_tof", "XfpObj_tof", 500,-70,-20,600,250,280);//PID1
		TH2D *ICSum_Obj = new TH2D ("ICSum_Obj", "ICSum_Obj",500,-70,-20,1000,50,750);//PID2



		TH2D *XfpObj_Obj_g = new TH2D ("XfpObj_Obj_g", "XfpObj_Obj",500,-120,-20,600,240,300);//PID1
	  TH2D *ICSum_Obj_g = new TH2D ("ICSum_Obj_g", "ICSum_Obj",500,-120,-20,1000,0,750);//PID2
	  TH2D *X_afp_g = new TH2D ("X_afp_g", "Xfp_vs_afp",500,-300,300,1000,-0.80,0.8);//PID3
		TH2D *XfpObj_Obj_ng = new TH2D ("XfpObj_Obj_ng", "XfpObj_Obj",500,-120,-20,600,240,300);//PID1
	  TH2D *ICSum_Obj_ng = new TH2D ("ICSum_Obj_ng", "ICSum_Obj",500,-120,-20,1000,0,750);//PID2
	  TH2D *X_afp_ng = new TH2D ("X_afp_ng", "Xfp_vs_afp",500,-300,300,1000,-0.80,0.8);//PID3


		TH1F* verZ = new TH1F("vertexZ","vertexZ",500,-200,1500);
		TH2F* vertexvsAngle = new TH2F("vertexvsAngle","vertexvsAngle",500,-200,1500, 180,0,180);
		TH2D *theta_range = new TH2D ("theta_range", "theta_range", 180, 0, 180, 250, 0, 600);
		TH2D *theta_kin = new TH2D ("theta_kin", "theta_kin", 180, 0, 180, 250, 0, 4);



		TH2D *dta_ata = new TH2D ("dta_ata", "dta_ata", 250, -0.25, 0.25, 100, -10, 10);//acceptance
		TH2D *x_y_crdc1 = new TH2D ("x_y_crdc1", "x_dy_crdc1", 300, -300, 300, 150, -150, 150);//positions crdc1
		TH2D *x_y_crdc2 = new TH2D ("x_y_crdc2", "x_dy_crdc2", 300, -300, 300, 150, -150, 150);//positions crdc2

		TH1D *hBr0;
		TH1D *hBr1;

		//-----
		Int_t ivt = 0,ivt_same=0;
		Double_t range_p1 = 0.,range_p2 =0.;
		Double_t eLoss_p1_reco = 0.0, eLoss_p2_reco = 0.0;
		Double_t epsilon_pp = -999;
		Double_t theta1=0., theta2=0., phi1=0., phi2=0., angle12=0.;
		Double_t mom1_norm_reco=0., mom2_norm_reco=0.;
		Double_t E_tot_he2=0., he2_mass_ex=0.;
		Double_t kin_He2=0., theta_He2=0.,kin_He2_same=0., theta_He2_same=0., phi_He2=0.;
		Double_t theta_cm=0., Ex4=0., Ex_reco_same=0.;
		Double_t lastX1=0.,lastX2=0.,lastY1=0.,lastY2=0.,lastZ1=0.,lastZ2=0., vertexX=0., vertexY=0., vertexZ=0.;
		Double_t ata=0., dta=0.;
		ULong64_t S800_timeStamp=0;
		Double_t S800_timeRf=0.,S800_x0=0.,S800_x1=0.,S800_y0=0.,S800_y1=0.,S800_E1up=0.,S800_E1down=0.,S800_tof=0.,
		S800_tofCorr=0.,S800_dE=0.,S800_dECorr=0.,S800_hodoSum=0.,S800_afp=0.,S800_bfp=0.,S800_ata=0.,S800_bta=0.,
		S800_yta=0.,S800_dta=0.,S800_thetaLab=0.,S800_phi=0.,S800_timeE1up=0.,S800_timeE1down=0.,S800_timeE1=0.,S800_timeXf=0.,S800_timeObj=0.;
		Double_t S800_XfObj_tof=0.,S800_ObjCorr=0., S800_Obj=0.;
		Double_t S800_ICSum=0.;
		Double_t MaxR1, MaxR2, MaxZ1, MaxZ2;
    Int_t CondMTDCObj = 0,CondMTDCXfObj = 0;
		Double_t beam_theta=0., beam_phi=0.;
		Double_t chargeTot1=0., chargeTot2=0.;
		Bool_t coincidence = kFALSE;

		TVector3 mom_proton1_reco, mom_proton2_reco,mom_He2_reco;
		std::vector<float> VerX, VerY, VerZ;
		std::vector<float> PlastX, PlastY, PlastZ;
		std::vector<float> Range;
		std::vector<float> ThetaLab;
		std::vector<float> PhiLab;
		std::vector<int> Nhits_track;
		std::vector<float> Charge;

		TTree *anatree = new TTree("anatree","new TTree");

		anatree->Branch("ivt",&ivt);
		anatree->Branch("VerX",&VerX);
		anatree->Branch("Very",&VerY);
		anatree->Branch("VerZ",&VerZ);
		anatree->Branch("PlastX",&PlastX);
		anatree->Branch("PlastY",&PlastY);
		anatree->Branch("PlastZ",&PlastZ);
		anatree->Branch("Range",&Range);
		anatree->Branch("ThetaLab",&ThetaLab);
		anatree->Branch("PhiLab",&PhiLab);
		anatree->Branch("Nhits_track",&Nhits_track);
		anatree->Branch("coincidence",&coincidence);
		anatree->Branch("Charge_track",&Charge);

		//----- S800
		anatree->Branch("S800_timeStamp",&S800_timeStamp,"S800_timeStamp/l");
		//anatree->Branch("S800_timeRf",&S800_timeRf);
		//anatree->Branch("S800_timeE1up",&S800_timeE1up);
		//anatree->Branch("S800_timeE1down",&S800_timeE1down);
		//anatree->Branch("S800_timeE1",&S800_timeE1);
		anatree->Branch("S800_timeXf",&S800_timeXf);
		anatree->Branch("S800_timeObj",&S800_timeObj);
		anatree->Branch("S800_tof",&S800_tof);
		anatree->Branch("S800_XfObj_tof",&S800_XfObj_tof);
		anatree->Branch("S800_ObjCorr",&S800_ObjCorr);
		anatree->Branch("S800_Obj",&S800_Obj);
		anatree->Branch("CondMTDCObj",&CondMTDCObj);
		anatree->Branch("CondMTDCXfObj",&CondMTDCXfObj);
		anatree->Branch("S800_ICSum",&S800_ICSum);

		anatree->Branch("S800_x0",&S800_x0);
		anatree->Branch("S800_x1",&S800_x1);
		anatree->Branch("S800_y0",&S800_y0);
		anatree->Branch("S800_y1",&S800_y1);
		anatree->Branch("S800_E1up",&S800_E1up);
		anatree->Branch("S800_E1down",&S800_E1down);
		anatree->Branch("S800_tof",&S800_tof);
		anatree->Branch("S800_tofCorr",&S800_tofCorr);
		anatree->Branch("S800_dE",&S800_dE);
		anatree->Branch("S800_dECorr",&S800_dECorr);
		anatree->Branch("S800_hodoSum",&S800_hodoSum);
		anatree->Branch("S800_afp",&S800_afp);
		anatree->Branch("S800_bfp",&S800_bfp);
		anatree->Branch("S800_ata",&S800_ata);
		anatree->Branch("S800_bta",&S800_bta);
		anatree->Branch("S800_yta",&S800_yta);
		anatree->Branch("S800_dta",&S800_dta);
		anatree->Branch("S800_thetaLab",&S800_thetaLab);
		anatree->Branch("S800_phi",&S800_phi);


		std::vector<ATTrack> trackCand;


		//----------------------- S800 -------------------------------------------------
		std::vector< std::string > mapList;
                mapList.push_back("invMap/invmap_14N/invmap_-05.inv");
                mapList.push_back("invMap/invmap_14N/invmap_-04.inv");
                mapList.push_back("invMap/invmap_14N/invmap_-03.inv");
                mapList.push_back("invMap/invmap_14N/invmap_-02.inv");
                mapList.push_back("invMap/invmap_14N/invmap_-01.inv");
                mapList.push_back("invMap/invmap_14N/invmap_00.inv");
                mapList.push_back("invMap/invmap_14N/invmap_01.inv");
                mapList.push_back("invMap/invmap_14N/invmap_02.inv");
                mapList.push_back("invMap/invmap_14N/invmap_03.inv");
                mapList.push_back("invMap/invmap_14N/invmap_04.inv");
                mapList.push_back("invMap/invmap_14N/invmap_05.inv");
                mapList.push_back("invMap/invmap_14N/invmap_06.inv");
                mapList.push_back("invMap/invmap_14N/invmap_07.inv");
                mapList.push_back("invMap/invmap_14N/invmap_08.inv");
                mapList.push_back("invMap/invmap_14N/invmap_09.inv");
                mapList.push_back("invMap/invmap_14N/invmap_10.inv");

		std::vector<Double_t> mapDist;
		for(int i=0;i<mapList.size();i++){
			mapDist.push_back(-0.5+0.1*(i));
		}

		//inv_map->SetDistPivotTarget(mapDist);
		//std::string mapFile="inv_map.inv";
		//TInverseMap *inv_map = new TInverseMap(mapFile.c_str());
		//TInverseMap *inv_map = new TInverseMap(mapList);
		TInverseMap *inv_map = new TInverseMap();
		inv_map->SetDistPivotTarget(mapDist);
		inv_map->ReadMultiMapFile(mapList);
		std::cout<<" mapDist "<<mapDist.size()<<" "<<mapDist.at(2)<<std::endl;


		/// --------------------- Event loop -------------------------------------------

		for(Int_t i=0;i<nEvents;i++){
			s800cal->Clear();
			trackCand.clear();
			bS800cal->GetEntry(i);
			reader.Next();

			S800_timeStamp=0;
			S800_timeRf=0.;S800_x0=0.;S800_x1=0.;S800_y0=0.;S800_y1=0.;S800_E1up=0.;S800_E1down=0.;S800_tof=0.;
			S800_tofCorr=0.;S800_dE=0.;S800_dECorr=0.;S800_hodoSum=0.;S800_afp=0.;S800_bfp=0.;S800_ata=0.;
			S800_bta=0.;S800_yta=0.;S800_dta=0.;S800_thetaLab=0.;S800_phi=0.;
			S800_timeE1up=0.;S800_timeE1down=0.;S800_timeE1=0.;S800_timeXf=0.;S800_timeObj=0.;
			S800_tof=0.;

			VerX.clear(); VerY.clear(); VerZ.clear();
			PlastX.clear(); PlastY.clear(); PlastZ.clear();
			Range.clear(); ThetaLab.clear(); PhiLab.clear(); Nhits_track.clear();
			Charge.clear();

			// ATRANSACN::ATRansac* fATRansac  = dynamic_cast<ATRANSACN::ATRansac*> (ransacArray->At(0));
		  // ATRansacMod* fATRansac  = dynamic_cast<ATRansacMod*> (ransacArray->At(0));
		  // ATMlesacMod* fATRansac  = dynamic_cast<ATMlesacMod*> (ransacArray->At(0));
		  ATLmedsMod* fATRansac  = dynamic_cast<ATLmedsMod*> (ransacArray->At(0));
			 if(fATRansac==nullptr){
			  std::cout<<" Null pointer fATRansac "<<"\n";
				continue;
				}

			trackCand = fATRansac->GetTrackCand();

			 //std::cout<<s800cal->GetTS()<<" "<<i<<"  "<<trackCand.size()<<"  "<<s800cal->GetIsInCut()<<std::endl;

			 //---------------------ungated


			 //--------------------











				//----------------------- S800 -------------------------------------------------

    CondMTDCObj = 0;
    CondMTDCXfObj = 0;
    vector<Float_t> S800_timeMTDCObj = s800cal->GetMultiHitTOF()->GetMTDCObj();
    vector<Float_t> S800_timeMTDCXf = s800cal->GetMultiHitTOF()->GetMTDCXf();
    Float_t S800_timeObjSelect=-999;
    Float_t S800_timeXfSelect=-999;
    S800_XfObj_tof = -999.;
    S800_ObjCorr = -999;
    S800_Obj = -999;
    Double_t ObjCorr1C1 = 70.; //70
    Double_t ObjCorr1C2 = 0.0085; //0.0085

    S800_ICSum = s800cal->GetIC()->GetSum();
    S800_x0 = s800cal->GetCRDC(0)->GetX();
    S800_x1 = s800cal->GetCRDC(1)->GetX();
    S800_y0 = s800cal->GetCRDC(0)->GetY();
    S800_y1 = s800cal->GetCRDC(1)->GetY();
    S800_afp = atan( (S800_x1-S800_x0)/1073. );
    S800_bfp = atan( (S800_y1-S800_y0)/1073. );

    for(int k=0; k<S800_timeMTDCXf.size(); k++){
    	if(S800_timeMTDCXf.at(k)>160 && S800_timeMTDCXf.at(k)<240) S800_timeXfSelect=S800_timeMTDCXf.at(k);
    }
    for(int k=0; k<S800_timeMTDCObj.size(); k++){
    	if(S800_timeMTDCObj.at(k)>-120 && S800_timeMTDCObj.at(k)<-20) S800_timeObjSelect=S800_timeMTDCObj.at(k);
    }


    if(S800_timeObjSelect!=-999){
	CondMTDCObj=1;
	S800_Obj = S800_timeObjSelect;
	S800_ObjCorr = S800_timeObjSelect + ObjCorr1C1*S800_afp + ObjCorr1C2*S800_x0;

    }
    if(S800_timeXfSelect!=-999 && S800_timeObjSelect!=-999) {
    	S800_XfObj_tof=S800_timeXfSelect-S800_timeObjSelect;
	CondMTDCXfObj=1;
    }



				S800_timeStamp = s800cal->GetTS();



			  S800_E1up = s800cal->GetSCINT(0)->GetDEup();
			  S800_E1down = s800cal->GetSCINT(0)->GetDEdown();
				S800_dE = sqrt( (corrGainE1up*S800_E1up) * (corrGainE1down* S800_E1down ) );
				S800_dECorr = S800_dE + afp_corr_dE*S800_afp + x0_corr_dE*fabs(S800_x0);
				for (Int_t j=0; j<32; j++) if (s800cal->GetHODOSCOPE(j)->GetEnergy()>=10 && s800cal->GetHODOSCOPE(j)->GetEnergy()<=4000) S800_hodoSum += s800cal->GetHODOSCOPE(j)->GetEnergy()*3000./coeff_hodo[j];


				// tof_dE->Fill(S800_tof,S800_dE);
				ICSum_Obj->Fill(S800_ObjCorr,S800_ICSum);
		    XfpObj_tof->Fill(S800_XfObj_tof,S800_ObjCorr);

			//dta_ata->Fill(S800_dta,S800_ata*180./TMath::Pi());//ata in deg
				x_y_crdc1->Fill(S800_x0,S800_y0);
				x_y_crdc2->Fill(S800_x1,S800_y1);

			//std::cout<<i<<" "<<" "<<S800_timeStamp<<" "<<S800_timeRf<<" "<<S800_y0<<std::endl;

			///no gate
			XfpObj_Obj_ng->Fill(S800_ObjCorr,S800_XfObj_tof);
			ICSum_Obj_ng->Fill(S800_ObjCorr,S800_ICSum);
			X_afp_ng->Fill(S800_x0,S800_afp);

			if(s800cal->GetIsInCut()==kTRUE){
			//gated
				XfpObj_Obj_g->Fill(S800_ObjCorr,S800_XfObj_tof);
				ICSum_Obj_g->Fill(S800_ObjCorr,S800_ICSum);
				X_afp_g->Fill(S800_x0,S800_afp);
			}

			coincidence = s800cal->GetIsInCut();

		   //if(trackCand.size()>0 && s800cal->GetIsInCut()==kTRUE) {
			if(trackCand.size()>0) {



				for(Int_t w=0;w<trackCand.size();w++){


					theta1=0.;  phi1=0.;  range_p1=0.;  eLoss_p1_reco=0.;  mom1_norm_reco=0.;  //reset variables
					theta_cm=0.; Ex4=0.;
			  	MaxR1=0.;  MaxZ1=0.;

					TVector3 vertexMean = trackCand.at(w).GetTrackVertex();
      		TVector3 lastPoint1 = trackCand.at(w).GetLastPoint();
      		MaxR1=sqrt(pow(lastPoint1.X(),2)+pow(lastPoint1.Y(),2));
      		MaxZ1=lastPoint1.Z();

					if(vertexMean.Z()<-100 || vertexMean.Z()>1200) continue;


					std::cout<<"Event:   "<<i<<std::endl;


					//------------------------------------------------------------------------------

					verZ->Fill(vertexMean.Z());


        	lastX1 = lastPoint1.X();
        	lastY1 = lastPoint1.Y();
        	lastZ1 = lastPoint1.Z();
        	vertexX = vertexMean.X();
        	vertexY = vertexMean.Y();
        	vertexZ = vertexMean.Z();

        	theta1 = trackCand.at(w).GetThetaPhi(vertexMean, lastPoint1,1).first;
        	phi1 = trackCand.at(w).GetThetaPhi(vertexMean, lastPoint1,1).second;

        	std::vector<Double_t> fitPar1 = trackCand.at(w).GetFitPar();

        	range_p1 = trackCand.at(w).GetLinearRange(vertexMean,lastPoint1);
        	eLoss_p1_reco = graphtable->Eval(range_p1);

					theta_range->Fill(theta1*180/3.1415,range_p1);
					theta_kin->Fill(theta1*180/3.1415,eLoss_p1_reco);

					vertexvsAngle->Fill(vertexMean.Z(),theta1*180/3.1415,range_p1);

					VerX.push_back(vertexMean.X());
					VerY.push_back(vertexMean.Y());
					VerZ.push_back(vertexMean.Z());
					PlastX.push_back(lastPoint1.X());
					PlastY.push_back(lastPoint1.Y());
					PlastZ.push_back(lastPoint1.Z());
					Range.push_back(range_p1);
					ThetaLab.push_back(theta1*180/3.1415);
					PhiLab.push_back(phi1*180/3.1415);

					std::vector<ATHit> *HitArray = trackCand.at(w).GetHitArray();
					Nhits_track.push_back(HitArray->size());
					double cargaInt = 0;
					for(auto j = HitArray->begin(); j != HitArray->end(); ++j){

				         double lacarga = j->GetCharge();
								 cargaInt +=lacarga;
							 }
					Charge.push_back(cargaInt);

				}//all tracks
		} //gated and  tracks.size()>0


		anatree->Fill();
	}// Event loop


	/// --------------------- End event loop ---------------------------------------
	outfile->cd();
	anatree->Write();

	/*
	//	tracks_z_r->Write ();
	//	tracks_x_y->Write ();
	theta_r_he2_reco->Write ();
	phi_r_he2_reco->Write ();
	kin_r_he2_reco->Write ();
	theta_kin_he2_reco->Write ();
	thetacm_he2_reco->Write ();
	Ex_reco->Write ();
	ex_he2_reco->Write ();
	thetacm_Ex_he2_reco->Write ();
	thetacm_epsilon_pp_reco->Write();
	epsilon_pp_reco->Write();

	XfpObj_tof->Write();
	ICSum_Obj->Write();
	dta_ata->Write();
	x_y_crdc1->Write();
	x_y_crdc2->Write();
	*/
	outfile->cd("nogate");
	XfpObj_Obj_ng->Write();
	ICSum_Obj_ng->Write();
	X_afp_ng->Write();
	outfile->cd("gated");
	XfpObj_Obj_g->Write();
	ICSum_Obj_g->Write();
	X_afp_g->Write();
	theta_range->Write();
	vertexvsAngle->Write();

	outfile->Close();

} //end main
