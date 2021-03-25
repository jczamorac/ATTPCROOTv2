

void random_angles(void)
{

double phi =0;
double r =0;



  TStopwatch timer;
  timer.Start();

  TString scriptfile = "e12014_pad_map_size.xml";
  TString dir = getenv("VMCWORKDIR");
  TString scriptdir = dir + "/scripts/"+ scriptfile;


  AtTpcMap* fAtMapPtr = new AtTpcMap();
  fAtMapPtr->GenerateATTPC();
  TH2Poly *fPadPlane = fAtMapPtr->GetATTPCPlane();
  fAtMapPtr->Dump();
  Bool_t MapIn = fAtMapPtr->ParseXMLMap(scriptdir);

  Float_t x=0;
  Float_t y=0;
  Float_t z=0;

  gSystem->Load("libATTPCReco.so");
 
  FairRunAna* run = new FairRunAna(); //Forcing a dummy run
 
std::vector<Float_t> PadCenterCoord;
std::vector<Int_t> vertexpads;


TH1D *phi_cal = new TH1D ("phi", "phi", 180, -190, 190);

for(int i=0; i<100000;i++){



phi = 2.0*3.14159*gRandom->Uniform(); 
r = 200.0*gRandom->Uniform() + 30.0; 

x = r*cos(phi);
y = r*sin(phi);
z =1;
Int_t bin=fPadPlane->Fill(x,y,z);
		//cout<<bin<<endl;
		if(bin<0 || bin> 10240) continue;
		
		PadCenterCoord = fAtMapPtr->CalcPadCenter(bin);
                  double valx = PadCenterCoord[0];
                  double valy = PadCenterCoord[1];
		PadCenterCoord.clear();

		double phi_in = 2.0*3.14159*gRandom->Uniform(); 
		int r_in = 30.0*gRandom->Uniform() ; 

		float x_in = r_in*cos(phi_in) ;
		float y_in = r_in*sin(phi_in) ;
		z =1;
		/*Int_t bin_in=fPadPlane->Fill(x_in,y_in,z);
		//cout<<bin<<endl;
		if(bin_in<0 || bin_in> 10240) continue;

		PadCenterCoord = fAtMapPtr->CalcPadCenter(bin_in);
                  double valx_in = PadCenterCoord[0];
                  double valy_in = PadCenterCoord[1];
		PadCenterCoord.clear();
		*/

		double Y = valy - y_in;
		double X = valx - gRandom->Gaus(valx,50);
		double angle_az = atan2(Y,X) * 180.0/3.1415;
		double Z = gRandom->Uniform(valx,1000);
		
		//if(X< 0.05){ 
		cout<<X<<"  "<<Y<<"  "<<angle_az<<endl;

		phi_cal->Fill(angle_az);

		TVector3 vp1(TMath::Sign(1,valx)*fabs(X),TMath::Sign(1,valy)*fabs(Y),Z);
		float phi_sim = vp1.Phi()*180./3.1415;

		phi_cal->Fill(phi_sim);
		//}
		
//cout<<phi<<endl;
}


		phi_cal->Draw("");
	    //fPadPlane->Draw("COL L0");
            //fPadPlane -> SetMinimum(1.0);
            //gStyle->SetOptStat(0);

   	    std::cout << std::endl << std::endl;
            std::cout << "Macro finished succesfully."  << std::endl << std::endl;
            // -----   Finish   -------------------------------------------------------
            timer.Stop();
            Double_t rtime = timer.RealTime();
            Double_t ctime = timer.CpuTime();
            cout << endl << endl;
            cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
            cout << endl;


}
