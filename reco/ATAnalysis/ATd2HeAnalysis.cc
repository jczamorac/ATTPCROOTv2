#include "ATd2HeAnalysis.hh"
#include "TMath.h"
#include <memory>


ClassImp(ATd2HeAnalysis)

ATd2HeAnalysis::ATd2HeAnalysis()
{


      THcm = 0.0;
      MissingMass = 0.0;
      theta_fit = 0.0;
      phi_fit = 0.0;
      ok = false;
      vertex = 0.0;
}

ATd2HeAnalysis::~ATd2HeAnalysis()
{
}

void ATd2HeAnalysis::SetThetaCM(Double_t value)    { THcm=value;}
void ATd2HeAnalysis::SetMissingMass(Double_t  value) { MissingMass=value;}

Double_t  ATd2HeAnalysis::GetThetaCM()    { return THcm;}
Double_t  ATd2HeAnalysis::GetMissingMass() { return MissingMass;}
Double_t  ATd2HeAnalysis::GetThetaFit() { return theta_fit;}
Double_t  ATd2HeAnalysis::GetPhiFit() { return phi_fit;}
bool         ATd2HeAnalysis::GetFitStatus() { return ok;}
Double_t  ATd2HeAnalysis::GetVertexFit() { return vertex;}

Double_t ATd2HeAnalysis::omega(Double_t x, Double_t  y, Double_t z) { return sqrt(x*x + y*y + z*z -2*x*y -2*y*z -2*x*z);}


Double_t ATd2HeAnalysis::EnergyFluctuation(Double_t energypoint){

        TRandom3* gRandome = new TRandom3();
        Double_t fanofactor = 0.2;
	Double_t ion_pot = 40e-6; //in MeV
	Int_t  Nelectrons = 0;
        Double_t eout ;

        Nelectrons = TMath::Floor(energypoint/ion_pot);
	Nelectrons = gRandome -> Gaus(Nelectrons, sqrt(Nelectrons*fanofactor));
	eout = Nelectrons*ion_pot;
        delete  gRandome;

        return eout;


}


void ATd2HeAnalysis::kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t K_proj, Double_t  thetalab, Double_t  K_eject){

 //in this definition: m1(projectile); m2(target); m3(ejectile); and m4(recoil);
        double Et1 = K_proj + m1;
        double Et2 = m2;
        double Et3 = K_eject + m3;
        double Et4  = Et1 + Et2 - Et3;
        double m4_ex, Ex, theta_cm;
        double s,t,u; //---Mandelstam variables



        s = pow(m1,2) + pow(m2,2) +2*m2*Et1;
        u = pow(m2,2) + pow(m3,2) - 2*m2*Et3;

        m4_ex = sqrt(  (cos(thetalab) * ATd2HeAnalysis::omega(s,pow(m1,2),pow(m2,2)) * ATd2HeAnalysis::omega(u,pow(m2,2),pow(m3,2)) - (s - pow(m1,2) - pow(m2,2))*(pow(m2,2) + pow(m3,2) - u) )/(2*pow(m2,2)) + s + u - pow(m2,2)  );
        Ex = m4_ex - m4;

        t =   pow(m2,2) + pow(m4_ex,2) - 2*m2*Et4;


        //for inverse kinematics Note: this angle corresponds to the recoil
        theta_cm = TMath::Pi() - acos( ( pow(s,2) +s*(2*t - pow(m1,2) - pow(m2,2) - pow(m3,2) - pow(m4_ex,2)) + (pow(m1,2) - pow(m2,2))*(pow(m3,2) - pow(m4_ex,2)) )/( ATd2HeAnalysis::omega(s,pow(m1,2),pow(m2,2))*ATd2HeAnalysis::omega(s,pow(m3,2),pow(m4_ex,2))) ) ;




	       THcm = theta_cm*TMath::RadToDeg();
         MissingMass =  Ex;


}


void ATd2HeAnalysis::FitTrack(TGraph2D *f2dtrack, TGraph *fpadtrack, std::vector<Double_t> *iniguess ){



	           ROOT::Fit::Fitter  fitter;
   		// make the functor objet
   		   SumDistance2 sdist(f2dtrack);

       	           ROOT::Math::Functor fcn(sdist,4);
   		   // set the function and the initial parameter values



	           Double_t pStart[4] = {iniguess->at(0),iniguess->at(1),iniguess->at(2),iniguess->at(3)};
   		   fitter.SetFCN(fcn,pStart);

		   fitter.Config().ParSettings(0).SetLimits(-10.0,10.0);
		   fitter.Config().ParSettings(1).SetLimits(-20,20);
		   fitter.Config().ParSettings(2).SetLimits(-10,1000);
		   fitter.Config().ParSettings(3).SetLimits(-20,20);

   		   // set step sizes different than default ones (0.3 times parameter values)
   		   for (int i = 0; i < 4; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);
   		   ok = fitter.FitFCN();

   		   if (!ok) {
      		        Error("line3Dfit","Line3D Fit failed");
   		   }
   		   const ROOT::Fit::FitResult & result = fitter.Result();

		   // get fit parameters
   		   const double * parFit = result.GetParams();

                    vertex = parFit[2];

                   theta_fit = fabs(atan(sqrt(parFit[1]*parFit[1] +1)/parFit[3] )*TMath::RadToDeg());
                   if(vertex>iniguess->at(5)) theta_fit = 180.0- theta_fit;


                   TFitResultPtr  res1 = fpadtrack->Fit("pol1","Sq");
                   phi_fit = atan2(signo(iniguess->at(4))*signo(res1->Value(1))*fabs(res1->Value(1)),signo(iniguess->at(4)))*TMath::RadToDeg();
                   if(phi_fit<0) phi_fit+=360;


}


TMatrixD ATd2HeAnalysis::CalcV(){

    TMatrixD Vval(16,16);
  	Vval.Zero();
  	TMatrixD Vsub1(4,4);
  	Vsub1.Zero();
    TMatrixD Vsub2(4,4);
  	Vsub2.Zero();
  	TMatrixD Vsub3(4,4);
  	Vsub3.Zero();
  	TMatrixD Vsub4(4,4);
  	Vsub4.Zero();

  	for(int i=0; i<4;i++){
  		for(int j=0; j<4;j++){
        if(i!=j) continue;
  			if(i==0){
  				 Vsub1[i][j] = 0.1;
           Vsub2[i][j] = 0.1;
  				 Vsub3[i][j] = 0.1;
  				 Vsub4[i][j] = 0.1;
  				}
  			if(i==1){
  				 Vsub1[i][j] = 0.1;
           Vsub2[i][j] = 0.1;
  				 Vsub3[i][j] = 0.1;
  				 Vsub4[i][j] = 0.1;
  				}
  			if(i==2){
  				 Vsub1[i][j] = 0.1;
           Vsub2[i][j] = 0.1;
  				 Vsub3[i][j] = 0.1;
  				 Vsub4[i][j] = 0.1;
  				}

  			if(i==3){
  				 Vsub1[i][j] = 0.1;
           Vsub2[i][j] = 0.1;
  				 Vsub3[i][j] = 0.1;
  				 Vsub4[i][j] = 0.1;
  				}
          /*if(i!=j){
            Vsub1[i][j] = -0.00001;
            Vsub2[i][j] = -0.00001;
   				  Vsub3[i][j] = -0.00001;
   				  Vsub4[i][j] = -0.00001;
          }
          */
  		}
  	}

  	Vval.SetSub(0,0,Vsub1);
    Vval.SetSub(4,4,Vsub2);
  	Vval.SetSub(8,8,Vsub3);
  	Vval.SetSub(12,12,Vsub4);

    //anticorrelation between the two protons to keep the relative energy approx. the same
    Vval[15][11] = -0.9*0.1; //Ener
    Vval[11][15] = -0.9*0.1;

    //Vval.Print();

  	return Vval;

}

TMatrixD ATd2HeAnalysis::CalcV_new(){

    TMatrixD Vval(12,12);
  	Vval.Zero();
  	TMatrixD Vsub1(4,4);
  	Vsub1.Zero();
    TMatrixD Vsub2(4,4);
  	Vsub2.Zero();
  	TMatrixD Vsub3(4,4);
  	Vsub3.Zero();

  	for(int i=0; i<3;i++){
  		for(int j=0; j<3;j++){
        if(i!=j) continue;
  			if(i==0){
  				 Vsub1[i][j] = 0.1;
           Vsub2[i][j] = 0.1;
  				 Vsub3[i][j] = 0.1;
  				}
  			if(i==1){
  				 Vsub1[i][j] = 0.1;
           Vsub2[i][j] = 0.1;
  				 Vsub3[i][j] = 0.1;
  				}
  			if(i==2){
  				 Vsub1[i][j] = 0.1;
           Vsub2[i][j] = 0.1;
  				 Vsub3[i][j] = 0.1;
  				}
  		}
  	}

  	Vval.SetSub(0,0,Vsub1);
    Vval.SetSub(4,4,Vsub2);
  	Vval.SetSub(8,8,Vsub3);

    //anticorrelation between the two protons to keep the relative energy approx. the same
    //Vval[15][11] = -0.9*0.1; //Ener
    //Vval[11][15] = -0.9*0.1;

    //Vval.Print();

  	return Vval;

}


TMatrixD ATd2HeAnalysis::CalcD(TMatrixD alpha0){

	TMatrixD Dval(4,16);
	Dval.Zero();
	TMatrixD Dsub1(4,4);
	Dsub1.Zero();
  TMatrixD Dsub2(4,4);
	Dsub2.Zero();
	TMatrixD Dsub3(4,4);
	Dsub3.Zero();
	TMatrixD Dsub4(4,4);
	Dsub4.Zero();

	for(int i=0; i<4;i++){
		for(int j=0; j<4;j++){
			if(i!=j) continue;
			if(i==0){
				 Dsub1[i][j] = 1;
         Dsub2[i][j] = -1;
				 Dsub3[i][j] = -1;
				 Dsub4[i][j] = -1;
				}
			if(i==1){
				 Dsub1[i][j] = 1;
         Dsub2[i][j] = -1;
				 Dsub3[i][j] = -1;
				 Dsub4[i][j] = -1;
				}
			if(i==2){
				 Dsub1[i][j] = 1;
         Dsub2[i][j] = -1;
				 Dsub3[i][j] = -1;
				 Dsub4[i][j] = -1;
				}

			if(i==3){
				 Dsub1[i][j] = 1;
         Dsub2[i][j] = -1;
				 Dsub3[i][j] = -1;
				 Dsub4[i][j] = -1;
				}
		}
	}

  /*
	TVector3 p1 (alpha0[0][0],alpha0[1][0],alpha0[2][0]);
	TVector3 p3 (alpha0[4][0],alpha0[5][0],alpha0[6][0]);
	TVector3 p4 (alpha0[8][0],alpha0[9][0],alpha0[10][0]);
	double theta3 = p1.Angle(p3);
	double m1 = Li6_mass;
	double m3 = Li6_mass;
	double m4 = C12_mass;
  */
	//Dsub1[4][0] = -2*p1.X(); Dsub1[4][1] = -2*p1.Y(); Dsub1[4][2] = -2*p1.Z(); Dsub1[4][3] = 2*alpha0[3][0];
	//Dsub3[4][0] = -2*p3.X(); Dsub3[4][1] = -2*p3.Y(); Dsub3[4][2] = -2*p3.Z(); Dsub3[4][3] = 2*alpha0[10][0];
	//Dsub4[4][0] = -2*p4.X(); Dsub4[4][1] = -2*p4.Y(); Dsub4[4][2] = -2*p4.Z(); Dsub4[4][3] = 2*alpha0[17][0];


	/*Dsub3[5][0] = 2*alpha0[7][0];
	Dsub3[5][1] = 2*alpha0[8][0];
	Dsub3[5][2] = 2*alpha0[8][0];
	Dsub3[5][3] = -2*alpha0[10][0];
	*/

	Dval.SetSub(0,0,Dsub1);
  Dval.SetSub(0,4,Dsub2);
	Dval.SetSub(0,8,Dsub3);
	Dval.SetSub(0,12,Dsub4);

	return Dval;


}


TMatrixD ATd2HeAnalysis::Calcd(TMatrixD alpha0){

	TMatrixD dval(4,1);
	dval.Zero();

	/*double m1 = Li6_mass;
	double m2 = C12_mass;
	double m3 = Li6_mass;
	double m4 = C12_mass;
	TVector3 p1 (alpha0[0][0],alpha0[1][0],alpha0[2][0]);
	TVector3 p3 (alpha0[4][0],alpha0[5][0],alpha0[6][0]);
	TVector3 p4 (alpha0[8][0],alpha0[9][0],alpha0[10][0]);
	double theta3 = p1.Angle(p3);
  */

  double mt = 2.01410177812 * 931.494; //target mass (2H)


	double h1 = alpha0[0][0] - alpha0[4][0] - alpha0[8][0] - alpha0[12][0]; //momentum conservation X
	double h2 = alpha0[1][0] - alpha0[5][0] - alpha0[9][0] - alpha0[13][0]; //momentum conservation Y
	double h3 = alpha0[2][0] - alpha0[6][0] - alpha0[10][0] - alpha0[14][0]; //momentum conservation Z
	double h4 = alpha0[3][0] +  mt - alpha0[7][0] - alpha0[11][0] - alpha0[15][0]; //Total Energy conservation
	//double h5 = sq(alpha0[3][0]) - p1.Mag2() - sq(m1) /*+ sq(alpha0[10][0]) - p3.Mag2() - sq(m3) + sq(alpha0[17][0]) - p4.Mag2() - sq(m4)*/  ;
	//double h6 = alpha0[7][0]*alpha0[7][0] + alpha0[8][0]*alpha0[8][0]  + alpha0[9][0]*alpha0[9][0] + sq(m3) - sq(alpha0[10][0]);

	dval[0][0] = h1;
	dval[1][0] = h2;
	dval[2][0] = h3;
	dval[3][0] = h4;
	//dval[4][0] = h5;
	//dval[5][0] = h6;

	return dval;

}

TMatrixD ATd2HeAnalysis::CalcD_new(TMatrixD alpha0){

	TMatrixD Dval(2,12);
	Dval.Zero();

  TVector3 Pb (alpha0[0][0],alpha0[1][0],alpha0[2][0]);
	TVector3 P1 (alpha0[4][0],alpha0[5][0],alpha0[6][0]);
	TVector3 P2 (alpha0[8][0],alpha0[9][0],alpha0[10][0]);
  Double_t Eb = alpha0[3][0];
  Double_t E1 = alpha0[7][0];
  Double_t E2 = alpha0[11][0];
  Double_t mt = 2.01410177812 * 931.494; //target mass (2H)


  Dval[0][0] = 0; Dval[0][1] = 0; Dval[0][2] = 0;  Dval[0][3] = 0; //dH1/dPbx, dH1/dPbx, dH1/dPbx, dH1/dEb
  Dval[0][4] = -2*(P1.X()+P2.X()); Dval[0][5] = -2*(P1.Y()+P2.Y()); Dval[0][6] = -2*(P1.Z()+P2.Z());  Dval[0][7] = 2*(E1+E2); //dH1/dP1x, dH1/dP1x, dH1/dP1x, dH1/dE1
  Dval[0][8] = -2*(P1.X()+P2.X()); Dval[0][9] = -2*(P1.Y()+P2.Y()); Dval[0][10] = -2*(P1.Z()+P2.Z());  Dval[0][11] = 2*(E1+E2); //dH1/dP2x, dH1/dP2x, dH1/dP2x, dH1/dE2

  Dval[1][0] = -2*(Pb.X()-P1.X()-P2.X()); Dval[1][1] = -2*(Pb.Y()-P1.Y()-P2.Y()); Dval[1][2] = -2*(Pb.Z()-P1.Z()-P2.Z());  Dval[1][3] = 2*(Eb+mt-E1-E2); //dH2/dPbx, dH2/dPbx, dH2/dPbx, dH2/dEb
  Dval[1][4] = 2*(Pb.X()-P1.X()-P2.X()); Dval[1][5] = 2*(Pb.Y()-P1.Y()-P2.Y()); Dval[1][6] = 2*(Pb.Z()-P1.Z()-P2.Z());  Dval[1][7] = -2*(Eb+mt-E1-E2); //dH2/dP1x, dH2/dP1x, dH2/dP1x, dH2/dE1
  Dval[1][8] = 2*(Pb.X()-P1.X()-P2.X()); Dval[1][9] = 2*(Pb.Y()-P1.Y()-P2.Y()); Dval[1][10] = 2*(Pb.Z()-P1.Z()-P2.Z());  Dval[1][11] = -2*(Eb+mt-E1-E2); //dH2/dP2x, dH2/dP2x, dH2/dP2x, dH2/dE2


	return Dval;


}


TMatrixD ATd2HeAnalysis::Calcd_new(TMatrixD alpha0, Double_t ex_2he, Double_t ex_eje){

	TMatrixD dval(2,1);
	dval.Zero();


  TVector3 Pb (alpha0[0][0],alpha0[1][0],alpha0[2][0]);
	TVector3 P1 (alpha0[4][0],alpha0[5][0],alpha0[6][0]);
	TVector3 P2 (alpha0[8][0],alpha0[9][0],alpha0[10][0]);
  Double_t Eb = alpha0[3][0];
  Double_t E1 = alpha0[7][0];
  Double_t E2 = alpha0[11][0];
  Double_t mt = 2.01410177812 * 931.494; //target mass (2H)
  Double_t m2he = 2.0*1.0078250322 * 931.494 + ex_2he; //2He mass
  Double_t meje = 14.00307400443 * 931.494 + ex_eje; //ejectile mass

  TVector3 mom2He = P1 + P2;
  TVector3 momEje = Pb - P1 - P2;

	double h1 = pow((E1+E2),2.0) - mom2He.Mag2() - pow(m2he,2.0); //2He MissingMass
	double h2 = pow((Eb+mt-E1-E2),2.0) - momEje.Mag2() - pow(meje,2.0); //Ejectile MissingMass

	dval[0][0] = h1;
	dval[1][0] = h2;

	return dval;

}


void ATd2HeAnalysis::KFit(std::vector<double> pin, std::vector<double>& pout){

    //Kinematic Fitting based on https://www.phys.ufl.edu/~avery/fitting.html
    //Input vector is composed by 4 four-momentum vectors of (projectile,ejectile,proton1,proton2)
    //Output vector contains the fitted coordinates in the same order

    TMatrixD alpha0(16,1);
  	alpha0.Zero();
  	TMatrixD alpha1(4,1);
  	alpha1.Zero();
    TMatrixD alpha2(4,1);
    alpha2.Zero();
  	TMatrixD alpha3(4,1);
  	alpha3.Zero();
  	TMatrixD alpha4(4,1);
  	alpha4.Zero();

  	alpha1[0][0] = pin[0]; alpha1[1][0] = pin[1]; alpha1[2][0] = pin[2]; alpha1[3][0] = pin[3]; //projectile (px,py,pz,E)
    alpha2[0][0] = pin[4]; alpha2[1][0] = pin[5]; alpha2[2][0] = pin[6]; alpha2[3][0] = pin[7]; //ejectile (px,py,pz,E)
  	alpha3[0][0] = pin[8]; alpha3[1][0] = pin[9]; alpha3[2][0] = pin[10]; alpha3[3][0] = pin[11]; //proton 1 (px,py,pz,E)
  	alpha4[0][0] = pin[12]; alpha4[1][0] = pin[13]; alpha4[2][0] = pin[14]; alpha4[3][0] = pin[15]; //proton 2 (px,py,pz,E)
  	alpha0.SetSub(0,0,alpha1);
    alpha0.SetSub(4,0,alpha2);
  	alpha0.SetSub(8,0,alpha3);
  	alpha0.SetSub(12,0,alpha4);

  	TMatrixD alpha(16,1);
  	TMatrixD delta_alpha(16,1);
  	delta_alpha.Zero();
  	alpha = alpha0;

  	TMatrixD chi2(1,1);
  	chi2.Zero();

  	TMatrixD Va = CalcV();


    double weighting=0.05;
	  int niter = 100;

		for(int w =0; w<niter; w++){


			TMatrixD D = CalcD(alpha);
			//D.Print();
			TMatrixD Dt = D;
			Dt.T();

			TMatrixD Vd = D*Va*Dt;
			//Vd.Print();
			Vd.Invert();
			//Vd.Print();

			// calculate Va
			Va = Va - Va*Dt*Vd*D*Va*weighting;

			//calculate labmda
			TMatrixD d = Calcd(alpha);
			TMatrixD temp1 = (d + D*delta_alpha);
			TMatrixD lambda = Vd*temp1;


			// calculate chi2
			TMatrixD lambdaT = lambda;
			chi2 = lambdaT.T()*temp1;

      //calculate alpha
			alpha = alpha0 - weighting*Va*Dt*lambda;
			delta_alpha = alpha - alpha0;
			alpha0 = alpha;

			if(fabs(chi2[0][0])<1.0) break;
			//if(fabs(chisqrd)<0.1) break;

		}

    for(int w=0;w<16;w++) pout.push_back(alpha[w][0]);


}

void ATd2HeAnalysis::KFit_new(std::vector<double> pin, std::vector<double>& pout){

    //Kinematic Fitting based on https://www.phys.ufl.edu/~avery/fitting.html
    //Input vector is composed by 3 four-momentum vectors of (projectile,proton1,proton2)
    //Output vector contains the fitted coordinates in the same order

    TMatrixD alpha0(12,1);
  	alpha0.Zero();
  	TMatrixD alpha1(4,1);
  	alpha1.Zero();
    TMatrixD alpha2(4,1);
    alpha2.Zero();
  	TMatrixD alpha3(4,1);
  	alpha3.Zero();

  	alpha1[0][0] = pin[0]; alpha1[1][0] = pin[1]; alpha1[2][0] = pin[2]; alpha1[3][0] = pin[3]; //projectile (px,py,pz,E)
    alpha2[0][0] = pin[4]; alpha2[1][0] = pin[5]; alpha2[2][0] = pin[6]; alpha2[3][0] = pin[7]; //proton 1 (px,py,pz,E)
  	alpha3[0][0] = pin[8]; alpha3[1][0] = pin[9]; alpha3[2][0] = pin[10]; alpha3[3][0] = pin[11]; //proton 2 (px,py,pz,E)

  	alpha0.SetSub(0,0,alpha1);
    alpha0.SetSub(4,0,alpha2);
  	alpha0.SetSub(8,0,alpha3);

    TVector3 Pb (alpha0[0][0],alpha0[1][0],alpha0[2][0]);
  	TVector3 P1 (alpha0[4][0],alpha0[5][0],alpha0[6][0]);
  	TVector3 P2 (alpha0[8][0],alpha0[9][0],alpha0[10][0]);
    Double_t Eb = alpha0[3][0];
    Double_t E1 = alpha0[7][0];
    Double_t E2 = alpha0[11][0];
    Double_t mt = 2.01410177812 * 931.494; //target mass (2H)
    Double_t m2he_gs = 2.0*1.0078250322 * 931.494 ; //2He mass
    Double_t meje_gs = 14.00307400443 * 931.494 ; //ejectile mass
    TVector3 mom2He = P1 + P2;
    TVector3 momEje = Pb - P1 - P2;
    Double_t ex_he2 = sqrt( pow(E1+E2,2.0) - mom2He.Mag2()) - m2he_gs;
  	//Double_t ex_eje = sqrt( pow(Eb+mt-E1-E2,2.0) - momEje.Mag2()) - meje_gs;

    //Double_t ex_he2 = 0;
  	Double_t ex_eje = 3.9;


  	TMatrixD alpha(12,1);
  	TMatrixD delta_alpha(12,1);
  	delta_alpha.Zero();
  	alpha = alpha0;

  	TMatrixD chi2(1,1);
  	chi2.Zero();

  	TMatrixD Va = CalcV_new();


    double weighting=0.5;
	  int niter = 100;

		for(int w =0; w<niter; w++){


			TMatrixD D = CalcD_new(alpha);
			//D.Print();
			TMatrixD Dt = D;
			Dt.T();

			TMatrixD Vd = D*Va*Dt;
			//Vd.Print();
			Vd.Invert();
			//Vd.Print();

			// calculate Va
			Va = Va - Va*Dt*Vd*D*Va*weighting;

			//calculate labmda
			TMatrixD d = Calcd_new(alpha,ex_he2,ex_eje);
			TMatrixD temp1 = (d + D*delta_alpha);
			TMatrixD lambda = Vd*temp1;


			// calculate chi2
			TMatrixD lambdaT = lambda;
			chi2 = lambdaT.T()*temp1;

      //calculate alpha
			alpha = alpha0 - weighting*Va*Dt*lambda;
			delta_alpha = alpha - alpha0;
			alpha0 = alpha;

      TVector3 Pb_new (alpha0[0][0],alpha0[1][0],alpha0[2][0]);
    	TVector3 P1_new (alpha0[4][0],alpha0[5][0],alpha0[6][0]);
    	TVector3 P2_new (alpha0[8][0],alpha0[9][0],alpha0[10][0]);
      Double_t Eb_new = alpha0[3][0];
      Double_t E1_new = alpha0[7][0];
      Double_t E2_new = alpha0[11][0];
      TVector3 mom2He_new = P1_new + P2_new;
      TVector3 momEje_new = Pb_new - P1_new - P2_new;
      Double_t ex_he2_new = sqrt( pow(E1_new+E2_new,2.0) - mom2He_new.Mag2()) - m2he_gs;
    	Double_t ex_eje_new = sqrt( pow(Eb_new+mt-E1_new-E2_new,2.0) - momEje_new.Mag2()) - meje_gs;
      ex_he2 = ex_he2_new;
      //ex_eje = ex_eje_new;

			if(fabs(chi2[0][0])<1.0) break;
			//if(fabs(chisqrd)<0.1) break;

		}

    for(int w=0;w<12;w++) pout.push_back(alpha[w][0]);


}
