#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

//---------------CONSTANTs initialization----------------------------
//--------------- Network Geometry ------------------------------------
#define I_TC    1     //  0 - No layer, 1 - Add Layer
#define I_RE    1     //  0 - No layer, 1 - Add Layer
#define I_CX    1     //  0 - No layer, 1 - Add Layer
#define I_IN    1     //  0 - No layer, 1 - Add Layer
#define I_GB    0     //  0 - No GABAb from IN to CX, 1 - Yes

#define M        50      
#define Mcx      100     
#define Min      25      
#define M1       1     
#define Mcx1     1     
#define Min1     1     
#define Max     ((M > Mcx) ? M : Mcx)
#define Max1    ((M1 > Mcx1) ? M1 : Mcx1)

//-------------Boundary conditions ---------------------------------------
#define BOUND      0   
#define BOUNDcx    0   
#define SELFING    0  
#define SELFINGcx  0   

//-------------- Define the connections between cells -----------------
#define MS_RE_RE 5
#define MS_RE_RE1 0
#define MS_RE_RE_MAX  ((MS_RE_RE > MS_RE_RE1) ? MS_RE_RE : MS_RE_RE1)
#define N_RE_RE  (2*MS_RE_RE+1)      //*(2*MS_RE_RE1+1)

#define MS_RE_TC 5 
#define MS_RE_TC1 0
#define MS_RE_TC_MAX ((MS_RE_TC > MS_RE_TC1) ? MS_RE_TC : MS_RE_TC1)
#define N_RE_TC  (2*MS_RE_TC+1)      //*(2*MS_RE_TC1+1)

#define MS_TC_RE 5 
#define MS_TC_RE1 0
#define MS_TC_RE_MAX ((MS_TC_RE > MS_TC_RE1) ? MS_TC_RE : MS_TC_RE1)
#define N_TC_RE  (2*MS_TC_RE+1)  

#define MS_CX_CX 5
#define MS_CX_CX1 0
#define MS_CX_CX_MAX  ((MS_CX_CX > MS_CX_CX1) ? MS_CX_CX : MS_CX_CX1)
#define N_CX_CX  (2*MS_CX_CX+1)         //*(2*MS_CX_CX1+1)

#define MS_CX_CX_NMDA 5
#define MS_CX_CX_NMDA1 0
#define MS_CX_CX_NMDA_MAX  ((MS_CX_CX_NMDA > MS_CX_CX_NMDA1) ? MS_CX_CX_NMDA : MS_CX_CX_NMDA1)
#define N_CX_CX_NMDA  (2*MS_CX_CX_NMDA+1)       //*(2*MS_CX_CX_NMDA1+1)

#define MS_CX_IN_NMDA 1  
#define MS_CX_IN_NMDA1 0
#define MS_CX_IN_NMDA_MAX  ((MS_CX_IN_NMDA > MS_CX_IN_NMDA1) ? MS_CX_IN_NMDA : MS_CX_IN_NMDA1)
#define N_CX_IN_NMDA (2*MS_CX_IN_NMDA+1)*Mcx/Min  //*(2*MS_CX_IN_NMDA1+1)*Mcx1/Min1

#define MS_CX_IN 1 
#define MS_CX_IN1 0
#define MS_CX_IN_MAX  ((MS_CX_IN > MS_CX_IN1) ? MS_CX_IN : MS_CX_IN1)
#define N_CX_IN  (2*MS_CX_IN+1)*Mcx/Min //*(2*MS_CX_IN1+1) * Mcx/Min * Mcx1/Min1

#define MS_IN_CX 5
#define MS_IN_CX1 0
#define MS_IN_CX_MAX  ((MS_IN_CX > MS_IN_CX1) ? MS_IN_CX : MS_IN_CX1)
#define N_IN_CX  ((2*MS_IN_CX+1)+4)*Min/Mcx  //* (2*MS_IN_CX1+1) *Min1/Mcx1

#define MS_TC_CX 10
#define MS_TC_CX1 0
#define MS_TC_CX_MAX  ((MS_TC_CX > MS_TC_CX1) ? MS_TC_CX : MS_TC_CX1)
#define N_TC_CX  (2*MS_TC_CX+1)  //*(2*MS_TC_CX1+1) * M1/Mcx1

#define MS_TC_IN 2
#define MS_TC_IN1 0
#define MS_TC_IN_MAX  ((MS_TC_IN > MS_TC_IN1) ? MS_TC_IN : MS_TC_IN1)
#define N_TC_IN  (2*MS_TC_IN+1)* M/Min //*(2*MS_TC_IN1+1) * M1/Min1

#define MS_CX_TC 5 
#define MS_CX_TC1 0
#define MS_CX_TC_MAX  ((MS_CX_TC > MS_CX_TC1) ? MS_CX_TC : MS_CX_TC1)
#define N_CX_TC  (2*MS_CX_TC+1)* Mcx/M //*(2*MS_CX_TC1+1) * Mcx1/M1

//------------Number of ODE for each cell -------------------------------
#define N_RE 7
#define N_TC 12 
#define N_GB 2

#define N_DEND   9
#define N_SOMA   4 //3
#define N_CX     (N_DEND + N_SOMA)
#define N_IN     N_CX   //4

#define N_EQ1  (N_RE*I_RE + N_TC*I_TC + N_RE_TC*N_GB*I_RE*I_TC)*M*M1
#define N_EQ2  (N_CX*I_CX + N_IN_CX*N_GB*I_CX*I_IN*I_GB)*Mcx*Mcx1
#define N_EQ3  (N_IN*I_IN)*Min*Min1
#define N_EQ   (N_EQ1 + N_EQ2 + N_EQ3)    //  Complete number of ODE

//++++++++++++++ Current DESCRIPTION ++++++++++++++++++++++++++++++++++++++
//---------------Low-threshold Ca2+ current (RE cell)---------------------
class IT_RE {
  static double Shift, Ca_0, Cels;
  double m_inf, tau_m, h_inf, tau_h, ratio, eca, Phi_m, Phi_h, eca0;
public:
  double iT, m0, h0, Qm, Qh;
  double G_Ca;
  IT_RE(double v) {
    G_Ca = 1.75;
    Qm = 5; //2.5; //3; 
    Qh = 3; //2.5; //5;
    Phi_m = pow(Qm,((Cels-24)/10));
    Phi_h = pow(Qh,((Cels-24)/10));
    m0 = 1/(1 + exp(-(v + Shift + 50)/7.4));
    h0 = 1/(1 + exp((v + Shift + 78)/5));
    eca0 = 1000*8.31441*(273.15 + Cels)/(2*96489);
    } 
  void calc(double m, double h, double &fm, double &fh, 
            double v, double cai, double x);
};

double IT_RE::Shift = 2, IT_RE::Ca_0 = 2, IT_RE::Cels = 36;

void IT_RE::calc(double m, double h, double &fm, double &fh,
                 double v, double cai, double x) {
  ratio = Ca_0/cai;
    if(ratio <= 0.)printf("\n LOG ERROR: RE: cai=%lf ratio=%lf",cai,ratio);
  eca = eca0 * log(ratio);
  iT = G_Ca*m*m*h*(v - eca);                              
  m_inf = 1/(1 + exp(-(v + 52)/7.4));
  tau_m = (3 + 1/(exp((v + 27)/10) + exp(-(v + 102)/15)))/Phi_m;
  h_inf = 1/(1 + exp((v + 80)/5));
  tau_h = (85 + 1/(exp((v + 48)/4) + exp(-(v + 407)/50)))/Phi_h;
  fm = -(m - m_inf)/tau_m;                                  
  fh = -(h - h_inf)/tau_h;                                  
}


//--------------fast Na and K current (RE and TC cells)------------------
class INaK {
  static double Cels;
  double Alpha1, Beta1, Alpha2, Beta2, Alpha3, Beta3, v2, v2K, Phi;
  double tau_m, m_inf, tau_h, h_inf, tau_n, n_inf;
public:
  static double E_Na, E_K;
  double iK, iNa, m0, h0, n0, mNa_TCRE, hNa_TCRE;
  double G_Na, G_K, Vtr, VtrK, Alpha_m, Beta_m, Alpha_h, Beta_h, tau_m_na,m_inf_na,tau_h_na,h_inf_na;
  INaK(double v) {
    G_K = 10;///////////////////////
    G_Na = 100;/////////////////////
    Vtr = -50;
    VtrK = -50;
    v2 = v - Vtr;
    v2K = v - VtrK;
    Phi = pow(3,((Cels-36)/10));
    Alpha1 = 0.32*(13 - v2)/(exp((13 - v2)/4) - 1);
    Beta1 = 0.28*(v2 - 40)/(exp((v2 - 40)/5) - 1);
    m0 = Alpha1/(Alpha1 + Beta1);

    Alpha2 = 0.128*exp((17 - v2)/18);
    Beta2 = 4/(exp((40 - v2)/5) + 1);
    h0 = Alpha2/(Alpha2 + Beta2);

    Alpha3 = 0.032*(15 - v2K)/(exp((15 - v2K)/5) - 1);
    Beta3 = 0.5*exp((10 - v2K)/40);
    n0 = Alpha3/(Alpha3 + Beta3);     } 
  void calc(double m, double h, double n, double &fm, double &fh, double &fn, 
            double v, double x);
};

double INaK::E_K = -95, INaK::E_Na = 50, INaK::Cels = 36; 

void INaK::calc(double m, double h, double n, double &fm, double &fh, double &fn,
                   double v, double x){
  v2 = v - Vtr;
  v2K = v - VtrK;
  iNa = G_Na*m*m*m*h*(v - E_Na);
  Alpha1 = 0.32*(13 - v2)/(exp((13 - v2)/4) - 1);
  Beta1 = 0.28*(v2 - 40)/(exp((v2 - 40)/5) - 1);
  tau_m = 1/(Alpha1 + Beta1) / Phi;
  m_inf = Alpha1/(Alpha1 + Beta1);

  Alpha2 = 0.128*exp((17 - v2)/18);
  Beta2 = 4/(exp((40 - v2)/5) + 1);
  tau_h = 1/(Alpha2 + Beta2) / Phi;
  h_inf = Alpha2/(Alpha2 + Beta2);
  
  Alpha_m=Alpha1;
  Beta_m=Beta1;
  tau_m_na=tau_m;
  m_inf_na=m_inf;
  tau_h_na=tau_h;
  h_inf_na=h_inf;
  Alpha_h=Alpha2;
  Beta_h=Beta2;

  fm = -(m - m_inf)/tau_m;                 
  fh = -(h - h_inf)/tau_h;       
  
  mNa_TCRE=m;   
  hNa_TCRE=h;       

  iK = G_K* n*n*n*n*(v - E_K);    
  Alpha3 = 0.032*(15 - v2K)/(exp((15 - v2K)/5) - 1);
  Beta3 = 0.5*exp((10 - v2K)/40);
  tau_n = 1/(Alpha3 + Beta3) / Phi;
  n_inf = Alpha3/(Alpha3 + Beta3);
  
  fn  = -(n - n_inf)/tau_n;                 
}

//------------------Ca-dynamics------------------------------------
class ICa {
  static double Ca_inf, K_T, K_d;
  double drive, drive0;                                 
public:
  double Taur, D, iCa, drive_save;
  ICa() {Taur = 5; D = 1.; //0.1;
         drive0 = 10.0/(2.*96489.); }
  void calc(double cai, double &fcai, double iT, double x);
};

double ICa::Ca_inf = 2.4e-4;
double ICa::K_T = 0.0001, ICa::K_d = 0.0001;

void ICa::calc(double cai, double &fcai, double iT, double x) {
  drive = -drive0 * iT / D;
  if(drive < 0) drive = 0;
  fcai = drive + (Ca_inf - cai)/Taur; // - K_T*cai/(cai + K_d);
  iCa = cai;
  drive_save=drive;
}

//------------------Low-theshold Ca2+ current (TC cell)-----------------
class IT_TC {
  static double Ca_0, Cels, Qm, Qh, Shift;
  double m_inf, tau_m, h_inf, tau_h, Phi_h, Phi_m, ratio, eca, eca0; 
public:
  double iT, m0, h0;
  double G_Ca;
  IT_TC(double v) {
     G_Ca = 2; //2;///////////////////////////////////
     Phi_m = pow(Qm,((Cels-24)/10));
     Phi_h = pow(Qh,((Cels-24)/10));
     m0 = 1 / (1+exp(-(v+59)/6.2));////////////////////////
     h0 = 1 / (1+exp((v+83)/4.0));
     eca0 = 1000*8.31441*(273.15 + Cels)/(2*96489);  } 
  void calc(double m, double h, double &fm, double &fh,  
            double v, double cai, double x);
};

double IT_TC::Shift = 2, IT_TC::Ca_0 = 2, IT_TC::Cels = 36;
double IT_TC::Qm = 3.55, IT_TC::Qh = 3; //2.8;

void IT_TC::calc(double m, double h, double &fm, double &fh,
                 double v, double cai, double x) {
  ratio = Ca_0/cai;
    if(ratio <= 0.)printf("\n LOG ERROR: RE: cai=%lf ratio=%lf",cai,ratio);
  eca = eca0 * log(ratio);
  iT = G_Ca*m*m*h*(v - eca); 

  m_inf = 1 / (1+exp(-(v+59)/6.2));//////////////////////////////////////
  h_inf = 1 / (1+exp((v+83)/4.)); /////////////////////////////////////

  tau_m = (1/(exp(-(v+131.6)/16.7)+exp((v+16.8)/18.2)) + 0.612) / Phi_m;
  tau_h = (30.8 + (211.4 + exp((v + Shift + 113.2)/5))/
           (1+exp((v + Shift + 84)/3.2))) / Phi_h;

  fm = -(1/tau_m)*(m - m_inf);                                
  fh = -(1/tau_h)*(h - h_inf);
}

//----------------- h-current (TC cell) -----------------------------------
class Ih_TC {
  static double E_h, Shift, Cels, k2, nca, nexp, taum;
  double h_inf, tau_s, alpha, beta, k3p, cc, Phi;
public:
  double ih, p10, o10, o20;
  double G_h, k1ca, ginc, cac, pc, k4;

  Ih_TC(double v, double cai) {
     G_h = 0.02; //////////////
     ginc = 1.5;  
     cac = 0.0015;
     pc = 0.01;
     k4 = 0.001;
     Phi = pow(3,((Cels-36)/10));
     h_inf = 1/(1 + exp((v + 75 - Shift)/5.5));
     tau_s = (taum + 1000 / (exp((v + 71.5 - Shift)/14.2) + 
                          exp(-(v + 89 - Shift)/11.6))) / Phi;
     alpha = h_inf/tau_s;
     beta = (1 - h_inf)/tau_s;
     p10 = 1/(1 + pow((cac/cai),nca));
     o10 = 1/(1 + beta/alpha + pow((p10/pc),nexp));
     o20 = pow((p10/pc),nexp) * o10;
  }
  void calc(double o1, double p1, double o2,  
                 double &fo1, double &fp1, double &fo2,  
                 double v, double cai, double x);
};

double Ih_TC::E_h = -40; 
double Ih_TC::Shift = 0, Ih_TC::Cels = 36;
double Ih_TC::k2 = 0.0004;
double Ih_TC::nca = 4, Ih_TC::nexp = 1, Ih_TC::taum = 20;

void Ih_TC::calc(double o1, double p1, double o2, double &fo1, double &fp1, 
                 double &fo2, double v, double cai, double x) {
  ih = G_h*(o1 + ginc * o2)*(v - E_h);
  h_inf = 1/(1 + exp((v + 75 - Shift)/5.5));
  tau_s = (taum + 1000 / (exp((v + 71.5 - Shift)/14.2) + 
                          exp(-(v + 89 - Shift)/11.6))) / Phi;
  alpha = h_inf/tau_s;
  beta = (1 - h_inf)/tau_s;
  k1ca = k2 * pow((cai/cac),nca);
  k3p = k4 * pow((p1/pc),nexp);
  fo1 = alpha * (1-o1-o2) - beta * o1; // + k4 * o2 - k3p * o1;
  fp1 = k1ca * (1-p1) - k2 * p1; // + k4 * o2 - k3p * o1;
  fo2 = k3p * o1 - k4 * o2;
}

//----------------------Potassium A-current (TC cell)------------------------
class IA_TC {
  static double E_K, Cels;     
  double m_inf, tau_m, h_inf, tau_h, Tad;                                 
public:
  double iA, m0, h0, G_A;
  IA_TC(double v) {
    G_A = 1; //2; //////////////////////////////////////////////
    Tad = pow(3,((Cels-23.5)/10));
    m0 = 1.0 / (1+exp(-(v+60)/8.5));
    h0 = 1.0/(1+exp((v+78)/6)); } 
  void calc(double m, double h, double &fm, double &fh, double v, double x);
};

double IA_TC::Cels = 36, IA_TC::E_K = -95;

void IA_TC::calc(double m, double h, double &fm, double &fh, double v, double x){
  iA = G_A*m*m*m*m*h*(v - E_K);
  tau_m = (1.0/( exp((v+35.82)/19.69)+exp(-(v+79.69)/12.7) ) +0.37) / Tad;
  m_inf = 1.0 / (1+exp(-(v+60)/8.5));
  tau_h = 1.0/((exp((v+46.05)/5)+exp(-(v+238.4)/37.45))) / Tad;
  if(v >= -63) 
      tau_h = 19.0/Tad;
  h_inf = 1.0/(1+exp((v+78)/6));
  fm = -(1/tau_m)*(m - m_inf);                                  
  fh = -(1/tau_h)*(h - h_inf);
}

//---------------------Hight-threshold Ca2+ current (CX cell)----------------
class IHVA_CX {
  static double Shift, Ca_0, Cels, Qm, Qh, E_Ca;
  double m_inf, tau_m, h_inf, tau_h, Phi_m, Phi_h, a, b, vm;
  double ratio, eca0, eca;
public:
  double iHVA, m0, h0, m_hva, h_hva;
  double G_HVA;
  IHVA_CX(double v) {
    G_HVA = 0.03;
    Phi_m = pow(Qm,((Cels-23)/10));
    Phi_h = pow(Qh,((Cels-23)/10));
    vm = v + Shift;
    a = 0.055*(-27 - vm)/(exp((-27-vm)/3.8) - 1);
    b = 0.94*exp((-75-vm)/17);
    m0 = a/(a+b);
    a = 0.000457*exp((-13-vm)/50);
    b = 0.0065/(exp((-vm-15)/28) + 1);
    h0 = a/(a+b);
    } 
  void calc(double m, double h, double &fm, double &fh, double v, double cai, double x);
};

double IHVA_CX::Shift = 0, IHVA_CX::Ca_0 = 2, IHVA_CX::E_Ca = 140; 
double IHVA_CX::Qm = 2.3, IHVA_CX::Qh = 2.3, IHVA_CX::Cels = 36;

void IHVA_CX::calc(double m, double h, double &fm, double &fh,
                 double v, double cai, double x) {
//------------ECa is fixed (=140mV) instead of using Nerst eq.-----------------

  iHVA = Phi_m * G_HVA * m*m*h * (v - E_Ca);
  vm = v + Shift;

  a = 0.055*(-27 - vm)/(exp((-27-vm)/3.8) - 1);
  b = 0.94*exp((-75-vm)/17);
  tau_m = (1/(a+b))/Phi_m;
  m_inf = a/(a+b);

  a = 0.000457*exp((-13-vm)/50);
  b = 0.0065/(exp((-vm-15)/28) + 1);
  tau_h = (1/(a+b))/Phi_h;
  h_inf = a/(a+b);

  fm = -(m - m_inf)/tau_m;                                  
  fh = -(h - h_inf)/tau_h;
  
  m_hva=m;
  h_hva=h;
}

//--------------Ca-dependent potassium current (CX cell)-----------------------
class IKCa_CX {
  static double E_KCa, Ra, Rb, Cels, Q, caix;     
  double m_inf, tau_m, Tad, a, b;                                 
public:
  double iKCa, m0;
  double G_KCa;
  IKCa_CX(double cai) {
    G_KCa = 0.3; 
    Tad = pow(Q,((Cels-23)/10));
    a = Ra * cai;  //------becouse caix = 1
    b = Rb;
    m0 = a/(a+b);
  }
  void calc(double m, double &fm, double v, double cai, double x);
};

double IKCa_CX::E_KCa = -90, IKCa_CX::Q = 2.3, IKCa_CX::caix = 1;
double IKCa_CX::Ra = 0.01, IKCa_CX::Rb = 0.02, IKCa_CX::Cels = 36;   

void IKCa_CX::calc(double m, double &fm, double v, double cai, double x){
  iKCa = Tad * G_KCa * m * (v - E_KCa);                         

//  a = Ra * pow(cai,caix);
  a = Ra * cai;  
  b = Rb;
  tau_m = (1/(a+b))/Tad;
  m_inf = a/(a+b);
  fm = -(1/tau_m)*(m - m_inf);                   
}

//--------------------Potassium M-current (CX cell)-------------------------
class IKm_CX {
  static double E_Km, Ra, Rb, Cels, Q, tha, qa;     
  double m_inf, tau_m, Tad, a, b;                                 
public:
  double iKm, m0;
  double G_Km;
  IKm_CX(double v) {
    G_Km = 0.01; 
    Tad = pow(Q,((Cels-23)/10));
    a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
    b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
    m0 = a/(a+b);
  }
  void calc(double m, double &fm, double v, double x);
};

double IKm_CX::E_Km = -90, IKm_CX::Q = 2.3;
double IKm_CX::tha = -30, IKm_CX::qa = 9;
double IKm_CX::Ra = 0.001, IKm_CX::Rb = 0.001, IKm_CX::Cels = 36;   

void IKm_CX::calc(double m, double &fm, double v, double x){
  iKm = Tad * G_Km * m * (v - E_Km);
  a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
  b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
  tau_m = (1/(a+b))/Tad;
  m_inf = a/(a+b);
  fm = -(1/tau_m)*(m - m_inf);                   
}

//--------------------Fast potassium current (CX cell)-------------------
class IKv_CX {
  static double Ra, Rb, Cels, Q, tha, qa;     
  double m_inf, tau_m, Tad, a, b;                                 
public:
  static double E_Kv;
  double iKv, g_Kv, m0;
  double G_Kv;
  IKv_CX(double v) {
    G_Kv = 150; 
    Tad = pow(Q,((Cels-23)/10));

    a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
    b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
    m0 = a/(a+b);
  }
  void calc(double m, double &fm, double v, double x);
};

double IKv_CX::E_Kv = -90, IKv_CX::Q = 2.3;
double IKv_CX::tha = 25 /*25*/, IKv_CX::qa = 9;
double IKv_CX::Ra = 0.02, IKv_CX::Rb = 0.002, IKv_CX::Cels = 36;   

void IKv_CX::calc(double m, double &fm, double v, double x){
  g_Kv = Tad * G_Kv * m;
  iKv = g_Kv * (v - E_Kv);                         
  a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
  b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
  tau_m = (1/(a+b))/Tad;
  m_inf = a/(a+b);
  fm = -(1/tau_m)*(m - m_inf);                   
}

//----------------Fast sodium current (CX cell)--------------------------
class INa_CX {
  static double Shift, Ca_0, Cels, Qm, Qh;
  static double tha, qa, Ra, Rb, thi1, thi2, qi, thinf, qinf, Rg, Rd; 
  double m_inf, tau_m, h_inf, tau_h, Phi_m, Phi_h, a, b, vm;
  double trap0(double v, double th, double a, double q) {
        if (fabs(v/th) > 1.0e-6) {
                return ( a * (v - th) / (1 - exp(-(v - th)/q)) );
        } else {
	  return (a * q ); }
        }
public:
  static double E_Na;
  double iNa, g_Na, m0, h0, m_Na, h_Na, v_m_E_Na;
  double G_Na;
  INa_CX(double v) {
    G_Na = 3000;
    Phi_m = pow(Qm,((Cels-23)/10));
    Phi_h = pow(Qh,((Cels-23)/10));
    vm = v + Shift;
    a = trap0(vm,tha,Ra,qa);
    b = trap0(-vm,-tha,Rb,qa);
    m0 = a/(a+b);

    a = trap0(vm,thi1,Rd,qi);
    b = trap0(-vm,-thi2,Rg,qi);
    h0 = 1/(1+exp((vm-thinf)/qinf));
    } 
  void calc(double m, double h, double &fm, double &fh, 
            double v, double x);
};

double INa_CX::Shift = -10, INa_CX::E_Na = 50; 
double INa_CX::Qm = 2.3, INa_CX::Qh = 2.3, INa_CX::Cels = 36;
double INa_CX::tha = -35, INa_CX::qa = 9;
double INa_CX::Ra = 0.182,INa_CX::Rb = 0.124; 
double INa_CX::thi1 = -50, INa_CX::thi2 = -75, INa_CX::qi = 5;
double INa_CX::thinf = -65, INa_CX::qinf = 6.2;
double INa_CX::Rg = 0.0091, INa_CX::Rd = 0.024; 

void INa_CX::calc(double m, double h, double &fm, double &fh,
                 double v, double x) {

  g_Na = Phi_m * G_Na * m*m*m*h;
  iNa = g_Na * (v - E_Na);
  v_m_E_Na = (v - E_Na);
  vm = v + Shift;

  a = trap0(vm,tha,Ra,qa);
  b = trap0(-vm,-tha,Rb,qa);
  tau_m = (1/(a+b))/Phi_m;
  m_inf = a/(a+b);

                //"h" inactivation 
  a = trap0(vm,thi1,Rd,qi);
  b = trap0(-vm,-thi2,Rg,qi);
  tau_h = (1/(a+b))/Phi_h;
  h_inf = 1/(1+exp((vm-thinf)/qinf));

  fm = -(m - m_inf)/tau_m;                                  
  fh = -(h - h_inf)/tau_h;
  m_Na = m;
  h_Na = h;
}

//-------------------Persist Sodium current (CX cell)-------------------
class INap_CX {
  double Tet, Sig, f, Cels, Q10, Phi_m, tau_m, m_inf;    
public:
  double m0, iNap, G_Nap, g_Nap;
  INap_CX(double v) {
    G_Nap = 2; 
    Tet = -42;
    Sig = 5;  
    f = 0.02;
    Cels = 36;
    Q10 = 2.7;
    Phi_m = pow(Q10,((Cels-22)/10));      
    tau_m = 0.8/Phi_m;
    m0 = f/(1 + exp(-(v - Tet)/Sig));
  }
  void calc(double, double&, double, double);
  double cond();
};

void INap_CX::calc(double m, double &fm, double v, double x){
  g_Nap = G_Nap * m;
  iNap = g_Nap * (v - INa_CX::E_Na);
  m_inf = f/(1 + exp(-(v - Tet)/Sig));
  fm = -(m - m_inf)/tau_m; 
}
double INap_CX::cond(){
  return(m_inf);
}


//===================Now we'll CREATE the CELLS==================================
//-------------------RE CELL-----------------------------------------------
class RE: public IT_RE, 
          public INaK, public ICa {
  static double Cai0, V0;
public:
  double G_kl, G_l, E_l, S_RE, iKl, iL, v_RE;

  RE() :IT_RE(V0), INaK(V0), ICa() {
        E_l = -70;
        G_l = 0.05;
        G_kl = 0.018; //0.012; //0.015;
        S_RE = 1.43e-4;
        }
  void init(double *y) {
        y[0] = V0;
	y[1] = Cai0;
        y[2] = IT_RE::m0;
        y[3] = IT_RE::h0;
        y[4] = INaK::m0;
        y[5] = INaK::h0;
        y[6] = INaK::n0;
        }
  void calc(double x, double *y, double *f); 
};  

double RE::V0 = -61, RE::Cai0 = 0.0001;

void RE::calc(double x, double *y, double *f){
    IT_RE::calc(y[2], y[3], f[2], f[3], y[0], y[1], x);
    INaK::calc(y[4], y[5], y[6], f[4], f[5], f[6], y[0], x);
    ICa::calc(y[1], f[1], iT, x);
    f[0] = -G_l * (y[0] - E_l) - iT - iNa - iK
                         - G_kl * (y[0] - INaK::E_K);
    iL = -G_l * (y[0] - E_l);
    iKl = - G_kl * (y[0] - INaK::E_K);
    v_RE = y[0];
}

//-------------------TC CELL-----------------------------------------------
class TC: public IT_TC, public Ih_TC, public INaK, public ICa, public IA_TC {
  static double G_l, Cai0, V0;
public:
  double G_kl, S_TC, E_l, DC, iL, iKl, v_TC;

  TC() :Ih_TC(V0,Cai0), IT_TC(V0), INaK(V0), ICa(), IA_TC(V0) {
        G_kl = 0.012;
        E_l = -70;
        DC = 0;
        INaK::G_Na = 90;
        INaK::G_K = 10;
        S_TC = 2.9e-4;
        }
  void init(double *y) {
        y[0] = V0;
	y[1] = Cai0;
        y[2] = IT_TC::m0;
        y[3] = IT_TC::h0;
        y[4] = Ih_TC::o10;
        y[5] = Ih_TC::p10;
        y[6] = Ih_TC::o20;
        y[7] = INaK::m0;
        y[8] = INaK::h0;
        y[9] = INaK::n0;
        y[10] = IA_TC::m0;
        y[11] = IA_TC::h0;
        }
  void calc(double x, double *y, double *f); 
};  

double TC::Cai0 = 0.0001;
double TC::G_l = 0.01;//0.05; //0.01;///
double TC::V0 = -68;

void TC::calc(double x, double *y, double *f){
    IT_TC::calc(y[2], y[3], f[2], f[3], y[0], y[1], x);
    Ih_TC::calc(y[4], y[5], y[6], f[4], f[5], f[6], y[0], y[1], x);
    INaK::calc(y[7], y[8], y[9], f[7], f[8], f[9], y[0], x); 
    IA_TC::calc(y[10], y[11], f[10], f[11], y[0], x);
    ICa::calc(y[1], f[1], iT, x);
    f[0] = -G_l * (y[0] - E_l) - iT - ih - iNa - iK - iA
              - G_kl * (y[0] - INaK::E_K) + DC;
    iL = -G_l * (y[0] - E_l);
    iKl = - G_kl * (y[0] - INaK::E_K);
    v_TC = y[0];
}

//-------------------CX CELL------------------------------------------------
//-------------------CX CELL (DENDRITE)-------------------------------------
class CX_DEND: public IHVA_CX, public IKCa_CX, public IKm_CX, 
               public INa_CX, public INap_CX, public ICa {
  static double G_l;
public:
  double iDEND, v_DEND, I_Stim1, E_l, G_kl, iKl, iL;
  CX_DEND(double V0, double Cai0) :IHVA_CX(V0), IKCa_CX(Cai0), IKm_CX(V0),
                                   INa_CX(V0), INap_CX(V0), ICa() { 
        E_l = -70; //-67;
        G_kl = 0; 
        INa_CX::G_Na = 0.8; //1.5; 
        I_Stim1 = 0; 
        ICa::Taur = 165; //200;
        INap_CX::G_Nap = 3.5; //2.5; //2; //3; //1.1; //1.1; //0.8;
        IKm_CX::G_Km = 0.01; //0.01;  
        IKCa_CX::G_KCa = 0.3; //0.3;
        IHVA_CX::G_HVA = 0.01; //0.02; //0.015; //0.03;
        }
  void init(double V0, double Cai0, double *y) {
        y[0] = V0;
	y[1] = Cai0;
        y[2] = IHVA_CX::m0;
        y[3] = IHVA_CX::h0;
        y[4] = IKCa_CX::m0;
        y[5] = IKm_CX::m0;
        y[6] = INa_CX::m0;
        y[7] = INa_CX::h0;
        y[8] = INap_CX::m0;
  }
  void calc(double x, double *y, double *f); 
};  
double CX_DEND::G_l = 1.0e3/30000;  //  mS/cm^2

void CX_DEND::calc(double x, double *y, double *f){
    IHVA_CX::calc(y[2], y[3], f[2], f[3], y[0], y[1], x);
    IKCa_CX::calc(y[4], f[4], y[0], y[1], x);
    IKm_CX::calc(y[5], f[5], y[0], x);
    INa_CX::calc(y[6], y[7], f[6], f[7], y[0], x);
    ICa::calc(y[1], f[1], iHVA, x);
    INap_CX::calc(y[8], f[8], y[0], x);
    iDEND =  -G_l * (y[0] - E_l) - iHVA - iKCa - iKm - iNa -iNap + I_Stim1
             -G_kl * (y[0] - INaK::E_K);
    iKl = -G_kl * (y[0] - INaK::E_K);
    iL = -G_l * (y[0] - E_l);
}

//-------------------CX CELL (SOMA)-------------------------------------
class CX_SOMA: public IKv_CX, public INa_CX, public INap_CX {
public:
  double v_SOMA, iSOMA, g1_SOMA, g2_SOMA, I_Stim2;
  CX_SOMA(double V0, double Cai0) :IKv_CX(V0), INa_CX(V0), INap_CX(V0){ 
        I_Stim2 = 0; 
        IKv_CX::G_Kv = 200; //150;
        INa_CX::G_Na = 3000; 
        INap_CX::G_Nap = 15;
        }
  void init(double V0, double Cai0, double *y) {
        v_SOMA = V0;
        y[0] = IKv_CX::m0;
        y[1] = INa_CX::m0;
        y[2] = INa_CX::h0;
        y[3] = INap_CX::m0;
        }
  void calc(double x, double *y, double *f); 
};  

void CX_SOMA::calc(double x, double *y, double *f){
    IKv_CX::calc(y[0], f[0], v_SOMA, x);
    INa_CX::calc(y[1], y[2], f[1], f[2], v_SOMA, x);
    INap_CX::calc(y[3], f[3], v_SOMA, x);
    g1_SOMA = g_Na + g_Kv + g_Nap;
    g2_SOMA = g_Na * INa_CX::E_Na + g_Kv * IKv_CX::E_Kv + 
              g_Nap * INa_CX::E_Na + 6.74172 + I_Stim2; 
    iSOMA =  - iNa - iKv - iNap;     
}

//------------CX CELL (connect DENDRITE and SOMA)---------------------------
class CX: public CX_DEND, public CX_SOMA {
  static double Cai0, V0, C;
public:
  double kappa, rho, S_CX_SOMA, S_CX_DEND, v_DEND, iGAP ; 
  CX() :CX_DEND(V0,Cai0), CX_SOMA(V0,Cai0) {
        kappa = 10.0e3;      // kOm: to get mS=1/kOm 
        rho = 165; //50; //200; //165;
        }
  void init(double *y) {
        CX_DEND::init(V0, Cai0, y);
	CX_SOMA::init(V0, Cai0, y+N_DEND);
        S_CX_SOMA = 1.0e-6;  // cm^2
        S_CX_DEND = S_CX_SOMA * rho;        
        }
  void calc(double x, double *y, double *f); 
};  

double CX::Cai0 = 0.0001, CX::V0 = -68;
double CX::C = 0.75;   // uF/cm^2

void CX::calc(double x, double *y, double *f){
    CX_SOMA::calc(x, y+N_DEND, f+N_DEND);
    v_SOMA = (y[0] + kappa * S_CX_SOMA * g2_SOMA) / 
                            (1 + kappa*S_CX_SOMA * g1_SOMA);
    CX_DEND::calc(x, y, f);
    f[0] = (1.0/C) * ( iDEND + 1.0/(kappa*S_CX_DEND) * (v_SOMA - y[0]) );
    v_DEND = y[0];
    iGAP = 1.0/(kappa*S_CX_DEND) * (v_SOMA - y[0]);
}

//---------SYNAPCES DESCRIPTION-------------------------------------------
//---------first order kinet model for GABA-A synapse---------------------
class Gaba_A {
  static double Cdur, Cmax, Deadtime, Prethresh;  
  double R, C, R0, R1;
  double lastrelease;
  double q, Rinf, Rtau;
  double exptable(double z)
    {
    if((z > -10) && (z < 10)) return( exp(z) );
    else return( 0 );
    }
public:
  double I, E_GABA, Alpha, Beta;
  Gaba_A() {
    E_GABA = -70; //-75; (-70 is more realistic?)
    R = 0, C = 0, R0 = 0, R1 = 0;
    Alpha = 10.5;
    Beta = 0.166;
    lastrelease = -100;
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
    Rtau = 1 / ((Alpha * Cmax) + Beta);
  }
  void calc(double g_GABA_A, double x, double y_pre, double y_post);
}; 
double Gaba_A::Cdur = 0.3, Gaba_A::Cmax = 0.5, Gaba_A::Deadtime = 1;
double Gaba_A::Prethresh = 0;
void Gaba_A::calc(double g_GABA_A, double x, double y_pre, 
                    double y_post) {

        q = ((x - lastrelease) - Cdur);         
        if (q > Deadtime) {
                if (y_post > Prethresh) {        
                        C = Cmax;                
                        R0 = R;
                        lastrelease = x;
                }
        } else if (q < 0) {                     
        } else if (C == Cmax) {                  
                R1 = R;
                C = 0.;
        }
        if (C > 0) {                            
           R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
        } else {                              
           R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
        }
       I = g_GABA_A * R * (y_pre - E_GABA);
}

//------second order kinet model (including G-proteins) for GABA-B synapse----
class Gaba_B {
  static double E_GABA, Cmax, Deadtime, Prethresh;   
  static double Kd, n; 
  double Gn, q;
public:
  double C, lastrelease;
  double I, r0, g0, Gn1, Cdur, K1, K2, K1K2, K3, K4, g_rec,r_rec, fr_rec,fg_rec;
  Gaba_B() {
    Cdur = 0.3; 
    K1 = 0.52; 
    K2 = 0.0013; 
    K3 = 0.098;
    K4 = 0.033;
    lastrelease = -10000000;
    C = 0, r0 = 0, g0 = 0;
  }
  void calc(double r, double g, double &fr, double &fg, 
            double g_GABA_B, double x, double y_pre, double y_post);
};
double Gaba_B::E_GABA = -95, Gaba_B::Cmax = 0.5, Gaba_B::Deadtime = 1;
double Gaba_B::Prethresh = 0;
double Gaba_B::Kd = 100, Gaba_B::n = 4; 

void Gaba_B::calc(double r, double g, double &fr, double &fg, 
                  double g_GABA_B, double x, double y_pre, double y_post) {
        g_rec=g;
        r_rec=r;
        Gn = pow(g,n); 
        Gn1 = Gn/(Gn + Kd);
        I = g_GABA_B * Gn1 * (y_pre - E_GABA);

        q = ((x - lastrelease) - Cdur);         
        if (q > Deadtime) {
                if (y_post > Prethresh) {        
                        C = Cmax;                
                        lastrelease = x; }
        } else if (q < 0) {                     
        } else if (C == Cmax) {   C = 0;  }
        fr = K1 * C * (1 - r) - r * K2;
        fg = K3 * r - K4 * g;
        fr_rec=fr;
        fg_rec=fg;
}

class GB: public Gaba_B {
public:
  GB() :Gaba_B() { }
  void init(double *y){
      lastrelease = -10000000;
      C = 0;
      y[0] = 0;
      y[1] = 0;
      }   
  void calc(double g_GABA_B, double x, double *y, double *f, 
                                             double y_pre, double y_post){
       Gaba_B::calc(y[0], y[1], f[0], f[1], g_GABA_B, x, y_pre, y_post); 
       } 
};  

//------------first order kiner model for AMPA synapse---------------------
class AMPA {
  static double E_AMPA;
  static double Cdur, Cmax, Deadtime, Prethresh, Cdel; 
  static double Alpha, Beta; 
  int s;
  double R, C, R0, R1;
  double lastrelease, lastspike;
  double q, Rinf, Rtau;
  double exptable(double z)
    {
    if((z > -10) && (z < 10)) return( exp(z) );
    else return( 0 );
    }
public:
  double I;
  AMPA() {
    R = 0, C = 0, R0 = 0, R1 = 0;
    lastrelease = -100;
    lastspike = -100;
    s = 1;
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
    Rtau = 1 / ((Alpha * Cmax) + Beta);
  }
  void calc(double g_AMPA, double x, double y_pre, double y_post);
};
double AMPA::E_AMPA = 0, AMPA::Cdur = 0.3, AMPA::Cmax = 0.5, AMPA::Deadtime = 1;
double AMPA::Cdel = 0;
double AMPA::Prethresh = 0, AMPA::Alpha = 0.94, AMPA::Beta = 0.18;
void AMPA::calc(double g_AMPA, double x, double y_pre, 
                    double y_post) {

        q = ((x - lastrelease) - Cdur); 
        
        if(q > Deadtime) {
               if(y_post > Prethresh) {
                  if( (x - lastspike) > (Cdel + Cdur) ){
                     lastspike = x;
                     s = 1; } }  //the flag that spike was but wasn't utilized yet

               if((s == 1) && ((x - lastspike) > Cdel)) {
                  s = 0;         //spike was utilized
                  C = Cmax;                
                  R0 = R;
                  lastrelease = x;  }
        } else if (q < 0) {                     
        } else if (C == Cmax) {                  
                R1 = R;
                C = 0.;
        }

        if (C > 0) {                            
           R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
        } else {                              
           R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
        }
       I = g_AMPA * R * (y_pre - E_AMPA);
}



//--first order kiner model for AMPA synapse WITH depression & spont releases------
class AMPA_D2 {
  static double E_AMPA;
  static double Cdur, Cmax, Deadtime, Prethresh, Cdel; 
  static double Alpha, Beta; 
  int s;
  double C, R0, R1;
  double lastrelease, lastrelease1, lastrandom, newrelease, Use, Tr;
  double q, q1, Rinf, Rtau;
  double Tau, Period, S, SS, factor;
  double exptable(double z)
    {
    if((z > -10) && (z < 10)) return( exp(z) );
    else return( 0 );
    }
public:
  double R, I, E, g1;
  AMPA_D2() {
    R = 0, C = 0, R0 = 0, R1 = 0;
    lastrelease = -10000;
    lastrelease1 = -10000;
    newrelease = 0;
    s = 1;
    g1=0.00006;
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
    Rtau = 1 / ((Alpha * Cmax) + Beta);
    E = 1;
    Use = 0.07; 
    Tr = 700;
    Period = 8000;
    Tau = 50; //50;
    factor = 1;
  }
  void iii(unsigned int seek) {srand(seek);}
  void calc(double, double, double, double, double);
};
double AMPA_D2::E_AMPA = 0, AMPA_D2::Cdur = 0.3;
double AMPA_D2::Cmax = 0.5, AMPA_D2::Deadtime = 1;
double AMPA_D2::Cdel = 0;
double AMPA_D2::Prethresh = 0, AMPA_D2::Alpha = 0.94, AMPA_D2::Beta = 0.18;
void AMPA_D2::calc(double g_AMPA, double g_AMPAmin, double x, double y_pre, 
                    double y_post) {

        q = ((x - lastrelease) - Cdur);
        q1 = ((x - lastrelease1) - Cdur);
        if(q > Deadtime) {
               if(y_post > Prethresh) {
                  g1 = g_AMPA;
                  factor = 1;
                  Use = 0.073; //0.08;
                  C = Cmax;                
                  R0 = R;
                  E = 1 - (1 - E*(1-Use)) * exptable(-q1/Tr);
                  lastrelease = x;
                  lastrelease1 = x;
                  }
                else if( ((x - lastrelease1) > 70.0) && 
                         ((x - lastrelease) > newrelease) ) {
                  SS = log((x - lastrelease1 +Tau)/Tau)/400;
                  S = rand()/(RAND_MAX + 1.0);
                  if(S < 0.000001) S = 0.000001;
                  newrelease = -(log(S))/SS;
                  g1 = g_AMPAmin;
                  factor = 2;
                  Use = 0; //0.01;
                  C = Cmax;                
                  R0 = R;
                  E = 1 - (1 - E*(1-Use)) * exptable(-q1/Tr);
                  lastrelease = x;
                  }
               }
        else if (q < 0) {                     
              } 
        else if (C == Cmax) {                  
                  R1 = R;
                  C = 0.;
        }

        if (C > 0) {                            
           R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
        } else {                              
           R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
        }
       I = g1 * R * E * (y_pre - E_AMPA);
}

//------------first order kiner model for NMDA synapse WITH depression------
class NMDA_D1 {
  static double E_NMDA;
  static double Cdur, Cmax, Deadtime, Prethresh, Cdel; 
  static double Alpha, Beta; 
  int s;
  double R, C, R0, R1;
  double lastrelease, lastspike, Use, Tr;
  double q, Rinf, Rtau, fn;
  double exptable(double z)
    {
    if((z > -10) && (z < 10)) return( exp(z) );
    else return( 0 );
    }
public:
  double I, E;
  NMDA_D1() {
    R = 0, C = 0, R0 = 0, R1 = 0;
    lastrelease = -100;
    lastspike = -100;
    s = 1;
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
    Rtau = 1 / ((Alpha * Cmax) + Beta);
    E = 1;
    Use = 0.0; 
    Tr = 750;
  }
  void iii(unsigned int seek) {srand(seek);}
  void calc(double g_NMDA, double x, double y_pre, double y_post);
};
double NMDA_D1::E_NMDA = 0, NMDA_D1::Cdur = 0.3, NMDA_D1::Cmax = 0.5;
double NMDA_D1::Deadtime = 1;
double NMDA_D1::Cdel = 0; 
double NMDA_D1::Prethresh = 0, NMDA_D1::Alpha = 1, NMDA_D1::Beta = 0.0067;

void NMDA_D1::calc(double g_NMDA, double x, double y_pre, double y_post) {

        q = ((x - lastrelease) - Cdur); 
        
        if(q > Deadtime) {
               if(y_post > Prethresh) {
                  if( (x - lastspike) > (Cdel + Cdur) ){
                     lastspike = x;
                     s = 1; } }  //the flag that spike was but wasn't utilized

               if((s == 1) && ((x - lastspike) > Cdel)) {
                  s = 0;         //spike was utilized
                  C = Cmax;                
                  R0 = R;
                  E = 1 - (1 - E*(1-Use)) * exptable(-q/Tr);
                  lastrelease = x;  }
        } else if (q < 0) {                     
        } else if (C == Cmax) {                  
                R1 = R;
                C = 0.;
        }

        if (C > 0) {                            
           R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
        } else {                              
           R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
        }
       fn = 1/(1+exp(-(y_pre - (-25))/12.5));
       I = g_NMDA * R * fn * E * (y_pre - E_NMDA);
}




//---first order kiner model for GABA-A synapse with DEPRESSION & spont IPSPs--
class Gaba_A_D2 {
  static double Cdur, Cmax, Deadtime, Prethresh; 
  static double Alpha, Beta;  
  double R, C, R0, R1;
  double lastrelease, lastrelease1, lastrandom, newrelease, Use, Tr;
  double q, q1, Rinf, Rtau;
  double Tau, Period, S, SS, factor;
  double exptable(double z)
    {
    if((z > -10) && (z < 10)) return( exp(z) );
    else return( 0 );
    }
public:
  double I, E_GABA, E;
  Gaba_A_D2() {
    E_GABA = -70; //-75; (-70 is more realistic?)
    R = 0, C = 0, R0 = 0, R1 = 0;
    lastrelease = -10000;
    lastrelease1 = -10000;
    newrelease = 0;
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
    Rtau = 1 / ((Alpha * Cmax) + Beta);
    E = 1;
    Use = 0; //0.07;
    Tr = 700;
    Period = 8000;
    Tau = 50; //50;
    factor = 1;
  }
  void iii(unsigned int seek) {srand(seek);}
  void calc(double g_GABA_A, double x, double y_pre, double y_post);
}; 
double Gaba_A_D2::Cdur = 0.3, Gaba_A_D2::Cmax = 0.5, Gaba_A_D2::Deadtime = 1;
double Gaba_A_D2::Prethresh = 0, Gaba_A_D2::Alpha = 10, Gaba_A_D2::Beta =0.25; 
void Gaba_A_D2::calc(double g_GABA_A, double x, double y_pre, 
                    double y_post) {

        q = ((x - lastrelease) - Cdur);
        q1 = ((x - lastrelease1) - Cdur);
        if (q > Deadtime) {
               if (y_post > Prethresh) {        
                  factor = 1;
                  Use = 0.07;
                  C = Cmax;                
                  R0 = R;
                  E = 1 - (1 - E*(1-Use)) * exptable(-q1/Tr);
                  lastrelease = x;
                  lastrelease1 = x;
               }
               else if( ((x - lastrelease1) > 70.0) && 
                        ((x - lastrelease) > newrelease) ){
                  SS = log((x - lastrelease1 +Tau)/Tau)/400;
                  S = rand()/(RAND_MAX + 1.0);
                  if(S < 0.000001) S = 0.000001;
                  newrelease = -(log(S))/SS;

                  factor = 10; //5;
                  Use = 0; //0.01;
                  C = Cmax;                
                  R0 = R;
                  E = 1 - (1 - E*(1-Use)) * exptable(-q1/Tr);
                  lastrelease = x;
                  }
               } 
        else if (q < 0) { } 
        else if (C == Cmax) {                  
                R1 = R;
                C = 0.;
        }
        if (C > 0) {                            
           R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
        } else {                              
           R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
        }
       I = (g_GABA_A/factor) * E * R * (y_pre - E_GABA);
}

//-----first order kiner model for AMPA synapse used for external stimulation----
class Extern_ampa {
  static double Cdur, Cmax, Deadtime, Prethresh; 
  double R, C, R0, R1;
  double lastrelease;
  double q, Rinf, Rtau;
  double TR, w, wom, RRR;
  double exptable(double z)
    {
    if((z > -10) && (z < 10)) return( exp(z) );
    else return( 0 );
    }
public:
  double g, Alpha, Beta;
  Extern_ampa() {
    Alpha = 0.94;
    Beta = 0.18;
    R = 0, C = 0, R0 = 0, R1 = 0;
    lastrelease = -100;
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
    Rtau = 1 / ((Alpha * Cmax) + Beta);
    TR = 1000, w=0.01, wom=0;
  }
  void iii(unsigned int seek) {srand(seek);}
  void calc(double g_Extern_ampa, double x);
};
double Extern_ampa::Cdur = 0.3, Extern_ampa::Cmax = 0.5, Extern_ampa::Deadtime = 1;
double Extern_ampa::Prethresh = 0;
void Extern_ampa::calc(double g_Extern_ampa, double x) {

        q = ((x - lastrelease) - Cdur);         
        if (q > Deadtime) {
                if ((x - lastrelease) > TR) {        
                        C = Cmax;                
                        R0 = R;
                        lastrelease = x;
                }
        } else if (q < 0) {                     
        } else if (C == Cmax) {                  
                R1 = R;
                C = 0.;
        }
        if (C > 0) {                            
           R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
        } else {                              
           R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
        }
       g = g_Extern_ampa * R;
}

//+++++++++++++++++++ MAIN PROGRAM +++++++++++++++++++++++++++++++++++++++++++

//----------external functions----------------------------------------------
void rk(unsigned, void (unsigned, double, double*, double*), double, double, 
                                double*, double*, double*, double*);
void fun(unsigned, double, double*, double*);
void solout (long, double, double, double*, unsigned, int*);

//----------external variables ---------------------------------------------
double g_gaba_a, g_gaba_a1, g_gaba_b, g_ampa;
double *g_ampa_cx_cx[Mcx][Mcx1], g_nmda_cx_cx, g_nmda_cx_in, g_ampa_cx_tc;
double g_ampa_cx_re, g_ampa_tc_cx, g_ampa_cx_cx_MINICE, g_ampa_cx_in_MINICE; 
double *g_ampa_cx_in[Min][Min1], *g_gaba_a_in_cx[Mcx][Mcx1];
double g_gaba_b_in_cx, g_ampa_tc_in; 
double g_ext_tc, g_ext_re, g_ext_cx, g_ext_in;
FILE *f1, *f2, *f3, *f4, *f6, *f7, *f8, *f9, *f10, *f11, *f12, *f13, *f14, *f15, *f16, *f17, *f18, *f19, *f20, *f21, *f22, *f23, *f24, *f25, *f26, *f27, *f28, *f29, *f30, *f31, *f32, *f33, *f34, *f35, *f36, *f37, *f38, *f39, *f40, *f41, *f42, *f43, *f44, *f45, *f46, *f47, *f48, *f49, *f50, *f51, *f52, *f53, *f54, *f55, *f56, *f57, *f58, *f59, *f60, *f61, *f62, *f63, *f64, *f65, *f66, *f67, *f68, *f69 , *f70 , *f71, *f72, *f73, *f74, *f75, *f76, *f77, *f78, *f79, *f80, *f81, *f82, *f83, *f84, *f85, *f86, *f87, *f88, *f89, *f90, *f91, *f92, *f93, *f94, *f95, *f96, *f97, *f98 ;

int no_re[M][M1][N_RE], no_tc[M][M1][N_TC], no_g[M][M1][N_RE_TC][N_GB];
int no_cx[Mcx][Mcx1][N_CX], no_in[Min][Min1][N_IN], no_gcx[Mcx][Mcx1][N_IN_CX][N_GB];

int C_RERE[M][M1][M][M1], C_RETC[M][M1][M][M1];
int C_TCRE[M][M1][M][M1];
int C_CXCX[Mcx][Mcx1][Mcx][Mcx1], C_CXIN[Mcx][Mcx1][Min][Min1];
int C_INCX[Min][Min1][Mcx][Mcx1];
int C_CXTC[Mcx][Mcx1][M][M1], C_CXRE[Mcx][Mcx1][M][M1];
int C_TCCX[M][M1][Mcx][Mcx1], C_TCIN[M][M1][Mcx][Mcx1];

int k_RERE[M][M1], k_RETC[M][M1];
int k_TCRE[M][M1];
int k_CXCX[Mcx][Mcx1], k_CXIN[Min][Min1];
int k_INCX[Mcx][Mcx1];
int k_CXTC[M][M1], k_CXRE[M][M1];
int k_TCCX[Mcx][Mcx1], k_TCIN[Mcx][Mcx1];

int k_REREmax=0, k_RETCmax=0;
int k_TCREmax=0;
int k_CXCXmax=0, k_CXINmax=0;
int k_INCXmax=0;
int k_CXTCmax=0, k_CXREmax=0;
int k_TCCXmax=0, k_TCINmax=0;

//----------external classes (beginning of initialization)------------------
Gaba_A       *g_a[M][M1];
Gaba_A       *g_a1[M][M1];
GB           *g_b[M][M1];
AMPA         *a_a[M][M1];
RE           re_cell[M][M1];
TC           tc_cell[M][M1];

AMPA_D2         *a_cx_cx[Mcx][Mcx1];
NMDA_D1         *nmda_cx_cx[Mcx][Mcx1];
NMDA_D1         *nmda_cx_in[Min][Min1];
AMPA_D2         *a_cx_in[Min][Min1];
Gaba_A_D2       *ga_in_cx[Mcx][Mcx1];
GB           *gb_in_cx[Mcx][Mcx1];

CX           cx_cell[Mcx][Mcx1];
CX           in_cell[Min][Min1];         

AMPA         *a_cx_tc[M][M1];
AMPA         *a_cx_re[M][M1];
AMPA         *a_tc_cx[Mcx][Mcx1];
AMPA         *a_tc_in[Min][Min1];

Extern_ampa  a_ext1, a_ext2, a_ext3, a_ext4;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int main(int argc,char **argv)
{
//---------allocate place for ALL variables and functions-----------
  double y_ini[N_EQ], f_ini[N_EQ];

//---------allocate place for TWO temporal arrays using by RK.C solver 
  double y1[N_EQ], y2[N_EQ];

//---------general parameters----------------------------------------------
  double t = 0, tmax, tmax1, t3D, ttime, TAU;
  double h = 0.02, red, t_ext, g_Extern_ampa1, R, r[2], av=0;
  int i, j, i1, j1, k, l, ii = 0, i_inp, i_gip = 0, ih, ni;
  double g_GABA_A, g_GABA_A1, g_GABA_B, g_AMPA, scale;
  double g_AMPA_CX_CX, g_NMDA_CX_CX, g_NMDA_CX_IN, g_AMPA_CX_TC, g_AMPA_CX_RE, g_AMPA_TC_CX, g_AMPA_CX_CX_MINICE, g_AMPA_CX_IN_MINICE;
  double g_AMPA_CX_IN, g_GABA_A_IN_CX, g_GABA_B_IN_CX, g_AMPA_TC_IN, g_Extern_ampa;

//---------solver parameters (DOPRI)-------------------------------------------
  int      res, iout = 2, itoler = 0; 
  double   atoler = 1.0E-5, rtoler = 1.0E-5; 

//----------arrays initialization----------------------------------------------
  for(i=0; i<N_EQ; i++){
    y_ini[i] = 0, f_ini[i] = 0; 
    y1[i] = 0, y2[i] = 0; }

//-------parameter initialization (from file)----------------------------------
int ret;

if (argc <= 1) {
  puts("Command parameters");
  puts("-----------------------");
  puts("Input File"); }

if (!(f1=fopen(argv[1],"r"))) {
   printf("%s doesn't exist\n",argv[1]);
   exit(0); }

  ret = fscanf(f1, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", 
       &tmax, &t3D, &ttime, &t_ext, &g_GABA_A, &g_GABA_A1, &g_GABA_B, &g_AMPA,
       &g_AMPA_CX_CX, &g_NMDA_CX_CX, &g_AMPA_CX_TC, &g_AMPA_CX_RE, 
       &g_AMPA_TC_CX, &g_AMPA_CX_IN, &g_NMDA_CX_IN, &g_GABA_A_IN_CX, 
       &g_GABA_B_IN_CX, 
       &g_AMPA_TC_IN, &g_Extern_ampa, &g_Extern_ampa1, &i_inp);
  printf("\n param: %lf %lf %lf %lf A=%lf A1=%lf B=%lf AM=%lf AM_CX_CX=%lf NMDA_CX_CX=%lf AM_CX_TC=%lf AM_CX_RE=%lf AM_TC_CX=%lf AM_CX_IN=%lf NMDA_CX_IN=%lf GA_IN_CX=%lf GB_IN_CX=%lf AM_TC_IN=%lf %lf %lf %d", 
       tmax, t3D, ttime, t_ext, g_GABA_A, g_GABA_A1, g_GABA_B, g_AMPA, 
       g_AMPA_CX_CX, g_NMDA_CX_CX, g_AMPA_CX_TC, g_AMPA_CX_RE, g_AMPA_TC_CX,
       g_AMPA_CX_IN, g_NMDA_CX_IN, g_GABA_A_IN_CX, g_GABA_B_IN_CX, 
       g_AMPA_TC_IN, g_Extern_ampa, g_Extern_ampa1, i_inp);
  
  fclose(f1);

  std::cout << "\n RE-RE_MAX=" << MS_RE_RE_MAX << " RE-TC_MAX=" << MS_RE_TC_MAX << 
                                          " TC-RE_MAX=" << MS_TC_RE_MAX;
printf("\n LLLLLLLLLLLLLLL, %lf",log(2.7));
//----------classes initialization (continue)----------------------------
  for(i=0; i < M; ++i)
  for(j=0; j < M1; ++j){
    g_a[i][j] = new Gaba_A[N_RE_RE];
    g_a1[i][j] = new Gaba_A[N_RE_TC];
    g_b[i][j] = new GB[N_RE_TC];
    a_a[i][j] = new AMPA[N_TC_RE]; 

    a_cx_tc[i][j] = new AMPA[N_CX_TC];
    a_cx_re[i][j] = new AMPA[N_CX_TC]; 
  }
  for(i=0; i < Mcx; ++i)
  for(j=0; j < Mcx1; ++j){ 
    a_cx_cx[i][j] = new AMPA_D2[N_CX_CX];
    nmda_cx_cx[i][j] = new NMDA_D1[N_CX_CX_NMDA];
    ga_in_cx[i][j] = new Gaba_A_D2[N_IN_CX];
    gb_in_cx[i][j] = new GB[N_IN_CX];

    a_tc_cx[i][j] = new AMPA[N_TC_CX];

    g_ampa_cx_cx[i][j] = new double[N_CX_CX];
    g_gaba_a_in_cx[i][j] = new double[N_IN_CX];
  }
  for(i=0; i < Min; ++i)
  for(j=0; j < Min1; ++j){ 
    nmda_cx_in[i][j] = new NMDA_D1[N_CX_IN_NMDA];
    a_cx_in[i][j] = new AMPA_D2[N_CX_IN];

    a_tc_in[i][j] = new AMPA[N_TC_IN];

    g_ampa_cx_in[i][j] = new double[N_CX_IN];
  }


   for(i=0; i < Min; ++i)
   for(j=0; j < Min1; ++j){
     in_cell[i][j].rho = 50;         
     in_cell[i][j].CX_DEND::G_Nap = 0.0;
     in_cell[i][j].CX_SOMA::G_Nap = 0.0;
     in_cell[i][j].CX_SOMA::G_Na = 2500;
   }

//----------creating the integer arrays containing the addresses--------
//----------of  ALL internal variables for ALL objects RE, TC ----------
//----------and GB classes (e.g., no_re[i][j][k] is the address --------
//----------of the variable y[k] for the object re_cell[i][j]) ---------
//----------NOTE: this is the relative addresses and you should ---------
//----------add the REAL address of the first element of the -----------
//----------original 1D array to use these arrays-----------------------
  for(i=0; i < M; ++i)
  for(j=0; j < M1; ++j)
  for(k=0; k < N_RE; ++k)
     no_re[i][j][k] = k + (j + i*M1) * N_RE;
  for(i=0; i < M; ++i)
  for(j=0; j < M1; ++j)
  for(k=0; k < N_TC; ++k)
     no_tc[i][j][k] = M*M1*N_RE*I_RE + k + (j + i*M1) * N_TC;  
  for(i=0; i < M; ++i)
  for(j=0; j < M1; ++j)
  for(k=0; k < N_RE_TC; ++k)
  for(l=0; l < N_GB; ++l)
     no_g[i][j][k][l] = M*M1*N_RE*I_RE + M*M1*N_TC*I_TC + l + 
                                   (k + (j + i*M1)*N_RE_TC) * N_GB;
  for(i=0; i < Mcx; ++i)
  for(j=0; j < Mcx1; ++j)
  for(k=0; k < N_CX; ++k)
     no_cx[i][j][k] = M*M1*N_RE*I_RE + M*M1*N_TC*I_TC + M*M1*N_RE_TC*N_GB*I_RE*I_TC +
                                    k + (j + i*Mcx1) * N_CX;
  for(i=0; i < Min; ++i)
  for(j=0; j < Min1; ++j)
  for(k=0; k < N_IN; ++k)
     no_in[i][j][k] = M*M1*N_RE*I_RE + M*M1*N_TC*I_TC + M*M1*N_RE_TC*N_GB*I_RE*I_TC +
                      Mcx*Mcx1*N_CX*I_CX + k + (j + i*Min1) * N_IN;
  for(i=0; i < Mcx; ++i)
  for(j=0; j < Mcx1; ++j)
  for(k=0; k < N_IN_CX; ++k)
  for(l=0; l < N_GB; ++l)
     no_gcx[i][j][k][l] = M*M1*N_RE*I_RE + M*M1*N_TC*I_TC + 
                          M*M1*N_RE_TC*N_GB*I_RE*I_TC +
                          Mcx*Mcx1*N_CX*I_CX + Min*Min1*N_IN*I_IN +
                          l + (k + (j + i*Mcx1)*N_IN_CX) * N_GB;
                                   
//---variable initialization (additional for standard constructor)----------
  for(i=0; i<M; i++)
    for(j=0; j<M1; j++){
      if(I_RE == 1) re_cell[i][j].init(y_ini+no_re[i][j][0]);
      if(I_TC == 1) tc_cell[i][j].init(y_ini+no_tc[i][j][0]);
      if(I_RE == 1 && I_TC == 1) for(k=0; k < N_RE_TC; k++){
                           g_b[i][j][k].init(y_ini+no_g[i][j][k][0]); }    
    }
  for(i=0; i<Mcx; i++)
    for(j=0; j<Mcx1; j++){
      if(I_CX == 1) cx_cell[i][j].init(y_ini+no_cx[i][j][0]);       
      if(I_CX == 1) for(k=0; k<N_CX_CX; k++) 
                       a_cx_cx[i][j][k].iii(1+i*(j+k));
      if(I_CX == 1 && I_IN == 1 && I_GB == 1) for(k=0; k < N_IN_CX; k++){
                       gb_in_cx[i][j][k].init(y_ini+no_gcx[i][j][k][0]); }
      if((I_CX == 1) && (I_IN == 1)) {
                 for(k=0; k<N_IN_CX; k++) ga_in_cx[i][j][k].iii(1+27+i*(j+k));
                 }
    }
  for(i=0; i<Min; i++)
    for(j=0; j<Min1; j++){
      if(I_IN == 1) in_cell[i][j].init(y_ini+no_in[i][j][0]);
      if((I_CX == 1) && (I_IN == 1)) {
                 for(k=0; k<N_CX_IN; k++) a_cx_in[i][j][k].iii(1+17+i*(j+k));
                 }
    }

//--------the scaling of the conductances-----------------------------------

if(I_RE == 1){
  g_gaba_a = 0; 
  g_ampa = 0; 
  g_ampa_cx_re = 0;
  g_ext_re = 0;
}
if(I_TC == 1){
  g_gaba_a1 = 0;
  g_gaba_b = 0;
  g_ampa_cx_tc = 0;
  g_ext_tc = 0;
}
if(I_CX == 1){
  for(i=0; i < Mcx; ++i)
   for(j=0; j < Mcx1; ++j)
    for(k=0; k < N_CX_CX; ++k)
      g_ampa_cx_cx[i][j][k] = 0;

  for(i=0; i < Mcx; ++i)
   for(j=0; j < Mcx1; ++j)
    for(k=0; k < N_IN_CX; ++k)
      g_gaba_a_in_cx[i][j][k] = 0;

  g_nmda_cx_cx = 0;
  g_ampa_tc_cx = 0;
  if(I_GB == 1) g_gaba_b_in_cx = 0;
  g_ext_cx = 0;
}

if(I_IN == 1){
  for(i=0; i < Min; ++i)
   for(j=0; j < Min1; ++j)
    for(k=0; k < N_CX_IN; ++k)
      g_ampa_cx_in[i][j][k] = 0;

  g_ampa_tc_in = 0;
  g_ext_in = 0;
}

//----------here we are changing some variables----------------------

  for(i=0; i < M; ++i)
  for(j=0; j < M1; ++j)
  for(k=0; k < N_RE_TC; ++k){  
     g_a[i][j][k].E_GABA = -70;  //GABA-A from RE to TC has more neg.revers. 
     g_a1[i][j][k].E_GABA = -83; //-85; //(J.Neur.1997,17(7),2348)
     g_b[i][j][k].K1 = 0.5; //0.09;
     g_b[i][j][k].K2 = 0.0012;
     g_b[i][j][k].K3 = 0.1; //0.18;
     g_b[i][j][k].K4 = 0.034;
     g_a[i][j][k].Alpha = 20;
     g_a[i][j][k].Beta = 0.162;
     g_a1[i][j][k].Alpha = 20;
     g_a1[i][j][k].Beta = 0.162;
  }
  for(i=0; i < M; ++i)
  for(j=0; j < M1; ++j){
     tc_cell[i][j].G_A = 0;
     tc_cell[i][j].ginc = 2;
     tc_cell[i][j].E_l = -70;
     tc_cell[i][j].G_Ca = 2.2; //2.7; //2.3; //2.;
     tc_cell[i][j].D = 2;
     tc_cell[i][j].pc = 0.007;
     tc_cell[i][j].k4 = 0.001;

     re_cell[i][j].E_l = -77; ///////////////////////////
     re_cell[i][j].G_Ca = 2.3; 
     tc_cell[i][j].INaK::Vtr = -40;
     tc_cell[i][j].INaK::VtrK = -25;
     tc_cell[i][j].G_K = 12;
     re_cell[i][j].G_kl = 0.005; //0.015; //0.005;
     tc_cell[i][j].G_h = 0.017;  
     tc_cell[i][j].G_kl = 0.03; //0.02; //0.0142; //0.025; //0.0142;
     tc_cell[i][j].DC = 0; //-0.05;
  }
  for(i=0; i < Mcx; ++i)
  for(j=0; j < Mcx1; ++j){
    g_AMPA_CX_CX_MINICE = 0.00006;
    g_AMPA_CX_IN_MINICE = 0.000025;
    cx_cell[i][j].G_kl = 0.0025;
    cx_cell[i][j].E_l = -68;
  }


//--------------open ALL files-------------------------------------
  f2 = fopen("dat", "w");
  f6 = fopen("graf_re", "w");
  f7 = fopen("time_re", "w");
  f8 = fopen("graf_tc", "w");
  f9 = fopen("time_tc", "w");
  f10 = fopen("graf_cx", "w");
  f11 = fopen("time_cx", "w");
  f12 = fopen("graf_in", "w");
  f13 = fopen("time_in", "w");
  f14 = fopen("ca_concentration_PY", "w");
  f15 = fopen("I_hva_PY","w");
  f16 = fopen("I_kca_PY","w");
  f17 = fopen("time_cx_dend", "w");
  f18 = fopen("I_hva_IN","w");
  f19 = fopen("I_kca_IN","w");
  f20 = fopen("I_na_PY_dend", "w");
  f21 = fopen("I_na_PY_soma", "w");
  f22 = fopen("I_na_IN_dend", "w");
  f23 = fopen("I_na_IN_soma", "w");
  f24 = fopen("I_nap_PY_dend", "w");
  f25 = fopen("I_nap_PY_soma", "w");
  f26 = fopen("I_nap_IN_dend", "w");
  f27 = fopen("I_nap_IN_soma", "w");
  f28 = fopen("I_k_PY", "w");
  f29 = fopen("I_k_IN", "w");
  f30 = fopen("I_kl_PY", "w");
  f31 = fopen("I_kl_IN", "w");
  f32 = fopen("I_km_PY", "w");
  f33 = fopen("I_km_IN", "w");
  f34 = fopen("I_l_PY", "w");
  f35 = fopen("I_l_IN", "w");
  f36 = fopen("Igap", "w");
  f37 = fopen("Igap_I", "w");
  f38 = fopen("ca_concentration_IN", "w");
  f39 = fopen("I_kl_TC", "w");
  f40 = fopen("I_kl_RE", "w");
  f41 = fopen("I_na_TC", "w");
  f42 = fopen("I_na_RE", "w");
  f43 = fopen("I_k_TC", "w");
  f44 = fopen("I_k_RE", "w");
  f45 = fopen("I_t_TC", "w");
  f46 = fopen("I_t_RE", "w");
  f47 = fopen("I_l_TC", "w");
  f48 = fopen("I_l_RE", "w");
  f49 = fopen("I_h_TC", "w");
  f50 = fopen("I_a_TC", "w");
  f51 = fopen("time_in_dend", "w");
  f52 = fopen("v_SOMA", "w");
  f53 = fopen("v_SOMA_IN", "w");
  f54 = fopen("m_na_PY", "w");
  f55 = fopen("h_na_PY", "w");
  f56 = fopen("m_na_IN", "w");
  f57 = fopen("h_na_IN", "w");
  f58 = fopen("v_RE", "w");
  f59 = fopen("v_TC", "w");
  f60 = fopen("a_cx_cx", "w");
  f61 = fopen("source", "w");
  f62 = fopen("target", "w");
  f63 = fopen("a_cx_in", "w");
  f64 = fopen("d_cx_in", "w");
  f65 = fopen("w_cx_in", "w");
  f66 = fopen("g1_cx_in", "w");
  f67 = fopen("g_Na_cx", "w");
  f68 = fopen("v_m_E_na_in", "w");
  f69 = fopen("E_na_in", "w");
  f70 = fopen("m_na_TC", "w");
  f71 = fopen("h_na_TC", "w");
  f72 = fopen("m_na_RE", "w");
  f73 = fopen("h_na_RE", "w");
  f74 = fopen("alpha_m_RE", "w");
  f75 = fopen("beta_m_RE", "w");
  f76 = fopen("alpha_h_RE", "w");
  f77 = fopen("beta_h_RE", "w");
  f78 = fopen("alpha_m_TC", "w");
  f79 = fopen("beta_m_TC", "w");
  f80 = fopen("alpha_h_TC", "w");
  f81 = fopen("beta_h_TC", "w");
  f82 = fopen("m_hva_cx", "w");
  f83 = fopen("h_hva_cx", "w");
  f84 = fopen("drive_cx", "w");
  f85 = fopen("a_cx_tc", "w");
  f86 = fopen("a_tc_cx", "w");
  f87 = fopen("a_nmda_cx_in", "w");
  f88 = fopen("a_gaba_in_cx", "w");
  f89 = fopen("a_gaba_a_re_tc", "w");
  f90 = fopen("a_gaba_b_re_tc", "w");
  f91 = fopen("g_rec_gabab", "w");
  f92 = fopen("r_rec_gabab", "w");
  f93 = fopen("tau_m_na_tc", "w");
  f94 = fopen("m_inf_na_tc", "w");
  f95 = fopen("tau_h_na_tc", "w");
  f96 = fopen("h_inf_na_tc", "w");
  f97 = fopen("fg_rec_gabab", "w");
  f98 = fopen("fr_rec_gabab", "w");

//----------Connection matrix-------------------------------

printf("\n Begin Connect Matrix");

for(i=0; i<M; i++)
  for(j=0; j<M1; j++){
     k_RERE[i][j]=0;
  for(i1=0; i1<M; i1++)
  for(j1=0; j1<M1; j1++) C_RERE[i][j][i1][j1]=0;
  }

for(i=0; i<M; i++)
  for(j=0; j<M1; j++){
     k_RETC[i][j]=0;
  for(i1=0; i1<M; i1++)
  for(j1=0; j1<M1; j1++) C_RETC[i][j][i1][j1]=0;
  }

for(i=0; i<M; i++)
  for(j=0; j<M1; j++){
     k_TCRE[i][j]=0;
  for(i1=0; i1<M; i1++)
  for(j1=0; j1<M1; j1++) C_TCRE[i][j][i1][j1]=0;
  }

for(i=0; i<Mcx; i++)
  for(j=0; j<Mcx1; j++){
     k_CXCX[i][j]=0;
  for(i1=0; i1<Mcx; i1++)
  for(j1=0; j1<Mcx1; j1++) C_CXCX[i][j][i1][j1]=0;
  }

for(i1=0; i1<Min; i1++)
  for(j1=0; j1<Min1; j1++){
     k_CXIN[i1][j1]=0;
  for(i=0; i<Mcx; i++)
  for(j=0; j<Mcx1; j++) C_CXIN[i][j][i1][j1]=0;
  }

for(i1=0; i1<Mcx; i1++)
  for(j1=0; j1<Mcx1; j1++){
     k_INCX[i1][j1]=0;
  for(i=0; i<Min; i++)
  for(j=0; j<Min1; j++) C_INCX[i][j][i1][j1]=0;
  }

for(i1=0; i1<M; i1++)
  for(j1=0; j1<M1; j1++){
     k_CXTC[i1][j1]=0;
  for(i=0; i<Mcx; i++)
  for(j=0; j<Mcx1; j++) C_CXTC[i][j][i1][j1]=0;
  }

for(i1=0; i1<M; i1++)
  for(j1=0; j1<M1; j1++){
     k_CXRE[i1][j1]=0;
  for(i=0; i<Mcx; i++)
  for(j=0; j<Mcx1; j++) C_CXRE[i][j][i1][j1]=0;
  }

for(i1=0; i1<Mcx; i1++)
  for(j1=0; j1<Mcx1; j1++){
     k_TCCX[i1][j1]=0;
  for(i=0; i<M; i++)
  for(j=0; j<M1; j++) C_TCCX[i][j][i1][j1]=0;
  }

for(i1=0; i1<Min; i1++)
  for(j1=0; j1<Min1; j1++){
     k_TCIN[i1][j1]=0;
  for(i=0; i<M; i++)
  for(j=0; j<M1; j++) C_TCIN[i][j][i1][j1]=0;
  }

printf("\n Connect Matrix1");

for(i1=0; i1<M; i1++)
for(j1=0; j1<M1; j1++)
for(i=0; i<M; i++)
for(j=0; j<M1; j++){
  scale = sqrt((double) ((i1-i)*(i1-i) + (j1-j)*(j1-j)));
  if(scale <= MS_RE_RE_MAX){
       C_RERE[i1][j1][i][j]=1;
       k_RERE[i][j]=k_RERE[i][j]+1;}
  if( (scale == 0) && (SELFING == 0)  ){
       C_RERE[i1][j1][i][j]=0;
       k_RERE[i][j]=k_RERE[i][j]-1;}       
  }

for(i1=0; i1<M; i1++)
for(j1=0; j1<M1; j1++)
for(i=0; i<M; i++)
for(j=0; j<M1; j++){
  scale = sqrt((double) ((i1-i)*(i1-i) + (j1-j)*(j1-j)));
  if(scale <= MS_TC_RE_MAX){
       C_TCRE[i1][j1][i][j]=1;
       k_TCRE[i][j]=k_TCRE[i][j]+1;}
  }

for(i1=0; i1<M; i1++)
for(j1=0; j1<M1; j1++)
for(i=0; i<M; i++)
for(j=0; j<M1; j++){
  scale = sqrt((double) ((i1-i)*(i1-i) + (j1-j)*(j1-j)));
  if(scale <= MS_RE_TC_MAX){
       C_RETC[i1][j1][i][j]=1;
       k_RETC[i][j]=k_RETC[i][j]+1;}
  }

for(i1=0; i1<Mcx; i1++)
for(j1=0; j1<Mcx1; j1++)
for(i=0; i<Mcx; i++)
for(j=0; j<Mcx1; j++){
  scale = sqrt((double) ((i1-i)*(i1-i) + (j1-j)*(j1-j)));
  if(scale <= MS_CX_CX_MAX){
       C_CXCX[i1][j1][i][j]=1;
       k_CXCX[i][j]=k_CXCX[i][j]+1;}
  if( (scale == 0) && (SELFINGcx == 0)  ){
       C_CXCX[i1][j1][i][j]=0;
       k_CXCX[i][j]=k_CXCX[i][j]-1;}
  }

for(i1=0; i1<Min; i1++)
for(j1=0; j1<Min1; j1++)
for(i=0; i<Mcx; i++)
for(j=0; j<Mcx1; j++){
  scale = sqrt((double) ((i1*Mcx/Min-i)*(i1*Mcx/Min-i) 
                       + (j1*Mcx1/Min1-j)*(j1*Mcx1/Min1-j)));
  if(scale <= MS_IN_CX_MAX){
       C_INCX[i1][j1][i][j]=1;
       k_INCX[i][j]=k_INCX[i][j]+1;}
  }

for(i1=0; i1<Mcx; i1++)
for(j1=0; j1<Mcx1; j1++)
for(i=0; i<Min; i++)
for(j=0; j<Min1; j++){
  scale = sqrt((double) ((i1*Min/Mcx-i)*(i1*Min/Mcx-i) 
                       + (j1*Min1/Mcx1-j)*(j1*Min1/Mcx1-j)));
  if(scale <= MS_CX_IN_MAX){
       C_CXIN[i1][j1][i][j]=1;
       k_CXIN[i][j]=k_CXIN[i][j]+1;//}
       printf("\n neuron CX %d and IN %d are connected", i1,i);}
  }

for(i1=0; i1<M; i1++)
for(j1=0; j1<M1; j1++)
for(i=0; i<Mcx; i++)
for(j=0; j<Mcx1; j++){
  scale = sqrt((double) ((i1*Mcx/M-i)*(i1*Mcx/M-i) 
                       + (j1*Mcx1/M1-j)*(j1*Mcx1/M1-j)));
  if(scale <= MS_TC_CX_MAX){
       C_TCCX[i1][j1][i][j]=1;
       k_TCCX[i][j]=k_TCCX[i][j]+1;}
  }

for(i1=0; i1<M; i1++)
for(j1=0; j1<M1; j1++)
for(i=0; i<Min; i++)
for(j=0; j<Min1; j++){
  scale = sqrt((double) ((i1*Min/M-i)*(i1*Min/M-i) 
                       + (j1*Min1/M1-j)*(j1*Min1/M1-j)));
  if(scale <= MS_TC_IN_MAX){
       C_TCIN[i1][j1][i][j]=1;
       k_TCIN[i][j]=k_TCIN[i][j]+1;}
}

for(i1=0; i1<Mcx; i1++)
for(j1=0; j1<Mcx1; j1++)
for(i=0; i<M; i++)
for(j=0; j<M1; j++){
  scale = sqrt((double) ((i1*M/Mcx-i)*(i1*M/Mcx-i) 
                       + (j1*M1/Mcx1-j)*(j1*M1/Mcx1-j)));
  if(scale <= MS_CX_TC_MAX){
       C_CXTC[i1][j1][i][j]=1;
       C_CXRE[i1][j1][i][j]=1;
       k_CXTC[i][j]=k_CXTC[i][j]+1;
       k_CXRE[i][j]=k_CXRE[i][j]+1;}
}
  
printf("\n Connect Matrix2");

k_REREmax=k_RERE[0][0];
for(i=0; i<M; i++){
for(j=0; j<M1; j++){
   if(k_RERE[i][j] > k_REREmax) k_REREmax=k_RERE[i][j];
   printf("k_re_re=%d , ", k_RERE[i][j]);}
}
printf("\n");

k_RETCmax=k_RETC[0][0];
for(i=0; i<M; i++){
for(j=0; j<M1; j++){
   if(k_RETC[i][j] > k_RETCmax) k_RETCmax=k_RETC[i][j];
   printf("k_re_tc=%d , ", k_RETC[i][j]);}
}
printf("\n");

k_TCREmax=k_TCRE[0][0];
for(i=0; i<M; i++){
for(j=0; j<M1; j++){
   if(k_TCRE[i][j] > k_TCREmax) k_TCREmax=k_TCRE[i][j];
   printf("k_tc_re=%d , ", k_TCRE[i][j]);}
}
printf("\n");

k_CXCXmax=k_CXCX[0][0];
for(i=0; i<Mcx; i++){
for(j=0; j<Mcx1; j++){
   if(k_CXCX[i][j] > k_CXCXmax) k_CXCXmax=k_CXCX[i][j];
   printf("k_cx_cx=%d , ", k_CXCX[i][j]);}
}
printf("\n");

k_CXINmax=k_CXIN[0][0];
for(i=0; i<Min; i++){
for(j=0; j<Min1; j++){
   if(k_CXIN[i][j] > k_CXINmax) k_CXINmax=k_CXIN[i][j];
   printf("k_cx_in=%d , ", k_CXIN[i][j]);}
}
printf("\n");
   

k_INCXmax=k_INCX[0][0];
for(i=0; i<Mcx; i++){
for(j=0; j<Mcx1; j++){
   if(k_INCX[i][j] > k_INCXmax) k_INCXmax=k_INCX[i][j];
   printf("k_in_cx=%d , ", k_INCX[i][j]);}
}
printf("\n");

k_CXTCmax=k_CXTC[0][0];
for(i=0; i<M; i++){
for(j=0; j<M1; j++){
   if(k_CXTC[i][j] > k_CXTCmax) k_CXTCmax=k_CXTC[i][j];
   printf("k_cx_tc=%d , ", k_CXTC[i][j]);}
}
printf("\n");

k_CXREmax=k_CXRE[0][0];
for(i=0; i<M; i++){
for(j=0; j<M1; j++){
   if(k_CXRE[i][j] > k_CXREmax) k_CXREmax=k_CXRE[i][j];
   printf("k_cx_re=%d , ", k_CXRE[i][j]);}
}
printf("\n");

k_TCCXmax=k_TCCX[0][0];
for(i=0; i<Mcx; i++){
for(j=0; j<Mcx1; j++){
   if(k_TCCX[i][j] > k_TCCXmax) k_TCCXmax=k_TCCX[i][j];
   printf("k_tc_cx=%d , ", k_TCCX[i][j]);}
}
printf("\n");

k_TCINmax=k_TCIN[0][0];
for(i=0; i<Min; i++){
for(j=0; j<Min1; j++){
   if(k_TCIN[i][j] > k_TCINmax) k_TCINmax=k_TCIN[i][j];
   printf("k_tc_in=%d , ", k_TCIN[i][j]);}
}
printf("\n");

printf("\n End Connect Matrix");

//----------------CALCULATION----------------------------------------
  printf("\n CALCULATION IN PROGRESS!!!: t= %lf: tmax= %lf", t,tmax);
  ih = (int) (1/h);
  TAU = h;

   while( t < tmax){ 
   tmax1 = t + TAU;
   ii = ii + 1;
   rk(N_EQ, fun, h, t, y_ini, f_ini, y1, y2);
   t = tmax1;

//--------the scaling of the conductances-----------------------------------
if( (t>50) && (t<(50+1.5*h)) ){
printf("\n CONDUCTANCE INITIALIZATION");

if(I_RE == 1){
  g_gaba_a = g_GABA_A / re_cell[0][0].S_RE; 
  g_ampa = g_AMPA / re_cell[0][0].S_RE; 
  g_ampa_cx_re = g_AMPA_CX_RE / re_cell[0][0].S_RE;
  g_ext_re = g_Extern_ampa/re_cell[0][0].S_RE;
}
if(I_TC == 1){
  g_gaba_a1 = g_GABA_A1 / tc_cell[0][0].S_TC;
  g_gaba_b = g_GABA_B / tc_cell[0][0].S_TC;
  g_ampa_cx_tc = g_AMPA_CX_TC / tc_cell[0][0].S_TC;
  g_ext_tc = g_Extern_ampa/tc_cell[0][0].S_TC;
}
if(I_CX == 1){
  g_ampa_cx_cx[0][0][0] = g_AMPA_CX_CX / cx_cell[0][0].S_CX_DEND;
  g_ampa_cx_cx_MINICE = g_AMPA_CX_CX_MINICE / cx_cell[0][0].S_CX_DEND;
  g_gaba_a_in_cx[0][0][0] = g_GABA_A_IN_CX / cx_cell[0][0].S_CX_DEND;
  g_nmda_cx_cx = g_NMDA_CX_CX / cx_cell[0][0].S_CX_DEND;
  g_ampa_tc_cx = g_AMPA_TC_CX / cx_cell[0][0].S_CX_DEND;
  if(I_GB == 1) g_gaba_b_in_cx = g_GABA_B_IN_CX / cx_cell[0][0].S_CX_DEND;
  g_ext_cx = g_Extern_ampa/cx_cell[0][0].S_CX_DEND/10; //2;
}

if(I_IN == 1){
  g_ampa_cx_in[0][0][0] = g_AMPA_CX_IN / in_cell[0][0].S_CX_DEND;  //S_IN;
  g_ampa_cx_in_MINICE = g_AMPA_CX_IN_MINICE / in_cell[0][0].S_CX_DEND;
  g_nmda_cx_in = g_NMDA_CX_IN / in_cell[0][0].S_CX_DEND;
  g_ampa_tc_in = g_AMPA_TC_IN / in_cell[0][0].S_CX_DEND;  //S_IN;
  g_ext_in = g_Extern_ampa/in_cell[0][0].S_CX_DEND/5; // S_IN /15; //3;
}

}
//-----------------------------------------------------------------------
   printf("\n param: g_ampa_cx_in=%.20lf", g_ampa_cx_in[0][0][0]);
   printf("\n param: g_ampa_cx_in_MINICE=%.20lf", g_ampa_cx_in_MINICE);
   printf("\n param: g_ampa_cx_cx=%.20lf", g_ampa_cx_cx[0][0][0]);
   printf("\n param: g_ampa_cx_cx_MINICE=%.20lf", g_ampa_cx_cx_MINICE);
   printf("\n param: g_gap_ex=%.20lf", 1.0/(10.0e3*cx_cell[0][0].S_CX_DEND));
   printf("\n param: g_gap_in=%.20lf", 1.0/(10.0e3*in_cell[0][0].S_CX_DEND));
   
   printf("\n param: g_gaba_a_re_re=%.20lf", g_gaba_a);
   printf("\n param : alpha gaba_a_re_re=%.2lf",g_a[25][0][5].Alpha);
   printf("\n param : beta gaba_a_re_re=%.2lf",g_a[25][0][5].Beta);
   printf("\n param : alpha gaba_a_re_tc=%.2lf",g_a1[25][0][5].Alpha);
   printf("\n param : beta gaba_a_re_tc=%.2lf",g_a1[25][0][5].Beta);
   printf("\n param: g_gaba_a_re_tc=%.20lf", g_gaba_a1);
   printf("\n param: g_gaba_b_re_tc=%.20lf", g_gaba_b);
   printf("\n");
   av=0;
   for(j = 0; j < Mcx1; ++j)
     for(i = 0; i < Mcx; ++i)
       av=av + cx_cell[i][j].v_SOMA;

   if((ii/(ih*2000))*(ih*2000) == ii) {

      printf("\n T= %lf ",t);
        for(j = 0; j < M1; ++j)
        for(i = 0; i < M; ++i){

        if(I_TC == 1){
          printf("\n TC: ");
          for(k = 0; k < N_TC; ++k) 
             printf("%lf ", y_ini[no_tc[i][j][k]]);}

        if(I_RE == 1){
          printf("\n RE: ");
          for(k = 0; k < N_RE; ++k)   
             printf("%lf ", y_ini[no_re[i][j][k]]);}
	}

        for(j = 0; j < Mcx1; ++j)
        for(i = 0; i < Mcx; ++i){

        if(I_CX == 1){
          printf("\n CX: ");
          for(k = 0; k < N_CX; ++k) 
             printf("%lf ", cx_cell[i][j].v_SOMA);}
	}
        for(j = 0; j < Min1; ++j)
        for(i = 0; i < Min; ++i){

        if(I_IN == 1){
          printf("\n IN: ");
          for(k = 0; k < N_IN; ++k)   
             printf("%lf ", y_ini[no_in[i][j][k]]);}
	}
   }



   if((t > ttime) && ((ii/(1))*(1) == ii)){
     if(I_RE == 1) fprintf(f7,"%lf ", t);
     if(I_TC == 1) fprintf(f9,"%lf ", t);
     if(I_CX == 1) {fprintf(f11,"%lf ", t);}
			    //fprintf(f14,"%lf ",t);}
     if(I_IN == 1) fprintf(f13,"%lf ", t);
     for (k=0;k<k_CXIN[12][0];k++){
     	if(I_IN == 1) fprintf(f63,"%lf ", a_cx_in[12][0][k].I);
     	fprintf(f87,"%lf ", nmda_cx_in[12][0][k].I);
     	}
     for (k=0;k<k_CXCX[50][0];k++){
     	if(I_CX == 1) fprintf(f60,"%lf ", a_cx_cx[50][0][k].I);
     	}
     for (k=0;k<k_CXTC[25][0];k++){
     	fprintf(f85,"%lf ", a_cx_tc[25][0][k].I);
     	}
     for (k=0;k<k_TCCX[50][0];k++){
     	fprintf(f86,"%lf ", a_tc_cx[50][0][k].I);
     	}
     for (k=0;k<k_RETC[25][0];k++){
     	fprintf(f89,"%lf ", g_a1[25][0][k].I);
     	fprintf(f90,"%lf ", g_b[25][0][k].I);
     	fprintf(f91,"%lf ", g_b[25][0][k].g_rec);
     	fprintf(f92,"%lf ", g_b[25][0][k].r_rec);
     	fprintf(f97,"%lf ", g_b[25][0][k].fr_rec);
     	fprintf(f98,"%lf ", g_b[25][0][k].fg_rec);
     	}
     for (k=0;k<k_INCX[50][0];k++){
     	fprintf(f88,"%lf ", ga_in_cx[50][0][k].I);
     	}
     for (k=0;k<k_CXIN[12][0];k++){
     	if(I_IN == 1) fprintf(f64,"%lf ", a_cx_in[12][0][k].E);
     	if(I_IN == 1) fprintf(f65,"%lf ", a_cx_in[12][0][k].R);
     	if(I_IN == 1) fprintf(f66,"%lf ", a_cx_in[12][0][k].g1);
     	}
     for(i = 0; i < M; ++i){
         if(I_RE == 1) {fprintf(f7,"%lf ", y_ini[no_re[i][M1/2][0]]);
			fprintf(f40,"%lf ", re_cell[i][M1/2].iKl);
			fprintf(f42,"%lf ", re_cell[i][M1/2].iNa);
			fprintf(f44,"%lf ", re_cell[i][M1/2].iK);
			fprintf(f46,"%lf ", re_cell[i][M1/2].iT);
			fprintf(f48,"%lf ", re_cell[i][M1/2].iL);
			fprintf(f58,"%lf ", re_cell[i][M1/2].v_RE); 
			fprintf(f72,"%lf ", re_cell[i][M1/2].mNa_TCRE);
			fprintf(f73,"%lf ", re_cell[i][M1/2].hNa_TCRE);
			fprintf(f74,"%lf ", re_cell[i][M1/2].Alpha_m);
			fprintf(f75,"%lf ", re_cell[i][M1/2].Beta_m);
			fprintf(f76,"%lf ", re_cell[i][M1/2].Alpha_h);
			fprintf(f77,"%lf ", re_cell[i][M1/2].Beta_h);
			}
         if(I_TC == 1) {fprintf(f9,"%lf ", y_ini[no_tc[i][M1/2][0]]);
			fprintf(f39,"%lf ", tc_cell[i][M1/2].iKl);
			fprintf(f41,"%lf ", tc_cell[i][M1/2].iNa);
			fprintf(f43,"%lf ", tc_cell[i][M1/2].iK);
			fprintf(f45,"%lf ", tc_cell[i][M1/2].iT);
			fprintf(f47,"%lf ", tc_cell[i][M1/2].iL);
			fprintf(f49,"%lf ", tc_cell[i][M1/2].ih);
			fprintf(f50,"%lf ", tc_cell[i][M1/2].iA);
			fprintf(f59,"%lf ", tc_cell[i][M1/2].v_TC);
			fprintf(f70,"%lf ", tc_cell[i][M1/2].mNa_TCRE);
			fprintf(f71,"%lf ", tc_cell[i][M1/2].hNa_TCRE);
			fprintf(f78,"%lf ", tc_cell[i][M1/2].Alpha_m);
			fprintf(f79,"%lf ", tc_cell[i][M1/2].Beta_m);
			fprintf(f80,"%lf ", tc_cell[i][M1/2].Alpha_h);
			fprintf(f81,"%lf ", tc_cell[i][M1/2].Beta_h);
			fprintf(f93,"%lf ", tc_cell[i][M1/2].tau_m_na);
			fprintf(f94,"%lf ", tc_cell[i][M1/2].m_inf_na);
			fprintf(f95,"%lf ", tc_cell[i][M1/2].tau_h_na);
			fprintf(f96,"%lf ", tc_cell[i][M1/2].h_inf_na);
			}
     }
     for(i = 0; i < Mcx; ++i)
         if(I_CX == 1) {fprintf(f52,"%lf ", cx_cell[i][Mcx1/2].v_SOMA);    
			fprintf(f14,"%lf ", cx_cell[i][Mcx1/2].iCa);
			fprintf(f15,"%lf ", cx_cell[i][Mcx1/2].iHVA);
		        fprintf(f16,"%lf ", cx_cell[i][Mcx1/2].iKCa);
			fprintf(f17,"%lf ", cx_cell[i][Mcx1/2].v_DEND);
			fprintf(f20,"%lf ", cx_cell[i][Mcx1/2].CX_DEND::iNa);
			fprintf(f21,"%lf ", cx_cell[i][Mcx1/2].CX_SOMA::iNa);
			fprintf(f24,"%lf ", cx_cell[i][Mcx1/2].CX_DEND::iNap);
			fprintf(f25,"%lf ", cx_cell[i][Mcx1/2].CX_SOMA::iNap);
			fprintf(f28,"%lf ", cx_cell[i][Mcx1/2].iKv);
			fprintf(f30,"%lf ", cx_cell[i][Mcx1/2].iKl);
			fprintf(f32,"%lf ", cx_cell[i][Mcx1/2].iKm);
			fprintf(f34,"%lf ", cx_cell[i][Mcx1/2].iL);
			fprintf(f36,"%lf ", cx_cell[i][Mcx1/2].iGAP);
			fprintf(f54,"%lf ", cx_cell[i][Mcx1/2].CX_SOMA::m_Na);
			fprintf(f55,"%lf ", cx_cell[i][Mcx1/2].CX_SOMA::h_Na);
			fprintf(f82,"%lf ", cx_cell[i][Mcx1/2].m_hva);
			fprintf(f83,"%lf ", cx_cell[i][Mcx1/2].h_hva);
			fprintf(f84,"%lf ", cx_cell[i][Mcx1/2].drive_save);
			}
     for(i = 0; i < Min; ++i)
         if(I_IN == 1) {fprintf(f53,"%lf ", in_cell[i][Min1/2].v_SOMA);
			fprintf(f51,"%lf ", in_cell[i][Min1/2].v_DEND);
			fprintf(f18,"%lf ", in_cell[i][Min1/2].iHVA);
			fprintf(f19,"%lf ", in_cell[i][Min1/2].iKCa);
			fprintf(f22,"%lf ", in_cell[i][Min1/2].CX_DEND::iNa);
			fprintf(f23,"%lf ", in_cell[i][Min1/2].CX_SOMA::iNa);
			fprintf(f67,"%lf ", in_cell[i][Min1/2].CX_SOMA::g_Na);
			fprintf(f68,"%lf ", in_cell[i][Min1/2].CX_SOMA::v_m_E_Na);
			fprintf(f69,"%lf ", in_cell[i][Min1/2].CX_SOMA::E_Na);
			fprintf(f26,"%lf ", in_cell[i][Min1/2].CX_DEND::iNap);
			fprintf(f27,"%lf ", in_cell[i][Min1/2].CX_SOMA::iNap);
			fprintf(f29,"%lf ", in_cell[i][Min1/2].iKv);
			fprintf(f31,"%lf ", in_cell[i][Min1/2].iKl);
			fprintf(f33,"%lf ", in_cell[i][Min1/2].iKm);
			fprintf(f35,"%lf ", in_cell[i][Min1/2].iL);
			fprintf(f37,"%lf ", in_cell[i][Min1/2].iGAP);
			fprintf(f38,"%lf ", in_cell[i][Min1/2].iCa);
			fprintf(f56,"%lf ", in_cell[i][Min1/2].CX_SOMA::m_Na);
			fprintf(f57,"%lf ", in_cell[i][Min1/2].CX_SOMA::h_Na);
			}

	 if(I_CX == 1) {
			fprintf(f62,"%lf ", cx_cell[50][Mcx1/2].v_DEND);
			fprintf(f61,"%lf ", cx_cell[45][Mcx1/2].v_SOMA);
			}

			

     fprintf(f7,"\n");
     fprintf(f9,"\n");
     fprintf(f11,"\n");
     fprintf(f13,"\n");
     fprintf(f14,"\n");
     fprintf(f15,"\n");
     fprintf(f16,"\n");
     fprintf(f17,"\n");
     fprintf(f18,"\n");
     fprintf(f19,"\n");
     fprintf(f20,"\n");
     fprintf(f21,"\n");
     fprintf(f22,"\n");
     fprintf(f23,"\n");
     fprintf(f24,"\n");
     fprintf(f25,"\n");
     fprintf(f26,"\n");
     fprintf(f27,"\n");
     fprintf(f28,"\n");
     fprintf(f29,"\n");
     fprintf(f30,"\n");
     fprintf(f31,"\n");
     fprintf(f32,"\n");
     fprintf(f33,"\n");
     fprintf(f34,"\n");
     fprintf(f35,"\n");
     fprintf(f36,"\n");
     fprintf(f37,"\n");
     fprintf(f38,"\n");
     fprintf(f39,"\n");
     fprintf(f40,"\n");
     fprintf(f41,"\n");
     fprintf(f42,"\n");
     fprintf(f43,"\n");
     fprintf(f44,"\n");
     fprintf(f45,"\n");
     fprintf(f46,"\n");
     fprintf(f47,"\n");
     fprintf(f48,"\n");
     fprintf(f49,"\n");
     fprintf(f50,"\n");
     fprintf(f51,"\n");
     fprintf(f52,"\n");
     fprintf(f53,"\n");
     fprintf(f54,"\n");
     fprintf(f55,"\n");
     fprintf(f56,"\n");
     fprintf(f57,"\n");
     fprintf(f58,"\n");
     fprintf(f59,"\n");
     fprintf(f60,"\n");
     fprintf(f61,"\n");
     fprintf(f62,"\n");
     fprintf(f63,"\n");
     fprintf(f64,"\n");
     fprintf(f65,"\n");
     fprintf(f66,"\n");
     fprintf(f67,"\n");
     fprintf(f68,"\n");
     fprintf(f69,"\n");
     fprintf(f70,"\n");
     fprintf(f71,"\n");
     fprintf(f72,"\n");
     fprintf(f73,"\n");
     fprintf(f74,"\n");
     fprintf(f75,"\n");
     fprintf(f76,"\n");
     fprintf(f77,"\n");
     fprintf(f78,"\n");
     fprintf(f79,"\n");
     fprintf(f80,"\n");
     fprintf(f81,"\n");
     fprintf(f82,"\n");
     fprintf(f83,"\n");
     fprintf(f84,"\n");
     fprintf(f85,"\n");
     fprintf(f86,"\n");
     fprintf(f87,"\n");
     fprintf(f88,"\n");
     fprintf(f89,"\n");
     fprintf(f90,"\n");
     fprintf(f91,"\n");
     fprintf(f92,"\n");
     fprintf(f93,"\n");
     fprintf(f94,"\n");
     fprintf(f95,"\n");
     fprintf(f96,"\n");
     fprintf(f97,"\n");
     fprintf(f98,"\n");

 }

}
//--------------------END CALCULATION-------------------------------

//-----------------close ALL files-----------------------------------
  fclose(f2);
  fclose(f6);
  fclose(f7);
  fclose(f8);
  fclose(f9);
  fclose(f10);
  fclose(f11);
  fclose(f12);
  fclose(f13);
  fclose(f14);
  fclose(f15);
  fclose(f16);
  fclose(f17);
  fclose(f18);
  fclose(f19);
  fclose(f20);
  fclose(f21);
  fclose(f22);
  fclose(f23);
  fclose(f24);
  fclose(f25);
  fclose(f26);
  fclose(f27);
  fclose(f28);
  fclose(f29);
  fclose(f30);
  fclose(f31);
  fclose(f32);
  fclose(f33);
  fclose(f34);
  fclose(f35);
  fclose(f36);
  fclose(f37);
  fclose(f38);
  fclose(f39);
  fclose(f40);
  fclose(f41);
  fclose(f42);
  fclose(f43);
  fclose(f44);
  fclose(f45);
  fclose(f46);
  fclose(f47);
  fclose(f48);
  fclose(f49);
  fclose(f50);
  fclose(f51);
  fclose(f52);
  fclose(f53);
  fclose(f54);
  fclose(f55);
  fclose(f56);
  fclose(f57);
  fclose(f58);
  fclose(f59);
  fclose(f60);
  fclose(f61);
  fclose(f62);
  fclose(f63);
  fclose(f64);
  fclose(f65);
  fclose(f66);
  fclose(f67);
  fclose(f68);
  fclose(f69);
  fclose(f70);
  fclose(f71);
  fclose(f72);
  fclose(f73);
  fclose(f74);
  fclose(f75);
  fclose(f76);
  fclose(f77);
  fclose(f78);
  fclose(f79);
  fclose(f80);
  fclose(f81);
  fclose(f82);
  fclose(f83);
  fclose(f84);
  fclose(f85);
  fclose(f86);
  fclose(f87);
  fclose(f88);
  fclose(f89);
  fclose(f90);
  fclose(f91);
  fclose(f92);
  fclose(f93);
  fclose(f94);
  fclose(f95);
  fclose(f96);
  fclose(f97);
  fclose(f98);


  printf("\n"); 
}

//+++++++++++ Function to calculate the index of BOUNDARY elements +++++++++++
int b(unsigned btype, unsigned m, int ind) {
      if(btype == 0) {  //------flow boundary conditions:
        if( (ind >= 0) && (ind <= (m-1)) ) return( ind );
        if(ind < 0) return( -ind - 1);
        if(ind > (m-1)) return( 2*m - ind - 1); }

      else if(btype == 1) {  //------periodic boundary conditions: 
        if( (ind >= 0) && (ind <= (m-1)) ) return( ind );
        if(ind < 0) return( ind + m);
        if(ind > (m-1)) return (ind - m); }

      else if(btype == 9) {  //------NO boundary conditions:
        return( ind ); }     

      else if(btype == 7) {  //------NO boundary conditions:
        return( ind ); }  
      
      return( ind ); //-----ADDED BY MATHILDE
}



//+++++++++++ Function to calculate the right sides for ALL ODE +++++++++++++++
void fun(unsigned neq, double x, double *y_ini, double *f_ini){
double scale;
int i, j, k, kmax, k1, i1, j1, ii, jj, nst, kk[Max][Max1], kk1[Max][Max1];


//========here the MAIN loop to calculate intrinsic conductances===========
//--------(f_ini IS changed, y_ini IS NOT changed)-------------------------
for(i=0; i < M; ++i)
for(j=0; j < M1; ++j){
  if(I_RE == 1) re_cell[i][j].calc(x, y_ini+no_re[i][j][0], f_ini+no_re[i][j][0]);
  if(I_TC == 1) tc_cell[i][j].calc(x, y_ini+no_tc[i][j][0], f_ini+no_tc[i][j][0]); 
}

for(i=0; i < Mcx; ++i)
for(j=0; j < Mcx1; ++j)
  if(I_CX == 1) cx_cell[i][j].calc(x, y_ini+no_cx[i][j][0], f_ini+no_cx[i][j][0]);

for(i=0; i < Min; ++i)
for(j=0; j < Min1; ++j)
  if(I_IN == 1) in_cell[i][j].calc(x, y_ini+no_in[i][j][0], f_ini+no_in[i][j][0]);

//========here the MAIN loop to calculate synaptic conductances=============
//--------(f_ini IS changed, y_ini IS NOT changed) -------------------------

for(i = 0; i < M; ++i)
for(j = 0; j < M1; ++j){

//--------reciprocal GABA-A between RE cells---------------------------------
if(I_RE == 1){

  for(i1 = 0, k = 0; i1 < M; ++i1)
  for(j1 = 0; j1 < M1; ++j1){
    if(C_RERE[i1][j1][i][j] > 0){
              g_a[i][j][k].calc(g_gaba_a, x, y_ini[no_re[i][j][0]], 
                                             y_ini[no_re[i1][j1][0]]);
              ++k; 
  } }
  kmax = k;
  for(k = 0; k < kmax; ++k){
    if(kmax == 0) { }  //NO synapses
    else { 
       g_a[i][j][k].I = g_a[i][j][k].I / k_RERE[i][j]; 
       f_ini[no_re[i][j][0]] = f_ini[no_re[i][j][0]] - g_a[i][j][k].I; }
  }}
 
//---------GABA-A and GABA-B from RE to TC cells------------------------------
if(I_RE == 1 && I_TC == 1){

  for(i1 = 0, k = 0; i1 < M; ++i1)
  for(j1 = 0; j1 < M1; ++j1){
    if(C_RETC[i1][j1][i][j] > 0){

         g_a1[i][j][k].calc(g_gaba_a1, x, y_ini[no_tc[i][j][0]], 
                             y_ini[no_re[i1][j1][0]]);

         g_b[i][j][k].calc(g_gaba_b, x, y_ini+no_g[i][j][k][0],
                             f_ini+no_g[i][j][k][0], y_ini[no_tc[i][j][0]], 
                             y_ini[no_re[i1][j1][0]]);
	 ++k; 
  } }
  kmax = k;
  for(k = 0; k < kmax; ++k){
    if(kmax == 0) { }  //NO synapses
    else {
       g_a1[i][j][k].I = g_a1[i][j][k].I / k_RETC[i][j]; 
       g_b[i][j][k].I = g_b[i][j][k].I / k_RETC[i][j]; 
       f_ini[no_tc[i][j][0]] = f_ini[no_tc[i][j][0]] - g_a1[i][j][k].I;
       f_ini[no_tc[i][j][0]] = f_ini[no_tc[i][j][0]] - g_b[i][j][k].I; }
  }}

//-----------AMPA from TC to RE cells--------------------------------------
if(I_RE == 1 && I_TC == 1){

  for(i1 = 0, k = 0; i1 < M; ++i1)
  for(j1 = 0; j1 < M1; ++j1){
    if(C_TCRE[i1][j1][i][j] > 0){
        a_a[i][j][k].calc(g_ampa, x, y_ini[no_re[i][j][0]], 
                           y_ini[no_tc[i1][j1][0]]);
    ++k; 
  } }
  kmax = k;
  for(k = 0; k < kmax; ++k){
    if(kmax == 0) { }  //NO synapses
    else {
       a_a[i][j][k].I = a_a[i][j][k].I / k_TCRE[i][j];
       f_ini[no_re[i][j][0]] = f_ini[no_re[i][j][0]] - a_a[i][j][k].I; }
}}

//------------AMPA from CX to TC cells-----------------------------------
if(I_CX == 1 && I_TC == 1){

  for(i1 = 0, k = 0; i1 < Mcx; ++i1)
  for(j1 = 0; j1 < Mcx1; ++j1){
    if(C_CXTC[i1][j1][i][j] > 0){
        a_cx_tc[i][j][k].calc(g_ampa_cx_tc, x, y_ini[no_tc[i][j][0]], 
                           cx_cell[i1][j1].v_SOMA);
        ++k; 
  } }
  kmax = k;
  for(k = 0; k < kmax; ++k){
    if(kmax == 0) { }  //NO synapses
    else {
       a_cx_tc[i][j][k].I = a_cx_tc[i][j][k].I / k_CXTC[i][j]; 
       f_ini[no_tc[i][j][0]] = f_ini[no_tc[i][j][0]] - a_cx_tc[i][j][k].I; }
}}

//-------------AMPA from CX to RE cells---------------------------------------
if(I_CX == 1 && I_RE == 1){

  for(i1 = 0, k = 0; i1 < Mcx; ++i1)
  for(j1 = 0; j1 < Mcx1; ++j1){
    if(C_CXRE[i1][j1][i][j] > 0){
        a_cx_re[i][j][k].calc(g_ampa_cx_re, x, y_ini[no_re[i][j][0]], 
                           cx_cell[i1][j1].v_SOMA);
        ++k; //} // number of synapses of the cell i-j 
  } }
  kmax = k;
  for(k = 0; k < kmax; ++k){
    if(kmax == 0) { }  //NO synapses
    else {
       a_cx_re[i][j][k].I = a_cx_re[i][j][k].I / k_CXRE[i][j]; 
       f_ini[no_re[i][j][0]] = f_ini[no_re[i][j][0]] - a_cx_re[i][j][k].I; }
}}
}

for(i = 0; i < Mcx; ++i)
for(j = 0; j < Mcx1; ++j){

//-----------reciprocal AMPA between CX cells--------------------------------
if(I_CX == 1){

  for(i1 = -MS_CX_CX, k = 0; i1 < Mcx+MS_CX_CX; ++i1)
  for(j1 = -MS_CX_CX1; j1 < Mcx1+MS_CX_CX1; ++j1){

    if( (i1<0) || (i1>(Mcx-1)) || (j1<0) || (j1>(Mcx1-1)) ){
//       scale = sqrt((double) ((i1-i)*(i1-i) + (j1-j)*(j1-j)));
       scale = abs(i1-i);
       if(scale <= MS_CX_CX_MAX){
          ii=b(BOUNDcx,Mcx,i1); jj=b(BOUNDcx,Mcx1,j1);
          a_cx_cx[i][j][k].calc(g_ampa_cx_cx[0][0][0], 0.0, 
                         x, y_ini[no_cx[i][j][0]], 
                         cx_cell[ii][jj].v_SOMA);
          ++k;}}
    
    else if(C_CXCX[i1][j1][i][j] > 0){
        a_cx_cx[i][j][k].calc(g_ampa_cx_cx[0][0][0], g_ampa_cx_cx_MINICE, 
                         x, y_ini[no_cx[i][j][0]], 
                         cx_cell[i1][j1].v_SOMA);
        ++k; }
  }
  kmax = k;
  for(k = 0; k < kmax; ++k){
    if(kmax == 0) { }  //NO synapses
    else {
       a_cx_cx[i][j][k].I = a_cx_cx[i][j][k].I / kmax; //k_CXCX[i][j]; 
       f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] - a_cx_cx[i][j][k].I;}
}}

//-----------reciprocal NMDA between CX cells--------------------------------
if(I_CX == 1){

  for(i1 = -MS_CX_CX_NMDA, k = 0; i1 < Mcx+MS_CX_CX_NMDA; ++i1)
  for(j1 = -MS_CX_CX_NMDA1; j1 < Mcx1+MS_CX_CX_NMDA1; ++j1){

    if( (i1<0) || (i1>(Mcx-1)) || (j1<0) || (j1>(Mcx1-1)) ){
//       scale = sqrt((double) ((i1-i)*(i1-i) + (j1-j)*(j1-j)));
       scale = abs(i1-i);
       if(scale <= MS_CX_CX_NMDA_MAX){
          ii=b(BOUNDcx,Mcx,i1); jj=b(BOUNDcx,Mcx1,j1);
          nmda_cx_cx[i][j][k].calc(g_nmda_cx_cx, x, y_ini[no_cx[i][j][0]], 
                        cx_cell[ii][jj].v_SOMA);
          ++k; }}

    else if(C_CXCX[i1][j1][i][j] > 0){
              nmda_cx_cx[i][j][k].calc(g_nmda_cx_cx, x, y_ini[no_cx[i][j][0]], 
                        cx_cell[i1][j1].v_SOMA);
        ++k; }
  }
  kmax = k;
  for(k = 0; k < kmax; ++k){
    if(kmax == 0) { }  //NO synapses
    else {
       nmda_cx_cx[i][j][k].I = nmda_cx_cx[i][j][k].I / kmax; //k_CXCX[i][j];
       f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] - nmda_cx_cx[i][j][k].I;}
}}

//--------------GABA-A from IN to CX cells------------------------------------
if(I_CX == 1 && I_IN == 1){

  for(i1 = 0, k = 0; i1 < Min; ++i1)
  for(j1 = 0; j1 < Min1; ++j1){
    if(C_INCX[i1][j1][i][j] > 0){
       ga_in_cx[i][j][k].calc(g_gaba_a_in_cx[0][0][0],x,y_ini[no_cx[i][j][0]],
                        in_cell[i1][j1].v_SOMA);
        ++k;} 
    } 
  kmax = k;
  for(k = 0; k < kmax; ++k){
    if(kmax == 0) { }  //NO synapses
    else {
       ga_in_cx[i][j][k].I = ga_in_cx[i][j][k].I / k_INCX[i][j];
       f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] - ga_in_cx[i][j][k].I; }
}}

//--------------GABA-B from IN to CX cells------------------------------------
if(I_CX == 1 && I_IN == 1 && I_GB == 1){

  for(i1 = 0, k = 0; i1 < Min; ++i1)
  for(j1 = 0; j1 < Min1; ++j1){
    if(C_INCX[i1][j1][i][j] > 0){
        gb_in_cx[i][j][k].calc(g_gaba_b_in_cx, x, y_ini+no_gcx[i][j][k][0],
                          f_ini+no_gcx[i][j][k][0], y_ini[no_cx[i][j][0]], 
                        in_cell[i1][j1].v_SOMA);
        ++k;}                                       
  } 
  kmax = k;
  for(k = 0; k < kmax; ++k){
    if(kmax == 0) { }  //NO synapses
    else {
       gb_in_cx[i][j][k].I = gb_in_cx[i][j][k].I / k_INCX[i][j];
       f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] - gb_in_cx[i][j][k].I; }
}}

//--------------AMPA from TC to CX cells--------------------------------------
if(I_CX == 1 && I_TC == 1){
  
  for(i1 = 0, k = 0; i1 < M; ++i1)
  for(j1 = 0; j1 < M1; ++j1){
    if(C_TCCX[i1][j1][i][j] > 0){
        a_tc_cx[i][j][k].calc(g_ampa_tc_cx, x, y_ini[no_cx[i][j][0]], 
                           y_ini[no_tc[i1][j1][0]]);
        ++k; 
  } }
  kmax = k;
  for(k = 0; k < kmax; ++k){
    if(kmax == 0) { }  //NO synapses
    else {
       a_tc_cx[i][j][k].I = a_tc_cx[i][j][k].I / k_TCCX[i][j]; 
       f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] - a_tc_cx[i][j][k].I; }
}}
}

for(i = 0; i < Min; ++i)
for(j = 0; j < Min1; ++j){

//-------------AMPA from CX to IN cells--------------------------------------
if(I_CX == 1 && I_IN == 1){
 
  for(i1 = 0, k = 0; i1 < Mcx; ++i1)
  for(j1 = 0; j1 < Mcx1; ++j1){
    if(C_CXIN[i1][j1][i][j] > 0){
        a_cx_in[i][j][k].calc(g_ampa_cx_in[0][0][0], g_ampa_cx_in_MINICE, 
                         x, y_ini[no_in[i][j][0]],
                         cx_cell[i1][j1].v_SOMA);
        ++k; } 
  }
  kmax = k;
  for(k = 0; k < kmax; ++k){
    if(kmax == 0) { }  //NO synapses
    else {
       a_cx_in[i][j][k].I = a_cx_in[i][j][k].I / k_CXIN[i][j];
       f_ini[no_in[i][j][0]] = f_ini[no_in[i][j][0]] - a_cx_in[i][j][k].I; }
}}

//-----------NMDA from CX to IN--------------------------------
if(I_CX == 1){

  for(i1 = 0, k = 0; i1 < Mcx; ++i1)
  for(j1 = 0; j1 < Mcx1; ++j1){
    if(C_CXIN[i1][j1][i][j] > 0){
              nmda_cx_in[i][j][k].calc(g_nmda_cx_in, x, y_ini[no_in[i][j][0]], 
                        cx_cell[i1][j1].v_SOMA);
        ++k; }}
  kmax = k;
  for(k = 0; k < kmax; ++k){
    if(kmax == 0) { }  //NO synapses
    else {
       nmda_cx_in[i][j][k].I = nmda_cx_in[i][j][k].I / k_CXIN[i][j];
       f_ini[no_in[i][j][0]] = f_ini[no_in[i][j][0]] - nmda_cx_in[i][j][k].I;}
}}

//-------------AMPA from TC to IN cells--------------------------------------
if(I_IN == 1 && I_TC == 1){

  for(i1 = 0, k = 0; i1 < M; ++i1)
  for(j1 = 0; j1 < M1; ++j1){
    if(C_TCIN[i1][j1][i][j] > 0){
       a_tc_in[i][j][k].calc(g_ampa_tc_in, x, y_ini[no_in[i][j][0]], 
                          y_ini[no_tc[i1][j1][0]]);
       ++k;
  } }
  kmax = k;
  for(k = 0; k < kmax; ++k){
    if(kmax == 0) { }  //NO synapses
    else {
       a_tc_in[i][j][k].I = a_tc_in[i][j][k].I / k_TCIN[i][j];
       f_ini[no_in[i][j][0]] = f_ini[no_in[i][j][0]] - a_tc_in[i][j][k].I; }
}}

}

//=============END of MAIN loop==============================================




//------------------expernal stimulation---------------------------------
nst = 1;   //type of external stimulation

//-----------train of shocks----------------------------------------
if(nst == -1){ 
  if( (x > 2000) && (x < 3100) ) {
    a_ext1.calc(g_ext_tc, x);
    a_ext2.calc(g_ext_re, x);
    a_ext3.calc(g_ext_cx, x);
    a_ext4.calc(g_ext_in, x);
    for(i = 0; i < Min; ++i)
    for(j = 0; j < Min1; ++j){
      if(I_IN == 1) f_ini[no_in[i][j][0]] = f_ini[no_in[i][j][0]] -
           a_ext4.g * exp(-0.1*sqrt( (double)((i-Min/2)*(i-Min/2)+(j-Min1/2)*(j-Min1/2)) ) ) * 
                                                       y_ini[no_in[i][j][0]];
    }
    for(i = 0; i < Mcx; ++i)
    for(j = 0; j < Mcx1; ++j){
      if(I_CX == 1) f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] -
           a_ext3.g * exp(-0.1*sqrt( (double)((i-Mcx/2)*(i-Mcx/2)+(j-Mcx1/2)*(j-Mcx1/2)) ) ) * 
                                                       y_ini[no_cx[i][j][0]];
   }}
}
else if(nst == 0){ 
  if( (x > 500) && (x < 5100) ) {
    a_ext1.calc(g_ext_tc, x);
    a_ext2.calc(g_ext_re, x);
    a_ext3.calc(g_ext_cx, x);
    a_ext4.calc(g_ext_in, x);

    for(i = Mcx/2-1; i <= Mcx/2+1; ++i)
    for(j = 0; j < Mcx1; ++j){
      if(I_CX == 1) f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] -
           a_ext3.g * exp(-0.5*sqrt( (double)((i-Mcx/2)*(i-Mcx/2)+(j-Mcx1/2)*(j-Mcx1/2)) ) ) * 
                                                       y_ini[no_cx[i][j][0]];
    }
    for(i = 0; i < M; ++i)
    for(j = 0; j < M1; ++j){
      if(I_RE == 1) f_ini[no_re[i][j][0]] = f_ini[no_re[i][j][0]] -
           a_ext2.g * exp(-0*sqrt( (double)((i-M/2)*(i-M/2)+(j-M1/2)*(j-M1/2)) ) ) * 
                                                       y_ini[no_re[i][j][0]];   
      if(I_TC == 1) f_ini[no_tc[i][j][0]] = f_ini[no_tc[i][j][0]] - 
           a_ext1.g * exp(-0*sqrt( (double)((i-M/2)*(i-M/2)+(j-M1/2)*(j-M1/2)) ) ) * 
                                                       y_ini[no_tc[i][j][0]]; } 
  }
} 
//-----------ONE shock--------------------------------------------- 
else if(nst == 1) {
  if( ((x > 1500) && (x < 1550)) ) {
     a_ext1.calc(g_ext_tc, x);
     a_ext2.calc(g_ext_re, x);
       if(I_RE == 1) f_ini[no_re[0][0][0]] = f_ini[no_re[0][0][0]] //- g_ext.I;
                     - a_ext2.g * y_ini[no_re[0][0][0]];
  }
}
}















