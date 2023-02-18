#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <thread>
#include <memory>
#include <chrono>
#include <ctime>
#include <sstream>
#include <algorithm>
#include "cmath"
#include "CAENDigitizer.h"
#include "keyb.h"


//ROOT Headers
#include "TPaletteAxis.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TObject.h"
#include "TList.h"
#include "TNamed.h"
#include "THttpCallArg.h"
#include "THttpEngine.h"
#include "THttpWSHandler.h"
#include <mutex>
#include <TMemFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <TH2F.h>
#include <TF1.h>
#include <TProfile.h>
#include <TFrame.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <THttpServer.h>
#include "TTree.h"
#include "TString.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include <TH1D.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TList.h>
#include <TStyle.h>
#include <TH2.h>
#include <TLegend.h>



//#define INDIVIDUAL_TRIGGER_INPUTS
// The following define must be set to the actual number of connected boards
#define MAXNB   1
// NB: the following define MUST specify the ACTUAL max allowed number of board's channels
// it is needed for consistency inside the CAENDigitizer's functions used to allocate the memory
#define MaxNChannels 16

// The following define MUST specify the number of bits used for the energy calculation
#define MAXNBITS 14

/* include some useful functions from file Functions.c
you can find this file in the src directory */
#include "Functions.h"



#define CAEN_USE_DIGITIZERS
#define IGNORE_DPP_DEPRECATED
using std::string;
using std::stringstream;

//#define WFRMPLOT
#define MATRICES
#define WRITEFILE





int checkCommand() {
	int c = 0;
	if(!kbhit())
			return 0;

	c = getch();
	switch (c) {
		case 's': 
			return 9;
			break;
		case 'k':
			return 1;
			break;
		case 'q':
			return 2;
			break;
	}
	return 0;
}





/* ###########################################################################
*  Functions
*  ########################################################################### */

/* --------------------------------------------------------------------------------------------------------- */
/*! \fn      int ProgramDigitizer(int handle, DigitizerParams_t Params, CAEN_DGTZ_DPPParamsPHA_t DPPParams)
*   \brief   Program the registers of the digitizer with the relevant parameters
*   \return  0=success; -1=error */
/* --------------------------------------------------------------------------------------------------------- */
int ProgramDigitizer(int handle, DigitizerParams_t Params, CAEN_DGTZ_DPP_PHA_Params_t DPPParams)
{
    /* This function uses the CAENDigitizer API functions to perform the digitizer's initial configuration */
    int i, ret = 0;

    /* Reset the digitizer */
    ret |= CAEN_DGTZ_Reset(handle);

    if (ret) {
        printf("ERROR: can't reset the digitizer.\n");
        return -1;
    }
    ret |= CAEN_DGTZ_WriteRegister(handle, 0x8000, 0x01000114);  // Channel Control Reg (indiv trg, seq readout) ??

    /* Set the DPP acquisition mode*/
    ret |= CAEN_DGTZ_SetDPPAcquisitionMode(handle, Params.AcqMode, CAEN_DGTZ_DPP_SAVE_PARAM_EnergyAndTime);
    
    // Set the digitizer acquisition mode (CAEN_DGTZ_SW_CONTROLLED or CAEN_DGTZ_S_IN_CONTROLLED)
    ret |= CAEN_DGTZ_SetAcquisitionMode(handle, CAEN_DGTZ_SW_CONTROLLED);
    
    // Set the digitizer syncronization mode
    ret |= CAEN_DGTZ_SetRunSynchronizationMode(handle, CAEN_DGTZ_RUN_SYNC_TrgOutTrgInDaisyChain);
    
    // Set the number of samples for each waveform
    ret |= CAEN_DGTZ_SetRecordLength(handle, Params.RecordLength);

    // Set the I/O level (CAEN_DGTZ_IOLevel_NIM or CAEN_DGTZ_IOLevel_TTL)
    ret |= CAEN_DGTZ_SetIOLevel(handle, Params.IOlev);

    /* Set the digitizer's behaviour when an external trigger arrives:

    CAEN_DGTZ_TRGMODE_DISABLED: do nothing
    CAEN_DGTZ_TRGMODE_EXTOUT_ONLY: generate the Trigger Output signal
    CAEN_DGTZ_TRGMODE_ACQ_ONLY = generate acquisition trigger
    CAEN_DGTZ_TRGMODE_ACQ_AND_EXTOUT = generate both Trigger Output and acquisition trigger
    see CAENDigitizer user manual, chapter "Trigger configuration" for details */
    ret |= CAEN_DGTZ_SetExtTriggerInputMode(handle, CAEN_DGTZ_TRGMODE_ACQ_ONLY);

    // Set the enabled channels
    ret |= CAEN_DGTZ_SetChannelEnableMask(handle, Params.ChannelMask);

    // Set how many events to accumulate in the board memory before being available for readout
    ret |= CAEN_DGTZ_SetDPPEventAggregation(handle, Params.EventAggr, 0);
    

    // Set the DPP specific parameters for the channels in the given channelMask
    ret |= CAEN_DGTZ_SetDPPParameters(handle, Params.ChannelMask, &DPPParams);
    
    for(i=0; i<MaxNChannels; i++) {
        if (Params.ChannelMask & (1<<i)) {
            // Set a DC offset to the input signal to adapt it to digitizer's dynamic range
            ret |= CAEN_DGTZ_SetChannelDCOffset(handle, i, 0x8000);
            
            // Set the Pre-Trigger size (in samples)
            ret |= CAEN_DGTZ_SetDPPPreTriggerSize(handle, i, 1000);
            
            // Set the polarity for the given channel (CAEN_DGTZ_PulsePolarityPositive or CAEN_DGTZ_PulsePolarityNegative)
            ret |= CAEN_DGTZ_SetChannelPulsePolarity(handle, i, Params.PulsePolarity);
        }
    }

    /* Set the virtual probes settings
    DPP-PHA can save:
    2 analog waveforms:
        the first and the second can be specified with the  ANALOG_TRACE 1 and 2 parameters
        
    2 digital waveforms:
        the first can be specified with the DIGITAL_TRACE_1 parameter
        the second  is always the trigger

    CAEN_DGTZ_DPP_VIRTUALPROBE_SINGLE	 -> Save only the ANALOG_TRACE_1 waveform
    CAEN_DGTZ_DPP_VIRTUALPROBE_DUAL      -> Save also the waveform specified in  ANALOG_TRACE_2

    Virtual Probes 1 types:
    CAEN_DGTZ_DPP_VIRTUALPROBE_Input
    CAEN_DGTZ_DPP_VIRTUALPROBE_Delta
    CAEN_DGTZ_DPP_VIRTUALPROBE_Delta2
    CAEN_DGTZ_DPP_VIRTUALPROBE_Trapezoid
    
    Virtual Probes 2 types:
    CAEN_DGTZ_DPP_VIRTUALPROBE_Input
    CAEN_DGTZ_DPP_VIRTUALPROBE_Threshold
    CAEN_DGTZ_DPP_VIRTUALPROBE_TrapezoidReduced
    CAEN_DGTZ_DPP_VIRTUALPROBE_Baseline
    CAEN_DGTZ_DPP_VIRTUALPROBE_None

    Digital Probes types:
    CAEN_DGTZ_DPP_DIGITALPROBE_TRGWin
    CAEN_DGTZ_DPP_DIGITALPROBE_Armed
    CAEN_DGTZ_DPP_DIGITALPROBE_PkRun
    CAEN_DGTZ_DPP_DIGITALPROBE_PileUp
    CAEN_DGTZ_DPP_DIGITALPROBE_Peaking
    CAEN_DGTZ_DPP_DIGITALPROBE_CoincWin
    CAEN_DGTZ_DPP_DIGITALPROBE_BLFreeze
    CAEN_DGTZ_DPP_DIGITALPROBE_TRGHoldoff
	CAEN_DGTZ_DPP_DIGITALPROBE_TRGVal
	CAEN_DGTZ_DPP_DIGITALPROBE_ACQVeto
	CAEN_DGTZ_DPP_DIGITALPROBE_BFMVeto
	CAEN_DGTZ_DPP_DIGITALPROBE_ExtTRG
	CAEN_DGTZ_DPP_DIGITALPROBE_Busy
	CAEN_DGTZ_DPP_DIGITALPROBE_PrgVeto*/

    ret |= CAEN_DGTZ_SetDPP_VirtualProbe(handle, ANALOG_TRACE_1, CAEN_DGTZ_DPP_VIRTUALPROBE_Delta2);
    ret |= CAEN_DGTZ_SetDPP_VirtualProbe(handle, ANALOG_TRACE_2, CAEN_DGTZ_DPP_VIRTUALPROBE_Input);
    ret |= CAEN_DGTZ_SetDPP_VirtualProbe(handle, DIGITAL_TRACE_1, CAEN_DGTZ_DPP_DIGITALPROBE_Peaking);

    if (ret) {
        printf("Warning: errors found during the programming of the digitizer.\nSome settings may not be executed\n");
        return ret;
    } else {
        return 0;
    }
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++ PULSE PROCESSING +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#define beg_baseline_fit_length  20
#define beg_baseline_fit_length_baf2  40
#define end_baseline_fit_length 100
#define the_noise_level 30
float threshold	= 2600;

struct baseline_data{
  Bool_t good_beg,good_end,stable,ok;
  float value;
};


float abss(float xx){if(xx<0){return -xx;}return xx;}


			      

//NEGSIG: GET THRESHOLD X IN SAMPLES
Short_t get_thr_x(float wf[]){    
  for(Int_t ii=0;ii<400;ii++){                             //Max sample to look for the threshold
    if(wf[ii]<threshold&&wf[ii+1]<threshold){return ii;}
  }
  return -10; 
}


////NEGSIG: LED IN SAMPLES 
Float_t get_led(float wf[]){
    int x_f = 0;
    int	x_s = 0;
    float b,m;
    for(UShort_t indx=1;indx<1023;indx++){
	if(wf[indx]>threshold && wf[indx+1]<threshold){x_f = indx; x_s = indx+1;
	break;}
	}
    m = (wf[x_f]-wf[x_s])/(x_f-x_s);
    b = wf[x_f]-(m*x_f);
    return (threshold-b)/m; 
}



////NEGSIG: LED_OLS IN SAMPLES 
Float_t get_led2(float wf[]){
    int x0 = 0;	
    int x1 = 0;
    int x2 = 0;
    int x3 = 0;
    float b,m;
    for(UShort_t tind=1;tind<1023;tind++){
	if(wf[tind]>threshold && wf[tind+1]<threshold && wf[tind+2]<threshold && wf[tind-1]>threshold){x0 = tind; x1 = tind+1; x2 = tind-1; x3 = tind+2;
	break;}
	}
	float xavg = 0.25 * (x0+x1+x2+x3);
	float yavg = 0.25 * (wf[x0]+wf[x1]+wf[x2]+wf[x3]);
	float sumup = (x0-xavg)*(wf[x0]-yavg)+(x1-xavg)*(wf[x1]-yavg)+(x2-xavg)*(wf[x2]-yavg)+(x3-xavg)*(wf[x3]-yavg);
	float sumdwn = pow(2,x0-xavg)+pow(2,x1-xavg)+ pow(2,x2-xavg) + pow(2,x3-xavg);
        m  = sumup/sumdwn ; 
        b = yavg - m*xavg ;
    return (threshold-b)/m; 
}




////NEGSIG: CFD IN SAMPLES 
Float_t get_cft(float wf2[],UShort_t risetime,Short_t thr_x,Float_t fraction, float baseline){
  UShort_t points=600;
  float wf[1023];
  Float_t delayed[points],inv[points],sum[points];
  UShort_t delay=(UShort_t)((1.-fraction)*(1.*risetime)+.5);
  for(UShort_t ip=0;ip<1023;ip++){
    wf[ip] = -wf2[ip]+baseline;
  }
  for(UShort_t ii=0;ii<points;ii++){
    inv[ii]=-fraction*wf[thr_x-1+ii];
    delayed[ii]=wf[thr_x-1-delay+ii];
    sum[ii]=inv[ii]+delayed[ii];
  }
  UShort_t jb0=0,ja0=0;
  for(UShort_t ii=1;ii<points;ii++){if(sum[ii]>0){ja0=ii;jb0=ii-1;break;}}
//  printf("sum A: %f 	sum B: %f \n",sum[ja0],sum[jb0]);
  if((int)(sum[ja0]-sum[jb0]*100)==0){return -10;}
  Float_t zc=-sum[jb0]/(sum[ja0]-sum[jb0])+jb0;
  return (1.*thr_x-1+zc);
}

////NEGSIG: CFD IN SAMPLES TEST DAVIDE
//Float_t get_cft(float wf2[],UShort_t risetime,Short_t thr_x,Float_t fraction, float baseline){
//  UShort_t points=1000;
//  float wf[1023];
//  Float_t delayed[points],inv[points],sum[points];
//  UShort_t delay=(UShort_t)(0.5*(1.*risetime));
//  for(UShort_t ip=0;ip<1023;ip++){
//    wf[ip] = -wf2[ip]+baseline;
//  }
//  for(UShort_t ii=0;ii<points;ii++){
//    inv[ii]=-fraction*wf[ii];
//    delayed[ii]=wf[ii-delay];
//    sum[ii]=inv[ii]+delayed[ii];
//  }
//  UShort_t jb0=0,ja0=0;
//  for(UShort_t ii=1;ii<points;ii++){if(sum[ii]>0){ja0=ii;jb0=ii-1;break;}}
//  printf("sum A: %f 	sum B: %f \n",sum[ja0],sum[jb0]);
//  if((int)(sum[ja0]-sum[jb0]*100)==0){return -10;}
//  Float_t zc=sum[jb0]/(sum[ja0]-sum[jb0])-jb0;
//  return (zc);
//}



//NEGSIG: BASELINE CALC
baseline_data get_baseline(float wf[]){
  UShort_t samples = 1023;
  TGraph *gr = new TGraph(1023); 
  UShort_t min_beg=4096,min_end=4096,max_beg=0,max_end=0;
  for(UShort_t ii=0;ii<samples;ii++){
    gr->SetPoint(ii,ii,wf[ii]);
    if(ii<beg_baseline_fit_length){
      if(min_beg>wf[ii]){min_beg=wf[ii];}
      if(max_beg<wf[ii]){max_beg=wf[ii];}}
    else if(ii>samples-1-end_baseline_fit_length){
      if(min_end>wf[ii]){min_end=wf[ii];}
      if(max_end<wf[ii]){max_end=wf[ii];}}
  }
  baseline_data b;
  b.value=0.;
  b.stable=false;
  if(max_beg-min_beg>the_noise_level){b.good_beg=false;}else{b.good_beg=true;}
  if(max_end-min_end>the_noise_level){b.good_end=false;}else{b.good_end=true;}
  Float_t base_line_beg=0.,base_line_end=0.;
  if(b.good_beg){TF1 *fbase0 = new TF1("fbase0","pol0",0,beg_baseline_fit_length);gr->Fit("fbase0","QRNC");base_line_beg=fbase0->GetParameter(0);delete fbase0;}
  if(b.good_end){TF1 *fbase1 = new TF1("fbase1","pol0",samples-1-end_baseline_fit_length,samples-1);gr->Fit("fbase1","QRNC");
  base_line_end=fbase1->GetParameter(0);delete fbase1;}
  delete gr;
  if(b.good_beg&&b.good_end){if(abss(base_line_beg-base_line_end)>the_noise_level){b.stable=false;}else{b.stable=true;}}
  if(b.stable||b.good_beg){b.value=base_line_beg;}
  if(!b.good_beg&&b.good_end){b.value=base_line_end;}
  if(b.good_beg||b.good_end){b.ok=true;}else{b.ok=false;}
  return b;
}





   
/* ########################################################################### */
/* MAIN                                                                        */
/* ########################################################################### */
int main(int argc, char* argv[])
{	CAEN_DGTZ_ErrorCode ret;
	int	handle[2];
	CAEN_DGTZ_BoardInfo_t BoardInfo;
	CAEN_DGTZ_EventInfo_t eventInfo;
	//CAEN_DGTZ_UINT16_EVENT_t *Evt = NULL;
	CAEN_DGTZ_X742_EVENT_t *Evt = NULL;
	char *buffer = NULL;
	int MajorNumber;
	int i,b;
	int c = 0;
	int Nevent = 0;
	char * evtptr = NULL;
	uint32_t size,bsize;
	uint32_t numEvents;
	Short_t threshold_sampl[32];
	Float_t CFD_sampl[32];
	float bsln[32];
	i = sizeof(CAEN_DGTZ_TriggerMode_t);
	

	CAEN_DGTZ_DPP_PHA_Event_t       *Events[MaxNChannels];  // events buffer
	CAEN_DGTZ_DPP_PHA_Waveforms_t   *Waveform=NULL;     // waveforms buffer

	/* The following variables will store the digitizer configuration parameters */
	CAEN_DGTZ_DPP_PHA_Params_t DPPParams[MAXNB];
	DigitizerParams_t Params[MAXNB];

	/* Arrays for data analysis */
	uint64_t PrevTime[MAXNB][MaxNChannels];
	uint64_t ExtendedTT[MAXNB][MaxNChannels];
	uint32_t *EHisto[MAXNB][MaxNChannels]; // Energy Histograms 
	int ECnt[MAXNB][MaxNChannels];
	int TrgCnt[MAXNB][MaxNChannels];
	int PurCnt[MAXNB][MaxNChannels];


	    /* Other variables */
	    int ch, ev;
	    int Quit=0;
	    int AcqRun = 0;
	    uint32_t AllocatedSize, BufferSize;
	    int Nb=0;
	    int DoSaveWave[MAXNB][MaxNChannels];
	    int BitMask = 0;
	    uint64_t CurrentTime, PrevRateTime, ElapsedTime;
	    uint32_t NumEvents[MaxNChannels];
	    uint32_t temp;
	    memset(DoSaveWave, 0, MAXNB*MaxNChannels*sizeof(int));
	    for (i = 0; i < MAXNBITS; i++)
		BitMask |= 1<<i; /* Create a bit mask based on number of bits of the board */
	
	
 	FILE *fptr;
	fptr = fopen("/run/media/root/JYU23/TEST/test_plot.txt","w");   // TO EDIT BASED ON THE MACHINE !!!!!
	//fptr = fopen("/run/media/root/JYU23/DATA/test_plot.txt","w");   // TO EDIT BASED ON THE MACHINE !!!!!
	
	
	
	/*******ROOT INIT******/
	
	TApplication app("App", &argc, argv);
	float timearr[1024];

	for (int ipt = 0; ipt<1023; ipt ++) {
	timearr[ipt] = ipt*.2  ;        // sample in ns
	}



        #ifdef WFRMPLOT

	TCanvas *canvas = new TCanvas("Live Waveforms", "Live Waveforms",200,10,800,600); 
	canvas->SetGrid();
	TH1F *f = canvas->DrawFrame(0,0,220,5000,"Waveform Monitor; Time [ns]; ADC Channel");
	#endif


	#ifdef MATRICES
	//TApplication appm("App", &argc, argv);
	
	Double_t w = 800;
    	Double_t h = 600;
	TCanvas *canvas2 = new TCanvas("Matrices", "Matrices",w,h);
    	canvas2->SetWindowSize(w + 4, h + 28);
 	auto toftof = new TH2F("toftof","",200,28,34,200,28,34);
	toftof->GetYaxis()->SetTitle("ARM1 TOF [ns]");
	toftof->GetXaxis()->SetTitle("ARM2 TOF [ns]");
	toftof->SetStats(0);
	gStyle->SetPalette(1);
	gPad->Update();
	canvas2->Update();
	toftof->Draw("COLZ");
	canvas2->Update();
	#endif

        //auto server = new THttpServer();

	
	stringstream st;
    	string nomefile;
	time_t now = time(NULL);


    	st <<"file:/run/media/root/JYU23/TEST/"<< "TOSCA_"<< now <<".root";
        //st <<"file:/run/media/root/JYU23/DATA/"<< "TOSCA_"<< now <<".root";
            
    	st >> nomefile;
    	const char *nome = nomefile.c_str();

	TFile *rf = new TFile(nome, "recreate");
   	TTree *T = new TTree("T", "Treedist");
//	float DataTr0[1024];
//	float Ch0[1024];
//	float Ch1[1024];
//	T->Branch("Tr0",&DataTr0,"Tr0[1024]/F");
//	T->Branch("Ch0",&Ch0,"Ch0[1024]/F");
//	T->Branch("Ch1",&Ch1,"Ch1[1024]/F");
//	T->Branch("time",&timearr,"time[1024]/F");
	T->Branch("TCh1",&CFD_sampl[1],"TCh1/F");
	T->Branch("TCh0",&CFD_sampl[0],"TCh0/F");
	T->Branch("TCh16",&CFD_sampl[16],"TCh16/F");
	T->Branch("TCh17",&CFD_sampl[17],"TCh17/F");
	float TOF_CFD1;   //ARM1 sp1-st1  ch1-ch0
	float TOF_CFD2;     //ARM2 sp2-st2  ch17-ch16
	T->Branch("TOF_CFD1",&TOF_CFD1,"TOF_CFD1/F");
	T->Branch("TOF_CFD2",&TOF_CFD2,"TOF_CFD2/F");
	
						     

   



/* **************************************************************** Digitizer INIT Routine ********************************************************** */

    /* *************************************************************************************** */
    /* Set V1742 Parameters                                                                          */
    /* *************************************************************************************** */
    memset(&Params, 0, MAXNB * sizeof(DigitizerParams_t));
    memset(&DPPParams, 0, MAXNB * sizeof(CAEN_DGTZ_DPP_PHA_Params_t));
    for (b = 0; b < MAXNB; b++) {
        for (ch = 0; ch < MaxNChannels; ch++)
            EHisto[b][ch] = NULL; //set all histograms pointers to NULL (we will allocate them later)

        /****************************\
        * Communication Parameters   *
        \****************************/
        
       	/* The following is for VME boards connected using A4818 with PID number p, through an V2718 CONET-VME Bridge
        in this case you must set <LikType> = CAEN_DGTZ_USB_A4818_V2718, <LinkNum_OR_A4818_PID> = p, <ConetNode> = 0 
		and <VMEBaseAddress> = <0xXXXXXXXX> (address of the VME board) */
//        ret = CAEN_DGTZ_OpenDigitizer(CAEN_DGTZ_USB_A4818_V2718, 23365, 0, 0, &handle[1]);
	
	
        Params[b].IOlev = CAEN_DGTZ_IOLevel_NIM;
        /****************************\
        *  Acquisition parameters    *
        \****************************/
        Params[b].AcqMode = CAEN_DGTZ_DPP_ACQ_MODE_Mixed;          // CAEN_DGTZ_DPP_ACQ_MODE_List or CAEN_DGTZ_DPP_ACQ_MODE_Oscilloscope
        Params[b].RecordLength = 2000;                              // Num of samples of the waveforms (only for Oscilloscope mode)
        Params[b].ChannelMask = 0xFFFF;                             // Channel enable mask
        Params[b].EventAggr = 0;                                   // number of events in one aggregate (0=automatic)
        Params[b].PulsePolarity = CAEN_DGTZ_PulsePolarityPositive; // Pulse Polarity (this parameter can be individual)

        /****************************\
        *      DPP parameters        *
        \****************************/
        for(ch=0; ch<MaxNChannels; ch++) {
			DPPParams[b].thr[ch] = 100;   // Trigger Threshold (in LSB)
			DPPParams[b].k[ch] = 3000;     // Trapezoid Rise Time (ns) 
			DPPParams[b].m[ch] = 900;      // Trapezoid Flat Top  (ns) 
			DPPParams[b].M[ch] = 50000;      // Decay Time Constant (ns) 
			DPPParams[b].ftd[ch] = 500;    // Flat top delay (peaking time) (ns) 
			DPPParams[b].a[ch] = 4;       // Trigger Filter smoothing factor (number of samples to average for RC-CR2 filter) Options: 1; 2; 4; 8; 16; 32
			DPPParams[b].b[ch] = 200;     // Input Signal Rise time (ns) 
			DPPParams[b].trgho[ch] = 1200;  // Trigger Hold Off
			DPPParams[b].nsbl[ch] = 4;     //number of samples for baseline average calculation. Options: 1->16 samples; 2->64 samples; 3->256 samples; 4->1024 samples; 5->4096 samples; 6->16384 samples
			DPPParams[b].nspk[ch] = 0;     //Peak mean (number of samples to average for trapezoid height calculation). Options: 0-> 1 sample; 1->4 samples; 2->16 samples; 3->64 samples
			DPPParams[b].pkho[ch] = 2000;  //peak holdoff (ns)
			DPPParams[b].blho[ch] = 500;   //Baseline holdoff (ns)
			DPPParams[b].enf[ch] = 1.0; // Energy Normalization Factor
			DPPParams[b].decimation[ch] = 0;  //decimation (the input signal samples are averaged within this number of samples): 0 ->disabled; 1->2 samples; 2->4 samples; 3->8 samples
			DPPParams[b].dgain[ch] = 0;    //decimation gain. Options: 0->DigitalGain=1; 1->DigitalGain=2 (only with decimation >= 2samples); 2->DigitalGain=4 (only with decimation >= 4samples); 3->DigitalGain=8( only with decimation = 8samples).
			DPPParams[b].otrej[ch] = 0;
			DPPParams[b].trgwin[ch] = 0;  //Enable Rise time Discrimination. Options: 0->disabled; 1->enabled
			DPPParams[b].twwdt[ch] = 100;  //Rise Time Validation Window (ns)
        }
    }





































    for(b=0; b<MAXNB; b++){
        

        /* ************** Print SW Info ********************* */
        

	printf("\n ************************************************************************* \n");
    printf(" ************ Readout Software For CAEN V1742 via A4818  ***************** \n");
	printf(" ************                              Version: 1.01 ***************** \n");
    printf(" ************            Authors: D. Panico, A. Di Nitto ***************** \n");
	printf(" ************************************************************************* \n");
                                

	

	/* *************** A4818 INIT *********************** */

	int32_t PID = 23365;
        
	
    ret = CAEN_DGTZ_OpenDigitizer2(CAEN_DGTZ_USB_A4818,&PID,0,0,&handle[0]);     //V1742
	ret = CAEN_DGTZ_OpenDigitizer2(CAEN_DGTZ_USB_A4818,&PID,0,1,&handle[1]);     //V1725
	
        if(ret != CAEN_DGTZ_Success) {
            printf("\n!!!!!!!!!!!!!         Error: CANNOT OPEN DIGITIZER!     !!!!!!!!!!!!!!!!! \n");
            goto QuitProgram;
        }


        

        /* Once we have the handler to the digitizer, we use it to call the other functions */

        ret = CAEN_DGTZ_GetInfo(handle[0], &BoardInfo);
        printf("\n|-----> Connected to CAEN Digitizer Model %s, recognized as board %d\n", BoardInfo.ModelName, 0);
        printf("\t|-----> ROC FPGA Release is %s\n", BoardInfo.ROC_FirmwareRel);
        printf("\t|-----> AMC FPGA Release is %s\n", BoardInfo.AMC_FirmwareRel);
        
        
        ret = CAEN_DGTZ_GetInfo(handle[1], &BoardInfo);
        printf("\n|-----> Connected to CAEN Digitizer Model %s, recognized as board %d\n", BoardInfo.ModelName, 1);
        printf("\t|-----> ROC FPGA Release is %s\n", BoardInfo.ROC_FirmwareRel);
        printf("\t|-----> AMC FPGA Release is %s\n", BoardInfo.AMC_FirmwareRel);
	    
		
		
	/* *************************************************************************************** */
        /* Program the V1742 digitizer                                                             */
        /* *************************************************************************************** */		

        /* *************** SAMPLING PARAMETERS ******************** */

        ret = CAEN_DGTZ_Reset(handle[b]);                                               /* Reset Digitizer */
	ret = CAEN_DGTZ_SetDRS4SamplingFrequency(handle[b],CAEN_DGTZ_DRS4_5GHz);        /* Set the sampling frequency (Default = 5GHz) */
	ret = CAEN_DGTZ_GetInfo(handle[b], &BoardInfo);                                 /* Get Board Info */
	ret = CAEN_DGTZ_SetRecordLength(handle[b],1024);                                /* Set the lenght of each waveform (in samples) */
	ret = CAEN_DGTZ_SetGroupEnableMask(handle[b],15);                                /* Enable all groups  */
	//ret = CAEN_DGTZ_WriteRegister(handle[b],
	ret = CAEN_DGTZ_SetPostTriggerSize(handle[b],15);

	ret = CAEN_DGTZ_WriteRegister(handle[b],0x1098,0xF6C00);                        /* DC offset GR0 -> 0x6C00 */
	ret = CAEN_DGTZ_WriteRegister(handle[b],0x1198,0xF6C00);       			/* DC offset GR1 -> 0x6C00 */
	ret = CAEN_DGTZ_WriteRegister(handle[b],0x1298,0xF6C00);			/* DC offset GR2 -> 0x6C00 */
	ret = CAEN_DGTZ_WriteRegister(handle[b],0x1398,0xF6C00);			/* DC offset GR3 -> 0x6C00 */
	

        /* ************** TRIGGER PARAMETERS ********************** */

       //ret = CAEN_DGTZ_SetGroupTriggerThreshold(handle[b],0,32768);                    /* Set selfTrigger threshold */
//    ret = CAEN_DGTZ_SetFastTriggerMode(handle[b],CAEN_DGTZ_TRGMODE_ACQ_ONLY,15);   /* Set trigger on channel 0 to be ACQ_ONLY */
    ret = CAEN_DGTZ_SetSWTriggerMode(handle[b],CAEN_DGTZ_TRGMODE_DISABLED);         /* Set the behaviour when a SW tirgger arrives */

	ret = CAEN_DGTZ_SetFastTriggerMode(handle[b],CAEN_DGTZ_TRGMODE_ACQ_ONLY);       /* Enable the TRn as the local trigger  */
	ret = CAEN_DGTZ_SetFastTriggerDigitizing(handle[b], CAEN_DGTZ_ENABLE);          /* Enable the TRn input signal digitization */
	ret = CAEN_DGTZ_SetGroupFastTriggerDCOffset(handle[b],15, 32768);          	/* Set the TRn input signal DC Offset */
	ret = CAEN_DGTZ_SetGroupFastTriggerThreshold(handle[b],15, 20934);          	/* Set the TRn input signal threshold */
 
    ret = CAEN_DGTZ_SetOutputSignalMode(handle[b],CAEN_DGTZ_FASTTRG_ACCEPTED);
        
        
	ret = CAEN_DGTZ_LoadDRS4CorrectionData(handle[b], CAEN_DGTZ_DRS4_5GHz);
	ret = CAEN_DGTZ_EnableDRS4Correction(handle[b]);

	ret = CAEN_DGTZ_SetMaxNumEventsBLT(handle[b],1000);                              /* Set the max number of events to transfer in a sigle readout */
        ret = CAEN_DGTZ_SetAcquisitionMode(handle[b],CAEN_DGTZ_SW_CONTROLLED);          /* Set the acquisition mode */


        
	/* *************************************************************************************** */
        /* Program the V1725 digitizer (see function ProgramDigitizer)                             */
        /* *************************************************************************************** */
         for (b = 0; b < MAXNB; b++) {
              ret = ProgramDigitizer(handle[1], Params[b], DPPParams[b]);
              if (ret) {
              printf("Failed to program the V1725 digitizer\n");
            goto QuitProgram;
           }
          }

    /* WARNING: The mallocs MUST be done after the digitizer programming,
    because the following functions needs to know the digitizer configuration
    to allocate the right memory amount */
    /* Allocate memory for the readout buffer */
    ret = CAEN_DGTZ_MallocReadoutBuffer(handle[1], &buffer, &AllocatedSize);
    /* Allocate memory for the events */
    ret = CAEN_DGTZ_MallocDPPEvents(handle[1], (void **) &Events, &AllocatedSize); 
    /* Allocate memory for the waveforms */
    ret = CAEN_DGTZ_MallocDPPWaveforms(handle[1], (void **) &Waveform, &AllocatedSize); 
    if (ret) {
        printf("Can't allocate memory buffers\n");
        goto QuitProgram;    
    }



	
        if(ret != CAEN_DGTZ_Success) {
            printf("\n !!!!!!!!!!!!  Error during Digitizer Configuration: %d     !!!!!!!!!!!!! \n", ret);
	    Sleep(100);
            goto QuitProgram;
	    
        } else { 
	  
	     printf("\n|----->  Successful Connection \n");
	     Sleep(100);
	}

    }




/* ******************************************* Acquisition Control ***************************************************/


    printf("\n\n|-----> Press 'S' to start the acquisition |\n");
        printf("|-----> Press 'K' to stop  the acquisition |\n");
        printf("|-----> Press 'Q' to quit  the application |\n\n");

    while (1) {
		c = checkCommand();
		if (c == 9) break;
		if (c == 2) return -1;
		Sleep(100);
    }




    
    /*************************************************************** V1742 Memory Readout Buffer Allocation *****************************************************/
    printf("|-----> Reading |\n");

    ret = CAEN_DGTZ_MallocReadoutBuffer(handle[b],&buffer,&size);

 
     /*************************************************************** Start Acquisition ********************************************************************/
    
    for(b=0; b<MAXNB; b++){
            ret = CAEN_DGTZ_SWStartAcquisition(handle[0]);
            ret = CAEN_DGTZ_SWStartAcquisition(handle[1]);
    
    /************************************************************** Start acquisition loop ****************************************************************/
    
    
	while(1) {
        for(b=0; b<MAXNB; b++) {
		    

		    ret = CAEN_DGTZ_SendSWtrigger(handle[b]); /* ???? Send a SW Trigger ???? */
		    ret = CAEN_DGTZ_ReadData(handle[b],CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT,buffer,&bsize); /* Read the buffer from the digitizer */

		    ret = CAEN_DGTZ_GetNumEvents(handle[b],buffer,bsize,&numEvents);

		   
		    
	


             
		    for (i=0;i<numEvents;i++) {
               				 /* Get the Infos and pointer to the event */
			    ret = CAEN_DGTZ_GetEventInfo(handle[b],buffer,bsize,i,&eventInfo,&evtptr);

               				 /* Decode the event to get the data */

			    //*************************************
			    // Event Elaboration
			    //*************************************

			      ret = CAEN_DGTZ_AllocateEvent(handle[b],(void**) &Evt);
			      ret = CAEN_DGTZ_DecodeEvent(handle[b],evtptr,(void**) &Evt);
			      Nevent +=1 ;

			
				threshold_sampl[0] = get_thr_x((Evt->DataGroup[0]).DataChannel[0]);
				bsln[0] = get_baseline((Evt->DataGroup[0]).DataChannel[0]).value;
				CFD_sampl[0] = get_cft((Evt->DataGroup[0]).DataChannel[0], 20, threshold_sampl[0], 0.45, bsln[0]);

				threshold_sampl[1] = get_thr_x((Evt->DataGroup[0]).DataChannel[1]);
				bsln[1] = get_baseline((Evt->DataGroup[0]).DataChannel[1]).value;
				CFD_sampl[1] = get_cft((Evt->DataGroup[0]).DataChannel[1], 20, threshold_sampl[1], 0.45, bsln[1]);

				threshold_sampl[17] = get_thr_x((Evt->DataGroup[2]).DataChannel[1]);
				bsln[17] = get_baseline((Evt->DataGroup[2]).DataChannel[1]).value;
				CFD_sampl[17] = get_cft((Evt->DataGroup[2]).DataChannel[1], 20, threshold_sampl[17], 0.45, bsln[17]);

				threshold_sampl[16] = get_thr_x((Evt->DataGroup[2]).DataChannel[0]);
				bsln[16] = get_baseline((Evt->DataGroup[2]).DataChannel[0]).value;
				CFD_sampl[16] = get_cft((Evt->DataGroup[2]).DataChannel[0], 20, threshold_sampl[16], 0.45, bsln[16]);


		
					
			
				TOF_CFD1 = (CFD_sampl[1]-CFD_sampl[0])*.2;    // Arm1 (ch1-ch0)
				TOF_CFD2 = (CFD_sampl[17]-CFD_sampl[16])*.2;    // Arm2 (ch17-ch16)
				
				

				printf("---> EVENT #%d | TOF ARM1: %f | TOF ARM2: %f \n",Nevent,TOF_CFD1,TOF_CFD2);




                // Clear Histograms and counters
                for (b = 0; b < MAXNB; b++) {
                    for (ch = 0; ch < MaxNChannels; ch++) {
                        EHisto[b][ch] = (uint32_t *)malloc((1 << MAXNBITS) * sizeof(uint32_t));
                        memset(EHisto[b][ch], 0, (1 << MAXNBITS) * sizeof(uint32_t));
                        TrgCnt[b][ch] = 0;
                        ECnt[b][ch] = 0;
                        PrevTime[b][ch] = 0;
                        ExtendedTT[b][ch] = 0;
                        PurCnt[b][ch] = 0;
                    }
                }
                PrevRateTime = get_time();
                AcqRun = 0;
                PrintInterface();
                
                            
                AcqRun = 1;
                    
                    /* Read data from the boards */
                    for (b = 0; b < MAXNB; b++) {
                        /* Read data from the board */
                        ret = CAEN_DGTZ_ReadData(handle[1], CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT, buffer, &BufferSize);
                        if (ret) {
                            printf("Readout Error\n");
                            goto QuitProgram;
                        }
                        if (BufferSize == 0)
                            continue;

                        Nb += BufferSize;
                        //ret = DataConsistencyCheck((uint32_t *)buffer, BufferSize/4);
                        ret |= CAEN_DGTZ_GetDPPEvents(handle[1], buffer, BufferSize, (void**) &Events, &NumEvents);
                        if (ret) {
                            printf("Data Error: %d\n", ret);
                            goto QuitProgram;
                        }

                        /* Analyze data */
                        //for(b=0; b<MAXNB; b++) printf("%d now: %d\n", b, Params[b].ChannelMask);
                        for (ch = 0; ch < MaxNChannels; ch++) {
                            if (!(Params[b].ChannelMask & (1<<ch)))
                                continue;
                            
                            /* Update Histograms */
                            for (ev = 0; ev < NumEvents[ch]; ev++) {
                                TrgCnt[b][ch]++;
                                /* Time Tag */
                                if (Events[ch][ev].TimeTag < PrevTime[b][ch])
                                    ExtendedTT[b][ch]++;
                                PrevTime[b][ch] = Events[ch][ev].TimeTag;
                                /* Energy */
                                if (Events[ch][ev].Energy > 0) {
                                    // Fill the histograms
                                    EHisto[b][ch][(Events[ch][ev].Energy)&BitMask]++;
                                    ECnt[b][ch]++;
                                } else {  /* PileUp */
                                    PurCnt[b][ch]++;
                                }
                                
                                
                                /* Get Waveforms (only from 1st event in the buffer) */
                                if ((Params[b].AcqMode != CAEN_DGTZ_DPP_ACQ_MODE_List) && DoSaveWave[b][ch] && (ev == 0)) {
                                    int size;
                                    int16_t *WaveLine;
                                    uint8_t *DigitalWaveLine;
                                    CAEN_DGTZ_DecodeDPPWaveforms(handle[b], &Events[ch][ev], Waveform);

                                    // Use waveform data here...
                                    size = (int)(Waveform->Ns); // Number of samples
                                    WaveLine = Waveform->Trace1; // First trace (ANALOG_TRACE_1)
                                    SaveWaveform(b, ch, 1, size, WaveLine);

                                    WaveLine = Waveform->Trace2; // Second Trace ANALOG_TRACE_2 (if single trace mode, it is a sequence of zeroes)
                                    SaveWaveform(b, ch, 2, size, WaveLine);

                                    DigitalWaveLine = Waveform->DTrace1; // First Digital Trace (DIGITALPROBE1)
                                    SaveDigitalProbe(b, ch, 1, size, DigitalWaveLine);

                                    DigitalWaveLine = Waveform->DTrace2; // Second Digital Trace (for DPP-PHA it is ALWAYS Trigger)
                                    SaveDigitalProbe(b, ch, 2, size, DigitalWaveLine);
                                    DoSaveWave[b][ch] = 0;
                                    printf("Waveforms saved to 'Waveform_<board>_<channel>_<trace>.txt'\n");
                                    
                                } // loop to save waves
                            } // loop on events
                        } // loop on channels
                    } // loop on boards












































			      		/* LIVE PLOT */

					#ifdef WFRMPLOT

					
					
					auto g1 = new TGraph(1023,timearr,((Evt->DataGroup[0]).DataChannel[0]));
					g1->SetLineColor(4);
					g1->SetLineWidth(2);
					
					auto g2 = new TGraph(1023,timearr,((Evt->DataGroup[0]).DataChannel[1]));
					g2->SetLineWidth(2);
					g2->SetLineColor(0);

					auto g3 = new TGraph(1023,timearr,((Evt->DataGroup[1]).DataChannel[0]));
					g3->SetLineWidth(2);
					g3->SetLineColor(1);

					auto g4 = new TGraph(1023,timearr,((Evt->DataGroup[1]).DataChannel[1]));
					g4->SetLineWidth(2);
					g4->SetLineColor(3);
					
                                      

					g1->Draw("L");
					g2->Draw("same");
					g3->Draw("same");
					g4->Draw("same");
					auto legend = new TLegend();
					legend->AddEntry(g4,"ARM2: STOP");
					legend->AddEntry(g2,"ARM1: STOP");
					legend->AddEntry(g3,"ARM2: START");
					legend->AddEntry(g1,"ARM1: START");
					legend->Draw();
					//serv->Register("graph/subfolder",g1);	
					canvas->Update();
					//getchar();
					//usleep(1000000);
					delete g1;
					delete g2;
					delete g3;
					delete g4;
					delete legend;
					
					#endif






					/* Plot Matrices */

					#ifdef MATRICES

					toftof->Fill(TOF_CFD2,TOF_CFD1);
					toftof->Draw("COLZ");
					canvas2->Update();

					
					
					

					#endif






			    		
					/* Write data to file */
					#ifdef WRITEFILE
					
//					for (int dp = 0; dp < 1023; dp++) {     
//						      DataTr0[dp] = ((Evt->DataGroup[0]).DataChannel[8])[dp];
//							Ch0[dp] = ((Evt->DataGroup[0]).DataChannel[0])[dp];
//							Ch1[dp] = ((Evt->DataGroup[0]).DataChannel[1])[dp];     
//							}
//					
					T->Fill();
					

		 			      //for (int k = 0; k < 2; k++){
						//for (int j = 0; j < 9; j++) { 
						  //fprintf(fptr,"\nGroup %d | Ch %d : ",k,j);
						    //for (int i = 0; i < 1023; i++) {     
						      //fprintf(fptr," %f ", ((Evt->DataGroup[k]).DataChannel[j])[i]);     
							//}      
						      //}
						    //}


					#endif
			

			    ret = CAEN_DGTZ_FreeEvent(handle[b],(void**) &Evt);
                CAEN_DGTZ_FreeDPPEvents(handle[1], (void**) &Events);
                CAEN_DGTZ_FreeDPPWaveforms(handle[1],(void**) &Waveform);
			    //delete f1;
					

		   }



  
		    c = checkCommand();
		    if (c == 1) goto Continue;
		    if (c == 2) goto Continue;
        } // end of loop on boards
    } // end of readout loop

Continue:
    for(b=0; b<MAXNB; b++){
	printf("\n|-----> Board %d: Retrieved %d Events\n",b, Nevent);
	printf("\n|----->  Output in TOSCA_%ld.root \n",now);
	}
    goto QuitProgram;



/******************************************* Quit program routine *************************************************/
QuitProgram:

    // Free the buffers and close the digitizers
	ret = CAEN_DGTZ_SWStopAcquisition(handle[0]);
        ret = CAEN_DGTZ_SWStopAcquisition(handle[1]);
	ret = CAEN_DGTZ_FreeReadoutBuffer(&buffer);
	fclose(fptr);

// Close Root file
       T->Write();
       rf->Write();
       rf->Close();
        
         

    for(b=0; b<MAXNB; b++){
        ret = CAEN_DGTZ_CloseDigitizer(handle[b]);
        ret = CAEN_DGTZ_CloseDigitizer(handle[1]);
    }
    printf("\n|-----> Press ENTER to quit... \n");
    getchar();
    printf("\n|-----> Quitting Program... \n");
   
 }
}












