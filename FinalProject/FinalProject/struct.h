//related structs
#ifndef struct_h_
#define struct_h_


#pragma once
#define _CRT_SECURE_NO_DEPRECATE
#define MAXTALKERS 99
#define PI 3.1415926535897932384626
#define LARRSIZE 10
#define ARRSIZE (1<<LARRSIZE) //1024
#define REALSIZE ((ARRSIZE>>1) + 1) //513?
#define SR 11000  //sampling rate

#include<string>

using namespace std;

//to store the feature data
static float database[MAXTALKERS][27];

//to store the feature of uncertain segment
static float uncertain[MAXTALKERS][28]; 

//store the classified label
static int cLabel[MAXTALKERS];

//the position of the now classified segment in the segSet
static int segPos = 0;
//the position of the now empty uncertain position.
static int uncPos = 0;

// Location for the PDA file, you need to change that to your own.
static const std::string PATH = "D:\\Thesis\\wavDeal\\20secs_longSeg";

//The label for each talker for PDA file, in the file order.
static char* vocNamePDA[16] = { "m1", "fm1", "m2", "m3", "m4", "fm2", "fm3", "fm4", "m5", "fm5", "m6", "fm6", "m7", "fm7", "fm8", "fm9" };


//the wave data struct
typedef struct WAV_DATA{
	short WavData;
}wav_data;

//the voice segment struct
typedef struct VOC_SEG{
	std::string fileName;
	int start;
	int end;
	int channel;
	char* vocName;
}voice_seg;




//all the 6 classes
static string model_name[6] = { "ffs", "mms", "ffd", "mmd", "fmd", "mfd" };

// the arrays to store all the information for 2-peak-GaussianMixModel

//[model_index][peaks_num]
static double weights_2g[6][2], weights_3g[6][3], weights_4g[6][4], weights_5g[6][5], weights_6g[6][6];

//[model_index][peaks_num]
static double consts_2g[6][2], consts_3g[6][3], consts_4g[6][4], consts_5g[6][5], consts_6g[6][6];

//[model_index][peaks_num][feature_dimension]
static double mu_2g[6][2][15], mu_3g[6][3][15], mu_4g[6][4][15], mu_5g[6][5][15], mu_6g[6][6][15];

//[model_index][peaks_num][feature_dimension][feature_dimension]
static double S_inv_2g[6][2][15][15], S_inv_3g[6][3][15][15], S_inv_4g[6][4][15][15], S_inv_5g[6][5][15][15], S_inv_6g[6][6][15][15];





#endif