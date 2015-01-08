#ifndef functions_h_
#define functions_h_

#define _CRT_SECURE_NO_DEPRECATE
#pragma once;
#include<ctime>
#include <iostream>
#include<windows.h>
#include<string>
#include"struct.h";
#include "Brian_Model_Matrix.h";
#include "PCA_matrix.h"
#include <fstream>
#include <ctime>
#include<iomanip>

#include <stdio.h>
#include <io.h>

//to calculate the fft for a given segment, extern from fft.c
extern "C"void realfft(float *snd, float *real, float *imag, long log2blocksiz, long dir);


//build the complete path string for all the voice files, the test dataset
void  buildFilePath(std::string path[16][16]);

//to reorder a segment set in a random way
voice_seg* Reorder(voice_seg *segSet, int num);

//given the segment and wanted length, come out with a input with random start position.
int inSegRandom(int start, int end, int length);

//recalculate the true label for the reordered segments based on vocName
int *makeTrueLabel(voice_seg *segSet, int num);


// given n, return the maximum value of the set before the nth position.
//This is used to judge the coming new label number of the current segment
int getMaxValue(int *array, int n);


//return the nth 2-character-set of the number.
//we use 2 digits to represent a talker, when segments from same talker been labeled with different labels,
//we use a long string to store all the labels, each one is 2 digits. This function is used to get the nth label.
int getCharInt(int targ, int n);


//this function is used inside the getAccuracy function
int getMaxValueTemp(int *array, int *offset, int n);



//return the Aaccuracy, given the true label and classfied label
double getAccuracy(int *trueLabel, int *Label, int num);


//calculate the pval for classify function
//this function calculate the probability given the segment and a specific model
//the parameters are the number of peaks of GMM model, and the dimension of the features.
//if you modified any of the above parameter, you need to modify the function inside to get a right answer
double getpval(float x[29], const double mu[2][29], const double S_inv[2][29][29], const double consts[2], const double weights[2], int gauss);


//calculate the pval for classify function
//this function calculate the probability given the segment and a specific model
//the parameters are the number of peaks of GMM model, and the dimension of the features.
//if you modified any of the above parameter, you need to modify the function inside to get a right answer
double getpval_PCA_2g(float x[15], const double mu[2][15], const double S_inv[2][15][15], const double consts[2], const double weights[2], int gauss);
double getpval_PCA_3g(float x[15], const double mu[3][15], const double S_inv[3][15][15], const double consts[3], const double weights[3], int gauss);
double getpval_PCA_4g(float x[15], const double mu[4][15], const double S_inv[4][15][15], const double consts[4], const double weights[4], int gauss);
double getpval_PCA_5g(float x[15], const double mu[5][15], const double S_inv[5][15][15], const double consts[5], const double weights[5], int gauss);
double getpval_PCA_6g(float x[15], const double mu[6][15], const double S_inv[6][15][15], const double consts[6], const double weights[6], int gauss);

//the main classification function
//input: voice segment
//input: length
//input/output: database
//input: threshold
float classify(short soundData[], int length, double pvals[MAXTALKERS][6], int clear_database, double threshold);


//given a voice segment, to decide the number of voiced frames
unsigned count_vframes(short soundData[], int length, int clear_database);


//test the average number of frames given different length of segment.
void count_vframes_test(unsigned frame_count[20]);




//the function to read all model information from txt file
void readTxtData();


//to build the real database for model testing
void buildDatabase();


//given the dataSet, the function returns the 27-dimension feature extracted from the data
//float* getFeat(short soundData[], unsigned length);
void getFeat(short soundData[], unsigned length, float* features, unsigned &flag);

#endif

