#include "functions.h"

//build path string set for PDA file.
void  buildFilePath(std::string path[16][16])
{
	for (int i = 1; i <= 9; i++)
	{
		for (int j = 1; j <= 16; j++)
		{
			path[i - 1][j - 1] = PATH + "\\aPDAm0" + std::to_string(i) + "_0" + std::to_string(j) + "_close.wav";
		}
	}
	for (int i = 10; i <= 16; i++)
	{
		for (int j = 1; j <= 16; j++)
		{
			path[i - 1][j - 1] = PATH + "\\aPDAm" + std::to_string(i) + "_0" + std::to_string(j) + "_close.wav";
		}
	}
}

//to reorder a segment set in a random way
voice_seg* Reorder(voice_seg *segSet, int num)
{
	LARGE_INTEGER nStartCounter;
	voice_seg *orderSet = new voice_seg[num];
	int *usedLabel = new int[num];
	int i = 0;

	memset(usedLabel, 0, sizeof(usedLabel));

	while (i < num)
	{
		::QueryPerformanceCounter(&nStartCounter);
		::srand((unsigned)nStartCounter.LowPart);
		int temp = rand() % num;
		while (usedLabel[temp] == 1)
		{
			::QueryPerformanceCounter(&nStartCounter);
			::srand((unsigned)nStartCounter.LowPart);
			temp = rand() % num;
		}
		usedLabel[temp] = 1;
		orderSet[i] = segSet[temp];
		i++;
	}
	return orderSet;
}


//given the segment and wanted length, come out with a input with random start position.
int inSegRandom(int start, int end, int length)
{
	if (length >= (end - start)) return start;

	else
	{
		LARGE_INTEGER nStartCounter;
		::QueryPerformanceCounter(&nStartCounter);
		::srand((unsigned)nStartCounter.LowPart);
		return (rand() % (end - start - length) + start);
	}
}



//come out with the true Label for a given random segment Set.
//the label is based on the vocName string given by the voice segment struct
int *makeTrueLabel(voice_seg *segSet, int num)
{
	int *trueLabel = new int[num];
	int labelCnt = 0;
	int i = 0;
	while (i<num)
	{
		if (i == 0)
		{
			trueLabel[i] = 1;
			labelCnt++;
		}
		else
		{
			for (int j = 0; j < i; j++)
			{
				if (segSet[i].vocName == segSet[j].vocName)
				{
					trueLabel[i] = trueLabel[j];
					goto next;
				}
			}
			labelCnt++;
			trueLabel[i] = labelCnt;
		}
	next:		i++;
	}
	return trueLabel;
}



// given n, return the maximum value of the set before the nth position.
int getMaxValue(int *array, int n)
{
	if (n == 0) return array[0];

	int temp = 0;
	for (int i = 0; i < n; i++)
	{
		if (temp < array[i]) temp = array[i];
	}
	return temp;
}



//return the nth 2-character-set of the number.
int getCharInt(int targ, int n)
{
	int temp = pow(100, n + 1);
	return(floor((targ%temp) / pow(100, n)));
}



//this function is used inside the getAccuracy function
int getMaxValueTemp(int *array, int *offset, int n)
{
	int temp = 0;
	for (int i = 0; i < n; i++)
	{
		if (temp < array[i] + offset[i]) temp = array[i] + offset[i];
	}
	return temp;
}



//return the Aaccuracy, given the true label and classified label
//TODO: some errors exist in the function
double getAccuracy(int *trueLabel, int *Label, int num)
{
	int *temp = new int[num];
	for (int i = 0; i <num; i++)
	{
		temp[i] = trueLabel[i];
	}

	int accCount = 0;
	int *offset = new int[num]; // for return to original label.
	int ignCount = 0;

	for (int i = 0; i < num; i++) offset[i] = 0;
	for (int i = 0; i < num; i++)
	{
		int accLabel = 0;
		if (Label[i] == 0)
		{
			ignCount++;
			continue;
		}
		//if true label is a new class, but classified as not, need to update the later labels
		if (i != 0 && Label[i] <= getMaxValue(Label, i) && temp[i] + offset[i]>getMaxValueTemp(temp, offset, i))
		{
			for (int j = i + 1; j < num; j++)
			{
				if (temp[j] + offset[j] >getMaxValueTemp(temp, offset, i + 1) && (temp[j]<100)){
					temp[j] -= 1;
					offset[j] += 1;
				}
			}
		}

		//if classified as new label but actually not, need to update the later labels
		if (i != 0 && Label[i] > getMaxValue(Label, i) && temp[i] + offset[i] <= getMaxValueTemp(temp, offset, i))
		{
			for (int j = i + 1; j < num; j++)
			{
				if (temp[j] + offset[j] >getMaxValueTemp(temp, offset, i + 1) && (temp[j]<100)){
					temp[j] += 1;
					offset[j] -= 1;
				}
			}
		}

		for (int j = 0; j < 10; j++)
		{
			if (Label[i] == getCharInt(temp[i], j)){
				accLabel = 1;
				break;
			}
		}
		accCount += accLabel;
		if (!accLabel&&Label[i]>getMaxValue(Label, i))
		{
			for (int m = i + 1; m < num; m++)
			{
				if (temp[m] == temp[i]) temp[m] = 100 * temp[i] + Label[i];
			}
		}
	}
	return((accCount - ignCount) / (1.0*(num - ignCount)));
}




//calculate the pval for classify function
//this function calculate the probability given the segment and a specific model
//the parameters are the number of peaks of GMM model, and the dimension of the features.
//if you modified any of the above parameter, you need to modify the function inside to get a right answer
double getpval(float x[29], const double mu[2][29], const double S_inv[2][29][29], const double consts[2], const double weights[2], int gauss)
{
	int i, j;
	double arr1[29], arr2[29], pval = 0;

	for (i = 0; i < 29; i++)
	{
		arr1[i] = x[i] - mu[gauss][i];
		arr2[i] = 0;
	}

	for (i = 0; i < 29; i++)
	{
		for (j = 0; j < 29; j++)
		{
			arr2[i] += arr1[j] * S_inv[gauss][j][i];
		}
	}

	for (i = 0; i < 29; i++)
	{
		pval += arr2[i] * arr1[i];
	}

	pval = weights[gauss] * consts[gauss] * exp(-pval / 2);
	return pval;
}

//calculate the pval for classify function
//this function calculate the probability given the segment and a specific model
//the parameters are the number of peaks of GMM model, and the dimension of the features.
//if you modified any of the above parameter, you need to modify the function inside to get a right answer
double getpval_PCA_2g(float x[15], const double mu[2][15], const double S_inv[2][15][15], const double consts[2], const double weights[2], int gauss)
{

	int i, j;
	double arr1[15], arr2[15], pval = 0;

	for (i = 0; i < 15; i++)
	{
		arr1[i] = x[i] - mu[gauss][i];
		arr2[i] = 0;
	}

	for (i = 0; i < 15; i++)
	{
		for (j = 0; j < 15; j++)
		{
			arr2[i] += arr1[j] * S_inv[gauss][j][i];
		}
	}

	for (i = 0; i < 15; i++)
	{
		pval += arr2[i] * arr1[i];
	}

	pval = weights[gauss] * consts[gauss] * exp(-pval / 2);
	return pval;
}

double getpval_PCA_3g(float x[15], const double mu[3][15], const double S_inv[3][15][15], const double consts[3], const double weights[3], int gauss)
{

	int i, j;
	double arr1[15], arr2[15], pval = 0;

	for (i = 0; i < 15; i++)
	{
		arr1[i] = x[i] - mu[gauss][i];
		arr2[i] = 0;
	}

	for (i = 0; i < 15; i++)
	{
		for (j = 0; j < 15; j++)
		{
			arr2[i] += arr1[j] * S_inv[gauss][j][i];
		}
	}

	for (i = 0; i < 15; i++)
	{
		pval += arr2[i] * arr1[i];
	}

	pval = weights[gauss] * consts[gauss] * exp(-pval / 2);
	return pval;
}

double getpval_PCA_4g(float x[15], const double mu[4][15], const double S_inv[4][15][15], const double consts[4], const double weights[4], int gauss)
{

	int i, j;
	double arr1[15], arr2[15], pval = 0;

	for (i = 0; i < 15; i++)
	{
		arr1[i] = x[i] - mu[gauss][i];
		arr2[i] = 0;
	}

	for (i = 0; i < 15; i++)
	{
		for (j = 0; j < 15; j++)
		{
			arr2[i] += arr1[j] * S_inv[gauss][j][i];
		}
	}

	for (i = 0; i < 15; i++)
	{
		pval += arr2[i] * arr1[i];
	}

	pval = weights[gauss] * consts[gauss] * exp(-pval / 2);
	return pval;
}
double getpval_PCA_5g(float x[15], const double mu[5][15], const double S_inv[5][15][15], const double consts[5], const double weights[5], int gauss)
{

	int i, j;
	double arr1[15], arr2[15], pval = 0;

	for (i = 0; i < 15; i++)
	{
		arr1[i] = x[i] - mu[gauss][i];
		arr2[i] = 0;
	}

	for (i = 0; i < 15; i++)
	{
		for (j = 0; j < 15; j++)
		{
			arr2[i] += arr1[j] * S_inv[gauss][j][i];
		}
	}

	for (i = 0; i < 15; i++)
	{
		pval += arr2[i] * arr1[i];
	}

	pval = weights[gauss] * consts[gauss] * exp(-pval / 2);
	return pval;
}
double getpval_PCA_6g(float x[15], const double mu[6][15], const double S_inv[6][15][15], const double consts[6], const double weights[6], int gauss)
{

	int i, j;
	double arr1[15], arr2[15], pval = 0;

	for (i = 0; i < 15; i++)
	{
		arr1[i] = x[i] - mu[gauss][i];
		arr2[i] = 0;
	}

	for (i = 0; i < 15; i++)
	{
		for (j = 0; j < 15; j++)
		{
			arr2[i] += arr1[j] * S_inv[gauss][j][i];
		}
	}

	for (i = 0; i < 15; i++)
	{
		pval += arr2[i] * arr1[i];
	}

	pval = weights[gauss] * consts[gauss] * exp(-pval / 2);
	return pval;
}

//the main classification function
//input: voice segment
//input: length
//input/output: database
//input: threshold
float classify(short soundData[], int length, double pvals[MAXTALKERS][6], int clear_database, double threshold)
{
	float dbRecordWeight[MAXTALKERS];
	int i, j, k = 0;
	int startPos, cepsMaxInd, voicedCount, numBlocks = 0;
	int width, pvalsMaxInd;
	static int featCount[MAXTALKERS], talkerPos;

	float specMag[REALSIZE], sndSq[REALSIZE], melFreqs[27], binFreqs[27],
		melTriangSums[25], specFeats[25], features[27], compfeats[29];

	float specMagMax, cepsMax,
		sumSq, fund, rah1, rah2, sum = 0;

	float voicedConf, voicedThresh, **pitchData, *pitchMeans, *pitchStdevs,
		meanPitch, stdevPitch, add, melStart, melEnd, melSpacing;

	double pvalsMax = 0;
	_declspec(align(16)) float snd[ARRSIZE], real[REALSIZE], imag[REALSIZE];

	//if need to clear database
	if (clear_database)
	{
		memset(database, 0, sizeof(database));
	}

	//get number of blocks
	startPos = 0;
	while (startPos + 0.01*SR * 50 + ARRSIZE < length)
	{ 
		++numBlocks;
		startPos += 25 * 0.01*SR;
	}
	
	//each block given 1 row, for 50 slices pitch data.
	startPos = 0;
	pitchData = (float **)malloc(numBlocks*sizeof(float));
	for (i = 0; i < numBlocks; i++){
		pitchData[i] = (float *)malloc(50 * sizeof(float));
	}

	//mean and sdv for each block over 50 entries
	pitchMeans = (float *)malloc(numBlocks*sizeof(float));
	pitchStdevs = (float *)malloc(numBlocks*sizeof(float));

	//for each block
	for (k = 0; k < numBlocks; k++)
	{
		voicedCount = 0;
		//each block has 50 entries
		for (i = 0; i < 50; i++)
		{
			cepsMax = 0; cepsMaxInd = 0; fund = 0; rah1 = 0; rah2 = 0; sumSq = 0;
			
			//each entry has 1024 data.
			for (j = 0; j < ARRSIZE; j++)
			{
				snd[j] = (0.54 - 0.46*cos(2.f*PI*j / (ARRSIZE - 1)))*soundData[startPos + j];
			}

			//do fft to the data
			realfft(snd, real, imag, LARRSIZE, 1);

			for (j = 0; j < REALSIZE; ++j)
			{
				//get the magnitude
				specMag[j] = sqrt(real[j] * real[j] + imag[j] * imag[j]);

				//doing the ceptrum 
				real[j] = log10(specMag[j]);
				imag[j] = 0;

				//doing filtering?
				if (j <= 500 / (SR / 2.f)*REALSIZE)
					real[j] *= 0.5*(1 - cos(2 * PI*j / (1000 / (SR / 2.f)*REALSIZE - 1)));
				if (j >= 3500 / (SR / 2.f)*REALSIZE && j <= 4000 / (SR / 2.f)*REALSIZE)
					real[j] *= 0.5*(1 - cos(2 * PI*(j - 3500 / (SR / 2.f)*REALSIZE + 500 / (SR / 2.f)*REALSIZE) / (1000 / (SR / 2.f)*REALSIZE - 1)));
				if (j > 4000 / (SR / 2.f)*REALSIZE)
					real[j] = 0;
			}
			//doing the ceptrum
			realfft(snd, real, imag, LARRSIZE, -1);

			//doing the ceptrum
			for (j = 0; j < REALSIZE; j++)
			{
				if (snd[j] < 0) snd[j] = 0;
				if (j<SR / 300.f || j>SR / 60.f) snd[j] = 0;
				if (snd[j]>cepsMax)
				{
					cepsMax = snd[j];
					cepsMaxInd = j;
				}
				sndSq[j] = snd[j] * snd[j]; //sound square
				sumSq += sndSq[j];
			}

			//getting the fundamental wave, h1 and h2
			for (j = -5; j <= 5; j++)
			{
				fund += sndSq[cepsMaxInd + j];  //fundamental wave
				if (2 * cepsMaxInd + 5 < REALSIZE)
				{
					rah1 += sndSq[2 * cepsMaxInd + j]; //harmonics 1
				}
				if (3 * cepsMaxInd + 5 < REALSIZE)
				{
					//rah1 += sndSq[3 * cepsMaxInd + j];   wrong?
					rah2 += sndSq[3 * cepsMaxInd + j]; //harmonics 2
				}
			}

			//get the voice confidence and threshold
			voicedConf = (fund + rah1 + rah2) / sumSq;
			voicedThresh = .2 / (SR / 60.f - SR / 300.f)*(SR / 300.f - cepsMaxInd) + .5;

			//decide if the entry is voiced or not
			if (voicedConf >= voicedThresh)
			{
				++voicedCount;
				pitchData[k][i] = 1 / ((float)cepsMaxInd / SR);
			}
			else pitchData[k][i] = 0;

			//the step for next entry is 10ms
			startPos += 0.01*SR;
		}
		//rewind for 250ms, there is a overlap between adjacent blocks
		startPos -= 25 * 0.01*SR;
	}

	//calculate the mean of pitchdata for each block
	for (i = 0; i < numBlocks; i++)
	{
		voicedCount = 0; sum = 0;
		for (j = 0; j < 50; j++)
		{
			if (pitchData[i][j] != 0)
			{
				sum += pitchData[i][j];
				++voicedCount;
			}
		}
		if (voicedCount != 0) pitchMeans[i] = sum / voicedCount;
		else pitchMeans[i] = 0;
	}

	//calculate the stdv of pitchdata for each block
	for (i = 0; i < numBlocks; ++i)
	{
		voicedCount = 0; sum = 0;
		for (j = 0; j < 50; ++j)
		{
			if (pitchData[i][j] != 0)
			{
				sum += (pitchData[i][j] - pitchMeans[i])*(pitchData[i][j] - pitchMeans[i]);
				++voicedCount;
			}
		}
		if (voicedCount != 0) pitchStdevs[i] = sqrt(sum / voicedCount);
		else pitchStdevs[i] = 0;
	}

	//get the mean of pitchdata over all blocks
	sum = 0;
	voicedCount = 0;
	for (i = 0; i < numBlocks; i++)
	{
		if (pitchMeans[i] != 0)
		{
			sum += pitchMeans[i];
			++voicedCount;
		}
	}
	if (voicedCount != 0) meanPitch = sum / voicedCount;
	else 
	{
		return 0;
	}

	//get the stdv of pitchdata over all blocks
	sum = 0; voicedCount = 0;
	for (i = 0; i < numBlocks; ++i)
	{
		if (pitchStdevs[i] != 0)
		{
			sum += pitchStdevs[i];
			++voicedCount;
		}
	}
	stdevPitch = sum / voicedCount;

	free(pitchData); free(pitchMeans); free(pitchStdevs);

	// find the 25 spectral features used in the features vector
	melStart = 2595 * log10(1 + 500 / 700.f);  //mel-frequency
	melEnd = 2595 * log10(1 + 4500 / 700.f);
	melSpacing = (melEnd - melStart) / 26.f;
	startPos = 0;
	voicedCount = 0;

	//target frequency
	for (i = 0; i < 27; i++)
	{
		melFreqs[i] = melStart + i*melSpacing;
		binFreqs[i] = 700 * (pow(10, melFreqs[i] / 2595.f) - 1);
	}

	//initialization
	for (i = 0; i < 25; i++) specFeats[i] = 0;

	//we use at most first 150 blocks (WHY ?)
	while (voicedCount < 150)
	{
		if (startPos + 0.01*SR + ARRSIZE >= length) break;

		for (i = 0; i < 25; i++) melTriangSums[i] = 0;
		specMagMax = 0; cepsMax = 0; cepsMaxInd = 0; fund = 0; rah1 = 0; rah2 = 0; sumSq = 0;
		snd[0] = soundData[startPos];    //the original dataset

		// build the wave data
		for (i = 1; i < ARRSIZE; i++)
			snd[i] = soundData[startPos + i] - 0.96*soundData[startPos + i - 1];

		for (i = 0; i < ARRSIZE; i++)
			snd[i] *= (0.54 - 0.46*cos(2.f*PI*i / (ARRSIZE - 1)));

		realfft(snd, real, imag, LARRSIZE, 1);

		//the same as before
		for (i = 0; i < REALSIZE; ++i)
		{
			specMag[i] = sqrt(real[i] * real[i] + imag[i] * imag[i]);
			if (specMag[i] > specMagMax)
				specMagMax = specMag[i];

			real[i] = log10(specMag[i]);
			imag[i] = 0;

			if (i <= 500 / (SR / 2.f)*REALSIZE)
				real[i] *= 0.5*(1 - cos(2 * PI*i / (1000 / (SR / 2.f)*REALSIZE - 1)));
			if (i >= 3500 / (SR / 2.f)*REALSIZE && i <= 4000 / (SR / 2.f)*REALSIZE)
				real[i] *= 0.5*(1 - cos(2 * PI*(i - 3500 / (SR / 2.f)*REALSIZE + 500 / (SR / 2.f)*REALSIZE) / (1000 / (SR / 2.f)*REALSIZE - 1)));
			if (i > 4000 / (SR / 2.f)*REALSIZE)
				real[i] = 0;
		}

		for (i = 0; i < REALSIZE; ++i)
			specMag[i] /= specMagMax;

		realfft(snd, real, imag, LARRSIZE, -1);


		for (i = 0; i < REALSIZE; ++i)
		{
			if (snd[i] < 0)
				snd[i] = 0;
			if (i < SR / 300.f || i > SR / 60.f)
				snd[i] = 0;
			if (snd[i] > cepsMax)
			{
				cepsMax = snd[i];
				cepsMaxInd = i;
			}
			sndSq[i] = snd[i] * snd[i];
			sumSq += sndSq[i];
		}

		for (i = -5; i <= 5; ++i)
		{
			fund += sndSq[cepsMaxInd + i];
			if (2 * cepsMaxInd + 5 < REALSIZE)
				rah1 += sndSq[2 * cepsMaxInd + i];
			if (3 * cepsMaxInd + 5 < REALSIZE)
				rah1 += sndSq[3 * cepsMaxInd + i];
		}

		voicedConf = (fund + rah1 + rah2) / sumSq;
		voicedThresh = .2 / (SR / 60.f - SR / 300.f)*(SR / 300.f - cepsMaxInd) + .5;

		if (voicedConf >= voicedThresh)
		{
			++voicedCount;
			for (i = 0; i < 25; i++){
				width = ((int)(binFreqs[i + 2] / (SR / 2.f)*REALSIZE)) - ((int)(binFreqs[i] / (SR / 2.f)*REALSIZE));

				for (j = 0; j < width; j++){
					add = (1 - fabs((j - (width - 1) / 2.f) / ((width + 1) / 2.f)))* //fabs: float abs
						specMag[j + ((int)(binFreqs[i] / (SR / 2.f)*REALSIZE))];
					melTriangSums[i] += add*add;
				}
			}

			for (i = 0; i < 25; i++) specFeats[i] += log(melTriangSums[i]);

		}
		startPos += 0.01*SR;
	}

	//can be modified
	if (voicedCount < 50) return 0;  //at least 0.75 seconds

	//the feature array
	for (i = 0; i < 25; i++) specFeats[i] /= voicedCount;
	features[0] = meanPitch;
	features[1] = stdevPitch;
	for (i = 0; i < 25; i++) features[i + 2] = specFeats[i];

	if (database[0][0] == 0)
	{
		++featCount[0];
		for (i = 0; i < 27; i++) database[0][i] = features[i];
		//update weight
		dbRecordWeight[0] = 100;
		cLabel[segPos] = 1;
		segPos++;
		return 1;

	}
	else
	{
		for (i = 0; i < MAXTALKERS; i++)
		{
			if (database[i][0] != 0)
			{
				//29 dimensions of comparision feature vector
				compfeats[0] = database[i][0];
				compfeats[1] = database[i][1];
				compfeats[2] = features[0];
				compfeats[3] = features[1];

				for (j = 0; j < 25; j++)
				{
					compfeats[j + 4] = database[i][j + 2] - features[j + 2];
				}
				//getting p-values for each of 6 subclasses
				pvals[i][0] = getpval(compfeats, mu_ffsame, S_inv_ffsame, consts_ffsame, weights_ffsame, 0);
				pvals[i][1] = getpval(compfeats, mu_mmsame, S_inv_mmsame, consts_mmsame, weights_mmsame, 0);
				pvals[i][2] = getpval(compfeats, mu_ffdiff, S_inv_ffdiff, consts_ffdiff, weights_ffdiff, 0);
				pvals[i][3] = getpval(compfeats, mu_mmdiff, S_inv_mmdiff, consts_mmdiff, weights_mmdiff, 0);
				pvals[i][4] = getpval(compfeats, mu_fmdiff, S_inv_fmdiff, consts_fmdiff, weights_fmdiff, 0);
				pvals[i][5] = getpval(compfeats, mu_mfdiff, S_inv_mfdiff, consts_mfdiff, weights_mfdiff, 0);

				for (int j = 1; j < 2; j++)
				{
					pvals[i][0] += getpval(compfeats, mu_ffsame, S_inv_ffsame, consts_ffsame, weights_ffsame, j);
					pvals[i][1] += getpval(compfeats, mu_mmsame, S_inv_mmsame, consts_mmsame, weights_mmsame, j);
					pvals[i][2] += getpval(compfeats, mu_ffdiff, S_inv_ffdiff, consts_ffdiff, weights_ffdiff, j);
					pvals[i][3] += getpval(compfeats, mu_mmdiff, S_inv_mmdiff, consts_mmdiff, weights_mmdiff, j);
					pvals[i][4] += getpval(compfeats, mu_fmdiff, S_inv_fmdiff, consts_fmdiff, weights_fmdiff, j);
					pvals[i][5] += getpval(compfeats, mu_mfdiff, S_inv_mfdiff, consts_mfdiff, weights_mfdiff, j);

				}
				// converting into percentages?
				sum = pvals[i][0] + pvals[i][1] + pvals[i][2] + pvals[i][3] + pvals[i][4] + pvals[i][5];

				for (j = 0; j < 5; ++j) pvals[i][j] /= sum / 100;

			}
			else
			{
				talkerPos = i;   //the new talker label should be what.
				break;
			}
		}
		for (i = 0; i < talkerPos; ++i)
		{
			if (pvals[i][0] + pvals[i][1] > pvalsMax)
			{
				pvalsMax = pvals[i][0] + pvals[i][1];
				pvalsMaxInd = i;
			}
		}
			
		//find the max pval
		for (i = 0; i < talkerPos; ++i)
		{
			if (pvals[i][0] + pvals[i][1] > pvalsMax)
			{
				pvalsMax = pvals[i][0] + pvals[i][1];
				pvalsMaxInd = i;
			}
		}
		//if it is classified as a old class
		if (pvalsMax > threshold)
		{
			++featCount[pvalsMaxInd];
			for (i = 0; i < 27; ++i)
				//update the record and weight
				//database[pvalsMaxInd][i] = ((float)(featCount[pvalsMaxInd] - 1) / featCount[pvalsMaxInd])*database[pvalsMaxInd][i] + features[i] / featCount[pvalsMaxInd];
				database[pvalsMaxInd][i] = ((float)(dbRecordWeight[pvalsMaxInd] / (dbRecordWeight[pvalsMaxInd] + pvalsMax)))*database[pvalsMaxInd][i] + features[i] * pvalsMax / (dbRecordWeight[pvalsMaxInd] + pvalsMax);

			dbRecordWeight[pvalsMaxInd] += pvalsMax;
			return pvalsMaxInd + 1;

		}
		//if it is a new class
		else if (pvalsMax < threshold)
		{
			++featCount[talkerPos];
			for (i = 0; i < 27; ++i)
				database[talkerPos][i] = features[i];
			dbRecordWeight[talkerPos] = 100;
			return talkerPos + 1;
		}		
	}
}

//given a voice segment, to decide the number of voiced frames
unsigned count_vframes(short soundData[], int length, int clear_database)
{
	int i, j, k = 0;
	int startPos, cepsMaxInd, voicedCount, numBlocks = 0;
	int width, pvalsMaxInd;

	float specMag[REALSIZE], sndSq[REALSIZE], melFreqs[27], binFreqs[27],
		melTriangSums[25];

	float specMagMax, cepsMax,
		sumSq, fund, rah1, rah2, sum = 0;

	float voicedConf, voicedThresh, **pitchData, *pitchMeans, *pitchStdevs,
		meanPitch, stdevPitch, add, melStart, melEnd, melSpacing;

	_declspec(align(16)) float snd[ARRSIZE], real[REALSIZE], imag[REALSIZE];

	//if need to clear database
	if (clear_database)
	{
		memset(database, 0, sizeof(database));
	}

	//get number of blocks
	startPos = 0;
	while (startPos + 0.01*SR * 50 + ARRSIZE < length)
	{
		++numBlocks;
		startPos += 25 * 0.01*SR;
	}


	//each block given 1 row, for 50 slices pitch data.
	startPos = 0;
	pitchData = (float **)malloc(numBlocks*sizeof(float));
	for (i = 0; i < numBlocks; i++){
		pitchData[i] = (float *)malloc(50 * sizeof(float));
	}

	//mean and sdv for each block over 50 entries
	pitchMeans = (float *)malloc(numBlocks*sizeof(float));
	pitchStdevs = (float *)malloc(numBlocks*sizeof(float));

	voicedCount = 0;
	//for each block
	for (k = 0; k < numBlocks; k++)
	{
		//voicedCount = 0;
		//each block has 50 entries
		for (i = 0; i < 50; i++)
		{
			cepsMax = 0; cepsMaxInd = 0; fund = 0; rah1 = 0; rah2 = 0; sumSq = 0;

			//each entry has 1024 data.
			for (j = 0; j < ARRSIZE; j++)
			{
				snd[j] = (0.54 - 0.46*cos(2.f*PI*j / (ARRSIZE - 1)))*soundData[startPos + j];
			}

			//do fft to the data
			realfft(snd, real, imag, LARRSIZE, 1);

			for (j = 0; j < REALSIZE; ++j)
			{
				//get the magnitude
				specMag[j] = sqrt(real[j] * real[j] + imag[j] * imag[j]);

				//doing the ceptrum 
				real[j] = log10(specMag[j]);
				imag[j] = 0;

				//doing filtering?
				if (j <= 500 / (SR / 2.f)*REALSIZE)
					real[j] *= 0.5*(1 - cos(2 * PI*j / (1000 / (SR / 2.f)*REALSIZE - 1)));
				if (j >= 3500 / (SR / 2.f)*REALSIZE && j <= 4000 / (SR / 2.f)*REALSIZE)
					real[j] *= 0.5*(1 - cos(2 * PI*(j - 3500 / (SR / 2.f)*REALSIZE + 500 / (SR / 2.f)*REALSIZE) / (1000 / (SR / 2.f)*REALSIZE - 1)));
				if (j > 4000 / (SR / 2.f)*REALSIZE)
					real[j] = 0;
			}
			//doing the ceptrum
			realfft(snd, real, imag, LARRSIZE, -1);

			//doing the ceptrum
			for (j = 0; j < REALSIZE; j++)
			{
				if (snd[j] < 0) snd[j] = 0;
				if (j<SR / 300.f || j>SR / 60.f) snd[j] = 0;
				if (snd[j]>cepsMax)
				{
					cepsMax = snd[j];
					cepsMaxInd = j;
				}
				sndSq[j] = snd[j] * snd[j]; //sound square
				sumSq += sndSq[j];
			}

			//getting the fundamental wave, h1 and h2
			for (j = -5; j <= 5; j++)
			{
				fund += sndSq[cepsMaxInd + j];  //fundamental wave
				if (2 * cepsMaxInd + 5 < REALSIZE)
				{
					rah1 += sndSq[2 * cepsMaxInd + j]; //harmonics 1
				}
				if (3 * cepsMaxInd + 5 < REALSIZE)
				{
					//rah1 += sndSq[3 * cepsMaxInd + j];   wrong?
					rah2 += sndSq[3 * cepsMaxInd + j]; //harmonics 2
				}
			}

			//get the voice confidence and threshold
			voicedConf = (fund + rah1 + rah2) / sumSq;
			voicedThresh = .2 / (SR / 60.f - SR / 300.f)*(SR / 300.f - cepsMaxInd) + .5;

			//decide if the entry is voiced or not
			if (voicedConf >= voicedThresh)
			{
				++voicedCount;
				pitchData[k][i] = 1 / ((float)cepsMaxInd / SR);
			}
			else pitchData[k][i] = 0;

			//the step for next entry is 10ms
			startPos += 0.01*SR;
		}
		//rewind for 250ms, there is a overlap between adjacent blocks
		startPos -= 25 * 0.01*SR;
	}

	return voicedCount;

}



//test the average number of frames given different length of segment.
void count_vframes_test(unsigned frame_count[20])
{

	wav_data wavData;
	FILE *wavFile;
	std::string Path[16][16];
	buildFilePath(Path);

	voice_seg segClosePDA[256];
	for (int i = 0; i < 256; i++)
	{
		segClosePDA[i].fileName = Path[i / 16][i % 16];
		segClosePDA[i].start = SR * 1;
		segClosePDA[i].end = SR * 15;
		segClosePDA[i].channel = 1;
		segClosePDA[i].vocName = vocNamePDA[i / 16];
	}

	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			//cout << i * 16 + j << endl;
			for (int time = 0; time < 20; time++)
			{
				for (int loop = 0; loop < 10; loop++)
				{
					wavFile = fopen(segClosePDA[i * 16 + j].fileName.c_str(), "r");
					unsigned seg_length = 0.1+0.5*time*SR;
					unsigned pos = inSegRandom(SR * 1, SR * 15, seg_length);

					short *dataSet = new short[seg_length];

					fseek(wavFile, 44 + pos, SEEK_SET);

					fread(dataSet, sizeof(short), seg_length, wavFile);
					fclose(wavFile);

					unsigned frame_num = count_vframes(dataSet, seg_length, 1);
					frame_count[time] += frame_num;
				}
			}
		}
	}
	for (int time = 0; time < 20; time++)
	{
		frame_count[time] /= (10*2*16);
		//cout << "time: " << 0.1 + 0.5*time << "secs, frame_num: " << frame_count[time] << endl;
		//system("pause");
	}
	    

}


//given the dataSet, the function returns the 27-dimension feature extracted from the data
void getFeat(short soundData[], unsigned length, float* features, unsigned &flag)
{

	int i, j, k = 0;
	int startPos, cepsMaxInd, voicedCount, numBlocks = 0;
	int width, pvalsMaxInd;

	float specMag[REALSIZE], sndSq[REALSIZE], melFreqs[27], binFreqs[27],
		melTriangSums[25], specFeats[25], compfeats[29];

	float specMagMax, cepsMax,
		sumSq, fund, rah1, rah2, sum = 0;

	float voicedConf, voicedThresh, **pitchData, *pitchMeans, *pitchStdevs,
		meanPitch, stdevPitch, add, melStart, melEnd, melSpacing;

	double pvalsMax = 0;
	_declspec(align(16)) float snd[ARRSIZE], real[REALSIZE], imag[REALSIZE];

	//get number of blocks
	startPos = 0;
	while (startPos + 0.01*SR * 50 + ARRSIZE < length)
	{
		++numBlocks;
		startPos += 25 * 0.01*SR;
	}

	//each block given 1 row, for 50 slices pitch data.
	startPos = 0;
	pitchData = (float **)malloc(numBlocks*sizeof(float));
	for (i = 0; i < numBlocks; i++){
		pitchData[i] = (float *)malloc(50 * sizeof(float));
	}

	//mean and sdv for each block over 50 entries
	pitchMeans = (float *)malloc(numBlocks*sizeof(float));
	pitchStdevs = (float *)malloc(numBlocks*sizeof(float));

	//for each block
	for (k = 0; k < numBlocks; k++)
	{
		voicedCount = 0;
		//each block has 50 entries
		for (i = 0; i < 50; i++)
		{
			cepsMax = 0; cepsMaxInd = 0; fund = 0; rah1 = 0; rah2 = 0; sumSq = 0;

			//each entry has 1024 data.
			for (j = 0; j < ARRSIZE; j++)
			{
				snd[j] = (0.54 - 0.46*cos(2.f*PI*j / (ARRSIZE - 1)))*soundData[startPos + j];
			}

			//do fft to the data
			realfft(snd, real, imag, LARRSIZE, 1);

			for (j = 0; j < REALSIZE; ++j)
			{
				//get the magnitude
				specMag[j] = sqrt(real[j] * real[j] + imag[j] * imag[j]);

				//doing the ceptrum 
				real[j] = log10(specMag[j]);
				imag[j] = 0;

				//doing filtering?
				if (j <= 500 / (SR / 2.f)*REALSIZE)
					real[j] *= 0.5*(1 - cos(2 * PI*j / (1000 / (SR / 2.f)*REALSIZE - 1)));
				if (j >= 3500 / (SR / 2.f)*REALSIZE && j <= 4000 / (SR / 2.f)*REALSIZE)
					real[j] *= 0.5*(1 - cos(2 * PI*(j - 3500 / (SR / 2.f)*REALSIZE + 500 / (SR / 2.f)*REALSIZE) / (1000 / (SR / 2.f)*REALSIZE - 1)));
				if (j > 4000 / (SR / 2.f)*REALSIZE)
					real[j] = 0;
			}
			//doing the ceptrum
			realfft(snd, real, imag, LARRSIZE, -1);

			//doing the ceptrum
			for (j = 0; j < REALSIZE; j++)
			{
				if (snd[j] < 0) snd[j] = 0;
				if (j<SR / 300.f || j>SR / 60.f) snd[j] = 0;
				if (snd[j]>cepsMax)
				{
					cepsMax = snd[j];
					cepsMaxInd = j;
				}
				sndSq[j] = snd[j] * snd[j]; //sound square
				sumSq += sndSq[j];
			}

			//getting the fundamental wave, h1 and h2
			for (j = -5; j <= 5; j++)
			{
				fund += sndSq[cepsMaxInd + j];  //fundamental wave
				if (2 * cepsMaxInd + 5 < REALSIZE)
				{
					rah1 += sndSq[2 * cepsMaxInd + j]; //harmonics 1
				}
				if (3 * cepsMaxInd + 5 < REALSIZE)
				{
					//rah1 += sndSq[3 * cepsMaxInd + j];   wrong?
					rah2 += sndSq[3 * cepsMaxInd + j]; //harmonics 2
				}
			}

			//get the voice confidence and threshold
			voicedConf = (fund + rah1 + rah2) / sumSq;
			voicedThresh = .2 / (SR / 60.f - SR / 300.f)*(SR / 300.f - cepsMaxInd) + .5;

			//decide if the entry is voiced or not
			if (voicedConf >= voicedThresh)
			{
				++voicedCount;
				pitchData[k][i] = 1 / ((float)cepsMaxInd / SR);
			}
			else pitchData[k][i] = 0;

			//the step for next entry is 10ms
			startPos += 0.01*SR;
		}

		//rewind for 250ms, there is a overlap between adjacent blocks
		startPos -= 25 * 0.01*SR;
	}

	if (voicedCount == 0)
	{
		flag = 1;
		return;
	}
		

	//calculate the mean of pitchdata for each block
	for (i = 0; i < numBlocks; i++)
	{
		voicedCount = 0; sum = 0;
		for (j = 0; j < 50; j++)
		{
			if (pitchData[i][j] != 0)
			{
				sum += pitchData[i][j];
				++voicedCount;
			}
		}
		if (voicedCount != 0) pitchMeans[i] = sum / voicedCount;
		else pitchMeans[i] = 0;
	}

	//calculate the stdv of pitchdata for each block
	for (i = 0; i < numBlocks; ++i)
	{
		voicedCount = 0; sum = 0;
		for (j = 0; j < 50; ++j)
		{
			if (pitchData[i][j] != 0)
			{
				sum += (pitchData[i][j] - pitchMeans[i])*(pitchData[i][j] - pitchMeans[i]);
				++voicedCount;
			}
		}
		if (voicedCount != 0) pitchStdevs[i] = sqrt(sum / voicedCount);
		else pitchStdevs[i] = 0;
	}

	//get the mean of pitchdata over all blocks
	sum = 0;
	voicedCount = 0;
	for (i = 0; i < numBlocks; i++)
	{
		if (pitchMeans[i] != 0)
		{
			sum += pitchMeans[i];
			++voicedCount;
		}
	}
	if (voicedCount != 0) meanPitch = sum / voicedCount;
	else
	{
		//printf("No data found in this segment.\n");
		return;
	}

	//get the stdv of pitchdata over all blocks
	sum = 0; voicedCount = 0;
	for (i = 0; i < numBlocks; ++i)
	{
		if (pitchStdevs[i] != 0)
		{
			sum += pitchStdevs[i];
			++voicedCount;
		}
	}
	stdevPitch = sum / voicedCount;

	free(pitchData); free(pitchMeans); free(pitchStdevs);

	// find the 25 spectral features used in the features vector
	melStart = 2595 * log10(1 + 500 / 700.f);  //mel-frequency
	melEnd = 2595 * log10(1 + 4500 / 700.f);
	melSpacing = (melEnd - melStart) / 26.f;
	startPos = 0;
	voicedCount = 0;

	//target frequency
	for (i = 0; i < 27; i++)
	{
		melFreqs[i] = melStart + i*melSpacing;
		binFreqs[i] = 700 * (pow(10, melFreqs[i] / 2595.f) - 1);
	}

	//initialization
	for (i = 0; i < 25; i++) specFeats[i] = 0;

	//we use at most first 150 blocks (WHY ?)
	//while (voicedCount < 150)
	while (true)
	{
		if (startPos + 0.01*SR + ARRSIZE >= length) break;

		for (i = 0; i < 25; i++) melTriangSums[i] = 0;
		specMagMax = 0; cepsMax = 0; cepsMaxInd = 0; fund = 0; rah1 = 0; rah2 = 0; sumSq = 0;
		snd[0] = soundData[startPos];    //the original dataset

		// build the wave data
		for (i = 1; i < ARRSIZE; i++)
			snd[i] = soundData[startPos + i] - 0.96*soundData[startPos + i - 1];

		for (i = 0; i < ARRSIZE; i++)
			snd[i] *= (0.54 - 0.46*cos(2.f*PI*i / (ARRSIZE - 1)));

		realfft(snd, real, imag, LARRSIZE, 1);

		//the same as before
		for (i = 0; i < REALSIZE; ++i)
		{
			specMag[i] = sqrt(real[i] * real[i] + imag[i] * imag[i]);
			if (specMag[i] > specMagMax)
				specMagMax = specMag[i];

			real[i] = log10(specMag[i]);
			imag[i] = 0;

			if (i <= 500 / (SR / 2.f)*REALSIZE)
				real[i] *= 0.5*(1 - cos(2 * PI*i / (1000 / (SR / 2.f)*REALSIZE - 1)));
			if (i >= 3500 / (SR / 2.f)*REALSIZE && i <= 4000 / (SR / 2.f)*REALSIZE)
				real[i] *= 0.5*(1 - cos(2 * PI*(i - 3500 / (SR / 2.f)*REALSIZE + 500 / (SR / 2.f)*REALSIZE) / (1000 / (SR / 2.f)*REALSIZE - 1)));
			if (i > 4000 / (SR / 2.f)*REALSIZE)
				real[i] = 0;
		}

		for (i = 0; i < REALSIZE; ++i)
			specMag[i] /= specMagMax;

		realfft(snd, real, imag, LARRSIZE, -1);


		for (i = 0; i < REALSIZE; ++i)
		{
			if (snd[i] < 0)
				snd[i] = 0;
			if (i < SR / 300.f || i > SR / 60.f)
				snd[i] = 0;
			if (snd[i] > cepsMax)
			{
				cepsMax = snd[i];
				cepsMaxInd = i;
			}
			sndSq[i] = snd[i] * snd[i];
			sumSq += sndSq[i];
		}

		for (i = -5; i <= 5; ++i)
		{
			fund += sndSq[cepsMaxInd + i];
			if (2 * cepsMaxInd + 5 < REALSIZE)
				rah1 += sndSq[2 * cepsMaxInd + i];
			if (3 * cepsMaxInd + 5 < REALSIZE)
				rah1 += sndSq[3 * cepsMaxInd + i];
		}

		voicedConf = (fund + rah1 + rah2) / sumSq;
		voicedThresh = .2 / (SR / 60.f - SR / 300.f)*(SR / 300.f - cepsMaxInd) + .5;

		if (voicedConf >= voicedThresh)
		{
			++voicedCount;
			for (i = 0; i < 25; i++){
				width = ((int)(binFreqs[i + 2] / (SR / 2.f)*REALSIZE)) - ((int)(binFreqs[i] / (SR / 2.f)*REALSIZE));

				for (j = 0; j < width; j++){
					add = (1 - fabs((j - (width - 1) / 2.f) / ((width + 1) / 2.f)))* //fabs: float abs
						specMag[j + ((int)(binFreqs[i] / (SR / 2.f)*REALSIZE))];
					melTriangSums[i] += add*add;
				}
			}

			for (i = 0; i < 25; i++) specFeats[i] += log(melTriangSums[i]);  //Z[m]
			//for (i = 0; i < 25; i++) specFeats[i] += melTriangSums[i];

		}
		startPos += 0.01*SR;
	}


	//the feature array
	for (i = 0; i < 25; i++) specFeats[i] /= voicedCount;      //take the average of all the voiced frames.
	features[0] = meanPitch;
	features[1] = stdevPitch;
	for (i = 0; i < 25; i++) features[i + 2] = specFeats[i];
	//for (i = 0; i < 25; i++) features[i + 2] = log(specFeats[i]);

	//return features;
}