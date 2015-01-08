#include "model_test.h"
#include <iostream>


//initialize all the static array lists
unsigned model_test::eer_count[100 * 15 * 3] = { 0 };

unsigned model_test::roc_count[100 * 3] = { 0 };

unsigned model_test::del_err[(20 + 20) * 15] = { 0 };

double model_test::mu_rej[2 * 15] = { 0 };

double model_test::mu_det[2 * 15] = { 0 };

//given the range of the length, return a random length from the range.
double model_test::len_random(double lbound, double ubound)
{
	LARGE_INTEGER nStartCounter;
	::QueryPerformanceCounter(&nStartCounter);
	::srand((unsigned)nStartCounter.LowPart);
	return (rand() / (1.0 * 0x7fff)*(ubound - lbound) + lbound);
}


//given the level, return the randomly chosen segment length
unsigned model_test::get_length(unsigned seg_lv)
{
	double seg_len = 0;
	unsigned len = 0;

	//given the level, decide the length of the segment1
	switch (seg_lv)
	{
	case 1:
		seg_len = this->len_random(0.8, 1.1);
		break;
	case 2:
		seg_len = this->len_random(1.1, 1.6);
		break;
	case 3:
		seg_len = this->len_random(1.6, 2.6);
		break;
	case 4:
		seg_len = this->len_random(2.6, 4.6);
		break;
	case 5:
		seg_len = this->len_random(4.6, 7.6);
		break;
	default:
		seg_len = 0;
	}

	len = seg_len*SR;
	return len;

}

//given specific level of 2 segments, doing the test 
//order represent the test number from (1,1) to (5,5), 15 in total
void model_test::test_once(unsigned seg1_lv, unsigned seg2_lv, unsigned order)
{

	//get the length of the 2 segments, length= time*sample rate
	unsigned len1 = get_length(seg1_lv);
	unsigned len2 = get_length(seg2_lv);

	//cout << "len1,2: " << len1 << ' ' << len2 << endl;
	//system("pause");


	//file pointer to the 2 segments
	FILE *wavFile1;
	FILE *wavFile2;

	//read the test data for both segments
	wavFile1 = fopen(segment1_->fileName.c_str(), "r");
	wavFile2 = fopen(segment2_->fileName.c_str(), "r");

	

	//for each given threshold p(w|X), given 2 dataSets from the segments,
	//decide if it is FN or FP, update the relate static array list.
	for (int i = 0; i < 50; i++)
	{
		//calculate the randomly chosen offset, given the test length of the segment
		unsigned offset1 = inSegRandom(segment1_->start, segment1_->end, len1);
		unsigned offset2 = inSegRandom(segment2_->start, segment2_->end, len2);


		//cout << "offset1,2: " << offset1 <<' '<< offset2 << endl;
		//system("pause");

		//two array list to store the read dataset
		short *dataSet1 = new short[len1];
		short *dataSet2 = new short[len2];


		fseek(wavFile1, 44 + sizeof(short)*offset1, SEEK_SET);
		fread(dataSet1, sizeof(short), len1, wavFile1);

		fseek(wavFile2, 44 + sizeof(short)*offset2, SEEK_SET);
		fread(dataSet2, sizeof(short), len2, wavFile2);
		if (i == 0)
		{
			//TO DO: can be put into once
				solve_PCA(dataSet1, dataSet2, len1, len2, order);
		}
		
		solve_others_PCA(dataSet1, dataSet2, len1, len2, order);

		delete[] dataSet1;
		delete[] dataSet2;
	}


	

}

//given specific level of 2 segments, doing the test 
//order represent the test number from (1,1) to (5,5), 15 in total
void model_test::test_ROC_once(unsigned seg1_lv, unsigned seg2_lv, unsigned order)
{

	//get the length of the 2 segments, length= time*sample rate
	unsigned len1 = get_length(seg1_lv);
	unsigned len2 = get_length(seg2_lv);

	//cout << "len1,2: " << len1 << ' ' << len2 << endl;
	//system("pause");


	//file pointer to the 2 segments
	FILE *wavFile1;
	FILE *wavFile2;

	//read the test data for both segments
	wavFile1 = fopen(segment1_->fileName.c_str(), "r");
	wavFile2 = fopen(segment2_->fileName.c_str(), "r");



	//for each given threshold p(w|X), given 2 dataSets from the segments,
	//decide if it is FN or FP, update the relate static array list.
	for (int i = 0; i < 10; i++)
	{
		//calculate the randomly chosen offset, given the test length of the segment
		unsigned offset1 = inSegRandom(segment1_->start, segment1_->end, len1);
		unsigned offset2 = inSegRandom(segment2_->start, segment2_->end, len2);


		//cout << "offset1,2: " << offset1 <<' '<< offset2 << endl;
		//system("pause");

		//two array list to store the read dataset
		short *dataSet1 = new short[len1];
		short *dataSet2 = new short[len2];


		fseek(wavFile1, 44 + sizeof(short)*offset1, SEEK_SET);
		fread(dataSet1, sizeof(short), len1, wavFile1);

		fseek(wavFile2, 44 + sizeof(short)*offset2, SEEK_SET);
		fread(dataSet2, sizeof(short), len2, wavFile2);

			//TO DO: can be put into once
		    solve_roc(dataSet1, dataSet2, len1, len2);

		delete[] dataSet1;
		delete[] dataSet2;
	}




}

//doing general test among all combination of levels
void model_test::test_general()
{
	int count = 0;

	for (int i = 1; i < 6; i++)
	{
		for (int j = i; j < 6; j++)
		{
			test_once(i, j, ++count);
			cout << count << endl;
		}
	}
}


//doing ROC test among all combination of levels
void model_test::test_ROC()
{
	int count = 0;

	for (int i = 1; i < 6; i++)
	{
		for (int j = i; j < 6; j++)
		{
			test_ROC_once(i, j, ++count);
			cout << count << endl;
		}
	}
}


//to start a new test, need to reinitialize all the static varibles
bool model_test::init()
{
	memset(roc_count, 0, sizeof(roc_count));
	memset(eer_count, 0, sizeof(eer_count));
	memset(del_err, 0, sizeof(del_err));
	memset(mu_rej, 0, sizeof(mu_rej));
	memset(mu_det, 0, sizeof(mu_det));
	return true;
}

//for each given threshold p(w|X), given 2 dataSets from the segments,
//decide if it is FN or FP, update the relate static array list.
void model_test::solve(short sData1[], short sData2[], unsigned length1,
	unsigned length2, unsigned order)
{
	//the array to store features
	//float *seg1_feat, *seg2_feat;
	float *seg1_feat = new float[27];
	float *seg2_feat = new float[27];

	float combFeat[29];

	//used to see if any of the segment is not initialized
	unsigned flag = 0;

	//calculate feature for each single segment
	getFeat(sData1, length1, seg1_feat,flag);
	getFeat(sData2, length2, seg2_feat,flag);

	if (flag == 1)
	{
		cout << "seg1 length: " << length1 << endl;
		cout << "seg2 length: " << length2 << endl;
		cout << "feat0: " << seg1_feat[0] << ' ' << seg2_feat[0] << endl;
		cout << "One of the segments is not initialized." << endl;
		return;
	}

	//cout << "feat[0] :"<<seg1_feat[0] << ' ' << seg2_feat[0] << endl;
	//system("pause");

	//rebuild the new combination feature
	combFeat[0] = seg1_feat[0];
	combFeat[1] = seg1_feat[1];
	combFeat[2] = seg2_feat[0];
	combFeat[3] = seg2_feat[1];

	for (int i = 0; i < 25; i++)
	{
		combFeat[i + 4] = seg1_feat[i + 2] - seg2_feat[i + 2];
	}

	//the array pvals to store the pvalue( p(w|X)) calculated from each of the 6 classes(model)
	double pvals[6];

	pvals[0] = getpval(combFeat, mu_ffsame, S_inv_ffsame, consts_ffsame, weights_ffsame, 0);
	pvals[1] = getpval(combFeat, mu_mmsame, S_inv_mmsame, consts_mmsame, weights_mmsame, 0);
	pvals[2] = getpval(combFeat, mu_ffdiff, S_inv_ffdiff, consts_ffdiff, weights_ffdiff, 0);
	pvals[3] = getpval(combFeat, mu_mmdiff, S_inv_mmdiff, consts_mmdiff, weights_mmdiff, 0);
	pvals[4] = getpval(combFeat, mu_fmdiff, S_inv_fmdiff, consts_fmdiff, weights_fmdiff, 0);
	pvals[5] = getpval(combFeat, mu_mfdiff, S_inv_mfdiff, consts_mfdiff, weights_mfdiff, 0);

	for (int j = 1; j < 2; j++)
	{
		pvals[0] += getpval(combFeat, mu_ffsame, S_inv_ffsame, consts_ffsame, weights_ffsame, j);
		pvals[1] += getpval(combFeat, mu_mmsame, S_inv_mmsame, consts_mmsame, weights_mmsame, j);
		pvals[2] += getpval(combFeat, mu_ffdiff, S_inv_ffdiff, consts_ffdiff, weights_ffdiff, j);
		pvals[3] += getpval(combFeat, mu_mmdiff, S_inv_mmdiff, consts_mmdiff, weights_mmdiff, j);
		pvals[4] += getpval(combFeat, mu_fmdiff, S_inv_fmdiff, consts_fmdiff, weights_fmdiff, j);
		pvals[5] += getpval(combFeat, mu_mfdiff, S_inv_mfdiff, consts_mfdiff, weights_mfdiff, j);
	}

	//for (unsigned i = 0; i < 6; i++) cout << pvals[i] << ' ';
	//cout << endl;
	//system("pause");


	double psum = 0;
	for (int i = 0; i < 6; i++) psum += pvals[i];

	double psame = (pvals[0] + pvals[1]) / psum;

	//cout << "psame: "<<psame << endl;
	//system("pause");

	//to judge the combFeat with the given threshold
	for (double thresh = 0; thresh < 1; thresh += 0.01)
	{
		//the offset in the eer_count array
		unsigned eer_offset = (thresh / 0.01) * 3 + (order - 1) * 100 * 3;

		//cout << "eer_offset: "<<eer_offset << endl;
		//system("pause");

		//if it is a false positive
		if (psame>thresh && this->is_same_ == false) eer_count[eer_offset]++;


		//cout << "psame: "<<psame << endl;
		//cout << "thresh: " << thresh << endl;
		//system("pause");

		//if it is a false negative
		if (psame <= thresh && this->is_same_ == true)	eer_count[eer_offset + 1]++;

		//the total count ++
		eer_count[eer_offset + 2]++;
	}
	
		

	//cout <<"is_same: "<< this->is_same_ << endl;
	//system("pause");

	//deal with del_err, mu_det, mu_rej
	//if (psame < 0.5)
	//{
	//	mu_rej[0 + (order - 1)*2] += psame;
	//	mu_rej[1 + (order - 1)*2]++;
	//}
	//else
	//{
	//	mu_det[0 + (order - 1) * 2] += psame;
	//	mu_det[1 + (order - 1) * 2]++;
	//}

	//count the datapoint dropped into each of the bin
	//unsigned derr_offset = (psame >= 0.025) ? ((psame - 0.025) / 0.05) : 0;

	//cout << "psame: " << psame << endl;
	//cout << "derr_offset: "<<derr_offset << endl;
	//system("pause");

	//count the Ps[k]
	//if (this->is_same_ == true) del_err[derr_offset + (order - 1) * 40]++;

	//count the q[k]
	//del_err[derr_offset + 20 + (order - 1) * 40]++;


}

//for each given threshold p(w|X), given 2 dataSets from the segments,
//decide if it is FN or FP, update the relate static array list.
void model_test::solve_PCA(short sData1[], short sData2[], unsigned length1,
	unsigned length2, unsigned order)
{
	float *seg1_feat = new float[27];
	float *seg2_feat = new float[27];

	float combFeat[25];
	float transFeat[25];
	float PCA_Feat[15];


	//used to see if any of the segment is not initialized
	unsigned flag = 0;

	//calculate feature for each single segment
	getFeat(sData1, length1, seg1_feat, flag);
	getFeat(sData2, length2, seg2_feat, flag);

	if (flag == 1)
	{
		cout << "seg1 length: " << length1 << endl;
		cout << "seg2 length: " << length2 << endl;
		cout << "feat0: " << seg1_feat[0] << ' ' << seg2_feat[0] << endl;
		cout << "One of the segments is not initialized." << endl;
		return;
	}

	//cout << "feat[0] :"<<seg1_feat[0] << ' ' << seg2_feat[0] << endl;
	//system("pause");

	//rebuild the new combination feature
	//combFeat[0] = seg1_feat[0];
	//combFeat[1] = seg1_feat[1];
	//combFeat[2] = seg2_feat[0];
	//combFeat[3] = seg2_feat[1];

	for (int i = 0; i < 25; i++)
	{
		combFeat[i] = seg1_feat[i + 2] - seg2_feat[i + 2];
		combFeat[i] -= PCA_mean_final[i];
	}
	//transform into the new coordinate
	memset(transFeat, 0, sizeof(transFeat));
	for (int j = 0; j < 25; j++)
	{
		for (int k = 0; k < 25; k++)
			transFeat[j] += combFeat[k] * PCA_coef_final[k][j];
	}

	//only use the first 11 dimensions
	PCA_Feat[0] = seg1_feat[0];
	PCA_Feat[1] = seg1_feat[1];
	PCA_Feat[2] = seg2_feat[0];
	PCA_Feat[3] = seg2_feat[1];
	for (int j = 4; j < 15; j++)
		PCA_Feat[j] = transFeat[j - 4];


	//the array pvals to store the pvalue( p(w|X)) calculated from each of the 6 classes(model)
	double pvals[6];

	pvals[0] = getpval_PCA_3g(PCA_Feat, mu_6g[0], S_inv_3g[0], consts_3g[0], weights_3g[0], 0);
	pvals[1] = getpval_PCA_3g(PCA_Feat, mu_3g[1], S_inv_3g[1], consts_3g[1], weights_3g[1], 0);
	pvals[2] = getpval_PCA_3g(PCA_Feat, mu_3g[2], S_inv_3g[2], consts_3g[2], weights_3g[2], 0);
	pvals[3] = getpval_PCA_3g(PCA_Feat, mu_3g[3], S_inv_3g[3], consts_3g[3], weights_3g[3], 0);
	pvals[4] = getpval_PCA_3g(PCA_Feat, mu_3g[4], S_inv_3g[4], consts_3g[4], weights_3g[4], 0);
	pvals[5] = getpval_PCA_3g(PCA_Feat, mu_3g[5], S_inv_3g[5], consts_3g[5], weights_3g[5], 0);

	for (int j = 1; j < 3; j++)
	{
		pvals[0] += getpval_PCA_3g(PCA_Feat, mu_3g[0], S_inv_3g[0], consts_3g[0], weights_3g[0], j);
		pvals[1] += getpval_PCA_3g(PCA_Feat, mu_3g[1], S_inv_3g[1], consts_3g[1], weights_3g[1], j);
		pvals[2] += getpval_PCA_3g(PCA_Feat, mu_3g[2], S_inv_3g[2], consts_3g[2], weights_3g[2], j);
		pvals[3] += getpval_PCA_3g(PCA_Feat, mu_3g[3], S_inv_3g[3], consts_3g[3], weights_3g[3], j);
		pvals[4] += getpval_PCA_3g(PCA_Feat, mu_3g[4], S_inv_3g[4], consts_3g[4], weights_3g[3], j);
		pvals[5] += getpval_PCA_3g(PCA_Feat, mu_3g[5], S_inv_3g[5], consts_3g[5], weights_3g[4], j);
	}
	//for (unsigned i = 0; i < 6; i++) cout << pvals[i] << ' ';
	//cout << endl;
	//system("pause");


	double psum = 0;
	for (int i = 0; i < 6; i++) psum += pvals[i];

	double psame = (pvals[0] + pvals[1]) / psum;

	//cout << "psame: "<<psame << endl;
	//system("pause");

	//to judge the combFeat with the given threshold
	for (double thresh = 0; thresh < 1; thresh += 0.01)
	{
		//the offset in the eer_count array
		unsigned eer_offset = (thresh / 0.01) * 3 + (order - 1) * 100 * 3;

		//cout << "eer_offset: "<<eer_offset << endl;
		//system("pause");

		//if it is a false positive
		if (psame>thresh && this->is_same_ == false) eer_count[eer_offset]++;


		//cout << "psame: "<<psame << endl;
		//cout << "thresh: " << thresh << endl;
		//system("pause");

		//if it is a false negative
		if (psame <= thresh && this->is_same_ == true)	eer_count[eer_offset + 1]++;

		//the total count ++
		eer_count[eer_offset + 2]++;
	}



	//cout <<"is_same: "<< this->is_same_ << endl;
	//system("pause");

	//deal with del_err, mu_det, mu_rej
	//if (psame < 0.5)
	//{
	//	mu_rej[0 + (order - 1)*2] += psame;
	//	mu_rej[1 + (order - 1)*2]++;
	//}
	//else
	//{
	//	mu_det[0 + (order - 1) * 2] += psame;
	//	mu_det[1 + (order - 1) * 2]++;
	//}

	//count the datapoint dropped into each of the bin
	//unsigned derr_offset = (psame >= 0.025) ? ((psame - 0.025) / 0.05) : 0;

	//cout << "psame: " << psame << endl;
	//cout << "derr_offset: "<<derr_offset << endl;
	//system("pause");

	//count the Ps[k]
	//if (this->is_same_ == true) del_err[derr_offset + (order - 1) * 40]++;

	//count the q[k]
	//del_err[derr_offset + 20 + (order - 1) * 40]++;


}

void model_test::solve_others(short sData1[], short sData2[], unsigned length1,
	unsigned length2, unsigned order)
{
	//the array to store features
	//float *seg1_feat, *seg2_feat;
	float *seg1_feat = new float[27];
	float *seg2_feat = new float[27];

	float combFeat[29];

	//used to see if any of the segment is not initialized
	unsigned flag = 0;

	//calculate feature for each single segment
	getFeat(sData1, length1, seg1_feat, flag);
	getFeat(sData2, length2, seg2_feat, flag);

	if (flag == 1)
	{
		cout << "seg1 length: " << length1 << endl;
		cout << "seg2 length: " << length2 << endl;
		cout << "feat0: " << seg1_feat[0] << ' ' << seg2_feat[0] << endl;
		cout << "One of the segments is not initialized." << endl;
		return;
	}

	//cout << "feat[0] :"<<seg1_feat[0] << ' ' << seg2_feat[0] << endl;
	//system("pause");

	//rebuild the new combination feature
	combFeat[0] = seg1_feat[0];
	combFeat[1] = seg1_feat[1];
	combFeat[2] = seg2_feat[0];
	combFeat[3] = seg2_feat[1];

	for (int i = 0; i < 25; i++)
	{
		combFeat[i + 4] = seg1_feat[i + 2] - seg2_feat[i + 2];
	}

	//the array pvals to store the pvalue( p(w|X)) calculated from each of the 6 classes(model)
	double pvals[6];

	pvals[0] = getpval(combFeat, mu_ffsame, S_inv_ffsame, consts_ffsame, weights_ffsame, 0);
	pvals[1] = getpval(combFeat, mu_mmsame, S_inv_mmsame, consts_mmsame, weights_mmsame, 0);
	pvals[2] = getpval(combFeat, mu_ffdiff, S_inv_ffdiff, consts_ffdiff, weights_ffdiff, 0);
	pvals[3] = getpval(combFeat, mu_mmdiff, S_inv_mmdiff, consts_mmdiff, weights_mmdiff, 0);
	pvals[4] = getpval(combFeat, mu_fmdiff, S_inv_fmdiff, consts_fmdiff, weights_fmdiff, 0);
	pvals[5] = getpval(combFeat, mu_mfdiff, S_inv_mfdiff, consts_mfdiff, weights_mfdiff, 0);

	for (int j = 1; j < 2; j++)
	{
		pvals[0] += getpval(combFeat, mu_ffsame, S_inv_ffsame, consts_ffsame, weights_ffsame, j);
		pvals[1] += getpval(combFeat, mu_mmsame, S_inv_mmsame, consts_mmsame, weights_mmsame, j);
		pvals[2] += getpval(combFeat, mu_ffdiff, S_inv_ffdiff, consts_ffdiff, weights_ffdiff, j);
		pvals[3] += getpval(combFeat, mu_mmdiff, S_inv_mmdiff, consts_mmdiff, weights_mmdiff, j);
		pvals[4] += getpval(combFeat, mu_fmdiff, S_inv_fmdiff, consts_fmdiff, weights_fmdiff, j);
		pvals[5] += getpval(combFeat, mu_mfdiff, S_inv_mfdiff, consts_mfdiff, weights_mfdiff, j);
	}

	//for (unsigned i = 0; i < 6; i++) cout << pvals[i] << ' ';
	//cout << endl;
	//system("pause");


	double psum = 0;
	for (int i = 0; i < 6; i++) psum += pvals[i];

	double psame = (pvals[0] + pvals[1]) / psum;

	//cout << "psame: "<<psame << endl;
	//system("pause");

	//to judge the combFeat with the given threshold

	//the offset in the eer_count array
	//unsigned eer_offset = (thresh / 0.01) * 3 + (order - 1) * 100 * 3;

	//cout << "eer_offset: "<<eer_offset << endl;
	//system("pause");

	//if it is a false positive
	//if (psame>thresh && this->is_same_ == false) eer_count[eer_offset]++;


	//cout << "psame: "<<psame << endl;
	//cout << "thresh: " << thresh << endl;
	//system("pause");

	//if it is a false negative
	//if (psame <= thresh && this->is_same_ == true)	eer_count[eer_offset + 1]++;





	//the total count ++
	//eer_count[eer_offset + 2]++;


	//cout <<"is_same: "<< this->is_same_ << endl;
	//system("pause");

	//deal with del_err, mu_det, mu_rej
	if (psame < 0.5)
	{
		mu_rej[0 + (order - 1) * 2] += psame;
		mu_rej[1 + (order - 1) * 2]++;
	}
	else
	{
		mu_det[0 + (order - 1) * 2] += psame;
		mu_det[1 + (order - 1) * 2]++;
	}

	//count the datapoint dropped into each of the bin
	unsigned derr_offset = (psame >= 0.025) ? ((psame - 0.025) / 0.05) : 0;

	//cout << "psame: " << psame << endl;
	//cout << "derr_offset: "<<derr_offset << endl;
	//system("pause");

	//count the Ps[k]
	if (this->is_same_ == true) del_err[derr_offset + (order - 1) * 40]++;

	//count the q[k]
	del_err[derr_offset + 20 + (order - 1) * 40]++;


}


void model_test::solve_others_PCA(short sData1[], short sData2[], unsigned length1,
	unsigned length2, unsigned order)
{
	//the array to store features
	//float *seg1_feat, *seg2_feat;
	float *seg1_feat = new float[27];
	float *seg2_feat = new float[27];

	float combFeat[25];
	float transFeat[25];
	float PCA_Feat[15];


	//used to see if any of the segment is not initialized
	unsigned flag = 0;

	//calculate feature for each single segment
	getFeat(sData1, length1, seg1_feat, flag);
	getFeat(sData2, length2, seg2_feat, flag);

	if (flag == 1)
	{
		cout << "seg1 length: " << length1 << endl;
		cout << "seg2 length: " << length2 << endl;
		cout << "feat0: " << seg1_feat[0] << ' ' << seg2_feat[0] << endl;
		cout << "One of the segments is not initialized." << endl;
		return;
	}

	//cout << "feat[0] :"<<seg1_feat[0] << ' ' << seg2_feat[0] << endl;
	//system("pause");

	//rebuild the new combination feature
	//combFeat[0] = seg1_feat[0];
	//combFeat[1] = seg1_feat[1];
	//combFeat[2] = seg2_feat[0];
	//combFeat[3] = seg2_feat[1];

	for (int i = 0; i < 25; i++)
	{
		combFeat[i] = seg1_feat[i + 2] - seg2_feat[i + 2];
		combFeat[i] -= PCA_mean_final[i];
	}
	//transform into the new coordinate
	memset(transFeat, 0, sizeof(transFeat));
	for (int j = 0; j < 25; j++)
	{
		for (int k = 0; k < 25; k++)
			transFeat[j] += combFeat[k] * PCA_coef_final[k][j];
	}

	//only use the first 11 dimensions
	PCA_Feat[0] = seg1_feat[0];
	PCA_Feat[1] = seg1_feat[1];
	PCA_Feat[2] = seg2_feat[0];
	PCA_Feat[3] = seg2_feat[1];
	for (int j = 4; j < 15; j++)
		PCA_Feat[j] = transFeat[j - 4];


	//the array pvals to store the pvalue( p(w|X)) calculated from each of the 6 classes(model)
	double pvals[6];

	pvals[0] = getpval_PCA_3g(PCA_Feat, mu_6g[0], S_inv_3g[0], consts_3g[0], weights_3g[0], 0);
	pvals[1] = getpval_PCA_3g(PCA_Feat, mu_3g[1], S_inv_3g[1], consts_3g[1], weights_3g[1], 0);
	pvals[2] = getpval_PCA_3g(PCA_Feat, mu_3g[2], S_inv_3g[2], consts_3g[2], weights_3g[2], 0);
	pvals[3] = getpval_PCA_3g(PCA_Feat, mu_3g[3], S_inv_3g[3], consts_3g[3], weights_3g[3], 0);
	pvals[4] = getpval_PCA_3g(PCA_Feat, mu_3g[4], S_inv_3g[4], consts_3g[4], weights_3g[4], 0);
	pvals[5] = getpval_PCA_3g(PCA_Feat, mu_3g[5], S_inv_3g[5], consts_3g[5], weights_3g[5], 0);

	for (int j = 1; j < 3; j++)
	{
		pvals[0] += getpval_PCA_3g(PCA_Feat, mu_3g[0], S_inv_3g[0], consts_3g[0], weights_3g[0], j);
		pvals[1] += getpval_PCA_3g(PCA_Feat, mu_3g[1], S_inv_3g[1], consts_3g[1], weights_3g[1], j);
		pvals[2] += getpval_PCA_3g(PCA_Feat, mu_3g[2], S_inv_3g[2], consts_3g[2], weights_3g[2], j);
		pvals[3] += getpval_PCA_3g(PCA_Feat, mu_3g[3], S_inv_3g[3], consts_3g[3], weights_3g[3], j);
		pvals[4] += getpval_PCA_3g(PCA_Feat, mu_3g[4], S_inv_3g[4], consts_3g[4], weights_3g[3], j);
		pvals[5] += getpval_PCA_3g(PCA_Feat, mu_3g[5], S_inv_3g[5], consts_3g[5], weights_3g[4], j);
	}

	//for (unsigned i = 0; i < 6; i++) cout << pvals[i] << ' ';
	//cout << endl;
	//system("pause");


	double psum = 0;
	for (int i = 0; i < 6; i++) psum += pvals[i];

	double psame = (pvals[0] + pvals[1]) / psum;

	//cout << "psame: "<<psame << endl;
	//system("pause");

	//to judge the combFeat with the given threshold

	//the offset in the eer_count array
	//unsigned eer_offset = (thresh / 0.01) * 3 + (order - 1) * 100 * 3;

	//cout << "eer_offset: "<<eer_offset << endl;
	//system("pause");

	//if it is a false positive
	//if (psame>thresh && this->is_same_ == false) eer_count[eer_offset]++;


	//cout << "psame: "<<psame << endl;
	//cout << "thresh: " << thresh << endl;
	//system("pause");

	//if it is a false negative
	//if (psame <= thresh && this->is_same_ == true)	eer_count[eer_offset + 1]++;





	//the total count ++
	//eer_count[eer_offset + 2]++;


	//cout <<"is_same: "<< this->is_same_ << endl;
	//system("pause");

	//deal with del_err, mu_det, mu_rej
	if (psame < 0.5)
	{
		mu_rej[0 + (order - 1) * 2] += psame;
		mu_rej[1 + (order - 1) * 2]++;
	}
	else
	{
		mu_det[0 + (order - 1) * 2] += psame;
		mu_det[1 + (order - 1) * 2]++;
	}

	//count the datapoint dropped into each of the bin
	unsigned derr_offset = (psame >= 0.025) ? ((psame - 0.025) / 0.05) : 0;

	//cout << "psame: " << psame << endl;
	//cout << "derr_offset: "<<derr_offset << endl;
	//system("pause");

	//count the Ps[k]
	if (this->is_same_ == true) del_err[derr_offset + (order - 1) * 40]++;

	//count the q[k]
	del_err[derr_offset + 20 + (order - 1) * 40]++;


}

//for each given threshold p(w|X), given 2 dataSets from the segments,
//decide if it is FN or FP, update the relate static array list.
void model_test::solve_roc(short sData1[], short sData2[], unsigned length1,
	unsigned length2)
{
	//the array to store features
	//float *seg1_feat, *seg2_feat;
	float *seg1_feat = new float[27];
	float *seg2_feat = new float[27];

	float combFeat[25];
	float transFeat[25];
	float PCA_Feat[15];


	//used to see if any of the segment is not initialized
	unsigned flag = 0;

	//calculate feature for each single segment
	getFeat(sData1, length1, seg1_feat, flag);
	getFeat(sData2, length2, seg2_feat, flag);

	if (flag == 1)
	{
		cout << "seg1 length: " << length1 << endl;
		cout << "seg2 length: " << length2 << endl;
		cout << "feat0: " << seg1_feat[0] << ' ' << seg2_feat[0] << endl;
		cout << "One of the segments is not initialized." << endl;
		return;
	}

	//cout << "feat[0] :"<<seg1_feat[0] << ' ' << seg2_feat[0] << endl;
	//system("pause");

	//rebuild the new combination feature
	//combFeat[0] = seg1_feat[0];
	//combFeat[1] = seg1_feat[1];
	//combFeat[2] = seg2_feat[0];
	//combFeat[3] = seg2_feat[1];

	for (int i = 0; i < 25; i++)
	{
		combFeat[i] = seg1_feat[i + 2] - seg2_feat[i + 2];
		combFeat[i] -= PCA_mean_final[i];
	}
	//transform into the new coordinate
	memset(transFeat, 0, sizeof(transFeat));
	for (int j = 0; j < 25; j++)
	{
		for (int k = 0; k < 25; k++)
			transFeat[j] += combFeat[k] * PCA_coef_final[k][j];
	}

	//only use the first 11 dimensions
	PCA_Feat[0] = seg1_feat[0];
	PCA_Feat[1] = seg1_feat[1];
	PCA_Feat[2] = seg2_feat[0];
	PCA_Feat[3] = seg2_feat[1];
	for (int j = 4; j < 15; j++)
		PCA_Feat[j] = transFeat[j - 4];
	

	//the array pvals to store the pvalue( p(w|X)) calculated from each of the 6 classes(model)
	double pvals[6];

	pvals[0] = getpval_PCA_3g(PCA_Feat, mu_6g[0], S_inv_6g[0], consts_6g[0], weights_6g[0], 0);
	pvals[1] = getpval_PCA_3g(PCA_Feat, mu_6g[1], S_inv_6g[1], consts_6g[1], weights_6g[1], 0);
	pvals[2] = getpval_PCA_3g(PCA_Feat, mu_6g[2], S_inv_6g[2], consts_6g[2], weights_6g[2], 0);
	pvals[3] = getpval_PCA_3g(PCA_Feat, mu_6g[3], S_inv_6g[3], consts_6g[3], weights_6g[3], 0);
	pvals[4] = getpval_PCA_3g(PCA_Feat, mu_6g[4], S_inv_6g[4], consts_6g[4], weights_6g[4], 0);
	pvals[5] = getpval_PCA_3g(PCA_Feat, mu_6g[5], S_inv_6g[5], consts_6g[5], weights_6g[5], 0);

	for (int j = 1; j < 6; j++)
	{
		pvals[0] += getpval_PCA_3g(PCA_Feat, mu_6g[0], S_inv_6g[0], consts_6g[0], weights_6g[0], j);
		pvals[1] += getpval_PCA_3g(PCA_Feat, mu_6g[1], S_inv_6g[1], consts_6g[1], weights_6g[1], j);
		pvals[2] += getpval_PCA_3g(PCA_Feat, mu_6g[2], S_inv_6g[2], consts_6g[2], weights_6g[2], j);
		pvals[3] += getpval_PCA_3g(PCA_Feat, mu_6g[3], S_inv_6g[3], consts_6g[3], weights_6g[3], j);
		pvals[4] += getpval_PCA_3g(PCA_Feat, mu_6g[4], S_inv_6g[4], consts_6g[4], weights_6g[3], j);
		pvals[5] += getpval_PCA_3g(PCA_Feat, mu_6g[5], S_inv_6g[5], consts_6g[5], weights_6g[4], j);
	}

	//for (unsigned i = 0; i < 6; i++) cout << pvals[i] << ' ';
	//cout << endl;
	//system("pause");


	double psum = 0;
	for (int i = 0; i < 6; i++) psum += pvals[i];

	double psame = (pvals[0] + pvals[1]) / psum;

	//cout << "psame: "<<psame << endl;
	//system("pause");

	//to judge the combFeat with the given threshold
	for (double thresh = 0; thresh < 1; thresh += 0.01)
	{
		//the offset in the eer_count array
		unsigned roc_offset = (thresh / 0.01) * 3;

		//if it is a false positive
		if (psame>thresh && this->is_same_ == false) roc_count[roc_offset + 1]++;

		//cout << "psame: "<<psame << endl;
		//cout << "thresh: " << thresh << endl;
		//system("pause");

		//if it is a true positive
		if (psame > thresh && this->is_same_ == true)	roc_count[roc_offset]++;

		//the total count ++
		roc_count[roc_offset + 2]++;
	}

}

//to read model in the txt file into arrays.
//you can rewrite it into a more convinient way.:)
void model_test::readTxtData()
{
	fstream in;

	//read consts data from txt for 6 models for different numbers of peaks
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m2_consts.txt");
		for (int j = 0; j < 2; j++) in >> consts_2g[i][j];
		in.close();

	}
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m3_consts.txt");
		for (int j = 0; j < 2; j++)	in >> consts_3g[i][j];
		in.close();
	}
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m4_consts.txt");
		for (int j = 0; j < 2; j++) in >> consts_4g[i][j];
		in.close();
	}
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m5_consts.txt");
		for (int j = 0; j < 2; j++) in >> consts_5g[i][j];
		in.close();
	}
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m6_consts.txt");
		for (int j = 0; j < 2; j++) in >> consts_6g[i][j];
		in.close();
	}

	//read weights data from txt for 6 models for different numbers of peaks
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m2_pi.txt");
		for (int j = 0; j < 2; j++) in >> weights_2g[i][j];
		in.close();
	}
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m3_pi.txt");
		for (int j = 0; j < 2; j++)	in >> weights_3g[i][j];
		in.close();
	}
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m4_pi.txt");
		for (int j = 0; j < 2; j++) in >> weights_4g[i][j];
		in.close();
	}
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m5_pi.txt");
		for (int j = 0; j < 2; j++) in >> weights_5g[i][j];
		in.close();
	}
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m6_pi.txt");
		for (int j = 0; j < 2; j++) in >> weights_6g[i][j];
		in.close();
	}

	//read mu data from txt for 6 models for different numbers of peaks
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m2_miu.txt");
		for (int j = 0; j < 15; j++)
		for (int k = 0; k < 2; k++) in >> mu_2g[i][k][j];
		in.close();
	}

	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m3_miu.txt");
		for (int j = 0; j < 15; j++)
		for (int k = 0; k < 3; k++) in >> mu_3g[i][k][j];
		in.close();
	}
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m4_miu.txt");
		for (int j = 0; j < 15; j++)
		for (int k = 0; k < 4; k++) in >> mu_4g[i][k][j];
		in.close();
	}
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m5_miu.txt");
		for (int j = 0; j < 15; j++)
		for (int k = 0; k < 5; k++) in >> mu_5g[i][k][j];
		in.close();
	}
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m6_miu.txt");
		for (int j = 0; j < 15; j++)
		for (int k = 0; k < 6; k++) in >> mu_6g[i][k][j];
		in.close();
	}


	//read inverse sigma data from txt for 6 models for different numbers of peaks
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m2_sigma.txt");
		for (int l = 0; l < 2; l++)
		for (int j = 0; j < 15; j++)
		for (int k = 0; k < 15; k++) in >> S_inv_2g[i][l][k][j];
		in.close();
	}
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m3_sigma.txt");
		for (int l = 0; l < 3; l++)
		for (int j = 0; j < 15; j++)
		for (int k = 0; k < 15; k++) in >> S_inv_3g[i][l][k][j];
		in.close();
	}
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m4_sigma.txt");
		for (int l = 0; l < 4; l++)
		for (int j = 0; j < 15; j++)
		for (int k = 0; k < 15; k++) in >> S_inv_4g[i][l][k][j];
		in.close();
	}
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m5_sigma.txt");
		for (int l = 0; l < 5; l++)
		for (int j = 0; j < 15; j++)
		for (int k = 0; k < 15; k++) in >> S_inv_5g[i][l][k][j];
		in.close();
	}
	for (int i = 0; i < 6; i++)
	{
		in.open("D:\\Thesis\\TIMIT_PCA_features\\PCA_feature_final\\train_result\\txt\\" + model_name[i] + "_out_m6_sigma.txt");
		for (int l = 0; l < 6; l++)
		for (int j = 0; j < 15; j++)
		for (int k = 0; k < 15; k++) in >> S_inv_6g[i][l][k][j];
		in.close();
	}


}