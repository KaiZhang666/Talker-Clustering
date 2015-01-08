#ifndef model_test_h_
#define model_test_h_


#include "functions.h"
#include "PCA_matrix.h"
#include "struct.h"

class model_test
{
public:
	//constructor
	model_test(voice_seg* seg1, voice_seg* seg2, bool is_same) : segment1_(seg1), segment2_(seg2),
	           is_same_(is_same){}

	//destructor
	~model_test(){}

	//dim1: p(w|X) as threshold from 0 to 1, step 0.01
	//dim2: 15 test group, from (1,1) to (5,5)
	//dim3: miss_det_count(FN), and FP_count, and total count
	static unsigned eer_count[100 * 15 * 3];

	//dim1: p(w|X) as threshold from 0 to 1, step 0.01
	//dim2: TP, and FP_count, and total count
	static unsigned roc_count[100 * 3];

	//dim1: Ps[k] count, for 20 linear bins based on p(w|X)
	//dim2: q[k] count, for 20 linear bins based on p(w|X)
	static unsigned del_err[(20 + 20) * 15];

	//first column to add  p(w|X)<0.5, second to count n0
	static double mu_rej[2 * 15];

	//first column to add  p(w|X)>0.5, second to count n1
	static double mu_det[2 * 15];

	//to read model in the txt file into arrays.
	//you can rewrite it into a more convinient way.:)
	void readTxtData();

	//to start a new test, need to reinitialize all the static varibles
	static bool init();

	//given specific level of 2 segments, doing the test 
	//order represent the test number from (1,1) to (5,5), 15 in total
	void test_once(unsigned seg1_lv, unsigned seg2_lv, unsigned order);

	//doing ROC test among all combination of levels
	void test_ROC();

	//doing general test among all combination of levels
	void test_general();

	//given specific level of 2 segments, doing the test 
	//order represent the test number from (1,1) to (5,5), 15 in total
	void test_ROC_once(unsigned seg1_lv, unsigned seg2_lv, unsigned order);

	//given the range of the length, return a random length from the range.
	double len_random(double lbound, double ubound);

	//given the level, return the randomly chosen segment length
	unsigned get_length(unsigned seg_lv);

	//for each given threshold p(w|X), given 2 dataSets from the segments,
	//decide if it is FN or FP, update the relate static array list.
	void solve(short sData1[], short sData2[], unsigned length1,
		unsigned length2, unsigned order);

	//for each given threshold p(w|X), given 2 dataSets from the segments,
	//decide if it is FN or FP, update the relate static array list.
	void solve_PCA(short sData1[], short sData2[], unsigned length1,
		unsigned length2, unsigned order);

	//for each given threshold p(w|X), given 2 dataSets from the segments,
	//decide if it is FN or FP, update the relate static array list.
	void solve_others(short sData1[], short sData2[], unsigned length1,
		unsigned length2,unsigned order);

	//for each given threshold p(w|X), given 2 dataSets from the segments,
	//decide if it is FN or FP, update the relate static array list.
	void solve_others_PCA(short sData1[], short sData2[], unsigned length1,
		unsigned length2, unsigned order);

	//for each given threshold p(w|X), given 2 dataSets from the segments,
	//decide if it is FN or FP, update the relate static array list.
	void solve_roc(short sData1[], short sData2[], unsigned length1,
		unsigned length2);

private:
	voice_seg *segment1_;
	voice_seg *segment2_;
	bool is_same_;
};

#endif