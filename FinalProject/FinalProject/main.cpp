#define _CRT_SECURE_NO_DEPRECATE
#pragma once

#include <iostream>
#include <string>
#include "functions.h"
#include "model_test.h"

using namespace std;

int main()
{
	//voice_seg seg1;
	//voice_seg seg2;

	//model_test::init();

	//cout << model_test::eer_count[0] << endl;
	//system("pause");

	


	double pvalSet[MAXTALKERS][6];
	wav_data wavData;
	FILE *wavFile;
	std::string Path[16][16];
	buildFilePath(Path);

	voice_seg segClosePDA[256];
	for (unsigned i = 0; i < 256; i++)
	{
		segClosePDA[i].fileName = Path[i / 16][i % 16];
		segClosePDA[i].start = SR * 1;
		segClosePDA[i].end = SR * 8;
		segClosePDA[i].channel = 1;
		segClosePDA[i].vocName = vocNamePDA[i / 16];
	}
	string file_path = "D:\\Thesis\\FinalProject\\count_file_new\\";

	//*********************************ROC TESTS***************************************
	////the location to store the count data.
	//string file_path = "D:\\Thesis\\FinalProject\\count_file_new\\";
	//for (unsigned i = 1+112; i < 255; i += 16)
	//{
	//	for (unsigned j = 3 + i; j < 255; j += 16)
	//	{
	//		model_test::init();
	//		string roc_path = file_path + "roc_count_" + to_string(i / 16) + "_" + to_string(j / 16) + "_(" + to_string(i % 16) + "," + to_string(j % 16) + ")" + "_6g.txt";
	//		fstream roc_out(roc_path, ios::out);
	//		if (segClosePDA[i].vocName == segClosePDA[j].vocName)
	//		{
	//			model_test test(&segClosePDA[i], &segClosePDA[j], true);
	//			test.readTxtData();
	//			test.test_ROC();

	//			for (unsigned k = 0; k < 300; k++) roc_out << model_test::roc_count[k] << ' ';
	//			roc_out.close();
	//		}
	//		else
	//		{
	//			model_test test(&segClosePDA[i], &segClosePDA[j], false);
	//          test.readTxtData();
	//			test.test_ROC();

	//			for (unsigned k = 0; k < 300; k++) roc_out << model_test::roc_count[k] << ' ';
	//			roc_out.close();
	//		}
	//	}
	//}





	//*********************************ORIGINAL TESTS***************************************
	unsigned file_count = 0;
	for (unsigned i = 1+192; i < 255; i+=16)
	{
		for (unsigned j = 4+i; j < 256; j+=16)
		{
	        //unsigned file_count = 0;

	        model_test::init();
			//unsigned i = 33;
			//unsigned j = 7+i;
			string eer_path = file_path + "eer_count_" + to_string(i / 16) + "_" + to_string(j / 16) + "_(" + to_string(i % 16) + "," + to_string(j % 16) + ")" + ".txt";
			string del_path = file_path + "del_count_" + to_string(i / 16) + "_" + to_string(j / 16) + "_(" + to_string(i % 16) + "," + to_string(j % 16) + ")" + ".txt";
			string rej_path = file_path + "rej_count_" + to_string(i / 16) + "_" + to_string(j / 16) + "_(" + to_string(i % 16) + "," + to_string(j % 16) + ")" + ".txt";
			string det_path = file_path + "det_count_" + to_string(i / 16) + "_" + to_string(j / 16) + "_(" + to_string(i % 16) + "," + to_string(j % 16) + ")" + ".txt";

			fstream eer_out(eer_path, ios::out);
			fstream del_out(del_path, ios::out);
			fstream rej_out(rej_path, ios::out);
			fstream det_out(det_path, ios::out);

			//if 2 segments from the same talker
			if (segClosePDA[i].vocName == segClosePDA[j].vocName)
			{
				//set all the related array to 0.
				//model_test::init();

				model_test test(&segClosePDA[i], &segClosePDA[j], true);
				test.readTxtData();
				test.test_general();

				//system("pause");
				//cout << "eer_count: " << endl;
				for (unsigned k = 0; k < 4500; k++) eer_out << model_test::eer_count[k] << ' ';				
				eer_out.close();

				//cout << endl << "del_err: " << endl;
				for (unsigned k = 0; k < 40 * 15; k++) del_out << model_test::del_err[k] << ' ';
				del_out.close();

				//cout << endl << "mu_rej: " << endl;
				for (unsigned k = 0; k < 30; k++) rej_out << model_test::mu_rej[k] << ' ';
				rej_out.close();

				//cout << endl << "mu_det: " << endl;
				for (unsigned k = 0; k < 30; k++) det_out << model_test::mu_det[k] << ' ';
				det_out.close();

				//system("pause");
			}
			//if not from the same talker
			else
			{
				//set all the related array to 0.
				//model_test::init();

				model_test test(&segClosePDA[i], &segClosePDA[j], false);
				test.readTxtData();
				test.test_general();

				//cout << "eer_count: " << endl;
				for (unsigned k = 0; k < 100 * 15 * 3; k++) eer_out << model_test::eer_count[k] << ' ';
				eer_out.close();

				//cout << endl << "del_err: " << endl;
				for (unsigned k = 0; k < 40 * 15; k++) del_out << model_test::del_err[k] << ' ';
				del_out.close();

				//cout << endl << "mu_rej: " << endl;
				for (unsigned k = 0; k < 30; k++) rej_out << model_test::mu_rej[k] << ' ';
				rej_out.close();

				//cout << endl << "mu_det: " << endl;
				for (unsigned k = 0; k < 30; k++) det_out << model_test::mu_det[k] << ' ';
				det_out.close();

				//system("pause");
			}

		}
	}

	return 0;
}