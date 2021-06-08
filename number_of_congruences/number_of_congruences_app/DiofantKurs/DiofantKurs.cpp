// DiofantKurs.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <src\arageli\arageli.hpp>
#include <cmath>
#include <fstream>
#include <omp.h>
#include <string>
#include <thread>
#include <future>

#pragma comment (lib, "arageli64r.lib")
using namespace std;
using namespace Arageli;


class answer{
public:
	Arageli::vector<big_int> resG;
	Arageli::vector<big_int> resB;
	big_int numberOfClasses = 0; // количество классов с решений
};

bool CheckFullSquare(big_int numb);				// проверка числа на полный квадрат
bool CheckSolutionExist(big_int D, big_int N);	// проверка на наличие решения

string DivWithRem(big_int numerator, big_int denom);

big_int powmod(big_int a, big_int p);
answer MainSolution(big_int D, big_int N, double& countComparison);

void ModuloComparision(big_int start, big_int end, big_int modd, big_int QQ0, Arageli::vector<big_int>* p0Vec);


int main()
{
	setlocale(0, "RUS");

	cout << "	|	Программа считает отношение суммарного количества решений сравнения P^2=D(N) для N < Nmax	|\n	|		к суммарному количеству примитивных решений уравения x^2-Dy^2=N для N < Nmax.		|" << endl;
	cout << "	|_______________________________________________________________________________________________________|" << endl;

	int isPrime = 0;
	cout << "Считать для простых(1) или для всех(2)?" << endl;
	cin >> isPrime;
	while (isPrime != 1 && isPrime != 2) {
		cout << "Введите (1) или (2): ";
		cin >> isPrime;
	};

	big_int D;
	big_int Dmin;
	big_int Dmax;
	big_int Z;
	cout << "Введите D_min:	";
	cin >> Dmin;
	D = Dmin;
	cout << "Введите D_max:	";
	cin >> Dmax;
	cout << "Введите максимальное N:	";
	cin >> Z;
	//double threshold = 999;
	/*double threshold = 0.0;
	cout << "Введите минимальную границу среднего примитивных решений\n (0 - отдельно записывать данные для всех D) :	";
	cin >> threshold;*/

	int threshold_n = Z;

	int count1_33 = 0;
	int count2 = 0;
	int count4 = 0;

	double min1_33 = 1.25;
	double max1_33 = 1.45;
	double min2 = 1.75;
	double max2 = 2.25;
	double min4 = 3.5;
	double max4 = 4.5;

	short int isChangeData = 0;
	while (isChangeData != 1 && isChangeData != 2) {
		cout << "Хотите ввести границы округления самостоятельно? (1) - Да. (2) - Нет: ";
		cin >> isChangeData;
	};

	if (isChangeData == 1) {
		cout << "	Нижняя граница для 1.33:	"; cin >> min1_33;
		cout << "	Верхняя граница для 1.33:	"; cin >> max1_33;
		cout << "	Нижняя граница для 2:		"; cin >> min2;
		cout << "	Верхняя граница для 2:		"; cin >> max2;
		cout << "	Нижняя граница для 4:		"; cin >> min4;
		cout << "	Верхняя граница для 4:		"; cin >> max4;
	}

	// cout << "введите n, при котором проверять эту границу:	";
	// cin >> threshold_n;

	big_int begin = 1;
	int isNegativeNumber = 1;
	cout << "Учитывать отрицательные N?\n (1) - Нет: 1 <= N < Nmax\n (2) - Да: -Nmax < N < Nmax" << endl;
	do {
		cout << "Введите (1) или (2): ";
		cin >> isNegativeNumber;
		
		if (isNegativeNumber == 2) {
			begin = -Z;
		}
	} while (isNegativeNumber != 1 && isNegativeNumber != 2);

	int isOnlyPrim = 1;
	cout << "Искать только примитивные решения?\n (1) - Да.\n (2) - нет." << endl;
	do {
		cout << "Введите (1) или (2): ";
		cin >> isOnlyPrim;
	} while (isOnlyPrim != 1 && isOnlyPrim != 2);


	stringstream ss;
	ss<< Dmin;
	string strD = ss.str();
	ss.str("");
	ss << Dmax;
	string strDmax = ss.str();
	ss.str("");
	ss << Z;
	string strN = ss.str();

	ofstream outFile_N_means; // файл, куда записываются значение N и среднее число классов (всех подряд. записывается при N == threshold_n)
	ofstream outFile_countComparison; // файл, куда записываются значение N, число сравнени	и число примитивных классов
	outFile_N_means.open("D_Means_values" + strD + "-" + strDmax + "N" + strN + ".txt");
	if (isOnlyPrim == 1) {
		outFile_N_means << "Знач-е_D	Среднее_прим-х	ОтношениеЧислоСравнений/числоПримРешений	ПослеОкругления" << endl;
	}
	else {
		outFile_N_means << "Знач-е_D	Среднее_прим-х	Среднее_всех	ОтношениеЧислоСравнений/числоПримРешений	ПослеОкругления" << endl;
	}

	switch (isPrime)
	{
	case 1:
		cout << "Считаем для простых..." << endl;
		cout << " СЧПримР	-	Среднее число прим-х решений\n СЧВсехР	-	Среднее число всех решений\n ОтнЧСр/ЧПримР	-	Отношение: ЧислоСравнений/числоПримРешений" << endl;
		while (D < Dmax) {

			cout << "********************* D = " << D << endl;
			//ofstream file2;
			//ofstream file3;

			stringstream ss;
			ss << D;
			string strD = ss.str();


			outFile_countComparison.open("D_countComp" + strD + "_" + "N" + strN + ".txt");
			outFile_countComparison << "Знач-е_D	Знач-е_N	ЧислоПримРеш.	ЧислоСравнений" << endl;

			//file2.open(strD + "_NumberOfClasses.txt");
			//string f2_str = (strD + "_NumberOfClasses.txt");
			//if (isOnlyPrim == 1) {
			//	file2 << "Зн-е_D	Знач_N	Число_примит-х	ОтношениеЧислоСравнений/числоПримРешений" << endl;
			//}
			//else {
			//	file2 << "Зн-е_D	Знач_N	Число_примит-х	ЧислоВсехКлассов	ОтношениеЧислоСравнений/числоПримРешений" << endl;
			//}
			////file2 << "Зн-е_D	Знач_N	Число_примит-х	ЧислоВсехКлассов	ОтношениеЧислоСравнений/числоПримРешений" << endl;
			//file3.open(strD + "_MeanNumberOfDecisions.txt");
			//string f3_str = (strD + "_MeanNumberOfDecisions.txt");
			//if (isOnlyPrim == 1) {
			//	file3 << "D	Значение N	Средняя_число_прим-х_решений" << endl;
			//}
			//else {
			//	file3 << "D	Значение N	Средняя_число_прим-х_решений	Среднее_число_всех_решений" << endl;
			//}
			big_int countN = 1;	// всего уравнений прошло
			big_int countWS = 0; // число уравнений с решениями

			if (CheckFullSquare(D)) {
				cout << "Число является полным квадратом!" << endl;
				continue;
			}

			////////// для x^2 - Dy^2 = N /////////
			/*big_int begin = 1;*/

			cout << "Значение N	СЧПримР			СЧВсехР			ОтнЧСр/ЧПримР" << endl;
			big_int primAnsSum = 0;
			big_int allAnsSum = 0;

			double countComparison = 0; // сумма количетсва решений сравнения P^2 = D (N)
			double sumCountComparison = 0; // сумма количетсва решений сравнения P^2 = D (N)
			double totalSumPrimSolutions = 0;	// сумма числа примитивных решений
			double ratioComparisonONPrimSum = 0;	// отношение  countComparison / totalSumPrimSolutions
			double countComparison2 = 0;
			double fullRoundRatio = 0;		// округление до 1.33, 2, 4

			string str_meanPrimCount = "";
			string str_meanAllCount = "";

			// информация для подсчета доли числа уравнений с решениями по отношению ко всем
			//cout << "Введите уже подсчитанное число уравнений с решениями. Пиши значение предыдущего шага (если с начала, то - 0): ";
			//cin >> countWS;
			//cout << "Введите сколько всего уравнений прошло (если с начала, то - 1): ";
			//cin >> countN;
			for (big_int N = begin; N <= Z; N++) {

				if (N == 0) {
					continue;
				}

				answer ans;
				answer resAns;
				big_int numberAllClasses = 0;	// количество всех классов с решениями

				if (isOnlyPrim == 2) {
					for (big_int n = intsqrt(Arageli::abs(N)); n > 1; n--) {		// чтобы искали максимальный квадрат сразу
						if (N % square(n) == 0) {
							big_int newN = N / square(n);

							ans = MainSolution(D, newN, countComparison2);
							for (big_int i = 0; i < ans.resG.size(); i++) {
								resAns.resG.push_back(ans.resG[i] * n);
								resAns.resB.push_back(ans.resB[i] * n);
							}
							//resAns.numberOfClasses = ans.numberOfClasses;
							numberAllClasses += ans.numberOfClasses;
							//break;
						}
					}
				}
				ans = MainSolution(D, N, countComparison);
				sumCountComparison += countComparison;
				for (big_int i = 0; i < ans.resG.size(); i++) {
					resAns.resG.push_back(ans.resG[i]);
					resAns.resB.push_back(ans.resB[i]);
				}
				resAns.numberOfClasses = ans.numberOfClasses;
				numberAllClasses += ans.numberOfClasses;
				if (resAns.resG.size() == 0) {
					//cout << "Нет решений!" << endl;
				}
				else {
					countWS++;
				}
				if (numberAllClasses != 0) {

					totalSumPrimSolutions += resAns.numberOfClasses;
					if (sumCountComparison == 0) {
						ratioComparisonONPrimSum = 0;
					}
					else {
						ratioComparisonONPrimSum = sumCountComparison / totalSumPrimSolutions;
					}

					outFile_countComparison << D << "	" << N << "	" << resAns.numberOfClasses << "	" << countComparison << endl;
					/*if (isOnlyPrim == 1) {
						file2 << D << "	" << N << "	" << resAns.numberOfClasses << "	" << ratioComparisonONPrimSum << endl;
					}
					else {
						file2 << D << "	" << N << "	" << resAns.numberOfClasses << "	" << numberAllClasses << "	" << ratioComparisonONPrimSum << endl;
					}*/

					primAnsSum += resAns.numberOfClasses;
					allAnsSum += numberAllClasses;

					str_meanPrimCount = DivWithRem(primAnsSum, countN);
					str_meanAllCount = DivWithRem(allAnsSum, countN);

					/*if (isOnlyPrim == 1) {
						file3 << D << "	" << N << "	" << str_meanPrimCount << endl;
					}
					else {
						file3 << D << "	" << N << "	" << str_meanPrimCount << "	" << str_meanAllCount << endl;
					}*/

					if (N % 100 == 0) {
						cout << N << "	" << setw(20) << str_meanPrimCount << "	" << setw(20) << str_meanAllCount << "		" << round(ratioComparisonONPrimSum * 100) / 100 << endl;
					}
				}
				if (N == threshold_n) {
					//if (N == 100) {

					double mean = 0.0;
					std::stringstream ss(str_meanPrimCount);
					ss >> mean;

					double roundRatio = round(ratioComparisonONPrimSum * 100) / 100;
					fullRoundRatio = roundRatio;
					if (roundRatio >= min1_33 && roundRatio <= max1_33) {
						fullRoundRatio = 1.33;
						count1_33++;
					}
					else if (roundRatio >= min2 && roundRatio <= max2) {
						fullRoundRatio = 2;
						count2++;
					}
					else if (roundRatio >= min4 && roundRatio <= max4) {
						fullRoundRatio = 4;
						count4++;
					}
					if (isOnlyPrim == 1) {
						outFile_N_means << D << "	" << mean << "		" << roundRatio << "		"  << fullRoundRatio << endl;
					} else {
						outFile_N_means << D << "	" << mean << "		" << str_meanAllCount << "		" << roundRatio << "		" << fullRoundRatio << endl;
					}
					//outFile_N_means << D << "	" << mean << "		" << str_meanAllCount << "		" << round(ratioComparisonONPrimSum * 100) / 100 << endl;

					//double mean = atof(str_meanPrimCount.c_str());
					//if (mean < threshold) {
					//	file2.close();
					//	file3.close();
					//	remove(f2_str.c_str());
					//	remove(f3_str.c_str());
					//	break;
					//}
				}
				countN++;
			}
			//file2.close();
			//file3.close();

			D = next_prime(D);
			outFile_countComparison.close();
		}
		break;
	case 2:
		cout << "Считаем для всех..." << endl;
		cout << " СЧПримР	-	Среднее число прим-х решений\n СЧВсехР	-	Среднее число всех решений\n ОтнЧСр/ЧПримР	-	Отношение: ЧислоСравнений/числоПримРешений" << endl;
		while (D < Dmax) {
			cout << "********************* D = " << D << endl;
			if (CheckFullSquare(D)) {
				cout << "Число является полным квадратом!" << endl;
				D++;
				continue;
			}

			//ofstream file2;
			//ofstream file3;

			stringstream ss;
			ss << D;
			string strD = ss.str();

			//file2.open(strD + "_NumberOfClasses.txt");
			//string f2_str = (strD + "_NumberOfClasses.txt");
			//if (isOnlyPrim == 1) {
			//	file2 << "Зн-е_D	Знач_N	Число_примит-х	ОтношениеЧислоСравнений/числоПримРешений" << endl;
			//}
			//else {
			//	file2 << "Зн-е_D	Знач_N	Число_примит-х	ЧислоВсехКлассов	ОтношениеЧислоСравнений/числоПримРешений" << endl;
			//}
			////file2 << "Зн-е_D	Знач_N	Число_примит-х	ЧислоВсехКлассов	ОтношениеЧислоСравнений/числоПримРешений" << endl;
			//file3.open(strD + "_MeanNumberOfDecisions.txt");
			//string f3_str = (strD + "_MeanNumberOfDecisions.txt");
			//if (isOnlyPrim == 1) {
			//	file3 << "D	Значение N	Средняя_число_прим-х_решений" << endl;
			//}
			//else {
			//	file3 << "D	Значение N	Средняя_число_прим-х_решений	Среднее_число_всех_решений" << endl;
			//}

			big_int countN = 1;	// всего уравнений прошло
			big_int countWS = 0; // число уравнений с решениями


			////////// для x^2 - Dy^2 = N /////////
			/*big_int begin = 1;*/
			cout << "Значение N	СЧПримР			СЧВсехР			ОтнЧСр/ЧПримР" << endl;
			big_int primAnsSum = 0;
			big_int allAnsSum = 0;
			double countComparison = 0; // количетсво решений сравнения P^2 = D (N)
			double sumCountComparison = 0; // сумма количетсва решений сравнения P^2 = D (N)
			double totalSumPrimSolutions = 0;	// сумма числа примитивных решений
			double ratioComparisonONPrimSum = 0;	// отношение  countComparison / totalSumPrimSolutions
			double countComparison2 = 0;
			double fullRoundRatio = 0;		// округление до 1.33, 2, 4

			string str_meanPrimCount = "";
			string str_meanAllCount = "";

			// информация для подсчета доли числа уравнений с решениями по отношению ко всем
			//cout << "Введите уже подсчитанное число уравнений с решениями. Пиши значение предыдущего шага (если с начала, то - 0): ";
			//cin >> countWS;
			//cout << "Введите сколько всего уравнений прошло (если с начала, то - 1): ";
			//cin >> countN;
			for (big_int N = begin; N <= Z; N++) {
				if (N == 0) {
					continue;
				}

				answer ans;
				answer resAns;
				big_int numberAllClasses = 0;	// количество всех классов с решениями

				if (isOnlyPrim == 2) {
					for (big_int n = intsqrt(Arageli::abs(N)); n > 1; n--) {		// чтобы искали максимальный квадрат сразу
						if (N % square(n) == 0) {
							big_int newN = N / square(n);

							ans = MainSolution(D, newN, countComparison2);
							// cout << newN << " countComparison:	" << countComparison << endl;
							for (big_int i = 0; i < ans.resG.size(); i++) {
								resAns.resG.push_back(ans.resG[i] * n);
								resAns.resB.push_back(ans.resB[i] * n);
							}
							//resAns.numberOfClasses = ans.numberOfClasses;
							numberAllClasses += ans.numberOfClasses;
							//break;
						}
					}
				}
				ans = MainSolution(D, N, countComparison);
				sumCountComparison += countComparison;
				// cout << N << " _countComparison:	" << countComparison << endl;
				for (big_int i = 0; i < ans.resG.size(); i++) {
					resAns.resG.push_back(ans.resG[i]);
					resAns.resB.push_back(ans.resB[i]);
				}
				resAns.numberOfClasses = ans.numberOfClasses;
				numberAllClasses += ans.numberOfClasses;
				if (resAns.resG.size() == 0) {
					//cout << "Нет решений!" << endl;
				}
				else {
					countWS++;
				}
				if (numberAllClasses != 0) {

					totalSumPrimSolutions += resAns.numberOfClasses;
					// cout << "totalSumPrimSolutions: " << totalSumPrimSolutions << endl;
					if (sumCountComparison == 0) {
						ratioComparisonONPrimSum = 0;
					} else {
						ratioComparisonONPrimSum = sumCountComparison / totalSumPrimSolutions;
					}
					// cout << "ratioComparisonONPrimSum: " << ratioComparisonONPrimSum << endl << endl;

					/*if (isOnlyPrim == 1) {
						file2 << D << "	" << N << "	" << resAns.numberOfClasses << "	" << ratioComparisonONPrimSum << endl;
					}
					else {
						file2 << D << "	" << N << "	" << resAns.numberOfClasses << "	" << numberAllClasses << "	" << ratioComparisonONPrimSum << endl;
					}*/

					primAnsSum += resAns.numberOfClasses;
					allAnsSum += numberAllClasses;


					str_meanPrimCount = DivWithRem(primAnsSum, countN);
					str_meanAllCount = DivWithRem(allAnsSum, countN);

					/*if (isOnlyPrim == 1) {
						file3 << D << "	" << N << "	" << str_meanPrimCount << endl;
					}
					else {
						file3 << D << "	" << N << "	" << str_meanPrimCount << "	" << str_meanAllCount << endl;
					}*/
					if (N % 100 == 0) {
						cout  << N << "	" << setw(20) << str_meanPrimCount << "	" << setw(20) << str_meanAllCount << "		" << round(ratioComparisonONPrimSum * 100) / 100  << endl;
					}
				}
				if (N == threshold_n) {
					double mean = 0.0;
					std::stringstream ss(str_meanPrimCount);
					ss >> mean;

					double roundRatio = round(ratioComparisonONPrimSum * 100) / 100;
					fullRoundRatio = roundRatio;
					if (roundRatio >= min1_33 && roundRatio <= max1_33) {
						fullRoundRatio = 1.33;
						count1_33++;
					}
					else if (roundRatio >= min2 && roundRatio <= max2) {
						fullRoundRatio = 2;
						count2++;
					}
					else if (roundRatio >= min4 && roundRatio <= max4) {
						fullRoundRatio = 4;
						count4++;
					}

					if (isOnlyPrim == 1) {
						outFile_N_means << D << "	" << mean << "		" << roundRatio << "	" << fullRoundRatio << endl;
					}
					else {
						outFile_N_means << D << "	" << mean << "		" << str_meanAllCount << "		" << roundRatio << "	" << fullRoundRatio << endl;
					}
					//outFile_N_means << D << "	" << mean << "		" << str_meanAllCount  << "		" << round(ratioComparisonONPrimSum * 100) / 100 << endl;

					//double mean = atof(str_meanPrimCount.c_str());
					//if (mean < threshold) {
					//	file2.close();
					//	file3.close();
					//	remove(f2_str.c_str());
					//	remove(f3_str.c_str());
					//	break;
					//}
				}
				countN++;
			}
			//file2.close();
			//file3.close();

			D++;
		}
		break;
	default:
		return 0;
		break;
	}


	outFile_N_means << endl << "На промежутки D = " << Dmin << " - " << Dmax << " после округления получилось:" << endl;
	outFile_N_means << "Отношение	количество" << endl;
	outFile_N_means << "	1.33	" << count1_33 << endl;
	outFile_N_means << "	2	" << count2 << endl;
	outFile_N_means << "	4	" << count4 << endl;

	outFile_N_means << endl << "Данные при введенных границах:" << endl;
	outFile_N_means << min1_33 << " <= 1.33 <= " << max1_33 << endl;
	outFile_N_means << min2 << " <= 2 <= " << max2 << endl;
	outFile_N_means << min4 << " <= 4 <= " << max4 << endl;

	outFile_N_means.close();

	system("pause");
    return 0;
}

//	проверка на полный квадрат
bool CheckFullSquare(big_int numb) {
	big_int d = intsqrt(numb);

	if (square(d) == numb) {
		//cout << "Число является полным квадратом!" << endl;
		return true;
	}
	return false;
}

// проверка на отсутсвие решения
bool CheckSolutionExist(big_int D, big_int N) {
	// 1.
	if (N == -3) {	// если N = -3, а D = (d^2+2), то нет решения
		if (CheckFullSquare(D - 2) == true) {
			//cout << "!Нет решения!\n[Код проверки: 1.1]\n" << endl;
			return false;
		}
	}
	// 2.
	if (D == pow(N, 2) - 1) {
		if (CheckFullSquare(N) == true) {
			return true;
		}
		else {
			//cout << "!Нет решения!\n[Код проверки: 1.2]\n" << endl;
			return false;
		}
	}
	return true;
}

// Выполняет бинарное возведение в квадрат по модулю a ^ 2 (mod p)
big_int powmod(big_int a, big_int p)
{
	big_int b = 2;
	big_int res = 1;
	while (b)
		if (b & (big_int)1)
			res = big_int(res * 1ll * a % p), --b;
		else
			a = big_int(a * 1ll * a % p), b >>= 1;
	return res;
}

// деление двух целых с вещественным результатом (но в строке)
string DivWithRem(big_int numerator, big_int denom) {
	string res;
	if ((numerator < 0 && denom >= 0)||(numerator >= 0 && denom < 0)) {
		res += "-";
		numerator = Arageli::abs(numerator);
		denom = Arageli::abs(denom);
	}

	big_int integ = numerator / denom;	// целая часть
	big_int remain = numerator % denom;	// остаток
	big_int temp = remain;

	stringstream ss;
	ss << integ;
	res += ss.str();

	if (remain != 0) {
		res += ".";
	}
	ss.str("");

	int count = 1;
	while (remain != 0) {
		temp = remain * 10;
		remain = temp % denom;
		temp = temp / denom;

		ss << temp;
		res += ss.str();
		ss.str("");
		count++;
		if (count > 6) {
			break;
		}
	}
	return res;
}

void ModuloComparision(big_int start, big_int end, big_int modd, big_int QQ0, Arageli::vector<big_int>* p0Vec) {
	for (big_int i = start; i < end; i++) {
		if (powmod(i, QQ0) == modd) {
			if (i != 0) {
				p0Vec->push_back(-i);
			}
			p0Vec->push_back(i);
		}
	}
}

answer MainSolution(big_int D, big_int N, double& countComparison) {
	answer res;
	// int countComparison = 0;

	Arageli::vector<big_int> p0Vec;

	big_int l = 0;
	big_int QQ0 = N;

	// Проверка 1.
	if (CheckSolutionExist(D, N) == false) {
		//cout << "!Нет решения!\n[Код проверки: 1]\n" << endl;
		return res;
	}
	big_int modd = Arageli::mod(D, QQ0);
	big_int end = Arageli::abs(QQ0) / 2;

	ModuloComparision(0, end, modd, QQ0, &p0Vec);

	big_int i = Arageli::abs(QQ0) / 2;	// дополнение чтобы левая граница была строгая, а правая нестрогая ( -Q0/2 < x <= Q02/ )
	if (Arageli::mod(square(i), QQ0) == modd) {	//	если N = нечетное, то левая граница тоже входит
		if (QQ0 % 2 != 0) {				
			if (-i != i) {
				p0Vec.push_back(-i);
			}
		}
		p0Vec.push_back(i);
	}
	
	// Проверка 2.
	if (p0Vec.is_empty()) {
		//cout << "Нет подходящих P0. Нет решения.\n[Код проверки: 2]\n";
		return res;
	}

	countComparison = 0;
	for (big_int i : p0Vec) {
		//cout << "ВЫБРАНО P0 = " << i << endl;
		//big_int PP0 = i;
		countComparison++;

		Arageli::vector<big_int> vecA;
		Arageli::vector<big_int> vecB;
		Arageli::vector<big_int> vecP;
		Arageli::vector<big_int> vecQ;

		big_int QQi = QQ0;
		big_int PPi = i;

		big_int qqi = 0;// , qq0 = 0;
		big_int int_part_sqrtD = 0;

		big_int AAi2 = 0;
		big_int AAi1 = 1;
		big_int BBi2 = 1;
		big_int BBi1 = 0;

		big_int Gt1 = 0;		// ответ при m = 0;
		big_int Bt1 = 0;

		big_int begin_P;
		big_int begin_Q;

		big_int AAi = 0, BBi = 0, GGi = 0;
		big_int count = 1;
		big_int start_index = 0;
		big_int t = 0;
		bool period_begin = false;

		while (true) {

			if (QQi > 0) {
				qqi = prquot((PPi + intsqrt(D)), QQi);
			}
			else if (QQi < 0) {
				big_int p;
				big_int q;
				prdivide((PPi + 1 + intsqrt(D)), QQi, p, q);
				if (q == 0) {
					qqi = p;
				}
				else {
					qqi = p - 1;
				}
			}

			int_part_sqrtD = intsqrt(D);

			PPi = qqi * QQi - PPi;
			QQi = (D - square(PPi)) / QQi;

			AAi = qqi*AAi1 + AAi2;
			BBi = qqi*BBi1 + BBi2;

			vecA.push_back(AAi);
			vecB.push_back(BBi);

			AAi2 = AAi1;
			BBi2 = BBi1;
			AAi1 = AAi;
			BBi1 = BBi;

			vecP.push_back(PPi);
			vecQ.push_back(QQi);

			if (QQi == 1 && t == 0) {	// если t = 1 и это встретилось впервые, запомним решение
				t = count;
				Gt1 = QQ0*AAi1 - i*BBi1;	// ответ при m = 0;
				Bt1 = BBi1;
			}

			if ((period_begin == true) && (begin_P == PPi && begin_Q == QQi)) {
				l = count - start_index;
				//cout << "КОНЕЦ ПЕРИОДА, Длина периода = " << l << endl;
				break;
			}

			if (QQi >= 0 && period_begin == false) {
				if (((QQi - PPi) <= int_part_sqrtD) && (PPi <= int_part_sqrtD) && ((PPi + QQi) > int_part_sqrtD)) {
					period_begin = true;
					//cout << "Начало периода:  P" << count << " = " << PPi << ", Q" << count << " = " << QQi << endl;
					begin_P = PPi;
					begin_Q = QQi;

					start_index = count;
				}
			}

			count++;
		}
		//t = count - 1;
		//cout << "Длина периода l = " << l << "; t = " << t << endl;
		if (t == 0) {
			//cout << "Нет решения! Нет Qi = 1." << endl;
			continue;
		}
		if ((l % 2 == 0) && (t % 2 != 0)) {	// Проверка
			//cout << "Нет решения! Длина периода четна, а t - нечетно." << endl;
			continue;
		}
		big_int m = 0;
		while ((l*m + t - 1) % 2 == 0) {
			m++;
		}
		//cout << "m= " << m << endl;

		// если l*m + t - 1 > t, то программа не досчитывает A7/B7 (закончила на Al/Bl)
		big_int r = l*m + t - 1;
		//cout << "r = " << r << endl;
		//cout << count - 1 << " - count -1" << endl;
		if (m == 0) {
			res.resG.push_back(Gt1);
			res.resB.push_back(Bt1);
			//cout << "Ответ: G" << (t - 1) << "= " << Gt1 << "; B" << (t - 1) << " = " << Bt1 << endl << endl;
		}
		else if (r >= t) {
			//else if (r >= count-1) {
			if (r < count - 1) {
				//cout << "(Не продолжаем.)" << endl;
				big_int resG = QQ0*vecA[l*m + t - 1] - i*vecB[l*m + t - 1];

				res.resG.push_back(resG);
				res.resB.push_back(vecB[l*m + t - 1]);
				//cout << "Ответ: G" << (l*m + t - 1) << "= " << resG << "; B" << (l*m + t - 1) << " = " << vecB[l*m + t - 1] << endl << endl;
			}
			else {
				//cout << "Продолжаем:" << endl;
				for (big_int i = count - 1/*t*/; i <= r;i++) {
					if (QQi > 0) {
						qqi = prquot((PPi + Arageli::intsqrt(D)), QQi);
					}
					else if (QQi < 0) {
						big_int p;
						big_int q;
						prdivide((PPi + 1 + Arageli::intsqrt(D)), QQi, p, q);
						if (q == 0) {
							qqi = p;
						}
						else {
							qqi = p - 1;
						}
					}

					AAi = qqi*AAi1 + AAi2;
					BBi = qqi*BBi1 + BBi2;

					PPi = qqi*QQi - PPi;
					QQi = (D - square(PPi)) / QQi;

					vecA.push_back(AAi);
					vecB.push_back(BBi);

					AAi2 = AAi1;
					BBi2 = BBi1;
					AAi1 = AAi;
					BBi1 = BBi;
				}

				big_int resG = QQ0*vecA[l*m + t - 1] - i*vecB[l*m + t - 1];

				res.resG.push_back(resG);
				res.resB.push_back(vecB[l*m + t - 1]);
				//cout << "Ответ: G" << (l*m + t - 1) << "= " << resG << "; B" << (l*m + t - 1) << " = " << vecB[l*m + t - 1] << endl << endl;
			}
		}
		res.numberOfClasses++;
	}
	return res;
}


