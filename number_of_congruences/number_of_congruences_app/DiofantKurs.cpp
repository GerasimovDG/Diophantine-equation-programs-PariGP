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
answer MainSolution(big_int D, big_int N);

void Task(big_int start, big_int end, big_int modd, big_int QQ0, Arageli::vector<big_int>* p0Vec);


int main()
{
	setlocale(0, "RUS");



	//ofstream file1; // создаем файл для записи значений
	//ofstream file2;
	//ofstream file3;
	//file1.open(to_string(1) + "ShareOfDecisions.txt");
	//file1 << "Зн-е_D	Знач_N	С_решениями	Всего	Доля	Число_примит-х	ЧислоВсехКлассов" << endl;
	//file2.open("NumberOfClasses.txt");
	//file2 << "Зн-е_D	Знач_N	Число_примит-х	ЧислоВсехКлассов" << endl;
	//file3.open("MeanNumberOfDecisions.txt");
	//file3 << "D	Значение N	Средняя_число_прим-х_решений	Среднее_число_всех_решений" << endl;
	//int countN = 1;	// всего уравнений прошло
	//int countWS = 0; // число уравнений с решениями
	//double share = 0.0;

	cout << "	[.Подсчет начинается с 1.]" << endl;
	big_int D;
	cout << "Введите число D:	";
	cin >> D;
	big_int Z;
	cout << "Введите максимальное N:	";
	cin >> Z;

	while (true) {

		cout << "********************* D = " << D << endl;
		ofstream file1; // создаем файл для записи значений
		ofstream file2;
		ofstream file3;

		stringstream ss;
		ss << D;
		string strD = ss.str();

		file1.open(strD + "_ShareOfDecisions.txt");
		//file1 << "Зн-е_D	Знач_N	С_решениями	Всего	Доля	Число_примит-х	ЧислоВсехКлассов" << endl;
		file1 << "Зн-е_D	Знач_N	С_решениями	Всего	Число_примит-х	ЧислоВсехКлассов" << endl;
		file2.open(strD + "_NumberOfClasses.txt");
		file2 << "Зн-е_D	Знач_N	Число_примит-х	ЧислоВсехКлассов" << endl;
		file3.open(strD + "_MeanNumberOfDecisions.txt");
		file3 << "D	Значение N	Средняя_число_прим-х_решений	Среднее_число_всех_решений" << endl;
		long int countN = 1;	// всего уравнений прошло
		long int countWS = 0; // число уравнений с решениями
		//big_float share = 0.0;




		if (CheckFullSquare(D)) {
			cout << "Число является полным квадратом!" << endl;
			//exit(0);
			continue;
		}

		////////// для x^2 - Dy^2 = N /////////

		big_int begin = 1;
		//cout << "Введите число, с которого начать подсчет(если сначала, то - 1): ";
		//cin >> begin;
		//int Z;
		//big_int N;
		//cout << "Введите максимальное N:	";
		////cin >> N;
		//cin >> Z;

		cout << "Значение N	Среднее_число_прим-х_решений	Среднее_число_всех_решений" << endl;
		big_int primAnsSum = 0;
		big_int allAnsSum = 0;
		//big_float meanPrimCount = 0.0;
		//big_float meanAllCount = 0.0;

		string str_meanPrimCount = "";
		string str_meanAllCount = "";

		// информация для подсчета доли числа уравнений с решениями по отношению ко всем
		//cout << "Введите уже подсчитанное число уравнений с решениями. Пиши значение предыдущего шага (если с начала, то - 0): ";
		//cin >> countWS;
		//cout << "Введите сколько всего уравнений прошло (если с начала, то - 1): ";
		//cin >> countN;
		for (big_int N = begin; N <= Z; N++) {
			//if (N % 100 == 0) {
			//	cout << N << endl;
			//}
			answer ans;
			answer resAns;
			big_int numberAllClasses = 0;	// количество всех классов с решениями


			//for (int n = 1; n <= Arageli::sqrt(Arageli::abs(N)); n++) {	// если делится на квадрат
			//for (int n = Arageli::sqrt(Arageli::abs(N)); n >= 1; n--) {		// чтобы искали максимальный квадрат сразу
			for (big_int n = intsqrt(Arageli::abs(N)); n >= 1; n--) {
				//long double newN = (long double)N / square(n);

				//if (newN == (int)newN) {
				if (N % square(n) == 0) {
					//cout << "N = " << newN << endl;
					big_int newN = N / square(n);
					if (newN == 1) {
						resAns.resG.push_back(Arageli::sqrt(N));
						resAns.resB.push_back(0);
						numberAllClasses = 1;
					}
					else {
						// поиск решения
						ans = MainSolution(D, newN);
						// умножаем ответ на n и добавляем в окончательный ответ
						for (big_int i = 0; i < ans.resG.size(); i++) {
							resAns.resG.push_back(ans.resG[i] * n);
							resAns.resB.push_back(ans.resB[i] * n);
						}
						numberAllClasses = ans.numberOfClasses;
						if (newN == N) {
							resAns.numberOfClasses = ans.numberOfClasses;
						}
					}
					if (newN != N) {
						//cout << "N= " << N << endl;
						ans = MainSolution(D, N);
						for (big_int i = 0; i < ans.resG.size(); i++) {
							resAns.resG.push_back(ans.resG[i]);
							resAns.resB.push_back(ans.resB[i]);
						}
						resAns.numberOfClasses = ans.numberOfClasses;
						numberAllClasses += ans.numberOfClasses;
					}
					break;
				}
			}
			if (resAns.resG.size() == 0) {
				//cout << "Нет решений!" << endl;
			}
			else {
				countWS++;
			}
			//share = (double)countWS / countN;	// доля уравнений с решениями
			//share = big_float(countWS) / big_float(countN);
			file1 << D << "	" << N << "	" << countWS << "	" << countN << "	" /*<< share << "	"*/ << resAns.numberOfClasses << "	" << numberAllClasses << endl;
			if (numberAllClasses != 0) {
				file2 << D << "	" << N << "	" << resAns.numberOfClasses << "	" << numberAllClasses << endl;

				primAnsSum += resAns.numberOfClasses;
				allAnsSum += numberAllClasses;

				//meanPrimCount = big_float(big_float(primAnsSum) / (big_float)N);
				//meanPrimCount = long double(primAnsSum / (long double)countN);
				//meanAllCount = big_float(big_float(allAnsSum) / (big_float)N);
				str_meanPrimCount = DivWithRem(primAnsSum, countN);
				str_meanAllCount = DivWithRem(allAnsSum, countN);

				file3 << D << "	" << N << "	" << str_meanPrimCount << "	" << str_meanAllCount << endl;
				//cout << setw(10) << left << N << "	" << setw(20) << str_meanPrimCount << "	" << setw(20) << str_meanAllCount << endl;
				if (N % 100 == 0) {
					cout << setw(10) << left << N << "	" << setw(20) << str_meanPrimCount << "	" << setw(20) << str_meanAllCount << endl;
				}
			}


			countN++;
		}
		////////////////////////////////////////new
		//Arageli::vector<big_int> vecResG;
		//Arageli::vector<big_int> vecResB;
		//
		//for (int i = 0; i < resAns.resG.size();i++) {
		//	cout << "(" << resAns.resG[i] << "; " << resAns.resB[i] << ")" << endl;
		//}
		//
		//Arageli::vector<int> p0Vec;
		//int l = 0;
		//big_int QQ0 = N;
		//
		//// Проверка 1.
		//if (CheckSolutionExist(D, N) == false) {
		//	//cout << "!Нет решения!\n[Код проверки: 1]\n" << endl;
		//	return 0;
		//}
		//
		//// находим подходящие P0 из интервала
		//for (big_int i = -Arageli::abs(QQ0) / 2; i <= Arageli::abs(QQ0) / 2;i++) {
		//	if (square(i) % QQ0 == D % QQ0) {
		//		p0Vec.push_back(i);
		//	}
		//}
		//// Проверка 2.
		//if (p0Vec.is_empty()) {
		//	cout << "Нет подходящих P0. Нет решения.\n[Код проверки: 2]\n";
		//	system("pause");
		//	return 0;
		//}
		//
		//cout << "Подходящие P0: " << endl;
		//for (big_int i : p0Vec) {
		//	cout << i << "-> " << endl;
		//}
		//
		//for (big_int i : p0Vec) {
		//	cout << "ВЫБРАНО P0 = "<< i << endl;
		//	big_int PP0 = i;
		//
		//	Arageli::vector<big_int> vecA;
		//	Arageli::vector<big_int> vecB;
		//	Arageli::vector<big_int> vecP;
		//	Arageli::vector<big_int> vecQ;
		//
		//	int QQi = QQ0;
		//	int PPi = PP0;
		//
		//	big_int qqi = 0, qq0 = 0;
		//	big_int PPi1 = 0;
		//	big_int QQi1 = 0;
		//	big_int AAi2 = 0;
		//	big_int AAi1 = 1;
		//	big_int BBi2 = 1;
		//	big_int BBi1 = 0;
		//
		//	big_int Gt1 = 0;		// ответ при m = 0;
		//	big_int Bt1 = 0;
		//
		//	//qq0 = Arageli::floor((i + Arageli::sqrt(long double(D))) / QQ0);
		//	//cout << "qq0 = " << qq0;
		//
		//	big_int AAi = 0, BBi = 0, GGi = 0;
		//	int count = 1;
		//	int t = 0;
		//	bool stopFlag = true;
		//	while (stopFlag == true) {
		//		qqi = Arageli::floor((PPi + Arageli::sqrt(long double(D))) / QQi);
		//		cout << "q" << count-1 << " = " << qqi ;
		//
		//		AAi = qqi*AAi1 + AAi2;
		//		BBi = qqi*BBi1 + BBi2;
		//
		//		vecA.push_back(AAi);
		//		vecB.push_back(BBi);
		//
		//		PPi1 = qqi*QQi - PPi;
		//		QQi1 = (D - square(PPi1)) / QQi;
		//
		//		//cout << "A" << count-1 << "/B" << count-1 << " = " << AAi << "/" << BBi << endl;
		//		cout << "; A" << count - 1 << "/B" << count - 1 << " = " << vecA[count - 1] << "/" << vecB[count - 1] << endl;
		//		cout << "P" << count << " = " << PPi1 << "; Q" << count << " = " << QQi1 << endl << endl;
		//
		//		PPi = PPi1;
		//		QQi = QQi1;
		//		AAi2 = AAi1;
		//		BBi2 = BBi1;
		//		AAi1 = AAi;
		//		BBi1 = BBi;
		//
		//		vecP.push_back(PPi);
		//		vecQ.push_back(QQi);
		//
		//		if (QQi == 1 && t == 0) {	// если t = 1 и это встретилось впервые, запомним решение
		//			t = count;
		//			Gt1 = QQ0*AAi1 - PP0*BBi1;	// ответ при m = 0;
		//			Bt1 = BBi1;
		//		}
		//		for (int i = 0; i < count-1; i++) {
		//			//cout << "i: " << i << "/ (Pi,Qi) = (" << PPi << ", " << QQi << ").  vecPi,vecQi = " << vecP[i] << ", " << vecQ[i] << endl;
		//			if (PPi == vecP[i] && QQi == vecQ[i]) {	//если значение повторилось, то нашли период
		//				//стоп
		//				l = count-1 - i;	// длина периода
		//				//cout << "i: "<< i << "/ (Pi,Qi) = (" << PPi << ", " << QQi << ").  vecPi,vecQi = " << vecP[i] << ", " << vecQ[i] << endl;
		//				stopFlag = false;
		//			}
		//		}
		//
		//		count++;
		//	}
		//	//t = count - 1;
		//	cout << "Длина периода l = " << l << "; t = " << t << endl;
		//	if (t == 0) {
		//		cout << "Нет решения! Нет Qi = 1." << endl;
		//		//system("pause");
		//		//return 0;
		//		continue;
		//	}
		//	if ((l % 2 == 0) && (t % 2 != 0)) {	// Проверка
		//		cout << "Нет решения! Длина периода четна, а t - нечетно." << endl;
		//		//system("pause");
		//		//return 0;
		//		continue;
		//	}
		//	int m = 0;
		//	while ((l*m + t - 1) % 2 == 0) {
		//		m++;
		//	}
		//	cout << "m= " << m << endl;
		//
		//	// если l*m + t - 1 > t, то программа не досчитывает A7/B7 (закончила на Al/Bl)
		//	int r = l*m + t - 1;
		//	cout << "r = " << r << endl;
		//	cout << count - 1 << " - count -1" << endl;
		//	if (m == 0) {
		//		vecResG.push_back(Gt1);
		//		vecResB.push_back(Bt1);
		//		cout << "Ответ: G" << (t - 1) << "= " << Gt1 << "; B" << (t - 1) << " = " << Bt1 << endl << endl;
		//	}
		//	else if (r >= t) {
		//	//else if (r >= count-1) {
		//		if (r < count - 1) {
		//			cout << "(Не продолжаем.)" << endl;
		//			big_int resG = QQ0*vecA[l*m + t - 1] - PP0*vecB[l*m + t - 1];
		//			vecResG.push_back(resG);
		//			vecResB.push_back(vecB[l*m + t - 1]);
		//			cout << "Ответ: G" << (l*m + t - 1) << "= " << resG << "; B" << (l*m + t - 1) << " = " << vecB[l*m + t - 1] << endl << endl;
		//		}
		//		else {
		//			cout << "Продолжаем:" << endl;
		//			for (int i = count - 1/*t*/; i <= r;i++) {
		//				//cout << "PPi = " << PPi << ", QQi = " << QQi << endl;
		//				qqi = Arageli::floor((PPi + Arageli::sqrt(long double(D))) / QQi);
		//				cout << "q" << i << " = " << qqi;
		//
		//				AAi = qqi*AAi1 + AAi2;
		//				BBi = qqi*BBi1 + BBi2;
		//
		//				PPi1 = qqi*QQi - PPi;
		//				QQi1 = (D - square(PPi1)) / QQi;
		//
		//				vecA.push_back(AAi);
		//				vecB.push_back(BBi);
		//				cout << "; A" << i << "/B" << i << " = " << vecA[i] << "/" << vecB[i] << endl;
		//				cout << "P" << i + 1 << " = " << PPi1 << "; Q" << i + 1 << " = " << QQi1 << endl << endl;
		//
		//				PPi = PPi1;
		//				QQi = QQi1;
		//				AAi2 = AAi1;
		//				BBi2 = BBi1;
		//				AAi1 = AAi;
		//				BBi1 = BBi;
		//			}
		//
		//			big_int resG = QQ0*vecA[l*m + t - 1] - PP0*vecB[l*m + t - 1];
		//			vecResG.push_back(resG);
		//			vecResB.push_back(vecB[l*m + t - 1]);
		//			cout << "Ответ: G" << (l*m + t - 1) << "= " << resG << "; B" << (l*m + t - 1) << " = " << vecB[l*m + t - 1] << endl << endl;
		//		}
		//	}
		//}
		//
		//if (vecResG.size() == 0) {
		//	cout << "Нет примитивных решений!" << endl;
		//}
		//else{
		//	cout << "Все примитивные решения: " << endl;
		//	for (int i = 0; i < vecResG.size();i++) {
		//		cout << "(" << vecResG[i] << "; " << vecResB[i] << ")" << endl;
		//	}
		//}
		//
		//resAns.resG = vecResG;
		//resAns.resB = vecResB;
		//for (int i = 0; i < resAns.resG.size();i++) {
		//	cout << "---------(" << resAns.resG[i] << "; " << resAns.resB[i] << ")-----------" << endl;
		//}


		file1.close();
		file2.close();
		file3.close();
		
		D = next_prime(D);
	}
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
			//system("pause");
			return false;
		}
	}
	// 2.
	if (D == pow(N, 2) - 1) {
		if (CheckFullSquare(N) == true) {
			//cout << "ОООООтвет: (" << Arageli::sqrt(N) << "; " << 0 << ").\n";
			//system("pause");
			//return false;
			return true;
		}
		else {
			//cout << "!Нет решения!\n[Код проверки: 1.2]\n" << endl;
			//system("pause");
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

//long int powmod1(long int a, long int p)
//{
//	int b = 2;
//	long int res = 1;
//	while (b)
//		if (b & (long int)1)
//			res = long int(res * 1ll * a % p), --b;
//		else
//			a = long int(a * 1ll * a % p), b >>= 1;
//	return res;
//}

void Task(big_int start, big_int end, big_int modd, big_int QQ0, Arageli::vector<big_int>* p0Vec) {
	for (big_int i = start; i < end; i++) {
		if (powmod(i, QQ0) == modd) {
			if (i != 0) {
				p0Vec->push_back(-i);
			}
			p0Vec->push_back(i);
		}
	}
}

string DivWithRem(big_int numerator, big_int denom) {
	big_int integ = numerator / denom;
	big_int remain = numerator % denom;
	big_int temp = remain;

	stringstream ss;
	ss << integ;
	string res = ss.str();
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
		//res += temp;
		count++;
		if (count > 6) {
			break;
		}
	}
	return res;
}

answer MainSolution(big_int D, big_int N) {
	answer res;

	Arageli::vector<big_int> p0Vec;
	Arageli::vector<big_int> p0Vec2;
	Arageli::vector<big_int> p0Vec3;
	Arageli::vector<big_int> p0Vec4;

	//Arageli::vector<big_int> p0Vec;
	big_int l = 0;
	big_int QQ0 = N;

	// Проверка 1.
	if (CheckSolutionExist(D, N) == false) {
		//cout << "!Нет решения!\n[Код проверки: 1]\n" << endl;
		return res;
		//exit(0);
	}
	
	// РАЗБИВАТЬ НА ПОТОКИ НЕТ СМЫСЛА. УСКОРЕНИЯ НЕТ
	//big_int modd = Arageli::mod(D, QQ0);
	//big_int drob = (Arageli::abs(QQ0) / 2);
	//big_int partdrob = drob / 4;

	//std::thread th(Task, 0, partdrob, modd, QQ0, &p0Vec);
	//th.join();
	//std::thread th2(Task, partdrob, 2 * partdrob, modd, QQ0, &p0Vec2);
	//th2.join();
	//std::thread th3(Task, 2 * partdrob, 3 * partdrob, modd, QQ0, &p0Vec3);
	//th3.join();
	//std::thread th4(Task, 3 * partdrob, drob, modd, QQ0, &p0Vec4);
	//th4.join();


	//std::copy(p0Vec2.begin(), p0Vec2.end(), std::back_inserter(p0Vec));
	//std::copy(p0Vec3.begin(), p0Vec3.end(), std::back_inserter(p0Vec));
	//std::copy(p0Vec4.begin(), p0Vec4.end(), std::back_inserter(p0Vec));


		big_int modd = Arageli::mod(D, QQ0);
		for (big_int i = 0; i < Arageli::abs(QQ0) / 2; i++) {
			//if (cmp(square(i) % QQ0, mod/* D % QQ0*/) == 0) {	// сравнение
			if (powmod(i, QQ0) == modd) {
				//if (Arageli::mod(square(i), QQ0) == modd) {
				if (i != 0) {
					p0Vec.push_back(-i);
				}
				p0Vec.push_back(i);
			}
		}
		big_int i = Arageli::abs(QQ0) / 2;	// дополнение чтобы левая граница была строгая, а правая нестрогая ( -Q0/2 < x <= Q02/ )
		//if (powmod1(i, QQ0) == modd) {
		//	//if (Arageli::mod(square(i), QQ0) == modd) {
		//	if (QQ0 % 2 != 0) {				//	если N = нечетное, то левая граница тоже входит
		//		p0Vec.push_back(-i);
		//	}
		//	p0Vec.push_back(i);
		//}
		if (Arageli::mod(square(i), QQ0) == modd) {
			if (QQ0 % 2 != 0) {				//	если N = нечетное, то левая граница тоже входит
				p0Vec.push_back(-i);
			}
			p0Vec.push_back(i);
		}


	
	// Проверка 2.
	if (p0Vec.is_empty()) {
		//cout << "Нет подходящих P0. Нет решения.\n[Код проверки: 2]\n";
		return res;
	}

	//cout << "Подходящие P0: " << endl;
	//for (big_int i : p0Vec) {
	//	cout << i << "-> " << endl;
	//}


	for (big_int i : p0Vec) {
		//cout << "ВЫБРАНО P0 = " << i << endl;
		//big_int PP0 = i;

		Arageli::vector<big_int> vecA;
		Arageli::vector<big_int> vecB;
		Arageli::vector<big_int> vecP;
		Arageli::vector<big_int> vecQ;

		big_int QQi = QQ0;
		big_int PPi = i;

		big_int qqi = 0;// , qq0 = 0;
		big_int int_part_sqrtD = 0;

		//big_int PPi1 = 0;
		//big_int QQi1 = 0;
		big_int AAi2 = 0;
		big_int AAi1 = 1;
		big_int BBi2 = 1;
		big_int BBi1 = 0;

		big_int Gt1 = 0;		// ответ при m = 0;
		big_int Bt1 = 0;

		//qq0 = Arageli::floor((i + Arageli::sqrt(long double(D))) / QQ0);
		//cout << "qq0 = " << qq0;

		big_int begin_P;
		big_int begin_Q;

		big_int AAi = 0, BBi = 0, GGi = 0;
		big_int count = 1;
		big_int start_index = 0;
		big_int t = 0;
		bool period_begin = false;

		while (true) {
			//qqi = Arageli::floor((PPi + Arageli::sqrt(long double(D))) / QQi);
			//cout << "q" << count - 1 << " = " << qqi;
			if (QQi > 0) {
				qqi = prquot((PPi + intsqrt(D)), QQi);
			}
			else if (QQi < 0) {
				big_int p;
				big_int q;
				prdivide((PPi + 1 + intsqrt(D)), QQi, p, q);
				if (q == 0) {
					//qqi = prquot((PPi + 1 + intsqrt(D)), QQi);
					qqi = p;
				}
				else {
					//qqi = prquot((PPi + 1 + intsqrt(D)), QQi) - 1;
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

			//PPi1 = qqi*QQi - PPi;
			//QQi1 = (D - square(PPi1)) / QQi;

			////cout << "A" << count-1 << "/B" << count-1 << " = " << AAi << "/" << BBi << endl;
			//cout << "; A" << count - 1 << "/B" << count - 1 << " = " << vecA[count - 1] << "/" << vecB[count - 1] << endl;
			//cout << "P" << count << " = " << PPi1 << "; Q" << count << " = " << QQi1 << endl << endl;

			//PPi = PPi1;
			//QQi = QQi1;
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
			//for (int i = 0; i < count - 1; i++) {
			//	if (PPi == vecP[i] && QQi == vecQ[i]) {	//если значение повторилось, то нашли период
			//		//стоп
			//		l = count - 1 - i;	// длина периода
			//		stopFlag = false;
			//	}
			//}
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
			//system("pause");
			//return 0;
			continue;
		}
		if ((l % 2 == 0) && (t % 2 != 0)) {	// Проверка
			//cout << "Нет решения! Длина периода четна, а t - нечетно." << endl;
			//system("pause");
			//return 0;
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
			//vecResG.push_back(Gt1);
			//vecResB.push_back(Bt1);

			res.resG.push_back(Gt1);
			res.resB.push_back(Bt1);
			//cout << "Ответ: G" << (t - 1) << "= " << Gt1 << "; B" << (t - 1) << " = " << Bt1 << endl << endl;
		}
		else if (r >= t) {
			//else if (r >= count-1) {
			if (r < count - 1) {
				//cout << "(Не продолжаем.)" << endl;
				big_int resG = QQ0*vecA[l*m + t - 1] - i*vecB[l*m + t - 1];
				//vecResG.push_back(resG);
				//vecResB.push_back(vecB[l*m + t - 1]);

				res.resG.push_back(resG);
				res.resB.push_back(vecB[l*m + t - 1]);
				//cout << "Ответ: G" << (l*m + t - 1) << "= " << resG << "; B" << (l*m + t - 1) << " = " << vecB[l*m + t - 1] << endl << endl;
			}
			else {
				//cout << "Продолжаем:" << endl;
				for (big_int i = count - 1/*t*/; i <= r;i++) {
					////cout << "PPi = " << PPi << ", QQi = " << QQi << endl;
					//qqi = Arageli::floor((PPi + Arageli::sqrt(long double(D))) / QQi);
					//cout << "q" << i << " = " << qqi;
					if (QQi > 0) {
						qqi = prquot((PPi + Arageli::intsqrt(D)), QQi);
					}
					else if (QQi < 0) {
						big_int p;
						big_int q;
						prdivide((PPi + 1 + Arageli::intsqrt(D)), QQi, p, q);
						if (q == 0) {
							//qqi = prquot((PPi + 1 + intsqrt(D)), QQi);
							qqi = p;
						}
						else {
							//qqi = prquot((PPi + 1 + intsqrt(D)), QQi) - 1;
							qqi = p - 1;
						}
					}


					AAi = qqi*AAi1 + AAi2;
					BBi = qqi*BBi1 + BBi2;

					PPi = qqi*QQi - PPi;
					QQi = (D - square(PPi)) / QQi;

					vecA.push_back(AAi);
					vecB.push_back(BBi);
					//cout << "; A" << i << "/B" << i << " = " << vecA[i] << "/" << vecB[i] << endl;
					//cout << "P" << i + 1 << " = " << PPi1 << "; Q" << i + 1 << " = " << QQi1 << endl << endl;

					//PPi = PPi1;
					//QQi = QQi1;
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


