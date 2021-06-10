checkSolutionExist = (D,N) -> {
	if (N == -3,
		\\ тут D-2 или D*D-2?
		if (issquare(D*D - 2),
			\\ print("No solution. [1]");
			return(0);
		);
	);
	if (D == (N*N) - 1,
		if (issquare(N),
			return(1),
				\\ print("No solution. [2]");
				return(0);
		);
	);
	return(1);
};


\\ division of two integers with a real result ( but in string)
divWithRem = (numeratr, denom) -> {
	my(res,integ,remain,temp,count);
	res = "";
	if(((numeratr < 0 && denom >= 0) || (numeratr >=0 && denom < 0)),
		res = Str(res, "-");
		numeratr = abs(numeratr);
		denom = abs(denom);
	);

	integ = numeratr \ denom;
	remain = numeratr % denom;
	temp  = remain;

	res = Str(res, integ);

	if(remain != 0,
		res = Str(res,".");
	);

	count = 1;
	while(remain != 0,
		temp = remain * 10;
		remain = temp % denom;
		temp = temp \ denom;

		res = Str(res,temp);

		count++;
		if(count > 6,
			break;
		);
	);
	return(res);
}

\\ *************************************************************
\\ *** all solutions of quadratic comparison x^2 = D (mod N) ***

chineseForArray = (mass, i = 1) -> {
	if(i == length(mass)-1,
		return(chinese(mass[i],mass[i+1])),
			return(chinese(mass[i], chineseForArray(mass, i+1)));
	);
}

changeRes = (mas, m, res) -> {
	my(i);
	for(i=1,m,
		if (mas[i] == 1,
			res[i] = -res[i];
		);
	);
	return(res);
}

setAllCombinations = (m, res) -> {
	my(n,i, binarMas, newRes);
	n = 2;
	binarMas = List();
	for(i= 1,m,
		listput(binarMas, 0);
	);
	i = 1;
	\\ print("All 2^s cases:");
	while(binarMas != 0,
		res[i] = changeRes(binarMas,m, res[i]);
		\\ print(i " = "res[i]);
		i++;
		binarMas = nextSet(binarMas, m, n);
		\\ print(i ": "binarMas);
	);
	newRes = List();
	listput(newRes, res[1]);
	for(i = 2,length(res),
		\\if(res[1] != res[i],
			listput(newRes, res[i]);
		\\);
	);
	newRes = List(Set(newRes));
	return(newRes);
}

nextSet = (a, m, n = 2) -> {
	my(j);
	j = m;
	while(j >= 1 && a[j] == n-1, j--);
	if(j < 1, return(0));
	if(a[j] >= n-1,
		j--;
	);
	a[j]++;
	if (j == m, return(a));
	for(k = j+1, m,
		a[k] = 0;
	);
	return(a);
}


convertBinaryToAnswer = (masX, a) -> {
	my(res,x);
	res = List();
	for(i = 1, #a,
		x = masX[i][a[i]+1];
		listput(res,x);
	);
	return(res);
}
lengthOfLongestArray = (arrays) -> {
	my(size);
	size = 0;
	for(i=1, #arrays,
		if( length(arrays[i]) > size,
			size = length(arrays[i]);
		);
	);
	return(size);

}
searchAllCombinations2 = (masX) -> {
	my(n,m,binaryArray,binaryArrays, res, masRes);
	n = lengthOfLongestArray(masX);
	m = length(masX);

	binaryArrays = List();
	binaryArray = List();
	for (i = 1, m, listput(binaryArray,0));

	while(binaryArray != 0,
		listput(binaryArrays, binaryArray);
		binaryArray = nextSet(binaryArray, m, n);
	);

	for(i = 1, length(binaryArrays),
		for( j =1, length(binaryArrays[i]),
			binaryArrays[i][j] = binaryArrays[i][j] % length(masX[j]);
		);
	);
	binaryArrays = List(Set(binaryArrays));

	res = List();
	for(i = 1,length(binaryArrays),
		masRes = convertBinaryToAnswer(masX, binaryArrays[i]);
		listput(res, masRes);
		\\print(i": "masRes);
	);
	return(res);

}

searchAllCombinations3 = (arrays) -> {
	my(ind,ind_const, k1, res, value, j, f);
	ind = List();
	ind_const = List();
	k1 = length(arrays);
	res = List();
	
	for(i = 1, k1,
		listput(ind, length(arrays[i]));
		listput(ind_const, length(arrays[i]));
	);

	while(ind[1] != 0,
		value = List();
		for(k=1, k1,
			listput(value, arrays[k][ind[k]]);
		);
		listput(res,value);
		j = 0;
		f = 1;

		while(f == 1,
			if (ind[k1 - j] != 1,
				ind[k1 - j]--;
				f = 0,
					if (k1 - j != 1,
						ind[k1-j] = ind_const[k1-j];
						j++,
							ind[1] = 0;
							break;
					);
			); 
		);
	);
	return(res);
}

quadraticComparisonForCoPrime = (D, N) -> {
	my(f,lengthF,masX,allPosSol,answers,allCombinations);

	f = factorint(abs(N));
	\\print("f = "f);
	lengthF = length(f[,1]);
	masX = List();
	for(i = 1,lengthF,
		listput(masX, List());
	);
	for(i = 1, lengthF,
		if(issquare(Mod(D, f[i,1]^f[i,2]), &x),
			if (f[i,1] == 2,
				if(f[i,2] == 1, listput(masX[i], x),
					if(f[i,2] == 2, 
					listput(masX[i], x);
					listput(masX[i], -x),
						if(f[i,2] > 2, 
							listput(masX[i], x);
							listput(masX[i], -x);
							listput(masX[i], x+2^(f[i,2]-1));
							listput(masX[i], -x-2^(f[i,2]-1));
						),
					);
				),
					listput(masX[i], x);
					listput(masX[i], -x);	
			),
				return(0);
		);
	);
	
	\\allCombinations = searchAllCombinations2(masX);
	allCombinations = searchAllCombinations3(masX);
	\\print("allCombinations =");
	\\for(i = 1,length(allCombinations),
	\\	print(allCombinations[i])
	\\);

	answers = List();
	for(j=1,length(allCombinations),
		if(length(allCombinations[1]) == 1,
			listput(answers, lift(allCombinations[j][1])),
			\\print1(lift(allCombinations[j][1])", "),
				listput(answers, lift(chineseForArray(allCombinations[j])));
				\\print1(lift(chineseForArray(allPosSol[j]))"_ ");
		);
	);

	\\print("answers:");
	\\for(i=1,#answers, print(answers[i]));
	
	return(answers);
}

quadraticComparisonForNotCoPrime3 = (D, N) -> {
	my(answer,fact,gcdDN, answers,n,f,r,pn,m,tmpRes, x, p, allCombinations, res);

	answers = List();
	fact = factorint(abs(N));
	\\print("f = "fact);
	for(i = 1, length(fact[,1]),
		answer = List();
		p = fact[i,1];
		n = fact[i,2];

		gcdDN = gcd(D,p^n);
		\\print("p = "p", n ="n", gcd ="gcdDN);

		if(gcdDN == 1,
			if(issquare(Mod(D, p^n), &x),
				if (p == 2,
					if(n == 1, listput(answer, x),
						if(n == 2, 
						listput(answer, x);
						listput(answer, -x),
							if(n > 2, 
								listput(answer, x);
								listput(answer, -x);
								listput(answer, x+2^(n-1));
								listput(answer, -x-2^(n-1));
							),
						);
					),
						listput(answer, x);
						listput(answer, -x);	
				),
					return(0);
			),
				f = factorint(gcdDN);
				r = f[1,2];
				\\print("r ="r);
				if (r == n,
					if (n % 2 == 0,
						m = n \ 2;
						if (n != 0 && m == 0, 
							m = n / 2.;
						);
						\\print("m =! "m);
						for(t = 0, p^m-1,
							x = p^m*t;

							listput(answer, Mod(x, p^n));
						),
							m = (n-1) \ 2;
							if ((n-1) != 0 && m == 0, 
								m = (n-1) / 2.;
							);
							\\	print("m == "m);
							\\print("n == "n);
							\\print("p == "p);
							for (t = 0, p^m-1,
								x = p^(m+1)*t;
								listput(answer,Mod(x, p^n));
							);
					),
						\\ 0 < r < n
						if (r % 2 == 1,
							\\print("No solutions");
							return(0),
								m = r \ 2;
								if (r != 0 && m == 0, 
									m = r / 2.;
								);
								\\	print("m = "m);
								tmpRes = List();
								tmpRes = quadraticComparisonForCoPrime(D \ p^r, p^(n-r));
								if(tmpRes == 0,
									\\print("No solutions");
									return(0);
								);
								\\print("tmpRes ="tmpRes);
								for(i = 1, length(tmpRes),
									for(t = 0,p^m-1,
										x = p^m * (tmpRes[i] + p^(n-r) * t);
										listput(answer, Mod(x,p^n));
									);
								);
						);
				);

		);
		listput(answers,answer);
		\\print("answer = "answer);
	);
		\\print("answers = "answers);

	\\allCombinations = searchAllCombinations2(answers);
	allCombinations = searchAllCombinations3(answers);
	\\print("allCombinations =");
	\\for(i = 1,length(allCombinations),
	\\	print(allCombinations[i])
	\\);

	res = List();
	for(j=1,length(allCombinations),
		if(length(allCombinations[1]) == 1,
			listput(res, lift(allCombinations[j][1])),
				listput(res, lift(chineseForArray(allCombinations[j])));
		);
	);

	return(res);
}

pellCalcFor = (D, N = 1) -> {
	my(d,a2,a1,ai,b2,b1,bi,p1,p0,pii,q1,q0,qi,count,l, x, y, xx, yy, res);

	res = vector(3);
	res[1] = List();
	res[2] = List();
	res[3] = 0;

	if (type(D) != "t_INT" || D < 1, 
		return(res);
	);

	d = sqrtint(D);
	if (sqr(d) == D,
		return(res);
	);
	a2 = 0; a1 = 1;
	b2 = 1; b1 = 0;
	p0 = 0; q0 = 1;
	q0 = ((p0+d)\q0);
	
	qi = q0;
	p1= 0; q1 = 1;

	count = 0;
	while(qi != 2*q0,
		qi = ((p1 + d) \ q1);

		ai = qi*a1 + a2;
		bi = qi*b1 + b2;

		pii = qi * q1 - p1;
		qii = (D - sqr(pii))\q1;

		p1 = pii;
		q1 = qii;

		a2 = a1;
		b2 = b1;
		a1 = ai;
		b1 = bi;

		count = count + 1;
	);
	l = count - 1;

	x = a2; y = b2;

	if(x*x - D*y*y == -1,
		xx = x*x+D*y*y;
		yy = 2*x*y;
	);

	if (N == x*x - D*y*y,
		listput(res[1], x);
		listput(res[2], y);
		res[3] = res[3] + 1,
			if (N == xx*xx - D*yy*yy,
				listput(res[1], xx);
				listput(res[2], yy);
				res[3] = res[3] + 1,));
	
	return(res);
};

\\ ******** main solution **************
mainSolution = (D,N) -> {
	my(l, QQ0, res, p0Vec, gcdDN, PPi, qqi, int_part_sqrtD,
	vecA, vecB, vecP, vecQ, AAi2, AAi1, BBi2, BBi1, Gt1, Bt1, 
	begin_P, begin_Q, GGi, AAi, BBi, count, start_index, t, period_begin, r, m);
	\\ print("N = "N);

	l = 0;
	QQ0 = N;

	\\ res[1] - resG;
	\\ res[2] - resB;
	\\ res[3] - numberOfClasses;
	res = vector(3);
	res[1] = List();
	res[2] = List();
	res[3] = 0;

	p0Vec = List();

	if (!checkSolutionExist(D,N), return(res));

	gcdDN = gcd(D,N);
	if (abs(N) == 1,
		listput(p0Vec, 0),
			if (gcdDN == 1,
				p0Vec = quadraticComparisonForCoPrime(D, N),
					p0Vec = quadraticComparisonForNotCoPrime3(D, N);
			);
	);

	\\print("p0Vec = "p0Vec);
	

	if(length(p0Vec) == 0,return(res));

	if (abs(N) == 1,
		if (N == 1,
			res = pellCalcFor(D),
				res = pellCalcFor(D , -1)),

			for(j=1,length(p0Vec),
				vecA = List(); vecB = List();
				vecP = List(); vecQ = List();

				QQi = QQ0;
				PPi = p0Vec[j];

				qqi = 0;
				int_part_sqrtD = 0;

				AAi2 = 0; AAi1 = 1;
				BBi2 = 1; BBi1 = 0;

				\\ answer if m = 0;
				Gt1 = 0; Bt1 = 0; 

				begin_P = 0; begin_Q = 0;

				AAi = 0; BBi = 0; GGi = 0;
				count = 1;
				start_index = 0;
				t = 0;
				\\ period_begin - boolean
				period_begin = 0;

				while(1,
					if (QQi > 0,
						qqi = ((PPi + sqrtint(D)) \ QQi),
							if(QQi < 0,
								pq = divrem(PPi + 1 + sqrtint(D), QQi);
								if (pq[2] == 0, 
									qqi = pq[1],
										qqi = pq[1] - 1;
								);
							);
					);
					int_part_sqrtD = sqrtint(D);

					PPi = qqi * QQi - PPi;
					QQi = (D - sqr(PPi)) \ QQi;

					AAi = qqi * AAi1 + AAi2;
					BBi = qqi * BBi1 + BBi2;

					listput(vecA, AAi);
					listput(vecB, BBi);

					AAi2 = AAi1; BBi2 = BBi1;
					AAi1 = AAi; BBi1 = BBi;

					listput(vecP, PPi);
					listput(vecQ, QQi);

					if (QQi == 1 && t == 0,
						t = count;
						Gt1 = QQ0 * AAi1 - p0Vec[j] * BBi1;
						\\ Gt1 = QQ0 * AAi1 - j * BBi1;
						Bt1 = BBi1;
					);

					if ((period_begin == 1) && (begin_P == PPi && begin_Q == QQi),
						l = count - start_index;
						break;
					);

					if (QQi >= 0 && period_begin == 0, 
						if (((QQi - PPi) <= int_part_sqrtD) && (PPi <= int_part_sqrtD) && ((PPi + QQi) > int_part_sqrtD),
							period_begin = 1;
							begin_P = PPi;
							begin_Q = QQi;
							start_index = count;
						);
					);
					count = count + 1;
				);

				if (t == 0, next);
				if ((l % 2 == 0) && (t % 2 != 0), next);

				m = 0;
				while (((l*m + t - 1) % 2 == 0), m = m + 1);

				r = l*m + t - 1;

				if (m == 0,
					listput(res[1], Gt1);
					listput(res[2], Bt1),
						if (r >= t,
							if (r < count - 1,
								resG = QQ0 * vecA[l*m + t] - p0Vec[j] * vecB[l*m + t];
								\\ resG = QQ0 * vecA[l*m + t] - j * vecB[l*m + t];
								listput(res[1], resG);
								listput(res[2], vecB[l*m + t]),
									for(i = count-1, r,
										if (QQi > 0,
											qqi = (PPi + sqrtint(D)) \ QQi,
											if (QQi < 0,
												pq = divrem(PPi + 1 + sqrtint(D), QQi);
												if (pq[2] == 0, 
													qqi = pq[1],
														qqi = pq[1] - 1;
												);
											);
										);

										AAi = qqi * AAi1 + AAi2;
										BBi = qqi * BBi1 + BBi2;

										PPi = qqi * QQi - PPi;
										QQi = (D - sqr(PPi)) \ QQi;

										listput(vecA, AAi);
										listput(vecB, BBi);

										AAi2 = AAi1;
										BBi2 = BBi1;
										AAi1 = AAi;
										BBi1 = BBi;	
									);
								resG = QQ0 * vecA[l * m + t] - p0Vec[j] * vecB[l*m + t];
								\\ resG = QQ0 * vecA[l* m + t] - j * vecB[l*m + t];

								listput(res[1], resG);
								listput(res[2], vecB[l*m+t]);
							);
							
						);

				);
				res[3] = res[3] + 1;
			);
	);
	return(res);
};


calculationOnlyPrimitiveSolution = (D, Nmin, Nmax, threashold, threshold_n, outFile = -1) -> {
	my(countN, countWS, primAnsSum, str_meanPrimCount, resAns, mean);

		print("****************** D = "D);

		\\ total number of completed equation
		countN = 1;
		\\ number of equation with solution
		countWS = 0;

		if(issquare(D),
			print("ERROR - D should not be a perfect square");
			return(0);
		);

		if (outFile == -1,
			outputPrimiviteInfo(D,Nmin, Nmax);
			return();
		);

		print("Value D 	mean_number_primitive_solution");
		primAnsSum = 0;
		str_meanPrimCount = "";

		for(N = Nmin, Nmax,
			if(N == 0, next);

			resAns = vector(3);
			resAns[1] = List();
			resAns[2] = List();
			resAns[3] = 0;

			resAns = mainSolution(D, N);

			if(length(resAns[1]) == 0,
				,
					countWS++;
			);

			if(resAns[3] != 0,
				primAnsSum += resAns[3];
				str_meanPrimCount = divWithRem(primAnsSum, countN);
			);
			if(N % 100 == 0,
				printf("%3d	%s\n", N, str_meanPrimCount);
			);

			if (N == threshold_n,
				mean = primAnsSum / countN;

				filewrite(outFile, Str(D, "	", str_meanPrimCount));
				fileflush(outFile);

				if (mean >= threashold,
					outputPrimiviteInfo(D,Nmin, Nmax);
					break;
				);
				break;
			);
			countN++;

		);
}


calculationAllSolution = (D, Nmin, Nmax, threashold, threshold_n, outFile = -1) -> {
	my(countN, countWS, primAnsSum, allAnsSum, str_meanPrimCount, str_meanAllCount,
	 resAns, ans, mean, numberAllClasses, newN);
		print("****************** D = "D);

		\\ total number of completed equation
		countN = 1;
		\\ number of equation with solution
		countWS = 0;

		if(issquare(D),
			print("ERROR - D should not be a perfect square");
			return(0);
		);
		if (outFile == -1,
			outputAllInfo(D,Nmin, Nmax);
			return();
		);

		print("Value D 	mean_number_primitive_solution");
		primAnsSum = 0;
		allAnsSum = 0;
		str_meanPrimCount = "";
		str_meanAllCount = "";

		for(N = Nmin, Nmax,
			if(N == 0, next);

			resAns = vector(3);
			resAns[1] = List();
			resAns[2] = List();
			resAns[3] = 0;

			ans = vector(3);
			ans[1] = List();
			ans[2] = List();
			ans[3] = 0;

			numberAllClasses = 0;

			forstep(n = sqrtint(abs(N)),2,-1,
				if (N % sqr(n) == 0,
					\\newN = N / sqr(n);
					newN = N \ sqr(n);
					ans = mainSolution(D, newN);
					for(i=1,length(ans[1]),
						listput(resAns[1], ans[1][i] * n);
						listput(resAns[2], ans[2][i] * n);
					);
					numberAllClasses += ans[3];
				);
			);

			ans = mainSolution(D, N);
			for(i = 1, length(ans[1]),
				listput(resAns[1], ans[1][i]);
				listput(resAns[2], ans[2][i]);
			);
			resAns[3] = ans[3];
			numberAllClasses += ans[3];


			if(length(resAns[1]) == 0,
				,
					countWS++;
			);

			if(numberAllClasses != 0,
				primAnsSum += resAns[3];
				allAnsSum += numberAllClasses;

				str_meanPrimCount = divWithRem(primAnsSum, countN);
				str_meanAllCount = divWithRem(allAnsSum, countN);
			);
			if(N % 100 == 0,
				printf("%3d	%s	%s\n", N, str_meanPrimCount, str_meanAllCount);
			);
			if (N == threshold_n,
				mean = primAnsSum / countN;

				filewrite(outFile, Str(D, "	", str_meanPrimCount, "	", str_meanAllCount));
				fileflush(outFile);

				if (mean >= threashold,
					outputAllInfo(D, Nmin, Nmax);
					break;
				);
				break;
			);
			countN++;

		);
}


outputPrimiviteInfo = (D, Nmin, Nmax) -> {
	my(f1,f2,f3, file1,file2,file3,countN, countWS, primAnsSum, str_meanPrimCount, resAns);
	\\f1 = Str(D, "p_ShareOfDecisions.txt");
	f2 = Str(D, "p_NumberOfClasses.txt");
	f3 = Str(D, "p_MeanNumberOfDecisions.txt");

	\\file1 = fileopen(f1, "w");
	file2 = fileopen(f2, "w");
	file3 = fileopen(f3, "w");

	\\filewrite(file1, "Value_D	Value_N	withSolution	All	numb_of_primitive");
	filewrite(file2, "Value_D	Value_N	numb_of_primitive");
	filewrite(file3, "D	Value_N	mean_numb_primitive_sol");

	countN = 1;
	countWS = 0;

	primAnsSum = 0;
	str_meanPrimCount = "";

	for(N = Nmin, Nmax,
		if(N == 0, next);

		resAns = vector(3);
		resAns[1] = List();
		resAns[2] = List();
		resAns[3] = 0;

		resAns = mainSolution(D, N);

		if(length(resAns[1]) == 0,
			,
				countWS++;
		);

		\\filewrite(file1, Str(D, "	", N, "	", countWS, "	", countN, "	", resAns[3]));

		if(resAns[3] != 0,
			filewrite(file2, Str(D, "	", N, "	", resAns[3]));

			primAnsSum += resAns[3];
			str_meanPrimCount = divWithRem(primAnsSum, countN);

			filewrite(file3, Str(D, "	", N, "	", str_meanPrimCount));
		);
		if(N % 100 == 0,
			printf("%3d	%s\n", N, str_meanPrimCount);
		);
		countN++;
	);
	\\fileclose(file1);
	fileclose(file2);
	fileclose(file3);

	\\print("Data is written to files:\n	"f1 "\n	"f2 "\n	"f3);
	print("Data is written to files:\n	"f2 "\n	"f3);
}

outputAllInfo = (D, Nmin, Nmax) -> {
	my(f1,f2,f3, file1,file2,file3,countN, countWS, primAnsSum, allAnsSum,  str_meanPrimCount, str_meanAllCount, resAns, ans, numberAllClasses, newN);
	\\f1 = Str(D, "_ShareOfDecisions.txt");
	f2 = Str(D, "_NumberOfClasses.txt");
	f3 = Str(D, "_MeanNumberOfDecisions.txt");

	\\file1 = fileopen(f1, "w");
	file2 = fileopen(f2, "w");
	file3 = fileopen(f3, "w");

	\\filewrite(file1, "Value_D	Value_N	withSolution	All	numb_of_primitive");
	filewrite(file2, "Value_D	Value_N	numb_of_primitive");
	filewrite(file3, "D	Value_N	mean_numb_primitive_sol");

	countN = 1;
	countWS = 0;

	primAnsSum = 0;
	str_meanPrimCount = "";

	for(N = Nmin, Nmax,
			if(N == 0, next);
			resAns = vector(3);
			resAns[1] = List();
			resAns[2] = List();
			resAns[3] = 0;

			ans = vector(3);
			ans[1] = List();
			ans[2] = List();
			ans[3] = 0;

			numberAllClasses = 0;

			forstep(n = sqrtint(abs(N)),2,-1,
				if (N % sqr(n) == 0,
					newN = N / sqr(n);
					ans = mainSolution(D, newN);
					for(i=1,length(ans[1]),
						listput(resAns[1], ans[1][i] * n);
						listput(resAns[2], ans[2][i] * n);
					);
					numberAllClasses += ans[3];
				);
			);

			ans = mainSolution(D, N);
			for(i = 1, length(ans[1]),
				listput(resAns[1], ans[1][i]);
				listput(resAns[2], ans[2][i]);
			);
			resAns[3] = ans[3];
			numberAllClasses += ans[3];


			if(length(resAns[1]) == 0,
				,
					countWS++;
			);

			\\filewrite(file1, Str(D, "	", N, "	", countWS, "	", countN, "	", resAns[3], "	", numberAllClasses));

			if(numberAllClasses != 0,
				filewrite(file2, Str(D,"	",N,"	", resAns[3],"	", numberAllClasses));

				primAnsSum += resAns[3];
				allAnsSum += numberAllClasses;

				str_meanPrimCount = divWithRem(primAnsSum, countN);
				str_meanAllCount = divWithRem(allAnsSum, countN);

				filewrite(file3, Str(D, "	",N,"	", str_meanPrimCount,"	", str_meanAllCount));
			);

			if(N % 100 == 0,
				printf("%3d	%s	%s\n", N, str_meanPrimCount, str_meanAllCount);
			);
			countN++;

		);

	\\fileclose(file1);
	fileclose(file2);
	fileclose(file3);

	\\print("Data is written to files:\n	"f1 "\n	"f2 "\n	"f3);
	print("Data is written to files:\n	"f2 "\n	"f3);
}


isOnlyPrimitiveSolution = () -> {
	my(isYes);
	while(isYes != "Y" && isYes != "y" && isYes != "N" && isYes != "n",
		print1("Calculate ONLY PRIMITIVE solution? Y - yes, N - no: ");
		isYes = Str(input());
	);
	if(isYes == "Y" || isYes == "y",
		return(1),
			return(0);
	);
};

number_of_solution_loop = (startD, endD, minN = 1, maxN = 1000, threshold = 0, thresholdN = 100) -> {
	my(outFileName, outfile);
	if(issquare(startD) == 1,
		print("ERROR! D should not be a perfect square!");
		return(-1);
	);
	default(timer, 1);
	if (isOnlyPrimitiveSolution() == 1,

		outFileName = Str("D_Means_values_", startD, "_", endD, "p.txt");
		outfile = fileopen(outFileName, "w");
		filewrite(outfile, "Value_D	mean_of_primvitive");
		forprime(D = startD, endD,
			calculationOnlyPrimitiveSolution(D, minN, maxN, threshold, thresholdN, outfile);
		),
			outFileName = Str("D_Means_values_", startD, "_", endD, ".txt");
			outfile = fileopen(outFileName, "w");
			filewrite(outfile, "Value_D	mean_of_primvitive	mean_of_all");

			forprime(D = startD, endD,
				calculationAllSolution(D, minN, maxN, threshold, thresholdN, outfile);
			);
	);

	fileclose(outfile);
	print("General data is written to file - "outFileName);
	print("******************end******************");
	print("time = "gettime()/1000.0" sec");
};


number_of_solution = (D, minN = 1, maxN = 1000) -> {
	if(issquare(D) == 1,
		print("ERROR! D should not be a perfect square!");
		return(-1);
	);
	threshold = 0;
	thresholdN = maxN;
	default(timer, 1);
	if (isOnlyPrimitiveSolution() == 1,
		calculationOnlyPrimitiveSolution(D, minN, maxN, threshold, thresholdN),
			calculationAllSolution(D, minN, maxN, threshold, thresholdN);
	);

	print("******************end******************");
	print("time = "gettime()/1000.0" sec");
};


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\




printProgress = (percentage) -> {
	my(val, lpad, rpad, PBSTR);
	val = percentage * 100;
	lpad = percentage * 60;
	rpad = 60 - lpad;
	PBSTR  = "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||";
	printf("%s%3d%% [%.*s%*s]", Strchr(13), val, lpad, PBSTR, rpad, "");
}


\\ only return period D
period = (D) -> {
	my(u,v,j,r,s);
	
	if (type(D) != "t_INT" || D < 2, return(-1));
	u = sqrtint(D); v = D-u^2;
	if (!v, return(0));
	s = v;
	r = u; j = 0;
	until (u == r && v == s,
		u = (r+u)\v * v - u;
		v = (D-u^2)\v; j++;
	); 
	return(j);
};

periodPrime = (D) -> {
	return(period(D));
}

periodAll = (D) -> {
	if(issquare(D),
		print("ERROR! D should not be a perfect square!");
		return(0);
	);
	return(period(D));
}


periodOfPrime = (startD, maxD, outFile) -> {
	my(count, len, endD, prevPersent, persent, denom);
	count = 0;
	len = 0;
	endD = 0;
	prevPersent = 0;
	\\ чтобы не учитывал последнее число, а считали его отдельно
	maxD = precprime(precprime(maxD)-1);
	denom = maxD - startD;

	forprime(D = startD, maxD,
		endD = D;
		len = periodPrime(D);

		filewrite(outFile, Str(D, "	", len));

		persent = (D - startD)*100 \ denom;
		if (prevPersent != persent,
			fileflush(outFile);
			printProgress(prevPersent / 100.0);
		);
		prevPersent = persent;
		count++;
	);
	\\ тут как-раз считает последнее отдельно, чтобы в файле не было переноса на последнюю строку
	len = periodPrime(nextprime(endD+1));
	filewrite1(outFile, Str(nextprime(endD+1), "	", len));
	printProgress(1);
}



periodOfAll = (startD, maxD, outFile) -> {
	my(count, len, endD, denom, prevPersent, persent);
	count = 0;
	len = 0;
	endD = 0;
	prevPersent = 0;
	denom = maxD - startD;
	printf("0%%");
	for(D = startD, maxD-2,
		endD = D;
		len = periodAll(D);

		if(len == 0,
			next;
		);

		filewrite(outFile, Str(D, "	", len));

		persent = (D - startD)*100 \ denom;
		if (prevPersent != persent,
			fileflush(outFile);
			printProgress(persent / 100.0);
		);
		prevPersent = persent;
		count = count + 1;
	);
	endD = endD + 1;
	len = periodAll(endD);
	if(len != 0,
		filewrite1(outFile, Str(endD, "	", len));
	);
	printProgress(1);
}


period_loop_prime = (D, maxD) -> {
	my(startD, fileName, file);
	default(timer, 1);
	\\ startTime = gettime();
	startD = D;
	fileName = Str("Period_length_prime" D "_" maxD ".txt");

	if(!isprime(D),
		D = nextprime(D);
	);

	file = fileopen(fileName, "w");
	filewrite(file, "Value_D	period_length");

	periodOfPrime(startD, maxD, file);

	fileclose(file);

	print("\nData is written to file - "fileName);
	\\ time = gettime() - startTime;
	print("time = "gettime()/1000.0" sec");

	print("******************end*******************");
}

period_loop_all = (D, maxD) -> {
	my(startD, fileName, file);
	default(timer, 1);
	startD = D;
	fileName = Str("Period_length_all" D "_" maxD ".txt");

	file = fileopen(fileName, "w");
	filewrite(file, "Value_D	period_length");

	periodOfAll(startD, maxD, file);

	fileclose(file);

	print("\nData is written to file - "fileName);
	print("time = "gettime()/1000.0" sec");

	print("******************end*******************");
}

substr(s,n = 1,m = #s) = {
	my(t);
	if (m > #s, m = #s, if (m < 1, m = #s+m));
	if (abs(n) > #s,
		n = 1,
			if (n < 0, n = #s+n+1);
	);
	if (n > m && ((n > 0 && m > 0) || (n < 0 && m < 0)), t = n; n = m; m = t);
	return(concat(Vec(s)[n..m]));
}

str_split_by_char(str, symbol) = {
	my(v);

	v = Vec(str);
	for (i = 1, #v,
		if (v[i] == symbol,
			return([substr(str, 1, i-1), substr(str, i+1, #str)]);
		);
	);
	return([str, ""]);
}

number_of_d_4k_plus_i(fileName, fileSize) = {
	my(n, outFileName0, outFileName1, outFileName2, outFileName3,
		out0, out1, out2, out3, firstStr, secondStr, firstInt, secondInt, count);

	n = fileopen(fileName, "r");

	outFileName0 = concat(fileName, "_4k.txt");
	outFileName1 = concat(fileName, "_4k+1.txt");
	outFileName2 = concat(fileName, "_4k+2.txt");
	outFileName3 = concat(fileName, "_4k+3.txt");

	out0 = fileopen(outFileName0, "w");
	out1 = fileopen(outFileName1, "w");
	out2 = fileopen(outFileName2, "w");
	out3 = fileopen(outFileName3, "w");

	filewrite(out0, "l=4k	numberOfD");
	filewrite(out1, "l=4k+1	numberOfD");
	filewrite(out2, "l=4k+2	numberOfD");
	filewrite(out3, "l=4k+3	numberOfD");

	print("\nCalculating number of D with period length l=4k,4k+1,4k+2,4k+3...");	

	count = 0;
	while(line = filereadstr(n),
		print1(Strchr(13));
		[firstStr, secondStr] = str_split_by_char(line, "\t");
		iferr(firstInt = eval(firstStr); secondInt = eval(secondStr), E, next);

		if (type(firstInt) == "t_INT" && type(secondInt) == "t_INT",
			count++;
			if (count % 10 == 0,
				if (fileSize,
					printProgress(count / fileSize * 1.0),
						print1(firstStr" - " secondStr));
			);
			if (firstInt % 4 == 0 && firstInt != 0,
				filewrite(out0, firstInt"	"secondInt);
				next;
			);
			if ((firstInt - 1) % 4 == 0,
				filewrite(out1, firstInt"	"secondInt);
				next;
			);
			if ((firstInt - 2) % 4 == 0,
				filewrite(out2, firstInt"	"secondInt);
				next;
			);
			if ((firstInt - 3) % 4 == 0,
				filewrite(out3, firstInt"	"secondInt);
				next;
			);
		);
	);

	printProgress(1);
	fileclose(n);
	fileclose(out0);
	fileclose(out1);
	fileclose(out2);
	fileclose(out3);

	print("\nData is writen to the files:\n "outFileName0"\n "outFileName1"\n "outFileName2"\n "outFileName3);
}


number_of_d_with_period_length(fileName) = {
	my(d, period_length, n, out, map, intD, intPeriodLength, count, fileSize);

	n = fileopen(fileName, "r");

	outFileName = concat(fileName, "_number_of_d_with_period_length.txt");
	out = fileopen(outFileName, "w");
	filewrite(out, "period_length	numberOfD");

	map = Map();

	count = 0;
	while(line = filereadstr(n),
		print1(Strchr(13));
		[d, periodLength] = str_split_by_char(line, "\t");
		iferr(intD = eval(d); intPeriodLength = eval(periodLength), E, next);

		if (type(intD) == "t_INT" && type(intPeriodLength) == "t_INT",
			if(mapisdefined(map, intPeriodLength),
				mapput(map, intPeriodLength, mapget(map, intPeriodLength)+1);
				if ((count < 1000) || (count > 1000 && (count % 500) == 0),
					print1("d = "intD"	period_length = "intPeriodLength);
				),
					mapput(map, intPeriodLength, 1);
					if ((count < 1000) || (count > 1000 && (count % 500) == 0),
						print1("d = "intD"	period_length = "intPeriodLength);
					);
			);
			count++;
		);
	);

	print("\nWriting to a file...");
	matMap = Mat(map);
	fileSize = #matMap[,1];
	for( i = 1, fileSize,
		filewrite(out, matMap[i,1]"	"matMap[i,2]);
		printProgress(i / fileSize);
	);
	printProgress(1);

	fileclose(n);
	fileclose(out);

	print("\nData is writen to the file - "outFileName);
	number_of_d_4k_plus_i(outFileName, fileSize);
}

addhelp(number_of_solution, " * number_of_solution(d, min_n, max_n) - Finds number of solution of the generalized Pell equation x^2-Dy^2=N with fixed D & changing N from min_n to max_n. Outputs the average number of solutions, the number of equations with solutions & the total number of equation.\n\n  	Input:	- d		- the value D of equation x^2-Dy^2=N.\n 		- min_n		- {1 by default} minimum value n on interval min_n <= N <= max_n (without N == 0) \n		- max_n		- {1000 by default} maximum value n on interval min_n <= N <= max_n (without N == 0)");
addhelp(number_of_solution_loop, " * number_of_solution_loop(start_d, end_d, min_n, max_n, threshold, threshold_n) - Finds number of solution of the generalized Pell equation on the interval start_d <= D <= end_d with fixed D & changing N from min_n to max_n. Outputs the average number of solutions, the number of equations with solutions & the total number of equation.\n\n  	Input:	- start_d	- the value of d that starts the calculation.\n 		- end_d		- the value of d at which the calculation ends.\n 		- min_n		- minimum value n on interval min_n <= N <= max_n (without N == 0) \n		- max_n		- maximum value n on interval min_n <= N <= max_n (without N == 0)\n		- threshold	- minimum bound of the mean of primitive solutions \n 		- threshold_n	- value at which the threshold is checked");
addhelp(period_loop_prime, " * period_loop_prime(start_d, end_d) - Finds period length for prime D on the interval start_d <= D <= end_d.\n 	Input:	- start_d	- positive integer\n 		- end_d 		- positive integer. \n 	The results are written to a file - Period_length_prime<start_d>_<end_d>.txt");
addhelp(period_loop_all, " * period_loop_all(start_d, end_d) - Finds period length for all D on the interval start_d <= D <= end_d.\n 	Input:	- start_d 	- positive integer\n 		- end_d 		- positive integer \n 	The results are written to a file - Period_length_all<start_d>_<end_d>.txt");
addhelp(number_of_d_with_period_length, "\n * number_of_d_with_period_length(fileName) - Calculates the number of d < D with a fixed period length l.\n 	Input: fileName - string name of the output file of the function period_loop_prime(...)\n 			or period_loop_all(...).\n That is, you need to first run the function 'period_loop_prime(...)' (or 'period_loop_all(...)'), and pass the resulting file to the function 'number_of_d_with_period_length(...)'.\n E.g: number_of_d_with_period_length(\"Period_length_all100_110.txt\")");


print("		*************************************************************************	");
print("		*	Number of solution of the generalized Pell equation  x^2-Dy^2=N	*\n		*		on the interval	startD <= D <= endD			*");
print("		*		with fixed D & changing N from 1 to maxN		*\n		*************************************************************************	");
print("		*	Search a period length on the interval startD <= D <= endD 	*\n		* 			for prime & all numbers				*");
print("		*************************************************************************	");
print("This package includes the function:\n 	* number_of_solution(d, min_n, max_n)	- Finds number of solution of the generalized\n 				 Pell equation with fixed D & changing N from min_n to max_n.");
print1(" 	* number_of_solution_loop(start_d, end_d, min_n, max_n, threshold, threshold_n)	- Finds number of solution\n 				of the generalized Pell equation on the interval start_d <= D <= end_d with\n 				fixed D & changing N from min_n to max_n.");



print1("\n 	* period_loop_prime(startD, endD)	- Finds period length for prime D on the interval startD <= D <= endD.\n	* period_loop_all(startD, endD)	- Finds period length for all D on the interval startD <= D <= endD.");

print("\n 	* number_of_d_with_period_length(fileName) - Calculates the number of d < D with a fixed period length l.\n 			Input: fileName - string name of the output file of the function period_loop_prime(...)\n 			or period_loop_all(...).");
print("\n 	* number_of_d_4k_plus_i(fileName) - Divides the source file into 4 parts, by the value of \n 					the first column value = 4k+i, i=[0,3].\n 			Input: File with a table of the form: value - count of values.")

print("______________________________________________________________________________________________________\n");

print("\n ---- For more information about the functions, use help: ? <method_name>");