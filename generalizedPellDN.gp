\\ !!!!!! REMOVED ALL UNNECESSARY FROM VERSION 1.2.2
\\ and using method searchAllCombinations3 instead searchAllCombinations2
\\ asks for confirmation of writing to the file
isWriteInFile = () -> {
	my(isYes);
	while(isYes != "Y" && isYes != "y" && isYes != "N" && isYes != "n",
		print1("Write data to the file? Y - yes, N - no: ");
		isYes = Str(input());
	);
	if(isYes == "Y" || isYes == "y",
		return(1),
			return(0);
	);
	
}

\\ *************************************************************
\\ *** all solutions of quadratic congruence x^2 = D (mod N) ***

chineseForArray = (mass, i = 1) -> {
	if(i == length(mass)-1,
		return(chinese(mass[i],mass[i+1])),
			return(chinese(mass[i], chineseForArray(mass, i+1)));
	);
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

printt = (a, size) -> {

	for(i = 1, size,
		print1(a[i]" ");
	);
	print("");
}

formatToCorrectView = (masX, a) -> {
	for(i = 1,#a,
		a[i] = a[i] % #masX[i];
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


quadraticCongruenceForCoPrime = (D, N, isPrint = 0) -> {
	my(f,lengthF,masX,allPosSol,answers,allCombinations);

	f = factorint(abs(N));
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

	if(isPrint == 1,
		print1(" Congruence:  ");
		for(i=1,#answers, print1(answers[i]", "));
		print("");
	);
	return(answers);
}

quadraticCongruenceForNotCoPrime3 = (D, N, isPrint = 0) -> {
	my(answer, answers,n,f,r,pn,m,tmpRes, x, p);
	\\print("quadraticCongruenceForNotCoPrime3");

	answers = List();
	fact = factorint(abs(N));
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
							print("No solutions");
							return(0),
								m = r \ 2;
								if (r != 0 && m == 0, 
									m = r / 2.;
								);
								\\	print("m = "m);
								tmpRes = List();
								tmpRes = quadraticCongruenceForCoPrime(D \ p^r, p^(n-r));
								if(tmpRes == 0,
									print("No solutions");
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
	\\	print("answers = "answers);
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

	if (isPrint == 1,
		print1(" Congruences: ");
		for(i=1,#res, print1(res[i]", "));
		print("");
	);

	return(res);
}

\\ ***************************************************************
\\ SEARCHING PERIOD OF NUMBER AND SOLUTION PELL EQUATION x^2-Dy^2=+-1

\\ only return period D
period = (D) -> {
	default(timer, 1);
	if (type(D) != "t_INT" || D < 1, 
		print("  ***   The number must be positive integer. (greater than 0)");
		return(-1);
	);


	my(d,p1,p0,pii,q1,q0,qi,count,l);

	d = sqrtint(D);
	if (sqr(d) == D,
		print("  ***   The number is a perfect square!");
		return(-1);
	);
	p0 = 0; q0 = 1;
	q0 = ((p0+d)\q0);
	
	qi = q0;
	p1= 0; q1 = 1;

	count = 0;
	while(qi != 2*q0,
		qi = ((p1 + d) \ q1);

		pii = qi * q1 - p1;
		qii = (D - sqr(pii))\q1;

		p1 = pii;
		q1 = qii;

		count = count + 1;
	);
	l = count - 1;
	print("time = "gettime()/1000.0" sec");
	return(l);
};


\\ with answer X Y
pell = (D) -> {
	if (type(D) != "t_INT" || D < 1, 
		print("  ***   The number must be positive integer. (greater than 0)");
		return(-1);
	);

	my(d,a2,a1,ai,b2,b1,bi,p1,p0,pii,q1,q0,qi,count,l, x, y, xx, yy);

	d = sqrtint(D);
	if (sqr(d) == D,
		print("  ***   The number is a perfect square!");
		return(-1);
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

	print("Period length = "l);

	x = a2; y = b2;
	print("Solution x^2-Dy^2="x*x - D*y*y":\n x = "x "\n y = "y);
	if(x*x - D*y*y == -1,
		xx = x*x+D*y*y;
		yy = 2*x*y;
		print1("Solution x^2-Dy^2=1:\n x = "xx "\n y = "yy"\n");
		write(fileName, "Solution x^2-Dy^2=1:\n x = "xx "\n y = "yy"\n");
	);

	if(isWriteInFile() == 1,
		fileName = Str("pell_D="D".txt");
		write(fileName, "D = ", D);
		write(fileName, "Period length = ", l);
	
		write(fileName, "Solution x^2-Dy^2=-1:\n x = "x "\n y = "y);

		if(x*x - D*y*y == -1,
			write(fileName, "Solution x^2-Dy^2=1:\n x = "xx "\n y = "yy"\n");
		);
		print("Data is written to file - "fileName);
	);
	print("******************end*******************");
};

\\ with answer X Y
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

	\\print("Period length = "l);

	x = a2; y = b2;
	\\print("Solution x^2-Dy^2="x*x - D*y*y":\n x = "x "\n y = "y);

	if(x*x - D*y*y == -1,
		xx = x*x+D*y*y;
		yy = 2*x*y;
		\\print1("Solution x^2-Dy^2=1:\n x = "xx "\n y = "yy"\n");
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


\\ *************** GENERALIZED PELL ********************
checkFullSquare = (D) -> {
	my(d);
	d = sqrtint(abs(D));
	if (sqr(d) == abs(D),
		\\ print("The number is a perfect square!");
		return(1);
	);
	return(0);
};

checkSolutionExist = (D,N) -> {
	if (N == -3,
		if (checkFullSquare(D - 2),
			print("No solution. [1]");
			return(0);
		);
	);
	if (D == (N*N) - 1,
		if (checkFullSquare(N),
			return(1),
				print("No solution. [2]");
				return(0);
		);
	);
	return(1);
};

checkEquivalence = (res, D, N) -> {
	my(count, i, j);
	count = 0;
	for(i = 1, length(res[1]),
		for(j = i + 1, length(res[2]),
			if(((res[1][i] * res[1][j] - D * res[2][i] * res[2][j]) % N == 0) && ((res[2][i]*res[1][j] - res[1][i]*res[2][j]) % N == 0),
				count = count + 1;
				print("! (" res[1][i]"; " res[2][i] ") ~~~ ( " res[1][j] "; " res[2][j] ")");
			);
		);
	);
	print(" Equivalence solution: " count);
}


allCongruences = (D, N) -> {
	my(p0Vec,gcdDN);
	p0Vec = List();
	gcdDN = gcd(D,N);
	if (N == 1,
		listput(p0Vec, 0),
			if (gcdDN == 1,
				p0Vec = quadraticCongruenceForCoPrime(D, N),
					p0Vec = quadraticCongruenceForNotCoPrime3(D, N);
			);
	);
	return(p0Vec);
}

allCongruencesTime = (D,N) -> {
	my(all);
	default(timer, 1);
	all = allCongruences(D,N);
	print("time = "gettime()/1000.0" sec");
	return(all);
}


mainSolution = (D,N) -> {
	my(l, QQ0, res, p0Vec, countRatio, modd, PPi, qqi, int_part_sqrtD,
	vecA, vecB, vecP, vecQ, AAi2, AAi1, BBi2, BBi1, Gt1, Bt1, 
	begin_P, begin_Q, GGi, AAi, BBi, count, start_index, t, period_begin, r, m);
	\\print("-- N = "N" --------");

	l = 0;
	\\ QQ0 = N;

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
				p0Vec = quadraticCongruenceForCoPrime(D, N),
					p0Vec = quadraticCongruenceForNotCoPrime3(D, N);
			);
	);
	
	if(length(p0Vec) == 0,return(res));

	if (abs(N) == 1,
		if (N == 1,
			res = pellCalcFor(D),
				res = pellCalcFor(D , -1)),

	\\if(N == -1, 
	\\	res = pellCalcFor(D , -1); print(res),
			for(j=1,length(p0Vec),
				vecA = List(); vecB = List();
				vecP = List(); vecQ = List();

				QQi = N;
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
						Gt1 = N * AAi1 - p0Vec[j] * BBi1;
						\\ Gt1 = N * AAi1 - j * BBi1;
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
								resG = N * vecA[l*m + t] - p0Vec[j] * vecB[l*m + t];
								\\ resG = N * vecA[l*m + t] - j * vecB[l*m + t];
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
								resG = N * vecA[l * m + t] - p0Vec[j] * vecB[l*m + t];
								\\ resG = N * vecA[l* m + t] - j * vecB[l*m + t];

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

generalizedPell = (D,N) -> {
	if (type(D) != "t_INT" || D < 1, 
		print("  ***   The number D must be positive integer. (greater than 0)");
		return(-1);
	);

	if (checkFullSquare(D),
		print("  ***   The number is a perfect square!");return(-1)
	);
	if (N == 0,
		print("  ***   The number N must not be zero!");return(-1));

	default(timer, 1);

	my(numberAllClasses, ans, resAns);

	numberAllClasses = 0;
	ans = vector(3);
	resAns = vector(3);
	resAns[1] = List();
	resAns[2] = List();

	forstep(n = sqrtint(abs(N)), 2, -1, 
		if(N % sqr(n) == 0,
			newN = N \ sqr(n);
			ans = mainSolution(D, newN);
			\\print("ans = " ans);

			for(i = 1,length(ans[1]),
				listput(resAns[1], ans[1][i] * n);
				listput(resAns[2], ans[2][i] * n);
			);
			numberAllClasses = numberAllClasses + ans[3];

		);
	);

	ans = mainSolution(D, N);
	\\print("ans = " ans);
	for(i = 1,length(ans[1]), 
		listput(resAns[1], ans[1][i]);
		listput(resAns[2], ans[2][i]);
	);

	resAns[3] = ans[3];
	numberAllClasses = numberAllClasses + ans[3];
	print("time = "gettime()/1000.0" sec");

	if((length(resAns[1]) == 0), 
		print("No solution!");
		return(-1),
			checkEquivalence(resAns, D, N);
			print("_______________________________________");
			print(" Number of primitive solution = " resAns[3]);
			print(" Number of all solution = " numberAllClasses);

			noPrim = numberAllClasses - resAns[3];
			print(" Primitive solutions:");
			for(i = noPrim + 1, length(resAns[1]),
				print(" 	(" resAns[1][i] "; " resAns[2][i] ")");
			);
			if(noPrim != 0,
				print(" Other solutions: ");
				for(i = 1, noPrim,
					print(" 	(" resAns[1][i] "; " resAns[2][i] ") ");
				);
			);
			print("_______________________________________");
			if(isWriteInFile() == 1,
				fileName = Str("generalizedPellD"D"_N"N".txt");
				write(fileName, "x^2-", D, "y^2=",N);
				write(fileName, " Number of primitive solution = ", resAns[3]);
				write(fileName, " Number of all solution = ", numberAllClasses);
				write(fileName, "Primitive solutions:");
				for(i = noPrim + 1, length(resAns[1]),
					write(fileName, " 	(" resAns[1][i] "; " resAns[2][i] ")");
				);
				if(noPrim != 0,
					write(fileName, " Other solutions: ");
					for(i = 1, noPrim,
						write(fileName, " 	(" resAns[1][i] "; " resAns[2][i] ") ");
					);
				);
				print("Data is written to file - "fileName);
			);
			
	);
	print("******************end*******************");
}


addhelp(period, " * period(x) - Return the length of the period of the continued fraction for sqrt(x).\n 	Input: x - positive integer that is not a perfect square.\n 	Output: period of a positive integer x");
addhelp(pell, " * pell(D) - Finds period length of the continued fraction for sqrt(D) & solution x,y of the equation x^2+Dy^2=+-1. \n If a solution to equation x^2+Dy^2=-1 is found, it also returns equation x^2+Dy^2=1.\n It is possible to write the results to a file.\n 	Input: D - positive integer that is not a perfect square.");
addhelp(generalizedPell," * generalizedPell(D,N) - Finds a solution to the generalized Pell equation x^2-Dy^2=N.\n  Returns the number of primitive solution, the number of all solution, primitive solutions, non-primitive solutions.\n  It is possible to write the results to a file.\n 	Input:	D - positive integer that is not a perfect square.\n 	 	N - integer.");
addhelp(allCongruences, " * allCongruences(D,N) - Return all solutions of the congruence x^2 = D (mod N).\n 	Input:	D - positive integer.\n 		N - integer.")
addhelp(allCongruencesTime, " * allCongruencesTime(D,N) - like allCongruences(D,N), but also outputs the run time.");

print("		*	Pell equation of the form	*\n		*	x^2-Dy^2=1 & x^2-Dy^2=-1	*\n		* 		&			*")
print("		*	Generalized Pell equation	*\n		*		x^2-Dy^2=N 		*\n		*****************************************	")
print1("This package includes the functions:\n 	*	period(x)	-	Return the length of the period of the continued fraction for sqrt(x).\n	*	pell(_D)	-	Prints period length & solution x,y of the equation x^2+Dy^2=+-1\n 	*	generalizedPell(D,N)	-	Finds a solution to the generalized Pell equation x^2-Dy^2=N.\n	*	allCongruences(D,N) - Return all solutions of the congruence x^2 = D (mod N).\n	*	allCongruencesTime(D,N) - like allCongruences(D,N), but also outputs the run time.");

print("\n---- For more information about the methods, use help: ? <method_name>");


