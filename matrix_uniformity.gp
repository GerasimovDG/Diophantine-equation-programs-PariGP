printProgress = (percentage, d) -> {
	my(val, lpad, rpad, PBSTR);
	val = percentage * 100;
	lpad = percentage * 60;
	rpad = 60 - lpad;
	PBSTR  = "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||";
	if (d,
		printf("%s%3d%% [%.*s%*s] d = %s", Strchr(13), val, lpad, PBSTR, rpad, "", d),
			printf("%s%3d%% [%.*s%*s]", Strchr(13), val, lpad, PBSTR, rpad, "");
	);
}

checkFullSquare = (D) -> {
	my(d);
	d = sqrtint(abs(D));
	if (sqr(d) == abs(D),
		\\ print("The number is a perfect square!");
		return(1);
	);
	return(0);
};

confirm = (message) -> {
	my(isYes);
	while(isYes != "Y" && isYes != "y" && isYes != "N" && isYes != "n",
		print1(message" Y - yes, N - no: ");
		isYes = Str(input());
	);
	if(isYes == "Y" || isYes == "y",
		return(1),
			return(0);
	);
	
}

searchTetta = (a, b ,c ,n) -> {
	my(quadraticEquation);
	tettaList = List();

	for (tetta = 0, abs(n)-1,
		quadraticEquation = a* sqr(tetta) + b* tetta + c;
		if (quadraticEquation % abs(n) == 0,
			listput(tettaList, tetta);
		);
	);
	return(tettaList);
}

validationInputData = (a, b, c) -> {

	my(discriminant);
	discriminant = sqr(b) - 4*a*c;

	if (discriminant < 0,
		print("Error! Negative discriminant");
		return(0);
	);
	if (sqr(sqrtint(discriminant)) == discriminant,
		print("Error! Discriminant are complete square!");
		return(0);
	);
	return(1);
}

PQ_evenDelta = (D, P, Q, N, tetta, ~result, isMain) -> {
	my(vecA, vecB, vecP, vecQ, resMas, res, Qi, Pi, qi, int_part_sqrtD,
	Ai2, Ai1, Bi2, Bi1, Ai, Bi, begin_P, begin_Q, count, t, l,
	start_index, period_begin, second_period, p, q, i, j, mOne,
	X, y, Y, x);

	if (Q == 0,
		print("Oops, something went wrong. Qi == 0");
		return;
	);

	if ((D - sqr(P)) % Q != 0,
		D = D * sqr(Q);
		P = P * abs(Q);
		Q = Q * abs(Q);
	);

	vecA = List(); vecB = List(); vecP = List(); vecQ = List();
	resMas = List();

	Qi = Q;
	Pi = P;

	qi = 0;
	int_part_sqrtD = 0;

	Ai2 = 0; Bi2 = 1;
	Ai1 = 1; Bi1 = 0;

	Ai = 0; Bi = 0;
	count = 1;
	t = 0;
	l = 0;

	start_index = 0;
	period_begin = 0;
	second_period = 0;

	while(1,
		if (Qi == 0,
			print("Oops, something went wrong. Qi == 0.");
			return;
		);
		if (Qi > 0,
			qi = divrem((Pi + sqrtint(D)), Qi)[1],
				if ( Qi < 0,
					p = divrem((Pi + 1 + sqrtint(D)), Qi)[1];
					q = divrem((Pi + 1 + sqrtint(D)), Qi)[2];
					if (q == 0,
						qi = p,
							qi = p-1;
					);
				);
		);

		int_part_sqrtD = sqrtint(D);
		i = count - 1;

		mOne = -1;
		j = i;
		if (isMain != 1,
			j = i + 1;
		);

		if (Qi == mOne^j * N \ abs(N),
			X = Ai;
			y = Bi;

			Y = y;
			x = y * tetta + abs(N) * X;
			res = List([x, Y]);
			listput(~result, res);
		);

		Pi = qi * Qi - Pi;
		Qi = (D - sqr(Pi)) \ Qi;

		Ai = qi * Ai1 + Ai2;
		Bi = qi * Bi1 + Bi2;

		listput(vecA, Ai);
		listput(vecB, Bi);

		Ai2 = Ai1;
		Bi2 = Bi1;
		Ai1 = Ai;
		Bi1 = Bi;

		listput(vecP, Pi);
		listput(vecQ, Qi);

		if (Qi == 1 && t == 0,
			t = count;
		);

		if ((period_begin == 1) && (begin_P == Pi && begin_Q == Qi),
			\\ если длина периода нечетная, то смотрим еще второй период
			if (second_period == 1,
				break;
			);
			l = count - start_index;
			if (l % 2 != 0,
				second_period = 1,
					break;
			);
		);


		if (Qi >= 0 && period_begin == 0,
			if (((Qi - Pi) <= int_part_sqrtD) && (Pi <= int_part_sqrtD) && ((Pi + Qi) > int_part_sqrtD),
				period_begin = 1;
				begin_P = Pi;
				begin_Q = Qi;

				start_index = count;
			);
		);
		count++;
	);
}


PQ_oddDelta = (D, P, Q, N, tetta, ~result, isMain, A) -> {
	my(vecA, vecB, vecP, vecQ, resMas, res, Qi, Pi, qi, int_part_sqrtD,
	Ai2, Ai1, Bi2, Bi1, Ai, Bi, begin_P, begin_Q, count, t, l,
	start_index, period_begin, second_period, p, q, i, j, mOne,
	X, y, Y, x);

	if (Q == 0,
		print("Oops, something went wrong. Qi == 0");
		return;
	);

	if ((D - sqr(P)) % Q != 0,
		D = D * sqr(Q);
		P = P * abs(Q);
		Q = Q * abs(Q);
	);

	vecA = List(); vecB = List(); vecP = List(); vecQ = List();
	resMas = List();

	Qi = Q;
	Pi = P;

	qi = 0;
	int_part_sqrtD = 0;

	Ai2 = 0; Bi2 = 1;
	Ai1 = 1; Bi1 = 0;

	Ai = 0; Bi = 0;
	count = 1;
	t = 0;
	l = 0;

	start_index = 0;
	period_begin = 0;
	second_period = 0;

	while(1,
		if (Qi == 0,
			print("Oops, something went wrong. Qi == 0.");
			return;
		);
		if (Qi > 0,
			qi = divrem((Pi + sqrtint(D)), Qi)[1],
				if ( Qi < 0,
					p = divrem((Pi + 1 + sqrtint(D)), Qi)[1];
					q = divrem((Pi + 1 + sqrtint(D)), Qi)[2];
					if (q == 0,
						qi = p,
							qi = p-1;
					);
				);
		);

		int_part_sqrtD = sqrtint(D);
		i = count - 1;

		mOne = -1;
		j = i;
		if (isMain != 1,
			j = i + 1;
		);
		if (Qi == mOne^j * 2 * N \ abs(N),
			X = Ai;
			y = Bi;

			Y = y;
			x = y * tetta + abs(N) * X;

			res = List([x, Y]);
			listput(result, res);

			if (D == 5 && (A * N < 0) && isMain == 1,
				X = Ai - Ai2;
				y = Bi - Bi2;

				Y = y;
				x = y * tetta + abs(N) * X;

				res[1] = x;
				res[2] = Y;
				listput(result, res);
			);
		);

		Pi = qi * Qi - Pi;
		Qi = (D - sqr(Pi)) \ Qi;

		Ai = qi * Ai1 + Ai2;
		Bi = qi * Bi1 + Bi2;

		listput(vecA, Ai);
		listput(vecB, Bi);

		Ai2 = Ai1;
		Bi2 = Bi1;
		Ai1 = Ai;
		Bi1 = Bi;

		listput(vecP, Pi);
		listput(vecQ, Qi);

		if (Qi == 1 && t == 0,
			t = count;
		);

		if ((period_begin == 1) && (begin_P == Pi && begin_Q == Qi),
			\\ если длина периода нечетная, то смотрим еще второй период
			if (second_period == 1,
				break;
			);
			l = count - start_index;
			if (l % 2 != 0,
				second_period = 1,
					break;
			);
		);

		if (Qi >= 0 && period_begin == 0,
			if (((Qi - Pi) <= int_part_sqrtD) && (Pi <= int_part_sqrtD) && ((Pi + Qi) > int_part_sqrtD),
				period_begin = 1;
				begin_P = Pi;
				begin_Q = Qi;

				start_index = count;
			);
		);
		count++;
	);
}

choosingBestSolution = (result) -> {
	my(resultEnd, minY, X, solution);

	if (#result != 0, 
		resultEnd = List();

		for( i = 1, #result,
			if (!(result[i][1] == 0 && result[i][2] == 0),
				listput(resultEnd, result[i]);
			);
		);

		minY = resultEnd[1][2];	\\ .Y
		X = resultEnd[1][1];	\\ .X

		\\ print("Answers (of these, we choose the one with the smallest y):");

		for (i = 1, #resultEnd,
			if (abs(resultEnd[i][2]) < abs(minY),
				minY = resultEnd[i][2];
				X = resultEnd[i][1];
			);
			if (resultEnd[i][2] == minY && resultEnd[i][1] < X,
				minY = resultEnd[i][2];
				X = resultEnd[i][1];
			);
		);

		\\ solution = List([X, minY]);
		solution = [X, minY];
		return(solution);
	);
	\\return(List());
}

SearchSolution = (tetta, A, B, C ,N) -> {
	my(delta,n,R,S, res, solution);

	delta = sqr(B) - 4 * A * C;
	n = 2 * A * tetta + B;

	R = n \ 2;

	S = A * abs(N);
	result = List();
	if (delta % 2 == 0,
		PQ_evenDelta((delta \ 4), -R, S, N, tetta, ~result, 1);
		PQ_evenDelta((delta \ 4), R, -S, N, tetta, ~result, 0),
			if (delta == 5 && (A * N == -1),
				\\main
				res = List([(-B + 1) \ 2,1]);
				listput(result, res);

				\\Special
				res = List([(-B - 1) \ 2, 1]);
				listput(result, res);

				\\!main
				res = List([(-B - 2), 2]);
				listput(result, res),
					PQ_oddDelta(delta, -(2 * R + 1), 2 * S, N, tetta, ~result, 1, A);
					PQ_oddDelta(delta, 2*R+1, -2*S, N, tetta, ~result, 0, A);
			);
	);

	if (#result != 0,
		solution = choosingBestSolution(result);
	);
	return(solution);
}

transformationABC = (a, b, c, N) -> {

	my(firstA, firstB, firstC, max, x, y, Nq, setQ, factorize, 
		X, Y, Ni, u, v, y2, y1, alpha, gamma, betta, delta, A, B, C);

	firstA = a;
	firstB = b;
	firstC = c;

	Nq = 1;

	factorize = factor(abs(N))[,1];
	setQ = List(factorize);

	for( i = 1, #setQ,
		Nq *= setQ[i];
	);

	X = 0;
	Y = 0;
	for (i = 1, #setQ,
		if(gcd(firstA, setQ[i]) == 1,
			x = 1;
			y = 0,
				if (gcd(firstC, setQ[i]) == 1,
					x = 0;
					y = 1,
						if (gcd(firstA, setQ[i]) != 1 && gcd(firstC, setQ[i]) != 1,
							x = 1;
							y = 1;
						);
				);
		);
		Ni = Nq \ setQ[i];

		for (j = 1, setQ[i],
			if ((Ni * j) % setQ[i] == x % setQ[i],
				y1 = j;
				break;
			);
		);
		X += Ni * y1;

		for (j = 1, setQ[i],
			if (((Ni * j) % setQ[i] == y % setQ[i]),
				y2 = j;
				break;
			);
		);
		Y += Ni * y2;
	);

	X = X % Nq;
	Y = Y % Nq;

	alpha = X \ gcd(X,Y);
	gamma = Y \ gcd(X,Y);

	A = a * sqr(alpha) + b * alpha * gamma + c * sqr(gamma);

	uvd = gcdext(alpha, -gamma);
	delta = uvd[1];
	betta = uvd[2];

	B = 2 * a * alpha * betta + b * alpha * delta + b * betta * gamma + 2 * c * gamma * delta;
	C = a * sqr(betta) + b * betta * delta + c * sqr(delta);

	a = A;
	b = B;
	c = C;

	if (gcd(a, gcd(b,c)) != 1,
		return;
		\\ return(List([a,b,c, alpha, betta, gamma, delta]));
	);
	return(List([a,b,c, alpha, betta, gamma, delta]));
}

mainSolution = (A, B, C ,N) -> {
	my(alpha, betta, gamma, deltaa, gcd_flag, result, final_res, solution,
		D, Delta, Q, P, tettaMas, x, y, newN);

	alpha = 0; betta = 0; gamma = 0; deltaa = 0;
	gcd_flag = 0;
	if (gcd(A, N) > 1,
		[A,B,C, alpha, betta, gamma, deltaa] = transformationABC(A, B, C, N);

		\\print("-------------\n|	alpha	betta	|\n|	gamma	delta	|:\n ------------");
		\\print("|	"alpha"	"betta"	|	\n|	"gamma"	"deltaa"	|");
		\\print("-------------");

		gcd_flag = 1;

		\\print("new A,B,C: "A" "B" " C);
	);

	D = sqr(B) - 4*A*C;
	Delta = D \ 4;
	Q = A * abs(N);

	final_res = List();
	tettaMas = searchTetta(A,B,C,N);
	if (#tettaMas == 0,
		\\ print("	No suitable theta.");
		return(final_res);
	);

	for(i = 1, #tettaMas,
		solution = SearchSolution(tettaMas[i], A, B, C, N);
		if (solution && !(solution[1] == 0 && solution[2] == 0),
			listput(final_res, solution);
		);
	);

	if (gcd_flag == 1,
		for( i = 1, #final_res,
			x = final_res[i][1];
			y = final_res[i][2];
			final_res[i][1] = alpha * x + betta * y;
			final_res[i][2] = gamma * x + deltaa * y;
		);
	);
	return(final_res);
}

ABC_N = (A,B,C,N) -> {
	my(final_res, add_res, final_add_res);
	if (!validationInputData(A,B,C),
		print("The input data did not pass the validation check!");
		return(0);
	);

	final_add_res = List();
	forstep( n = sqrtint(abs(N)), 1, -1,
		if (N % sqr(n) == 0,
			newN = N \ sqr(n);
			if (newN != N,
				add_res = mainSolution(A,B,C, newN);
				for (i = 1, #add_res,
					add_res[i][1] = add_res[i][1] * n;
					add_res[i][2] = add_res[i][2] * n;

					listput(final_add_res, add_res[i]);
				),
					final_res = mainSolution(A,B,C,N);
			);
		);
	);

	if (#final_res == 0,
		\\print("No primitive solutions found."),
			\\ print("primitive solutions:");
			\\for (i = 1, #final_res,
			\\	print ("(X,Y) = ("final_res[i][1]", "final_res[i][2]")");
			\\);
	);

	if (#final_add_res == 0,
		\\print("No non-primitive solutions found"),
			if (gcd_flag == 1,
				for( i = 1, #final_add_res,
					x = final_add_res[i][1];
					y = final_add_res[i][2];
					final_add_res[i][1] = alpha * x + betta * y;
					final_add_res[i][2] = gamma * x + deltaa * y;
				);
			);
			\\print("non-primitive solutions: ");
			\\for ( i = 1, #final_add_res,
			\\	print("(X,Y) = ("final_add_res[i][1]", "final_add_res[i][2]")");
			\\);
	);
	return([final_res, final_add_res]);
}



\\\\\\\\\\\\ matrix conformity \\\\\\\\\\\\\\\\

is_size_of_matrix_2x2(matrix) = {
	if (#matrix == 2 && #matrix[1,] == 2 && #matrix[2,] == 2,
		return(1));
	return(0);
}

is_matrix_uniformity(matrA, matrB, isReturnX = 0) = {
	local(matrixS, basis, sol_1, sol, sol_minus1, A, B, C, X, res, x,y);

	if (!is_size_of_matrix_2x2(matrA) || !is_size_of_matrix_2x2(matrB),
		print("	!! Sorry, Wrong matrix size! 2x2 matrices are required!");
		return;
	);
	if ((matdet(matrA) != matdet(matrB)) || (matrA[1,1]+matrA[2,2] != matrB[1,1]+matrB[2,2]),
		print(" !! Sorry, Traces of matrices or determinants of matrices NOT EQUAL!");
		return;
	);

	matrixS = [
		matrA[1,1]-matrB[1,1],	matrA[1,2],			-matrB[2,1],		0;
		matrA[2,1],			matrA[2,2]-matrB[1,1],	0,				-matrB[2,1];
		-matrB[1,2],		0,				matrA[1,1]-matrB[2,2],	matrA[1,2];
		0,				-matrB[1,2],		matrA[2,1],			matrA[2,2]-matrB[2,2]
	];

	basis = matkerint(matrixS);

	A = basis[1,1]*basis[4,1]-basis[2,1]*basis[3,1];
	B = basis[1,2]*basis[4,1]+basis[1,1]*basis[4,2]-basis[2,2]*basis[3,1]-basis[2,1]*basis[3,2];
	C = basis[1,2]*basis[4,2]-basis[2,2]*basis[3,2];

	sol = ABC_N(A,B,C,1);
	if (!sol || !#sol[1],
		sol = ABC_N(A,B,C,-1);
	);	
	if (sol && #sol[1],
		if (isReturnX,
			res = Vec(sol[1]);
			x = res[1][1]; y = res[1][2];

			X = [basis[1,1]*x + basis[1,2]*y, basis[3,1]*x+basis[3,2]*y; 
				basis[2,1]*x+basis[2,2]*y, basis[4,1]*x+basis[4,2]*y];
			\\print("A = "matrA, "	B = "matrB"	X = "X);

			return(X),
				return(1);
		);
		\\return(1);
	);
	return(0);
}


matrix_similarity_classes(d, isPrintMaxtices = 0) = {
	local(a, aMax,bc,b,c, divisors_bc, A, B,
		classesMatrix, iStart, iEnd,
		isUniformity, matrixX);

	aMax = floor((sqrt(2*d-1)-1) \ 2);

	A = [0, 1; d, 0];

	classesMatrix = List();
	listput(classesMatrix, A);

	for( a = 0, aMax,
		bc = d - a^2;
		divisors_bc = divisors(bc);

		if (a == 0,
			iStart = 2; iEnd = #divisors_bc,
				iStart = 1; iEnd = #divisors_bc;
		);
		for( i = iStart, iEnd,
			b = divisors_bc[i];
			c = bc \ b;

			B = [a,b;c,-a];

			for( j = 1, #classesMatrix,
				isUniformity = is_matrix_uniformity(classesMatrix[j], B);

				if (isUniformity,
					break;
				);
			);

			if (!isUniformity,
				listput(classesMatrix, B);
			);
		);
	);

	if (isPrintMaxtices, 
		return([#classesMatrix, classesMatrix]);
	);
	return([#classesMatrix]);
}

debug_is_matrix_uniformity(matrA, matrB) = {
	local(matrixS, basis, sol_1, sol, sol_minus1, A, B, C, X, res, x,y);

	print("	----[debug_is_matrix_uniformity - start]--------");
	print("	matrA = "matrA"\n	matrB = "matrB);

	if (!is_size_of_matrix_2x2(matrA) || !is_size_of_matrix_2x2(matrB),
		print("	!! Sorry, Wrong matrix size! 2x2 matrices are required!");
		return;
	);
	if ((matdet(matrA) != matdet(matrB)) || (matrA[1,1]+matrA[2,2] != matrB[1,1]+matrB[2,2]),
		print(" !! Sorry, Traces of matrices or determinants of matrices NOT EQUAL!");
		return;
	);

	matrixS = [
		matrA[1,1]-matrB[1,1],	matrA[1,2],			-matrB[2,1],		0;
		matrA[2,1],			matrA[2,2]-matrB[1,1],	0,				-matrB[2,1];
		-matrB[1,2],		0,				matrA[1,1]-matrB[2,2],	matrA[1,2];
		0,				-matrB[1,2],		matrA[2,1],			matrA[2,2]-matrB[2,2]
	];

	print("	S = "matrixS);

	basis = matkerint(matrixS);

	print("	matkerint = "basis);

	A = basis[1,1]*basis[4,1]-basis[2,1]*basis[3,1];
	B = basis[1,2]*basis[4,1]+basis[1,1]*basis[4,2]-basis[2,2]*basis[3,1]-basis[2,1]*basis[3,2];
	C = basis[1,2]*basis[4,2]-basis[2,2]*basis[3,2];

	print("	A = "A);
	print("	B = "B);
	print("	C = "C);

	sol = ABC_N(A,B,C,1);
	print("	Ax^2+Bxy+Cy^2=1		[prim_sol, non_prim_sol] -" sol);
	if (!sol || !#sol[1],
		sol = ABC_N(A,B,C,-1);
		print("	Ax^2+Bxy+Cy^2=-1	[prim_sol, non_prim_sol] -" sol);
	);	
	
	if (sol && #sol[1],
		res = Vec(sol[1]);

		x = res[1][1];
		y = res[1][2];

		X = [basis[1,1]*x + basis[1,2]*y, basis[3,1]*x+basis[3,2]*y; 
			basis[2,1]*x+basis[2,2]*y, basis[4,1]*x+basis[4,2]*y];
		print("	X = "X);
		print("	----[debug_is_matrix_uniformity - end]--------");
		return(1);
	);
	print("	----[debug_is_matrix_uniformity - end]--------");
	return(0);
}

debug_matrix_similarity_classes(d, isPrintMaxtices = 0) = {
	local(a, aMax,bc,b,c, divisors_bc, A, B,
		classesMatrix, iStart, iEnd,
		isUniformity);

	aMax = floor((sqrt(2*d-1)-1) \ 2);
	print("0 <= a <= "aMax);

	A = [0, 1; d, 0];

	classesMatrix = List();
	listput(classesMatrix, A);

	for( a = 0, aMax,

		print("a = "a);
		bc = d - a^2;
		print("	bc = "bc);
		divisors_bc = divisors(bc);
		print("	divisors("bc")="divisors_bc);


		if (a == 0,
			iStart = 2; iEnd = #divisors_bc,
				iStart = 1; iEnd = #divisors_bc;
		);
		for( i = iStart, iEnd,
			b = divisors_bc[i];
			c = bc \ b;


			B = [a,b;c,-a];

			print("   B = "B);

			for( j = 1, #classesMatrix,
				isUniformity = debug_is_matrix_uniformity(classesMatrix[j], B);

				print("   isUniformity: "classesMatrix[j]", "B" - " isUniformity);
				print("   --------------------------------------------");

				if (isUniformity,
					break;
				);
			);

			if (!isUniformity,
				listput(classesMatrix, B);
			);

			
		);
	);

	if (isPrintMaxtices, 
		print("RESULT:");
		for( i = 1, #classesMatrix, print(classesMatrix[i]));
		print("NUMBER of matrix similarity classes: " #classesMatrix);
		return([#classesMatrix, classesMatrix]);
	);
	
	print("NUMBER of matrix similarity classes: " #classesMatrix);
	return([#classesMatrix]);
}

matrix_similarity_classes_interval_d(dStart, dEnd) = {
	local(countOfClasses, classes, denom, persent, 
		outFileName, outFile, isClassesWriteToo, forPrime);

	forPrime = 0;
	if(confirm("Should it be calculated only for prime numbers?"),
		forPrime = 1;
	);

	isClassesWriteToo = 0;
	if(confirm("Should representatives of similarity classes also be written?"),
		isClassesWriteToo = 1;
	);

	if (forPrime,
		outFileName = Str("matrix_similarity_D="dStart"-"dEnd"_prime.txt"),
			outFileName = Str("matrix_similarity_D="dStart"-"dEnd".txt")
	);
	\\outFileName = Str("matrix_similarity_D="dStart"-"dEnd".txt");
	outFile = fileopen(outFileName, "w");
	filewrite(outFile, "D	countOfClasses");

	if (isClassesWriteToo,
		if (forPrime,
			outFileName2 = Str("matrix_similarity_D="dStart"-"dEnd"_prime_with_classes.txt"),
				outFileName2 = Str("matrix_similarity_D="dStart"-"dEnd"_with_classes.txt");
		);
		\\outFileName2 = Str("matrix_similarity_D="dStart"-"dEnd"_with_classes.txt");
		outFile2 = fileopen(outFileName2, "w");
		filewrite(outFile2, "D	countOfClasses	representatives_of_similarity_classes");
	);
	
	printf("0%%");
	denom = dEnd - dStart;

	if (isClassesWriteToo,
		mainSolutionWithClassesWrite(dStart, dEnd, outFile, outFile2, forPrime),
			mainSolutionWithoutClassesWrite(dStart, dEnd, outFile, forPrime);
	);

	printProgress(1);
	fileclose(outFile);
	print("\nData is written to file - "outFileName);
	if (isClassesWriteToo,
		fileclose(outFile2);
		print("Data is written to file - "outFileName2);
	);
}


mainSolutionWithoutClassesWrite(dStart, dEnd, outFile, forPrime) = {
	local(countOfClasses, classes, persent, denom);

	denom = dEnd - dStart;

	if (forPrime,
		forprime( d = dStart, dEnd,
			if (!checkFullSquare(d),
				[countOfClasses, classes] = matrix_similarity_classes(d, 1);

				persent = (d - dStart)*100 \ denom;
				printProgress(persent / 100.0, d);

				filewrite(outFile, d"	"countOfClasses);
				fileflush(outFile);
			);
		),
			for( d = dStart, dEnd,
				if (!checkFullSquare(d),
					[countOfClasses, classes] = matrix_similarity_classes(d, 1);

					persent = (d - dStart)*100 \ denom;
					printProgress(persent / 100.0, d);

					filewrite(outFile, d"	"countOfClasses);
					fileflush(outFile);
				);
			);
	);

}
mainSolutionWithClassesWrite(dStart, dEnd, outFile, outFile2, forPrime) = {
	local(countOfClasses, classes, persent);

	if (forPrime,
		forprime( d = dStart, dEnd,
			if (!checkFullSquare(d),
				[countOfClasses, classes] = matrix_similarity_classes(d, 1);

				persent = (d - dStart)*100 \ denom;
				printProgress(persent / 100.0, d);

				filewrite(outFile, d"	"countOfClasses);
				fileflush(outFile);

				\\filewrite(outFile2, d"	"countOfClasses"	"Vec(classes));
				filewrite(outFile2, d"	"countOfClasses"	"classes[1]);
				if(#classes > 1,
					for( j = 2, #classes,
						filewrite(outFile2, "		"classes[j]);
					);
				);
				filewrite(outFile2,"");
				fileflush(outFile2);
			);
		),
			for( d = dStart, dEnd,
				if (!checkFullSquare(d),
					[countOfClasses, classes] = matrix_similarity_classes(d, 1);

					persent = (d - dStart)*100 \ denom;
					printProgress(persent / 100.0, d);

					filewrite(outFile, d"	"countOfClasses);
					fileflush(outFile);

					\\filewrite(outFile2, d"	"countOfClasses"	"Vec(classes));
					filewrite(outFile2, d"	"countOfClasses"	"classes[1]);
					if(#classes > 1,
						for( j = 2, #classes,
							filewrite(outFile2, "		"classes[j]);
						);
					);
					filewrite(outFile2,"");
					fileflush(outFile2);
				);
			);
	);
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


mean_count_of_classes(fileName) = {
	my(line,d,countOfClasses, mean, n, out,
		intD, intCountOfClasses, sum_countOfClassses, count);

	n = fileopen(fileName, "r");
	out = fileopen(concat(fileName, "_average_count_of_classes.txt"), "w");
	filewrite(out, "D	CountOfClasses	Average");

	mean = 0.0;
	sum_countOfClassses = 0;
	count = 1;
	print("Start:");
	while(line = filereadstr(n),
		print1(Strchr(13));
		[d, countOfClasses] = str_split_by_char(line, "\t");
		iferr(intD = eval(d); intCountOfClasses = eval(countOfClasses), E, next);

		if (type(intD) == "t_INT" && type(intCountOfClasses) == "t_INT",
			sum_countOfClassses += intCountOfClasses;

			mean = sum_countOfClassses / count * 1.0;

			filewrite(out, d"	"countOfClasses"	"strprintf("%.4f", mean));
			fileflush(out);

			print1("D = "d"	Count = "countOfClasses"	Average = "strprintf("%.4f", mean));
			count++;
		);
	);

	fileclose(n);
	fileclose(out);
	print("\nData is writen to the file - "concat(fileName, "_average_count_of_classes.txt"));
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

	filewrite(out0, "value=4k	numberOfValue");
	filewrite(out1, "value=4k+1	numberOfValue");
	filewrite(out2, "value=4k+2	numberOfValue");
	filewrite(out3, "value=4k+3	numberOfValue");

	print("\nCalculating number of value with value=4k,4k+1,4k+2,4k+3...");	

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

number_of_d_with_classes(fileName) = {
	my(n, out, map);

	n = fileopen(fileName, "r");
	outFileName = concat(fileName, "_number_classes_at_d.txt");
	out = fileopen(outFileName, "w");
	filewrite(out, "countOfClasses	numberOfD");

	map = Map();
	while(line = filereadstr(n),
		print1(Strchr(13));
		[d, countOfClasses] = str_split_by_char(line, "\t");
		iferr(intD = eval(d); intCountOfClasses = eval(countOfClasses), E, next);

		if (type(intD) == "t_INT" && type(intCountOfClasses) == "t_INT",
			if(mapisdefined(map, intCountOfClasses),
				mapput(map, intCountOfClasses, mapget(map, intCountOfClasses)+1);
				print1("d = "intD"	countOfClasses = "intCountOfClasses),
					mapput(map, intCountOfClasses, 1);
					print1("d = "intD"	countOfClasses = "intCountOfClasses);
			);
		);
	);

	print("\nWriting to a file...");
	matMap = Mat(map);
	for( i = 1, #matMap[,1],
		filewrite(out, matMap[i,1]"	"matMap[i,2]);
	);

	fileclose(n);
	fileclose(out);
	print("Data is writen to the file - "outFileName);
}


addhelp(is_matrix_uniformity, " * is_matrix_uniformity(matrixA, matrixB, {returnX = 0}) - Returns 1 if the matrices are similarly, returns 0 else. \n	Input:  Required two 2x2 matrices.\n 		returnX - when =1, returns transform matrix X (AX=XB) instead of 1 if the matrices are similar. \n	Output: 1 - True (transform matrix if the input parameter returnX = 1)), 0 - False");
addhelp(matrix_similarity_classes, " * matrix_similarity_classes(d, {flag = 0}) - Returns the number of matrix similarity classes for a matrix of the form [0, 1; d, 0]\n 	Input:	d - for matrix [0,1;d,0].\n		flag - (default = 0) if 1 - printing representatives of similarity classes.\n 	Output: [count of classes] - if flag = 0\n 		[count of classes, list representatives of similarity classes] - if flag = 1");
addhelp(matrix_similarity_classes_interval_d, " * matrix_similarity_classes_interval_d(dStart, dEnd)	- Perform matrix_similarity_classes(d) method for d in range dStart <= d <= dEnd.\n	Output in file.");

addhelp(mean_count_of_classes, " * mean_count_of_classes(file_name)	- Calculates the average number of classes from a file with a data table of the form: D - number_of_classes. (You can get such a file using the \"matrix_similarity_classes_interval_d(dStart, dEnd)\" function.)");
addhelp(number_of_d_with_classes, " * number_of_d_with_classes(file_name)	- Calculates how many times a certain number of classes have been encountered in file with a data table of the form: D-number_of_classes. (You can get such a file using the \"matrix_similarity_classes_interval_d(dStart, dEnd)\" function.)");


addhelp(debug, "To switch to the detailed mode, add the word 'debug_'at the beginning of the function name.\n	  For example, 'debug_key_generation(37,19) inctead of 'key_generation(37,19)'.\n The debug version is available for the functions:\n 	* debug_is_matrix_uniformity(matrixA, matrixB)\n 	* debug_matrix_similarity_classes(d, {flag = 0})");


print("		*****************************************");
print("		*	Matrix similarity classes	*")
print("		*****************************************");
print(" This package includes the functions:");
print(" 	* is_matrix_uniformity(matrixA, matrixB, {returnX = 0})	- Returns 1 if the matrices are similarly, 0 else.");
print(" 	* matrix_similarity_classes(d, {flag = 0})	- Returns the number of matrix similarity classes \n 						  	for a matrix of the form [0, 1; d, 0].");
print("	* matrix_similarity_classes_interval_d(dStart, dEnd)	- Perform matrix_similarity_classes(d) method for d\n 								in range dStart <= d <= dEnd.");

print("\n 	* mean_count_of_classes(file_name)	- Calculates the average number of classes from a file\n 						  with a data table of the form: D-number_of_classes.");
print(" 	* number_of_d_with_classes(file_name)	- Calculates how many times a certain number of classes have been encountered\n						  in file with a data table of the form: D-number_of_classes. ");


print("\n Additional functions:\n 	*	ABC_N(A,B,C,N)			- Finding primitive solutions of the equation Ax2+Bxy+Cy2=N.");
print(" 	*	validationInputData(A,B,C)	- Checks the input data for correctness.");
print(" 							Return: 1 - True, 0 - False");
print(" 	*	searchTetta(A,B,C,N)		- Search all theta for a*theta^2 + b*theta + c = 0 (mod |N|).");
print(" 							Return theta list");
print(" 	*	transformationABC(A,B,C,N)	- Conversion to an equation with GCD(A, N) = 1.");
print("							Return list [newA,newB,newC,alpha,betta,gamma,delta]");

print("\n---- For more information about the functions, use help: ? <method_name>");
print("\n--@@@@@@@@-- For more information when executing the function, use the DEBUG mod.");
print("--@@@@@@@@-- Command for getting HELP on DEBUG mod: ? debug");

