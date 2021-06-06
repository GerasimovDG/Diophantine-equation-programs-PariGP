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
	print("D,P,Q: \n" D "\nP = "P "\n Q = "Q);

	if ((D - sqr(P)) % Q != 0,
		D = D * sqr(Q);
		P = P * abs(Q);
		Q = Q * abs(Q);
		print(" New D,P,Q: \nD = "D "\nP = "P "\nQ = "Q);
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

		print("q"count-1" = "qi", P"count - 1" = "Pi",  Q"count - 1" = "Qi);

		mOne = -1;
		j = i;
		if (isMain != 1,
			j = i + 1;
		);

		if (Qi == mOne^j * N \ abs(N),
			print("!!!!! Qi = "Qi", Ai = "Ai", Bi = "Bi);
			X = Ai;
			y = Bi;

			Y = y;
			x = y * tetta + abs(N) * X;
			res = List([x, Y]);
			print("!!!!!!!!! resX = "res[1]", resY = "res[2]);
			listput(~result, res);
		);

		Pi = qi * Qi - Pi;
		Qi = (D - sqr(Pi)) \ Qi;

		Ai = qi * Ai1 + Ai2;
		Bi = qi * Bi1 + Bi2;

		listput(vecA, Ai);
		listput(vecB, Bi);

		print(";     A"count-1"/B"count - 1" = "vecA[count]"/"vecB[count]);

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
				print("--------------------------------------------------------");
			);
		);
		count++;
	);
	print("Perion length l = "l);

	if (t != 0,
		print(" t = "t);
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
	print("D,P,Q:\nD = "D"\nP = "P"\nQ = "Q);

	if ((D - sqr(P)) % Q != 0,
		D = D * sqr(Q);
		P = P * abs(Q);
		Q = Q * abs(Q);
		print(" New D,P,Q: \nD = "D "\nP = "P "\nQ = "Q);

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

		print("q"count-1" = "qi", P"count - 1" = "Pi",  Q"count - 1" = "Qi);

		mOne = -1;
		j = i;
		if (isMain != 1,
			j = i + 1;
		);
		if (Qi == mOne^j * 2 * N \ abs(N),
			print("\n!!!!!!!!! Qi = "Qi", Ai = "Ai", Bi ="Bi);
			X = Ai;
			y = Bi;

			Y = y;
			x = y * tetta + abs(N) * X;

			res = List([x, Y]);
			print("!!!!! resX = "res[1]", resY = "res[2]);
			listput(result, res);

			if (D == 5 && (A * N < 0) && isMain == 1,
				X = Ai - Ai2;
				y = Bi - Bi2;

				Y = y;
				x = y * tetta + abs(N) * X;

				res[1] = x;
				res[2] = Y;

				print("Special test: resX = "res[1]", resY = "res[2]);
				listput(result, res);
			);
		);

		Pi = qi * Qi - Pi;
		Qi = (D - sqr(Pi)) \ Qi;

		Ai = qi * Ai1 + Ai2;
		Bi = qi * Bi1 + Bi2;

		listput(vecA, Ai);
		listput(vecB, Bi);

		print("A"count-1"/B"count-1" = "vecA[count]"/"vecB[count]);

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
				print("--------------------------------------------------------");
			);
		);
		count++;
	);

	print("Perion length l = "l);
	print(result " - result2");
	if (t != 0,
		print(" t = "t);
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

		print("Answers (of these, we choose the one with the smallest y):");

		for (i = 1, #resultEnd,
			print(resultEnd[i][1]", " resultEnd[i][2]);
			if (abs(resultEnd[i][2]) < abs(minY),
				minY = resultEnd[i][2];
				X = resultEnd[i][1];
			);
			if (resultEnd[i][2] == minY && resultEnd[i][1] < X,
				minY = resultEnd[i][2];
				X = resultEnd[i][1];
			);
		);

		solution = List([X, minY]);
		print("solution = "X" "minY);
		return(solution);
	);
	\\return(List());
}

SearchSolution = (tetta, A, B, C ,N) -> {
	my(delta,n,R,S, res, solution);

	delta = sqr(B) - 4 * A * C;
	n = 2 * A * tetta + B;
	print("n=2*a*tetta+b="n);

	R = n \ 2;

	print("R = floor(n\2) = "R);
	S = A * abs(N);
	print("S = A * |N| = "S);

	\\res = List();
	result = List();
	print(delta);
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
		if (!(solution[1] == 0 && solution[2] == 0),
			print(solution[1]" _ "solution[2]" - answer.");
		);
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
		print1(setQ[i]", ");
		Nq *= setQ[i];
	);
	print();

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
		print("xi: "x". yi: "y);

		print("Nq = "Nq);
		print("setQ[i] = "setQ[i]);
		Ni = Nq \ setQ[i];

		print("Ni = "Ni", q = "setQ[i]".");

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

	print("X -"X"\nY -"Y);

	alpha = X \ gcd(X,Y);
	gamma = Y \ gcd(X,Y);

	A = a * sqr(alpha) + b * alpha * gamma + c * sqr(gamma);

	print("A = "A". GCD(alpha, gamma) = "gcd(alpha, gamma)". GCD(A, N) = "gcd(A,N));

	uvd = gcdext(alpha, -gamma);
	delta = uvd[1];
	betta = uvd[2];

	B = 2 * a * alpha * betta + b * alpha * delta + b * betta * gamma + 2 * c * gamma * delta;
	C = a * sqr(betta) + b * betta * delta + c * sqr(delta);

	a = A;
	b = B;
	c = C;

	if (gcd(a, gcd(b,c)) != 1,
		print("!!!!!!!!!!!!!!!! GCD(A,B,C) != 1 !!!!!!!");
		return(List([a,b,c, alpha, betta, gamma, delta]));
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

		print("-------------\n|	alpha	betta	|\n|	gamma	delta	|:\n ------------");
		print("|	"alpha"	"betta"	|	\n|	"gamma"	"deltaa"	|");
		print("-------------");

		gcd_flag = 1;

		print("new A,B,C: "A" "B" " C);
	);

	D = sqr(B) - 4*A*C;
	Delta = D \ 4;
	Q = A * abs(N);

	final_res = List();
	tettaMas = searchTetta(A,B,C,N);
	if (#tettaMas == 0,
		print("There are no matching thetas.");
		return(final_res);
	);

	print1("matching thetas:");
	for(i = 1, #tettaMas,
		print1(tettaMas[i]", ");
	);
	print("");

	for(i = 1, #tettaMas,
		print("--------theta ="tettaMas[i]": ");
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
			print("new N: "newN);
			if (newN != N,
				add_res = mainSolution(A,B,C, newN);
				for (i = 1, #add_res,
					add_res[i][1] = add_res[i][1] * n;
					add_res[i][2] = add_res[i][2] * n;

					print(add_res[i][1]" -;- "add_res[i][2]);
					listput(final_add_res, add_res[i]);
				),
					final_res = mainSolution(A,B,C,N);
			);
		);
	);

	if (#final_res == 0,
		print("No primitive solutions found."),
			print("primitive solutions:");
			for (i = 1, #final_res,
				print ("(X,Y) = ("final_res[i][1]", "final_res[i][2]")");
			);
	);

	if (#final_add_res == 0,
		print("No non-primitive solutions found"),
			if (gcd_flag == 1,
				for( i = 1, #final_add_res,
					x = final_add_res[i][1];
					y = final_add_res[i][2];
					final_add_res[i][1] = alpha * x + betta * y;
					final_add_res[i][2] = gamma * x + deltaa * y;
				);
			);
			print("non-primitive solutions: ");
			for ( i = 1, #final_add_res,
				print("(X,Y) = ("final_add_res[i][1]", "final_add_res[i][2]")");
			);
	);
}


print("		*********************************************************************************");
print("		*	Finding primitive solutions of the diophantine equation ax2+bxy+cy2=n,	*\n		*		where b2-4ac > 0 and is not a perfect square.			*")
print("		*********************************************************************************");
print("This package includes the functions:\n 	*	ABC_N(A,B,C,N)	-	Finding primitive solutions of the equation Ax2+Bxy+Cy2=N.");
print("Additional functions:\n 	*	validationInputData(A,B,C)	-	Checks the input data for correctness. Return: 1 - True, 0 - False");
print(" 	*	searchTetta(A,B,C,N)		-	Search all theta for a*theta^2 + b*theta + c = 0 (mod |N|). Return theta list");
print(" 	*	transformationABC(A,B,C,N)	-	Conversion to an equation with GCD(A, N) = 1. Return list [newA,newB,newC,alpha,betta,gamma,delta]");