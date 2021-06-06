key_generation(pp, qq) = {
	local(n, N, e, d, p, q);

	\\p = prime(random() % (1000 - 100) + 100);
	\\q = prime(random() % (1000 - 100) + 100);

	if (pp && qq,
		p = pp;
		q = qq,
			p = randomprime(2^1024);
			q = randomprime(2^1024);
			while (p == q,
				\\q = prime(random() % (1000 - 100) + 100);
				q = randomprime(2^1024);
			);
	);

	n = p * q;
	N = lcm(p-1,q-1);

	for( i = 2, N, 
		if (gcd(i,N) == 1,
			e = i;
			break;
		);
	);

	d = lift(Mod(1/e, N));
	
	return([[e,n],[p,q,d]]);
}

\\ Scheme 2.
isBelongZnCheck(number, n) = {
	if (number > n || gcd(number, n) != 1,
		return(0)
	);
	return(1);
}

encryption(Mx, My, en) = {
	local(Z1, Z1_reverse, X, Y, a, D, C, M, ciphertext, e, n);

	e = en[1]; n = en[2];

	if (!isBelongZnCheck(Mx,n) || !isBelongZnCheck(My,n),
		print("\nERROR: Mx, My should belong Z*("n").\n");
		return;
	);

	Z1 = Mod(Mx*My, n);
	Z1_reverse = 1/Z1;

	Y = My;
	X = (Z1 + Z1_reverse)/2;

	a = (Z1_reverse - X)/Y;

	if (gcd(lift(a), n) != 1,
		print("\nERROR: gcd(a,n) != 1. Please, generate and give me another public key (3th- argument)\n");
		return;
	);

	D = a^2;
	M =  X-a*Y;
	C = M^e;

	ciphertext = [lift(C),lift(a)];
	return(ciphertext);
}

decryption(Ca, k_pr) = {
	local(M, M_reverse, Mx, My, X, Y, C, a, p, q, d, n);

	p = k_pr[1]; q = k_pr[2]; d=k_pr[3];
	n = p * q;

	C = Ca[1]; a = Ca[2];

	M = Mod(C, n)^d;
	M_reverse = 1/M;

	X = (M_reverse + M)/2;
	Y = (M_reverse - M)/(2*a);

	My = lift(Y);
	Mx = lift(Mod(M/Y, n));

	return([Mx,My]);
}

\\ Scheme 3.

save_encryption(Mx, My, en) = {
	local(Z1, Z1_reverse, X, Y, a, D, C C0, C1, r, binary_r, b, M,
		k, l, r_minus_MSBZ, f_r, ciphertext, e, n);

	e = en[1];
	n = en[2];

	if (!isBelongZnCheck(Mx,n) || !isBelongZnCheck(My,n),
		print("\nERROR: Mx, My should belong Z*("n").\n");
		return;
	);

	Z1 = Mod(Mx*My, n);
	Z1_reverse = 1/Z1;

	Y = My;
	X = (Z1 + Z1_reverse)/2;

	a = (Z1_reverse - X)/Y;

	if (gcd(a, n) != 1,
		print("\nERROR: gcd(a,n) != 1. Please, generate and give me another public key (3th- argument)\n");
		return;
	);

	D = a^2;
	M = X-a*Y;

	r = random(n);
	binary_r = binary(r);
	C0 = lift(Mod(r,n)^e);

	k = #binary_r;
	l = random() % (k-k\2) + k\2;

	f_r = f_MSBZ(binary_r, l, e, n);

	C1 = lift(f_r + M * C0);

	b = lift(a + r^2);

	return([[C0, C1, b], l]);
}

save_decryption(C0C1b, en, pqd, l) = {
	local(r, a, M, M_reverse, X ,Y, Z1, Mx, My,
		C0_reverse, binary_r, f_r, k, C0, C1, b, p, q, d, e, n);

	C0 = C0C1b[1]; C1 = C0C1b[2]; b = C0C1b[3];

	p = pqd[1]; q = pqd[2]; d = pqd[3];
	\\ n = p * q;
	n = en[2];
	e = en[1];

	r = Mod(C0, n)^d;
	a = b - r^2;

	C0_reverse = Mod(1/C0, n);

	binary_r = binary(lift(r));
	k = #binary_r;

	f_r = f_MSBZ(binary_r, l, e, n);

	M = C0_reverse * (C1 - f_r);
	M_reverse = 1/M;

	X = (M_reverse + M)/2;
	Y = (M_reverse - M)/(2*a);

	My = lift(Y);
	Mx = lift(Mod(M/Y, n));

	return([Mx, My]);
}

string_to_intASCII(text) = {
	local(vec, vecLength, m, t, charLength);

	vec = Vecsmall(text);
	vecLength = #vec;

	m = "";
	for( i = 1, #vec,
		charLength = #strexpand(vec[i]);
		if (charLength == 2,
			m = strjoin([m, "0", vec[i]]);
			next;
		);
		if (charLength == 1,
			m = strjoin([m, "00", vec[i]]);
			next;
		);
		m = strjoin([m, vec[i]]);
	);

	\\m = 0;
	\\for (i = 0, #vec-1,
	\\	t = vec[vecLength - i];
	\\	m = m + t*10^(3*i);
	\\);
	return(eval(m));
}

intASCII_to_string(number) = {
	local(vec, modVec, str);

	vec = digits(number, 1000);
	modVec = lift(Mod(vec, 256));

	for( i = 1, #modVec,
		if (modVec[i] < 1,
			modVec[i] = 64;
		);
	);

	str = strchr(modVec);
	return(str);
}

enc_dec_perform = (Mx, My, key_public, key_private, l) -> {
	my(cyphertext, l resMx, resMy);

	[cyphertext, l] = save_encryption(Mx, My, key_public);
	print("\nCyphertext:\n"cyphertext);
	print("\nCyphertext convert to string:\n"
		intASCII_to_string(cyphertext[1]),
		intASCII_to_string(cyphertext[2]),
		intASCII_to_string(cyphertext[3])
	);

	[resMx, resMy] = save_decryption(cyphertext, key_public, key_private, l);
	print("\nDecripted integer:\n Mx = "resMx"\n My = " resMy);
	print("\nDecripted text:\n Mx = \""strchr(resMx)"\"\n My = \""strchr(resMy)"\"");

	return([strchr(resMx), strchr(resMy)]);
}

crypto_text_char_by_char(text) = {
	local(Mx, My, key_public, key_private, asciiIntList,
			cyphertext, l, resText, resMxText, resMyText);

	asciiIntList = Vec(Vecsmall(text));
	print("#asciiIntList = "#asciiIntList);

	print("key generation...");
	[key_public, key_private] = key_generation();

	print("\nkey_public:\n"key_public);
	print("\nkey_private:\n"key_private);

	resText = "";
	forstep( i = 1, #asciiIntList-1, 2,
		Mx = asciiIntList[i];
		My = asciiIntList[i+1];

		[resMxText, resMyText] = enc_dec_perform(Mx, My, key_public, key_private);
		resText = strjoin([resText, resMxText, resMyText]);
	);
	\\ last character
	if (#asciiIntList % 2, 
		Mx = asciiIntList[#asciiIntList];
		[resMxText, resMyText] = enc_dec_perform(Mx, Mx, key_public, key_private);
		resText = strjoin([resText, resMxText]);
	);
	print("resText = "resText);
	return(resText);
}

crypto_text(text) = {
	local(Mx, My, resMx, resMy, key_public, key_private, MxText, MyText
		cyphertext, l,
		resDescriptedTextMx, resDescriptedTextMy, resDescriptedText);

	[MxText, MyText] = split_string_in_half(text);

	print("MxText: "MxText);
	print("MyText: "MyText);

	Mx = string_to_intASCII(MxText);
	My = string_to_intASCII(MyText);

	print("Mx = "Mx);
	print("My = "My);

	print("key generation...");
	[key_public, key_private] = key_generation();

	print("\nkey_public:\n"key_public);
	print("\nkey_private:\n"key_private);

	[cyphertext, l] = save_encryption(Mx, My, key_public);

	print("\nCyphertext:\n"cyphertext);

	print("\nCyphertext convert to string:\n"
		intASCII_to_string(cyphertext[1]),
	 	intASCII_to_string(cyphertext[2]),
	  	intASCII_to_string(cyphertext[3]));

	[resMx, resMy] = save_decryption(cyphertext, key_public, key_private, l);

	print("\nDecripted integer:\n Mx = "resMx"\n My = " resMy);

	resDescriptedTextMx = intASCII_to_string(resMx);
	resDescriptedTextMy = intASCII_to_string(resMy);

	resDescriptedText = strjoin([resDescriptedTextMx, resDescriptedTextMy]);
	print("\nDecripted text:\n"resDescriptedText);
	return(resDescriptedText);
}

split_string_on_parts(text, part_length) = {
	my(res);
	res = List();
	forstep(i = 1, #text, part_length,
		listput(res, ssubstr(text, i, part_length));
	);
	return(Vec(res));
}

crypto_text2(text) = {
	local(permissible_text_length, Mx, My, resText, resMxText, resMyText,
		key_public, key_private, list_text_parts);
	permissible_text_length = 205; \\ (308*2) \ 3

	print("key generation...");
	[key_public, key_private] = key_generation();

	if(#text <= permissible_text_length,
		return(crypto_text(text));
	);

	resText = "";
	list_text_parts = split_string_on_parts(text, permissible_text_length);
	forstep( i = 1, #list_text_parts-1, 2,
		Mx = string_to_intASCII(list_text_parts[i]);
		My = string_to_intASCII(list_text_parts[i+1]);

			[cyphertext, l] = save_encryption(Mx, My, key_public);
			print("\nCyphertext:\n"cyphertext);
			print("\nCyphertext convert to string:\n"
				intASCII_to_string(cyphertext[1]),
				intASCII_to_string(cyphertext[2]),
				intASCII_to_string(cyphertext[3])
			);

			[resMx, resMy] = save_decryption(cyphertext, key_public, key_private, l);
			print("\nDecripted integer:\n Mx = "resMx"\n My = " resMy);
			resMxText = intASCII_to_string(resMx);
			resMyText = intASCII_to_string(resMy);
			print("\nDecripted text:\n Mx = \""resMxText"\"\n My = \""resMyText"\"");
		resText = strjoin([resText, resMxText, resMyText]);
	);
	\\ last part
	if (#list_text_parts % 2, 
		Mx = string_to_intASCII(list_text_parts[#list_text_parts]);
		
			[cyphertext, l] = save_encryption(Mx, Mx, key_public);
			print("\nCyphertext:\n"cyphertext);
			print("\nCyphertext convert to string:\n"
				intASCII_to_string(cyphertext[1]),
				intASCII_to_string(cyphertext[2]),
				intASCII_to_string(cyphertext[3])
			);

			[resMx, resMy] = save_decryption(cyphertext, key_public, key_private, l);
			print("\nDecripted integer:\n Mx = "resMx);
			resMxText = intASCII_to_string(resMx);
			print("\nDecripted text:\n Mx = \""resMxText"\"");
		resText = strjoin([resText, resMxText]);
	);
	print("resText = "resText);
	return(resText);
}


ssubstr(str,s=1,n=0)={
	my(vt=Vecsmall(str),ve,vr,vtn=#str,n1);
	if(vtn==0,return(""));
	if(s<1||s>vtn,return(str));
	n1=vtn-s+1; if(n==0,n=n1); if(n>n1,n=n1);
	ve=vector(n,z,z-1+s); vr=vecextract(vt,ve); return(Strchr(vr));
}

split_string_in_half(str) = {
	local(leftHalf, rightHalf);

	if (#str < 2, return([str, ""]));
	leftHalf = ssubstr(str, 1, floor(#str/2));
	rightHalf = ssubstr(str, floor(#str/2)+1);
	return([leftHalf, rightHalf]);
}

\\\\\\\\\\\\\\\\\\\\\\\ CRYPTO TIME TEST \\\\\\\\\\\\\\\\\\\\\
crypto_text_char_by_char_test(text, key_public, key_private) = {
	local(Mx, My, asciiIntList, resMx, resMy,
			cyphertext, l, resText, resMxText, resMyText);

	asciiIntList = Vec(Vecsmall(text));
	resText = "";
	forstep( i = 1, #asciiIntList-1, 2,
		Mx = asciiIntList[i];
		My = asciiIntList[i+1];

		[cyphertext, l] = save_encryption(Mx, My, key_public);
		[resMx, resMy] = save_decryption(cyphertext, key_public, key_private, l);
		resText = strjoin([resText, strchr(resMx), strchr(resMy)]);
	);
	\\ last character
	if (#asciiIntList % 2, 
		Mx = asciiIntList[#asciiIntList];
		[cyphertext, l] = save_encryption(Mx, Mx, key_public);
		[resMx, resMy] = save_decryption(cyphertext, key_public, key_private, l);
		resText = strjoin([resText, strchr(resMx)]);
	);
	return(resText);
}
crypto_text_test(text, key_public, key_private) = {
	local(Mx, My, resMx, resMy, MxText, MyText, resDescriptedTextMx, resDescriptedTextMy,
		resDescriptedText, cyphertext, l);

	[MxText, MyText] = split_string_in_half(text);

	Mx = string_to_intASCII(MxText);
	My = string_to_intASCII(MyText);

	[cyphertext, l] = save_encryption(Mx, My, key_public);
	[resMx, resMy] = save_decryption(cyphertext, key_public, key_private, l);

	resDescriptedTextMx = intASCII_to_string(resMx);
	resDescriptedTextMy = intASCII_to_string(resMy);

	resDescriptedText = strjoin([resDescriptedTextMx, resDescriptedTextMy]);
	return(resDescriptedText);
}
crypto_text2_test(text, key_public, key_private) = {
	my(permissible_text_length, resText, Mx, My, cyphertext, l, resMx, resMy,
		resMxText,  resMyText, list_text_parts);

	permissible_text_length = 205;
	if(#text <= permissible_text_length,
		return(crypto_text_test(text, key_public, key_private));
	);

	resText = "";
	list_text_parts = split_string_on_parts(text, permissible_text_length);
	forstep( i = 1, #list_text_parts-1, 2,
		Mx = string_to_intASCII(list_text_parts[i]);
		My = string_to_intASCII(list_text_parts[i+1]);

		[cyphertext, l] = save_encryption(Mx, My, key_public);
		[resMx, resMy] = save_decryption(cyphertext, key_public, key_private, l);
		resMxText = intASCII_to_string(resMx);
		resMyText = intASCII_to_string(resMy);
		resText = strjoin([resText, resMxText, resMyText]);
	);
	\\ last part
	if (#list_text_parts % 2, 
		Mx = string_to_intASCII(list_text_parts[#list_text_parts]);
		
		[cyphertext, l] = save_encryption(Mx, Mx, key_public);
		[resMx, resMy] = save_decryption(cyphertext, key_public, key_private, l);
		resMxText = intASCII_to_string(resMx);
		resText = strjoin([resText, resMxText]);
	);
	return(resText);
}

crypto_time_test(text) = {
	local(key_public, key_private, resTextCbC, resText, resText2, timeStart);

	print1("\nKey generation...");
	[key_public, key_private] = key_generation();
	print1(Strchr(13));

	print("crypto_text_char_by_char TEST:");
	timeStart = getwalltime();
	resTextCbC = crypto_text_char_by_char_test(text, key_public, key_private);
	print(" Time: "getwalltime() - timeStart" ms.");
	print(" resTextCbC = "resTextCbC);
	print(" Is inputText equal to resultText?: "text == resTextCbC"\n");

	print("crypto_text2 TEST:");
	timeStart = getwalltime();
	resText2 = crypto_text2_test(text, key_public, key_private);
	print(" Time: "getwalltime() - timeStart" ms.");
	print(" resText = "resText2);
	print(" Is inputText equal to resultText?: "text == resText2"\n");


	print("crypto_text TEST:");
	timeStart = getwalltime();
	resText = crypto_text_test(text, key_public, key_private);
	print(" Time: "getwalltime() - timeStart" ms.");
	print(" resText = "resText);
	print(" Is inputText equal to resultText?: "text == resText"\n");
}


\\\\\\\\\\\\\\\\\\\\\\\ DEBUG_MODE \\\\\\\\\\\\\\\\\\\\\\
debug_key_generation(pp, qq) = {
	local(n, N, e, d, p, q);

	if (pp && qq,
		p = pp;
		q = qq,
			p = randomprime(10^307);
			q = randomprime(10^307);
			\\p = prime(random() % (1000 - 100) + 100);
			\\q = prime(random() % (1000 - 100) + 100);
			while (p == q,
				q = randomprime(10^307);
				\\q = prime(random() % (1000 - 100) + 100);
			);
	);

	print(" p = "p);
	print(" q = "q);

	if (!isprime(p) || !isprime(q),
		error("ERROR: Number q or q - is not prime\n");
	);

	n = p * q;
	N = lcm(p-1,q-1);

	print("	n = "n);
	print("	N = "N);

	for( i = 2, N, 
		if (gcd(i,N) == 1,
			e = i;
			break;
		);
	);

	d = lift(Mod(1/e, N));

	print("	e = "e);
	print("	d = "d);
	
	return([[e,n],[p,q,d]]);
}

\\ Scheme 2.

debug_encryption(Mx, My, en) = {
	local(Z1, Z1_reverse, X, Y, a, D, C, M, ciphertext, e, n);

	e = en[1]; n = en[2];

	if (!isBelongZnCheck(Mx,n) || !isBelongZnCheck(My,n),
		print("\nERROR: Mx, My should belong Z*("n").\n");
		return;
	);

	Z1 = Mod(Mx*My, n);
	print("	Z1 =		"Z1);

	Z1_reverse = 1/Z1;
	print("	Z1_reverse =	"Z1_reverse);

	Y = My;
	X = (Z1 + Z1_reverse)/2;

	print("	X =		"X"\n	Y =		"Y);

	a = (Z1_reverse - X)/Y;
	print("	a =		"a);

	if (gcd(lift(a), n) != 1,
		print("\nERROR: gcd(a,n) != 1. Please, generate and give me another public key (3th- argument)\n");
		return;
	);

	D = a^2;
	print("	D =		"D);

	M =  X-a*Y;
	print("	M =		"M);

	C = M^e;
	print("	C =		"C);

	ciphertext = [lift(C),lift(a)];
	return(ciphertext);
}

debug_decryption(Ca, k_pr) = {
	local(M, M_reverse, Mx, My, X, Y, C, a, p, q, d, n);

	p = k_pr[1]; q = k_pr[2]; d=k_pr[3];
	n = p * q;

	C = Ca[1]; a = Ca[2];

	M = Mod(C, n)^d;
	print("	M =		"M);
	M_reverse = 1/M;
	print("	M_reverse =	"M_reverse);

	X = (M_reverse + M)/2;
	Y = (M_reverse - M)/(2*a);

	print("	X =		"X"\n	Y =		"Y);

	My = lift(Y);
	Mx = lift(Mod(M/Y, n));

	print("	Mx =		"Mx"\n	My =		"My);

	return([Mx,My]);
}


\\ Scheme 3.

debug_save_encryption(Mx, My, en) = {
	local(Z1, Z1_reverse, X, Y, a, D, C C0, C1, r, binary_r, b, M,
		k, l, r_minus_MSBZ, f_r, ciphertext, e, n);

	e = en[1];
	n = en[2];

	if (!isBelongZnCheck(Mx,n) || !isBelongZnCheck(My,n),
		print("\nERROR: Mx, My should belong Z*("n").\n");
		return;
	);

	Z1 = Mod(Mx*My, n);
		print("	Z1		= "Z1);

	Z1_reverse = 1/Z1;
		print("	Z1_reverse	= "Z1_reverse);

	Y = My;
	X = (Z1 + Z1_reverse)/2;

		print("	X		= "X"\n	Y		= "Y);

	a = (Z1_reverse - X)/Y;
		print("	a		= "a);
	if (gcd(a, n) != 1,
		print("\nERROR: gcd(a,n) != 1. Please, generate and give me another public key (3th- argument)\n");
		return;
	);

	D = a^2;
		print("	D		= "D);

	M = X-a*Y;
		print("	M		= "M);

	r = random(n);
	binary_r = binary(r);
	C0 = lift(Mod(r,n)^e);

		print("	r		= "r);
		print("	C0		= "C0);

	k = #binary_r;
	l = random() % (k-k\2) + k\2;

		print("	k (length(binary_r)) = "k);
		print("	l		= "l);

	f_r = debug_f_MSBZ(binary_r, l, e, n);
		print("	f(r)		= " f_r);

	C1 = lift(f_r + M * C0);
		print("	C1		= " C1);

	b = lift(a + r^2);
		print("	b		= " b);

	return([[C0, C1, b], l]);
}

debug_save_decryption(C0C1b, en, pqd, l) = {
	local(r, a, M, M_reverse, X ,Y, Z1, Mx, My,
	C0_reverse, binary_r, f_r, k, C0, C1, b, p, q, d, e, n);

	C0 = C0C1b[1]; C1 = C0C1b[2]; b = C0C1b[3];

	p = pqd[1]; q = pqd[2]; d = pqd[3];
	\\ n = p * q;
	n = en[2];
	e = en[1];

	print(" p	= "p);
	print(" q	= "q);
	print(" d	= "d);
	print(" e	= "e);

	r = Mod(C0, n)^d;
		print("	r		= " r);

	a = b - r^2;
		print("	a		= " a);

	C0_reverse = Mod(1/C0, n);
		print("	C0^(-1)		= " C0_reverse);

	binary_r = binary(lift(r));
	k = #binary_r;

		print("	k (length(binary_r)) = "k);
		print("	l		= "l);

	f_r = debug_f_MSBZ(binary_r, l, e, n);
		print("	f(r)		= " f_r);

	M = C0_reverse * (C1 - f_r);
		print("	M		= " M);
	M_reverse = 1/M;
		print("	M^(-1)		= " M_reverse);

	X = (M_reverse + M)/2;
	Y = (M_reverse - M)/(2*a);

		print("	X		= " X"\n	Y		= "Y);

	My = lift(Y);
	Mx = lift(Mod(M/Y, n));

		print("	Mx		= "Mx"\n	My		= "My);

	return([Mx, My]);
}


\\ additional
to_decimal = (binary_num)-> {
	my(sum,l);
	sum=0;
	l=length(binary_num);
	forstep(i=l,1,-1,
		sum+=binary_num[i]*2^(l-i);
	);
	return(sum);
}

f_MSBZ = (binary_num, l, e, n)->{
	my(sum,len, r_minus_MSBZ);

	sum=0;
	len=length(binary_num);
	forstep(i=len,len-l,-1,
		sum+=binary_num[i]*2^(len-i);
	);

	r_minus_MSBZ = lift(Mod(sum^e, n));

	return(r_minus_MSBZ);
}

debug_f_MSBZ = (binary_num, l, e, n)->{
	my(sum,len, r_minus_MSBZ);

	print("****************** f_MSBZ start ***********");
	print(" binary_r = "binary_num);
	print(" l = "l"; e = "e"; n = "n);

	sum=0;
	len=length(binary_num);
	forstep(i=len,len-l,-1,
		sum+=binary_num[i]*2^(len-i);
		\\print(sum);
		\\print("i = "i);
		\\print("binary_num[i] ="binary_num[i]);
		\\print("degree = "len - i);
	);

	\\print("sum = "sum);

	r_minus_MSBZ = lift(Mod(sum^e, n));

	print(" f(r) = "r_minus_MSBZ);

	print("****************** f_MSBZ end ***********");

	return(r_minus_MSBZ);
}


addhelp(key_generation, " 	* key_generation() -> [[e,n],[p,q,d]]	- Generation cryptosystem key. \n			Output: [e,n]	- public key,\n				[p,q,d] - private key.\n");
addhelp(encryption, " 	* encryption(Mx, My, [e, n]) -> [C, a]	- Encryption of the message Mx and My. \n			Input:	Mx, My	- message, Mx,My belongs to Z*(n);\n				[e, n]	- public key.\n			Output:	[C, a] - ciphertext (encripted message).\n");
addhelp(decryption, " 	* decryption([C, a], [p, q, d]) -> [Mx, My]	- Decryption of the encripted message (C,a). \n			Input:	[C, a]	- ciphertext (encripted message);\n				[p, q, d] - private key from key_generation result.\n			Output:	[Mx, My] - decripted meassage.\n");
addhelp(save_encryption, " 	* save_encryption(Mx, My, [e, n]) -> [[C0, C1, b], l]   - Encryption of the message Mx and My. \n			Input:	Mx, My	- message, Mx,My belongs to Z*(n);\n				[e, n]	- public key.\n			Output:	[C0, C1, b] - ciphertext.\n				l	- large number for calculate f(r)\n");
addhelp(save_decryption, " 	* save_decryption([C0, C1, b], [e, n], [p, q, d], l) -> [Mx, My]   - Decryption of the encripted message (C0, C1, b).\n			Input:	[C0, C1, b]	- encripted message;\n				[e, n]		- public key.\n				[p, q, d]	- private key from key_generation result.\n				l		- large number for calculate f(r).\n			Output:	[Mx, My] - decripted meassage.\n");
addhelp(crypto_text, " 	* crypto_text(message) - [Not recommended for use]\n			Visualization of message encryption and decryption processes.\n			This method divides the message in half. First part - Mx, Second part - My\n			Input:	message - text meassage.");
addhelp(crypto_text2, " 	* crypto_text2(message) - Visualization of message encryption and decryption processes.\n			This method divides the message into several parts with the maximum allowed length.\n			Input:	message - text meassage.");
addhelp(crypto_text_char_by_char_test, " 	* crypto_text_char_by_char_test(message) - Visualization of message encryption and decryption processes.\n			This method performs character-by-character encryption.\n			Input:	message - text meassage.");
addhelp(crypto_time_test, "Compares the execution time of the following methods:\n		1) crypto_text_char_by_char(message);\n		2) crypto_text2(message);\n		3) crypto_text(message);");

addhelp(debug, "To switch to the detailed mode, add the word 'debug_'at the beginning of the function name.\n	  For example, 'debug_key_generation(37,19) instead of 'key_generation(37,19)'.\n Allowed to DEBUG mode methods:\n 	debug_key_generation();\n	debug_encryption(Mx, My, en)\n	debug_decryption(Ca, k_pr)\n	debug_save_encryption(Mx, My, en)\n	debug_save_decryption(C0C1b, en, pqd, l)\n	debug_f_MSBZ(binary_num, l, e, n)")

print("		*****************************************************************");
print("		*	A PUBLIC KEY CRYPTOSYSTEM BASED ON PELL EQUATION	*")
print("		*****************************************************************");
print("This package includes the functions:\n 	* key_generation()		-> [[e,n],[p,q,d]]	- Generation cryptosystem key.");
print("\n 	* encryption(Mx, My, [e, n])	-> [C, a]	- Encryption of the message Mx and My.");
print(" 	* decryption([C, a], [p, q, d])	-> [Mx, My]	- Decryption of the encripted message (C,a).");

print("\n 	* save_encryption(Mx, My, [e, n]) -> [[C0, C1, b], l]   - Encryption of the message Mx and My.");
print(" 	* save_decryption([C0, C1, b], [e, n], [p, q, d], l) -> [Mx, My]   - Decryption of the encripted message (C0, C1, b).");
print("\n 	* crypto_text(message) - [Not recommended for use] Visualization of message encryption and decryption processes.\n			(divides the message in half).");
print(" 	* crypto_text2(message) - Visualization of message encryption and decryption processes.\n			(divides the message into several parts).");
print(" 	* crypto_text_char_by_char_test(message) - Visualization of message encryption and decryption processes.\n			(character-by-character encryption).");
print(" 	* crypto_time_test(message) - Compares the execution time of different crypto_ methods.");

print("\n\n--@@@@@@@@-- For more information when executing the method, use the DEBUG mod.");
print("--@@@@@@@@-- Command for getting HELP on DEBUG mod: ? debug");
print("\n---- For more information about the methods, use help: ? <method_name>");
print("\n 	version: v2_mod - removed unnecessary conversions from Mod to normal arithmetic.\n 	Now most operations are calculated in modular arithmetic.")