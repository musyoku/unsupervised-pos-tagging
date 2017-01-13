#ifndef _util_
#define _util_

double min(double a, double b){
	return a > b ? b : a;
}

double factorial(double n) {
	if (n == 0){
		return 1;
	}
	return n * factorial(n - 1);
}

#endif