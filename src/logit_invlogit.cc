# include <cmath>

// Function logit
double logit (double x) {
	return std::log(x) - std::log(1-x);
}

// Function invlogit
double invlogit (double x) {
	if (x > 0) {
		return 1 / (1 + std::exp(-x));
	}
	else {
		return std::exp(x) / (1 + std::exp(x));
	}
}

// End