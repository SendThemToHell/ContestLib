inline int legandre_symbol(int a, int mod) {
	return pwm(a, (mod - 1) >> 1, mod);
}

inline int get_bad(int mod) {
	while (true) {
		int tr = rd() % (mod - 1) + 1;
		if (legandre_symbol(tr, mod) != 1) {
			return tr;
		}
	}
}

/*returns one of the solutions to x^2 = n (mod p), or -1, if there is no such x*/
int tonelli_shanks(int n, int mod) {
	if (n == 0) {
		return -1;
	}
	if (mod == 2) {
		return 1;
	}
	if (legandre_symbol(n, mod) != 1) {
		return -1;
	}
	int s = 0;
	int q = mod - 1;
	while (!(q & 1)) {
		++s;
		q >>= 1;
	}
	assert(s);

	if (s == 1) {
		return pwm(n, (mod + 1) >> 2, mod);
	}

	int r = pwm(n, (q + 1) >> 1, mod);
	int t = pwm(n, q, mod);
	int c = pwm(get_bad(mod), q, mod);

	while (t != 1) {
		int it;
		int tt = t;
		for (it = 0; tt != 1; ++it) {
			tt = int(ll(tt) * tt % mod);
		}
		int b = pwm(c, 1 << (s - it - 1), mod);
		r = int(ll(r) * b % mod);
		c = int(ll(b) * b % mod);
		t = int(ll(t) * c % mod);
		s = it;
	}
	return r;
}

