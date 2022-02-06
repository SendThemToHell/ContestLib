inline ll mul(ll a, ll b, ll n) {
#ifndef LOCAL
	int128_t x = a;
	x *= b;
	return x % n;
#endif // !LOCAL
	ll ret = 0;
	for (; b; b >>= 1) {
		if (b & 1) {
			ret += a;
			if (ret >= n) ret -= n;
		}
		a <<= 1;
		if (a >= n) a -= n;
	}
	return ret;
}

inline ll f(ll x, ll n) {
	// return (x * x + 1) % n; USE THIS FOR INTEGER n !!!!

	ll x2 = mul(x, x, n);
	return x2 + 1 == n ? 0 : x2 + 1;
}

ll rho(ll n, int iter = 500, int iter2 = 5) { // iter ~ sqrt(min_d), iter2 ~ log(prob_of_fail), (500, 5) works for n <= 1e9
	for (int i = 0; i < iter2; ++i) { // in general, failing prob is ~ (1 - e^{-L})^{iter2}, where L = (iter / n^{1/4})^2 / 2
		ll x2i = rd() % n;
		ll xi = x2i;

		for (int i = 0; i < iter; ++i) {
			x2i = f(f(x2i, n), n);
			xi = f(xi, n);
			ll d = gcd(n, abs(xi - x2i));
			if (d > 1 && d < n) {
				assert(n % d == 0);
				return d;
			}
		}
	}
	return -1;
}