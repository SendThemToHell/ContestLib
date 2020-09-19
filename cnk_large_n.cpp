/*
O(p) precalc
O(log_p(n)) query
*/

const int P = 1009;
int fct[P];
int ifct[P];
int invvl[P];

pair<ll, int> calc_fct(ll n) {
	int ans = 1;
	ll cnt = 0;
	while (n > 1) {
		ll nxtn = n / P;
		cnt += nxtn;
		if (nxtn & 1) {
			ans = P - ans;
		}
		ans = ll(ans) * fct[n % P] % P;
		n = nxtn;
	}
	return{ cnt, ans };
}

int cnk(ll n, ll k) {
	auto p1 = calc_fct(n);
	auto p2 = calc_fct(k);
	auto p3 = calc_fct(n - k);
	p2.x += p3.x;
	p2.y = ll(p2.y) * p3.y % P;

	if (p1.x == p2.x) {
		return int(ll(p1.y) * invvl[p2.y] % P);
	}
	return 0;
}

void prec() {
	fct[0] = 1;
	for (int i = 1; i < P; ++i) {
		fct[i] = fct[i - 1] * ll(i) % P;
	}
	ifct[P - 1] = P - 1;
	for (int i = P - 1; i >= 0; --i) {
		ifct[i - 1] = ll(ifct[i]) * i % P;
	}
	for (int i = 1; i < P; ++i) {
		invvl[i] = ll(fct[i - 1]) * ifct[i] % P;
	}
}

int main() {
	prec();

	ll n, k;
	cin >> n >> k;

	cout << cnk(n, k) << endl;
}