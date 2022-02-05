const int L = 21;
const int N = (1 << L);
const int W0 = pwm(3, (FFT_MOD - 1) / (1 << L), FFT_MOD);
const int WI0 = inv(W0, FFT_MOD);

int a[N], b[N];
int rev[L + 1][N];

void fft(int *a, bool inv, int LG) { // LG <= L - 1
	for (int i = 0; i < (1 << LG); ++i) {
		if (i < rev[LG][i]) {
			swap(a[i], a[rev[LG][i]]);
		}
	}

	for (int lg = 0; lg < LG; ++lg) {
		int bl = (1 << lg);
		int ml = inv ? pwm(WI0, 1 << (L - 1 - lg), FFT_MOD) : pwm(W0, 1 << (L - 1 - lg), FFT_MOD);
		for (int pos = 0; pos < (1 << LG); pos += bl) {
			int ww = 1;
			for (int ptr = 0; ptr < bl; ++ptr, ++pos) {
				int hui = ww * ll(a[pos + bl]) % FFT_MOD;
				int tmp = a[pos] + hui >= FFT_MOD ? a[pos] + hui - FFT_MOD : a[pos] + hui;
				a[pos + bl] = a[pos] - hui < 0 ? a[pos] - hui + FFT_MOD : a[pos] - hui;
				a[pos] = tmp;

				ww = ll(ww) * ml % FFT_MOD;
			}
		}
	}


	if (inv) {
		int in = ::inv(1 << LG, FFT_MOD);
		for (int i = 0; i < (1 << LG); ++i) {
			a[i] = ll(a[i]) * in % FFT_MOD;
		}
	}
}

vector<int> mult(const vector<int>& aa, const vector<int>& bb, int lg) {
	if (lg <= 8) {
		vector<int> ans(sz(aa) + sz(bb) - 1);
		for (int i = 0; i < sz(aa); ++i) {
			for (int j = 0; j < sz(bb); ++j) {
				ans[i + j] = (ans[i + j] + ll(aa[i]) * bb[j]) % FFT_MOD;
			}
		}
		while (ans.back() == 0) {
			ans.pop_back();
		}
		if (sz(ans) == 0) {
			ans.push_back(0);
		}

		return ans;
	}

	memset(a, 0, (1 << lg) * sizeof(ll));
	memset(b, 0, (1 << lg) * sizeof(ll));

	for (int i = 0; i < sz(aa); ++i) {
		a[i] = aa[i];
	}
	for (int i = 0; i < sz(bb); ++i) {
		b[i] = bb[i];
	}

	fft(a, false, lg);
	fft(b, false, lg);

	for (int i = 0; i < (1 << lg); ++i) {
		a[i] = a[i] * ll(b[i]) % FFT_MOD;
	}

	fft(a, true, lg);
	vector<int> ans(a, a + sz(aa) + sz(bb));

	while (ans.back() == 0) {
		ans.pop_back();
	}
	if (sz(ans) == 0) {
		ans.push_back(0);
	}

	return ans;
}

void prec() {
	for (int lg = 1; lg <= L; ++lg) {
		for (int i = 1; i < (1 << lg); ++i) {
			rev[lg][i] = (rev[lg][i >> 1] >> 1) | ((i & 1) << (lg - 1));
		}
	}
}