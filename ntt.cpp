namespace NTT {
	const int L = 12;
	const int N = (1 << L);
	const int G = 3; // For FFT_MOD = 119 * 2^23 + 1
	const int W0 = pwm(G, 119 * (1 << (23 - L)), FFT_MOD); // W0^{2^L} = 1, W0^{2^{L - 1}} != 1
	const int WI = inv(W0, FFT_MOD);
	const int IN = inv(N, FFT_MOD);

	int a[N], b[N];
	int rev[N];

	void ntt(int *a, bool inv) {
		for (int i = 0; i < N; ++i) {
			if (i < rev[i]) {
				swap(a[i], a[rev[i]]);
			}
		}

		for (int lg = 0; lg < L; ++lg) {
			int bl = (1 << lg);
			int w = (inv ? pwm(WI, 1 << (L - lg - 1), FFT_MOD) : pwm(W0, 1 << (L - lg - 1), FFT_MOD));
			assert(pwm(w, bl, FFT_MOD) != 1);
			assert(pwm(w, 2 * bl, FFT_MOD) == 1);
			for (int pos = 0; pos < N; pos += bl) {
				int ww = 1;
				for (int ptr = 0; ptr < bl; ++ptr, ++pos) {
					int tmp = (a[pos] + ll(ww) * a[pos + bl] % FFT_MOD) % FFT_MOD;
					a[pos + bl] = (a[pos] + FFT_MOD - ll(ww) * a[pos + bl] % FFT_MOD) % FFT_MOD;
					a[pos] = tmp;

					ww = ll(ww) * w % FFT_MOD;
				}
			}
		}

		if (inv) {
			for (int i = 0; i < N; ++i) {
				a[i] = (ll(a[i]) * IN) % FFT_MOD;
			}
		}
	}

	vector<int> mul(int* aa, int* bb, int n, int m) {
		memset(a, 0, sizeof a);
		memset(b, 0, sizeof b);
		for (int i = 0; i < n; ++i) {
			a[i] = aa[i];
		}
		for (int i = 0; i < m; ++i) {
			b[i] = bb[i];
		}

		ntt(a, false);
		ntt(b, false);

		for (int i = 0; i < N; ++i) {
			a[i] = ll(a[i]) * b[i] % FFT_MOD;
		}

		ntt(a, true);
		vector<int> ans;
		for (int i = 0; i < m + n - 1; ++i) {
			ans.push_back(a[i]);
		}
		return ans;
	}


	void prec() {
		rev[0] = 0;
		for (int i = 1; i < N; ++i) {
			rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (L - 1));
		}
	}
};