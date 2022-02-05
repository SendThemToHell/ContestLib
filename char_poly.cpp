const int N = 505;
int mt[N][N];

void sw_rw(int i, int j, int n) {
	for (int k = 0; k < n; ++k) {
		swap(mt[i][k], mt[j][k]);
	}
}

void sw_cl(int i, int j, int n) {
	for (int k = 0; k < n; ++k) {
		swap(mt[k][i], mt[k][j]);
	}
}

void make_upper_hess(int n) {
	for (int j = 0; j < n - 2; ++j) {
		int i;
		for (i = j + 2; i < n; ++i) {
			if (mt[i][j] != 0) {
				break;
			}
		}
		if (i != n) {
			if (mt[j + 1][j] == 0) {
				sw_rw(i, j + 1, n);
				sw_cl(i, j + 1, n);
			}
			for (int k = j + 2; k < n; ++k) {
				int u = int(ll(mt[k][j]) * inv(mt[j + 1][j], FFT_MOD) % FFT_MOD);

				for (int l = 0; l < n; ++l) {
					mt[k][l] = (mt[k][l] + mt[j + 1][l] * ll(FFT_MOD - u) % FFT_MOD) % FFT_MOD;
				}
				for (int l = 0; l < n; ++l) {
					mt[l][j + 1] = (mt[l][j + 1] + mt[l][k] * ll(u) % FFT_MOD) % FFT_MOD;
				}
			}
		}
	}
}

vector<int> poly[N];

vector<int> calc_poly(int n) {
	poly[0].push_back(1);
	for (int k = 0; k < n; ++k) {
		poly[k + 1].push_back(0);
		poly[k + 1].insert(poly[k + 1].end(), all(poly[k]));
		for (int i = 0; i < sz(poly[k]); ++i) {
			poly[k + 1][i] = (poly[k + 1][i] + FFT_MOD - ll(mt[k][k]) * poly[k][i] % FFT_MOD) % FFT_MOD;
		}

		int t = 1;
		for (int i = 1; i <= k; ++i) {
			t = ll(t) * mt[k - i + 1][k - i] % FFT_MOD;
			for (int j = 0; j < sz(poly[k - i]); ++j) {
				poly[k + 1][j] = (poly[k + 1][j] + FFT_MOD - ll(t) * mt[k - i][k] % FFT_MOD * poly[k - i][j] % FFT_MOD) % FFT_MOD;
			}
		}
	}
	return poly[n];
}

/*
run make_upper_hess, then calc_poly
*/