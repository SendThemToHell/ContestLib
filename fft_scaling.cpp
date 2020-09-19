struct cm
{
	ld re, im;

	cm(ld x = 0, ld y = 0) : re(x), im(y) {}

	cm operator*(const cm &oth) {
		return cm(re * oth.re - im * oth.im, re * oth.im + oth.re * im);
	}

	cm operator*(ld vl) {
		return cm(re * vl, im * vl);
	}

	cm operator+(const cm &oth) {
		return cm(re + oth.re, im + oth.im);
	}

	cm operator-(const cm &oth) {
		return cm(re - oth.re, im - oth.im);
	}
};

const int L = 20;
const int N = (1 << L);

cm a[N], b[N];
cm w[N], wi[N];
int rev[L + 1][N];

void fft(cm *a, bool inv, int LG) {
	for (int i = 0; i < (1 << LG); ++i) {
		if (i < rev[LG][i]) {
			swap(a[i], a[rev[LG][i]]);
		}
	}

	for (int lg = 0; lg < LG; ++lg) {
		int bl = (1 << lg);
		for (int pos = 0; pos < (1 << LG); pos += bl) {
			cm* ww = inv ? wi + bl - 1 : w + bl - 1;
			for (int ptr = 0; ptr < bl; ++ptr, ++pos) {
				cm tmp = a[pos] + *ww * a[pos + bl];
				a[pos + bl] = a[pos] - *ww * a[pos + bl];
				a[pos] = tmp;

				++ww;
			}
		}
	}

	if (inv) {
		for (int i = 0; i < (1 << LG); ++i) {
			a[i].re /= (1 << LG);
		}
	}
}

/*
Precision-safe version, 5 fft per multiplication

cm c[N], d[N];

vector<ll> mult(const vector<ll>& aa, const vector<ll>& bb, int lg) {
	memset(a, 0, (1 << lg) * sizeof(cm));
	memset(c, 0, (1 << lg) * sizeof(cm));
	for (int i = 0; i < sz(aa); ++i) {
		a[i].re = aa[i] & 1023;
	}
	for (int i = 0; i < sz(bb); ++i) {
		a[i].im = bb[i] & 1023;
	}
	for (int i = 0; i < sz(aa); ++i) {
		c[i].re = aa[i] >> 10;
	}
	for (int i = 0; i < sz(bb); ++i) {
		c[i].im = bb[i] >> 10;
	}

	fft(a, false, lg);
	fft(c, false, lg);

	b[0] = { a[0].im, 0 };
	for (int i = 1; i < (1 << lg); ++i) {
		b[i] = { (a[i].im + a[(1 << lg) - i].im) / 2, (a[(1 << lg) - i].re - a[i].re) / 2 };
	}
	a[0] = { a[0].re, 0 };
	for (int i = 1; i <= (1 << (lg - 1)); ++i) {
		cm tmp = a[i];
		a[i] = { (a[i].re + a[(1 << lg) - i].re) / 2, (a[i].im - a[(1 << lg) - i].im) / 2 };
		a[(1 << lg) - i] = { (tmp.re + a[(1 << lg) - i].re) / 2, (a[(1 << lg) - i].im - tmp.im) / 2 };
	}

	d[0] = { c[0].im, 0 };
	for (int i = 1; i < (1 << lg); ++i) {
		d[i] = { (c[i].im + c[(1 << lg) - i].im) / 2, (c[(1 << lg) - i].re - c[i].re) / 2 };
	}
	c[0] = { c[0].re, 0 };
	for (int i = 1; i <= (1 << (lg - 1)); ++i) {
		cm tmp = c[i];
		c[i] = { (c[i].re + c[(1 << lg) - i].re) / 2, (c[i].im - c[(1 << lg) - i].im) / 2 };
		c[(1 << lg) - i] = { (tmp.re + c[(1 << lg) - i].re) / 2, (c[(1 << lg) - i].im - tmp.im) / 2 };
	}

	for (int i = 0; i < (1 << lg); ++i) {
		cm a1, b1, c1, d1;
		a1 = a[i] * b[i];
		b1 = a[i] * d[i];
		c1 = c[i] * b[i];
		d1 = c[i] * d[i];

		a[i] = a1;
		b[i] = b1 + c1;
		c[i] = d1;
	}

	fft(a, true, lg);
	fft(b, true, lg);
	fft(c, true, lg);

	vector<ll> ans;
	for (int i = 0; i < sz(aa) + sz(bb); ++i) {
		ll vl1 = a[i].re + 0.5;
		ll vl2 = b[i].re + 0.5;
		ll vl3 = c[i].re + 0.5;

		ans.push_back(vl1 + vl2 * 1024LL + vl3 * 1024LL * 1024LL);
	}
	while (ans.back() == 0) {
		ans.pop_back();
	}
	if (sz(ans) == 0) {
		ans.push_back(0);
	}
	return ans;
}
*/

/*
Default version, 2 fft per multiplication
*/
vector<ll> mult(const vector<ll>& aa, const vector<ll>& bb, int lg) {
	memset(a, 0, (1 << lg) * sizeof(cm));
	for (int i = 0; i < sz(aa); ++i) {
		a[i].re = aa[i];
	}
	for (int i = 0; i < sz(bb); ++i) {
		a[i].im = bb[i];
	}

	fft(a, false, lg);

	b[0] = { a[0].im, 0 };
	for (int i = 1; i < (1 << lg); ++i) {
		b[i] = { (a[i].im + a[(1 << lg) - i].im) / 2, (a[(1 << lg) - i].re - a[i].re) / 2 };
	}
	a[0] = { a[0].re, 0 };
	for (int i = 1; i <= (1 << (lg - 1)); ++i) {
		cm tmp = a[i];
		a[i] = { (a[i].re + a[(1 << lg) - i].re) / 2, (a[i].im - a[(1 << lg) - i].im) / 2 };
		a[(1 << lg) - i] = { (tmp.re + a[(1 << lg) - i].re) / 2, (a[(1 << lg) - i].im - tmp.im) / 2 };
	}
	for (int i = 0; i < (1 << lg); ++i) {
		a[i] = a[i] * b[i];
	}

	fft(a, true, lg);
	vector<ll> ans;
	for (int i = 0; i < sz(aa) + sz(bb); ++i) {
		ans.push_back(ll(a[i].re + 0.5));
	}

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
	int ptr = 0;
	for (int lg = 0; lg < L; ++lg) {
		int bl = (1 << lg);
		for (int i = 0; i < bl; ++i) {
			w[ptr] = cm(cos(i * PI / bl), sin(i * PI / bl));
			wi[ptr] = cm(w[ptr].re, -w[ptr].im);
			++ptr;
		}
	}
}