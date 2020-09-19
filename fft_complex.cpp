struct cm
{
	ld re, im;

	cm(ld x = 0, ld y = 0) : re(x), im(y) {}

	cm operator*(const cm &oth) {
		return cm(re * oth.re - im * oth.im, re * oth.im + oth.re * im);
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
int rev[N];

void fft(cm *a, bool inv) {
	for (int i = 0; i < N; ++i) {
		if (i < rev[i]) {
			swap(a[i], a[rev[i]]);
		}
	}

	for (int lg = 0; lg < L; ++lg) {
		int bl = (1 << lg);
		for (int pos = 0; pos < N; pos += bl) {
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
		for (int i = 0; i < N; ++i) {
			a[i].re /= N;
		}
	}
}

vector<int> mult(const vector<int>& aa, const vector<int>& bb) {
	memset(a, 0, sizeof a);
	for (int i = 0; i < sz(aa); ++i) {
		a[i].re = aa[i];
	}
	for (int i = 0; i < sz(bb); ++i) {
		a[i].im = bb[i];
	}

	fft(a, false);

	b[0] = { a[0].im, 0 };
	for (int i = 1; i < N; ++i) {
		b[i] = { (a[i].im + a[N - i].im) / 2, (a[N - i].re - a[i].re) / 2 };
	}
	a[0] = { a[0].re, 0 };
	for (int i = 1; i <= N / 2; ++i) {
		cm tmp = a[i];
		a[i] = { (a[i].re + a[N - i].re) / 2, (a[i].im - a[N - i].im) / 2 };
		a[N - i] = { (tmp.re + a[N - i].re) / 2, (a[N - i].im - tmp.im) / 2 };
	}
	for (int i = 0; i < N; ++i) {
		a[i] = a[i] * b[i];
	}

	fft(a, true);
	vector<int> ans;
	for (int i = 0; i < sz(aa) + sz(bb); ++i) {
		ans.push_back(int(a[i].re + 0.5));
	}
	return ans;
}


void prec() {
	for (int i = 1; i < N; ++i) {
		rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (L - 1));
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

int main() {
	ld start = ld(clock());

	prec();

	ld endprec = ld(clock());
	int iter = 20;
	cerr << "precalc time = " << (endprec - start) / CLOCKS_PER_SEC << endl;
	for (int i = 0; i < iter; ++i) {
		vector<int> v1 = { 1, 2, 1 };
		vector<int> v2 = { 1, 3, 3, 1 };
		cout << mult(v1, v2) << endl;
	}
	cerr << "avg multiplication time = " << (ld(clock()) - endprec) / (iter * CLOCKS_PER_SEC) << endl;
	return 0;
}
