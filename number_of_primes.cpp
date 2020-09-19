int setBit(int n, int pos) { return n = n | (1 << pos); }
int resetBit(int n, int pos) { return n = n & ~(1 << pos); }
bool checkBit(ll n, ll pos) { return (bool)(n & (1LL << pos)); }


namespace pcf {
	///   Prime-Counting Function
	///   initialize once by calling init()
	///   Legendre(n) and Lehmer(n) returns the number of primes less than or equal to m
	///   Lehmer(n) is faster

#define MAXN 1000010 /// initial sieve limit
#define MAX_PRIMES 1000010 /// max size of the prime array for sieve
#define PHI_N 100000 ///
#define PHI_K 100

	unsigned int ar[(MAXN >> 6) + 5] = { 0 };
	int len = 0; /// total number of primes generated by sieve
	int primes[MAX_PRIMES];
	int counter[MAXN]; /// counter[m] --> number of primes <= i
	int phi_dp[PHI_N][PHI_K]; /// precal of phi(n,k)

	bool isComp[MAXN];  // ara[i] is true if i is composite

	void Sieve(int N, bool* isComp, int* primes, int& len) {
		int  i, j, sq = sqrt(N);
		isComp[1] = true;
		for (i = 4; i <= N; i += 2) isComp[i] = true;
		for (i = 3; i <= sq; i += 2) {
			if (!isComp[i]) {
				for (j = i * i; j <= N; j += i + i) isComp[j] = 1;
			}
		}
		for (i = 1; i <= N; i++) {
			if (!isComp[i]) primes[len++] = i;
			counter[i] = len;
		}
	}

	void Sieve(int n) {
		Sieve(n, isComp, primes, len);
	}

	void init() {
		Sieve(MAXN - 1);

		/// precalculation of phi upto size (PHI_N,PHI_K)
		int k, n, res;
		for (n = 0; n < PHI_N; n++) phi_dp[n][0] = n;
		for (k = 1; k < PHI_K; k++) {
			for (n = 0; n < PHI_N; n++) {
				phi_dp[n][k] = phi_dp[n][k - 1] - phi_dp[n / primes[k - 1]][k - 1];
			}
		}
	}

	/// returns number of integers less or equal n which are
	/// not divisible by any of the first k primes
	/// recurrence --> phi( n , k ) = phi( n , k-1 ) - phi( n / p_k , k-1)
	long long phi(long long n, int k) {
		if (n < PHI_N && k < PHI_K) return phi_dp[n][k];
		if (k == 1) return ((++n) >> 1);
		if (primes[k - 1] >= n) return 1;
		return phi(n, k - 1) - phi(n / primes[k - 1], k - 1);
	}


	long long Legendre(long long n) {
		if (n < MAXN) return counter[n];

		int lim = sqrt(n) + 1;
		int k = upper_bound(primes, primes + len, lim) - primes;
		return phi(n, k) + (k - 1);
	}


	long long Lehmer(long long n) {
		if (n < MAXN) return counter[n];

		long long w, res = 0;
		int i, j, a, b, c, lim;
		b = sqrt(n), c = Lehmer(cbrt(n)), a = Lehmer(sqrt(b)), b = Lehmer(b);
		res = phi(n, a) + (((b + a - 2) * (b - a + 1)) >> 1);

		for (i = a; i < b; i++) {
			w = n / primes[i];
			lim = Lehmer(sqrt(w)), res -= Lehmer(w);

			if (i <= c) {
				for (j = i; j < lim; j++) {
					res += j;
					res -= Lehmer(w / primes[j]);
				}
			}
		}
		return res;
	}
}
