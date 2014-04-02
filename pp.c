#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <libgen.h>

#define MAXN 26
// ppoly[i] is a primitive polynomial of degree i over GF(2)
const unsigned int ppoly[] = {-1, -1, 07, 013, 023, 045, 0103, 0211, 0435, 01021, 02011, 04005, 010123, 020033, 042103, 0100003, 0210013, 0400011, 01000201, 02000047, 04000011, 010000005, 020000003, 040000041, 0100000207, 0200000011};

unsigned int mat[MAXN][MAXN], A[MAXN][MAXN];
char v[1<<(MAXN-1)];

int gcd(int a, int b) {
	return b == 0 ? a : gcd(b, a%b);
}

void swap(unsigned int *a, unsigned int *b) {
	unsigned int t = *a;
	*a = *b;
	*b = t;
}

// A = A * B
void matrix_mul(int n, unsigned int A[][MAXN], unsigned int B[][MAXN]) {
	unsigned int T[MAXN][MAXN];
	int i, j, k;

	for (i=0; i<n; ++i) for (j=0; j<n; ++j) {
		T[i][j] = 0;
		for (k=0; k<n; ++k) {
			T[i][j] ^= A[i][k] * B[k][j];
		}
	}
	for (i=0; i<n; ++i) for (j=0; j<n; ++j) {
		A[i][j] = T[i][j];
	}
}

void human_readable(int deg, unsigned int poly) {
	int i;

	for (i=deg; i>0; --i) if (poly & 1 << i) {
		printf("x^%d + ", i);
	}
	printf("1\n");
}

unsigned int char_poly(int n, unsigned int mat[][MAXN]) {
	int i, j, k;
	unsigned int t[MAXN][MAXN], u[MAXN][MAXN];
	unsigned int ret;

	// Householder reduction to Hessenberg form
	for (i=0; i<n; ++i) for (j=0; j<n; ++j) {
		u[i][j] = mat[i][j];
	}
	for (i=1; i<n; ++i) {
		for (j=i; j<n; ++j) if (u[j][i-1]) break;
		if (j == n) continue;
		for (k=0; k<n; ++k) swap(u[j]+k, u[i]+k);
		for (k=0; k<n; ++k) swap(u[k]+j, u[k]+i);
		for (j=i+1; j<n; ++j) if (u[j][i-1]) {
			for (k=i-1; k<n; ++k) u[j][k] ^= u[i][k];
			for (k=0; k<n; ++k) u[k][i] ^= u[k][j];
		}
	}

	for (j=n-1; j>=0; --j) {
		for (i=0; i<=j; ++i) {
			for (k=n-j-2; k>=0; --k) {
				t[k+1][i] = (u[i][j] * t[k][j+1]) ^ (u[j+1][j] * t[k][i]);
			}
			t[0][i] = u[i][j];
		}
		for (k=0; k<n-j-1; ++k) {
			t[k][j] ^= t[k][j+1];
		}
	}
	ret = 1;
	for (i=0; i<n; ++i) {
		ret = (ret << 1) + t[i][0];
	}

	return ret;
}

// deg: the degree of primitive polynomial
// h: print in a human readable way or a bit compressed way
void pp(int deg, int h) {
	int i, j, n;
	unsigned int cur;

	// A is the companion matrix of ppoly[deg]
	cur = ppoly[deg];
	memset(A, 0, sizeof(A));
	for (i=0; i<deg; ++i) {
		A[i][deg-1] = !!(cur & 1 << i);
	}
	for (i=1; i<deg; ++i) {
		A[i][i-1] = 1;
	}

	for (i=0; i<deg; ++i) for (j=0; j<deg; ++j) {
		mat[i][j] = A[i][j];
	}

	n = (1 << deg) - 1;
	memset(v, 0, (n>>1)*sizeof(char));
	for (i=1; i<n; ++i) {
		if (!v[i] && gcd(i, n) == 1) {
			cur = char_poly(deg, mat);
			if (h) {
				human_readable(deg, cur);
			} else {
				printf("%u\n", cur);
			}
			// mark all conjungate roots
			for (j=i; !v[j]; j=(j<<1)%n) {
				v[j] = 1;
			}
		}
		matrix_mul(deg, mat, A);
	}
}

void usage(char name[]) {
	printf("calculate primitive polynomial over GF(2) of the given degree\n");
	printf("usage: %s degree [-H]\n\n", basename(name));
	printf("\tdegree: the degree of the primitive polynomial\n");
	printf("\t-H(uman-readable) print in a readable way\n");
}

int main(int argc, char *argv[]) {
	int deg, i, H;

	if (argc < 2 || argc > 3) {
		usage(argv[0]);
		return 1;
	}
	if (argc == 3) {
		if (strcmp(argv[2], "-H") != 0 && strcmp(argv[2], "-h") != 0) {
			usage(argv[0]);
			return 1;
		}
	}
	deg = 0;
	for (i=0; argv[1][i]; ++i) {
		assert('0' <= argv[1][i] && argv[1][i] <= '9');
		deg = deg * 10 + argv[1][i] - '0';
	}
	assert(2 <= deg && deg <= MAXN);
	H = argc == 3 && strcmp(argv[2], "-H") == 0;
	pp(deg, H);

	return 0;
}

