#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include "align.h"
#include "util.h"

void copy_until_space(char *dest, char **src)
{
	size_t i;
	for (i = 0; **src && !isspace(**src); i++) {
		dest[i] = **src;
		(*src)++;
	}

	dest[i] = '\0';
	(*src)++;  // skip the last space
}

char *escape(char *s)  // adapted from BWA
{
	char *p, *q;

	for (p = q = s; *p; ++p) {
		if (*p == '\\') {
			++p;
			if (*p == 't') *q++ = '\t';
			else if (*p == 'n') *q++ = '\n';
			else if (*p == 'r') *q++ = '\r';
			else if (*p == '\\') *q++ = '\\';
		} else *q++ = *p;
	}
	*q = '\0';

	return s;
}

bc_t encode_bc_default(const char *bc)
{
#define BC_ADD_BASE(x) (encoded_bc |= (x))

	bc_t encoded_bc = 0UL;
	char *base = (char *)&bc[BC_LEN-1];
	for (int i = 0; i < BC_LEN; i++) {
		encoded_bc <<= 2;
		switch (*base--) {
		case 'A': case 'a': BC_ADD_BASE(0UL); break;
		case 'C': case 'c': BC_ADD_BASE(1UL); break;
		case 'G': case 'g': BC_ADD_BASE(2UL); break;
		case 'T': case 't': BC_ADD_BASE(3UL); break;
		default: assert(0); break;
		}
	}

	return encoded_bc;

#undef BC_ADD_BASE
}

bc_t encode_bc_haptag(const char *bc)
{
#define CharToInt(c) ((c) - '0')
#define TwoCharToInt(s) ((10 * CharToInt(s[0])) + CharToInt(s[1])) 
#define PackHaplotagString(s) PackHaplotag(TwoCharToInt((s+1)),TwoCharToInt((s+7)),TwoCharToInt((s+4)),TwoCharToInt((s+10)))
#define PackHaplotag(a,b,c,d) ((((uint32_t)(a)) << 24) | (((uint32_t)(c)) << 16) | (((uint32_t)(b)) << 8) | (uint32_t)(d))
	return (bc_t)PackHaplotagString(bc);
}

bc_t encode_bc(const char *bc, const int is_haplotag)
{
	if (is_haplotag) return encode_bc_haptag(bc);
	else return encode_bc_default(bc);
}

void decode_bc_default(bc_t bc, char *out)
{
	for (int i = 0; i < BC_LEN; i++) {
		out[i] = "ACGT"[bc & 0x3];
		bc >>= 2;
	}
}

void decode_bc_haptag(bc_t bc, char *out)
{
	sprintf(out, "A%02uC%02uB%02uD%02u", (bc >> 24) & 127, (bc >> 16) & 127, (bc >> 8) & 127, bc & 127);
}

void decode_bc(bc_t bc, char *out, const int is_haplotag)
{
	if (is_haplotag) decode_bc_haptag(bc, out);
	else decode_bc_default(bc, out);
}

size_t count_lines(FILE *f)
{
	size_t lines = 0;
	while (!feof(f)) {
		if (fgetc(f) == '\n')
			++lines;
	}
	rewind(f);
	return lines + 1;
}

size_t trim_after_space(char *s)
{
	size_t i;
	for (i = 0; s[i] != '\0'; i++) {
		if (isspace(s[i])) {
			s[i] = '\0';
			break;
		}
	}
	return i;
}

uint32_t hash_ident(const char *ident)
{
	uint32_t h = 0;
	for (size_t i = 0; ident[i] != '\0'; i++) {
		h = 31*h + ident[i];
	}
	return h;
}

void normalize_log_probs(double *p, const size_t n)
{
#define EPSILON 1e-50

	if (n == 1) {
		p[0] = 1.0;
		return;
	}

	const double thresh = log(EPSILON) - log(n);
	double p_max = p[0];
	for (size_t i = 1; i < n; i++) {
		const double p_cur = p[i];
		if (p_cur > p_max)
			p_max = p_cur;
	}

	double total = 0;
	for (size_t i = 0; i < n; i++) {
		p[i] -= p_max;

		if (p[i] < thresh)
			p[i] = 0;
		else
			p[i] = exp(p[i]);

		total += p[i];
	}

	for (size_t i = 0; i < n; i++) {
		p[i] /= total;
	}

#undef EPSILON
}

void *safe_malloc(const size_t n)
{
	void *p = malloc(n);
	assert(n == 0 || p != NULL);
	return p;
}

void *safe_calloc(size_t num, size_t size)
{
	void *p = calloc(num, size);
	assert(p != NULL || num == 0 || size == 0);
	return p;
}

void *safe_realloc(void *p, size_t n)
{
	void *new_p = realloc(p, n);
	assert(!(new_p == NULL && n > 0));
	return new_p;
}

void serialize_uint64(FILE *out, const uint64_t x)
{
	/* assumes little-endian (or, at least, that the binary
	   barcode dict files will be generated by the same machine
	   as the one on which they are used) */
	assert(fwrite(&x, sizeof(uint64_t), 1, out));
}

void serialize_uint32(FILE *out, const uint32_t x)
{
	/* assumes little-endian (or, at least, that the binary
	   barcode dict files will be generated by the same machine
	   as the one on which they are used) */
	assert(fwrite(&x, sizeof(uint32_t), 1, out));
}

void serialize_uint8(FILE *out, const uint8_t x)
{
	assert(fwrite(&x, sizeof(uint8_t), 1, out));
}

uint64_t read_uint64(FILE *in)
{
	uint64_t x;
	assert(fread(&x, sizeof(uint64_t), 1, in) == 1);
	return x;
}

uint32_t read_uint32(FILE *in)
{
	uint32_t x;
	assert(fread(&x, sizeof(uint32_t), 1, in) == 1);
	return x;
}

uint8_t read_uint8(FILE *in)
{
	uint8_t x;
	assert(fread(&x, sizeof(uint8_t), 1, in) == 1);
	return x;
}
