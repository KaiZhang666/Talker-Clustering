#if 0 //To compile as a stand-alone test program, type "nmake fft.c"
fft.exe: fft.c; cl fft.c /Ox /G6 /MD /DFFTTEST /link /opt:nowin98
!if 0
#endif

/*****************************************************************************
10/10/2002. Ken Silverman's optimized FFT code for power-of-2 block sizes.
Requirements:
	CPU: Intel Pentium III, Pentium 4 (and above) or AMD AthlonXP (and above)
	OS: Microsoft Windows (95?)/98/ME/NT/2K/XP
	Compilers: Microsoft Visual C/C++ 6.0 (VC6) with Processor Pack installed /
				  Microsoft Visual Studio .Net (VC7)

FFT code optimized for PURE REAL input. Currently, this code is a stand-alone
	program. To compile, type "nmake fft.c". It is easy to make this an
	externally callable library. Just copy these function declarations into
	your code:

	//Real value optimized FFT. dir:1=forward fft,-1=inverse fft
extern void realfft (float *snd, float *real, float *imag, long log2blocksiz,
	long dir);

	//In-place complex-to-complex FFT. dir=1:forward fft, dir=-1:inverse fft
extern void compfft (float *real, float *imag, long log2blocksiz, long dir);

	//Fast Discrete Cosine Transform (type 2). dir:1=forward fft,-1=inverse fft
extern void fastdct (float *tim, float *frq, long log2blocksiz, long dir);

----------------------------------- NOTES: -----------------------------------
	  snd: the time domain signal
	 real: the    real   component of the frequency domain signal
	 imag: the imaginary component of the frequency domain signal
blocksiz: can be anything that fits in memory, but MUST be a power of 2

You must allocate (blocksiz/2+1) indices for the real/imag buffers. To see why,
study this example diagram for a 4-point FFT:

	?TIME DOMAIN: (blocksiz = 4) ?FREQUENCY DOMAIN: 

	                            real[0]          
	snd[0] snd[1] snd[2] snd[3] real[1] imag[1]  
	                            real[2]          

NOTE: imag[0] & imag[2] are always 0 after a forward transform (due to the
		symmetry of PURE REAL transform).

Here's a snippet of C code which shows how to use realfft():
	__declspec(align(16)) float snd[1024], real[513], imag[513];
	for(i=0;i<1024;i++) snd[i] = (float)(rand()-16384);
	realfft(snd,real,imag,10,1);
	for(i=0;i<=512;i++) printf("%f + %fj\n",real[i],imag[i]);

------------------------------------------------------------------------------
Use either MS VC6.0 or VC Studio .Net to compile. For VC6, you'll need to
install the Processor Pack.  There some issues with the installation of the
Processor Pack, so follow these instructions:

1. Go here: http://msdn.microsoft.com/vstudio/downloads/ppack/default.asp
	download the Processor Pack and start the install.
2. If it gives no error message, then ignore the rest of these instructions.
	If it gives an error message like "your version not supported", then follow
	these additional steps: BEFORE you click 'ok' to this error message, go to
	the command prompt and do this:
3. Find the temporary directory where it de-compressed the files. The file
	you're interested in is: "C2.DLL"
4. Replace the current C2.DLL in your Visual C's "bin" directory with this
	file.  Once you do this, your compiler is upgraded.

I'm still a DOS fan, so I prefer compiling on the command line.  To compile
the code as a stand-alone program inside the VC environment, be sure to select
"Win32 Console application" and make sure "FFTTEST" is defined. You can either
do this in the makefile with the /D option or as a #define at the top of the
program. To compile the code as a library, just don't define "FFTTEST".

------------------------------------------------------------------------------
Because of my extensive use of Pentium III assembly code (SSE), all arrays
passed to my FFT functions must be aligned on a 16-byte boundary. If you don't
do this, you're likely to get a crash inside my code. Alignment is not hard to
do. To align static arrays, just stick "__declspec(align(16))" in front of the
array declaration, like this:

	__declspec(align(16)) float snd[1024];

For dynamically allocated arrays, you can either use alignmalloc(numbytes,16)/
alignfree from my FFT.C source code, _aligned_malloc(numbytes,16)/
_aligned_free from the Processor Pack or align it yourself by allocating 15
extra bytes and using this formula:

	new_pointer = (((long)pointer+15)&~15);

I could make my library handle unaligned arrays, but this would require extra
steps which would slow it down.

------------------------------------------------------------------------------
Ken's official website: http://www.advsys.net/ken
*****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <memory.h>

//========================== ALIGNED MALLOC BEGINS ==========================
#define ALIGNMALLOCMAX 256
static void *mallocorig[ALIGNMALLOCMAX], *mallocuse[ALIGNMALLOCMAX];
static long malloccnt = 0;
void alignfree (void *v)
{
	long i;
	for(i=malloccnt-1;i>=0;i--)
		if (mallocuse[i] == v)
		{
			free(mallocorig[i]);
			malloccnt--;
			mallocorig[i] = mallocorig[malloccnt];
			mallocuse[i] = mallocuse[malloccnt];
			return;
		}
}
void *alignmalloc (long n, long align)
{
	void *ptr;
	if (malloccnt >= ALIGNMALLOCMAX) return(0);
	ptr = malloc(n+align-1); if (!ptr) return(0);
	mallocorig[malloccnt] = ptr;
	mallocuse[malloccnt] = (void *)((((long)ptr)+align-1)&~(align-1));
	//memset(mallocuse[malloccnt],205,n); //HACK for debugging
	malloccnt++;
	return(mallocuse[malloccnt-1]);
}
//=========================== ALIGNED MALLOC ENDS ===========================

#if (FFTTEST == 1)
static long lb0, lb1;
#endif

#define PI 3.141592653589793
#define SQRT2O2 0.7071067811865476
#define COSPIO8 0.9238795325112867
#define COS3PIO8 .3826834323650898

	//Note: If not using compfft, these array sizes can all be halved!
static long cbl = 0, *bitrev = 0, *bitrev2 = 0, *bitrevreal;
static float *sinlut = 0, *coslut, *sinlut2 = 0, *coslut2;
static float *radix4tab = 0, oneovern, twoovern;
static void (*fftfunc)(float *, float *, long);

	//  //How to do a printf inside _asm block:
	//static const char snoto[] = "%f %f %f %f\n";
	//   push eax ;----
	//   push edx ;----
	//   push ecx ;----
	//   mov ecx, 4
	//l0:shufps xmm0, xmm0, 0x93
	//   movss f, xmm0
	//   fld dword ptr f
	//   fstp qword ptr [esp-8]
	//   sub esp, 8
	//   dec ecx
	//   jnz short l0
	//   push offset snoto
	//   call printf
	//   add esp, 8*4+4
	//   pop ecx  ;----
	//   pop edx  ;----
	//   pop eax  ;----

#define FFTRADIX4ASM 1
static __declspec(align(16)) long neg01[4] = {0x80000000,0x80000000,0,0};
static __declspec(align(16)) long neg12[4] = {0,0x80000000,0x80000000,0};
static void fftradix4 (float *r, float *i, long n)
{
#if (FFTRADIX4ASM == 0)
	long i0, i1, i2, i3, pw, a, aa, ainc, msk, pw3;
	float f, g, r1, q1, r2, q2, r3, q3, *nr, *ni;

	pw = 1; a = 0; ainc = (n>>1); msk = ainc-1;
	for(i0=n-4;i0>=0;i0-=4)
	{
		nr = &r[i0]; ni = &i[i0];
		r1 = nr[1]; r2 = nr[2]; r3 = nr[3];
		q1 = ni[1]; q2 = ni[2]; q3 = ni[3];
		f = nr[0]-r2; g = q1-q3; nr[1] = f+g; nr[3] = f-g;
		f = ni[0]-q2; g = r1-r3; ni[1] = f-g; ni[3] = f+g;
		f = nr[0]+r2; g = r1+r3; nr[0] = f+g; nr[2] = f-g;
		f = ni[0]+q2; g = q1+q3; ni[0] = f+g; ni[2] = f-g;
	}
	while (ainc > 4)
	{
		pw <<= 2; ainc >>= 2; pw3 = ~(pw+pw+pw);
		i0 = n-1;
		do
		{
			a = ((a-ainc)&msk);
			i0 &= pw3; i1 = i0+pw; i2 = i1+pw; i3 = i2+pw;
			r1 = coslut[ a]*r[i1] + sinlut[ a]*i[i1];
			q1 = coslut[ a]*i[i1] - sinlut[ a]*r[i1]; aa = a+a;
			r2 = coslut[aa]*r[i2] + sinlut[aa]*i[i2];
			q2 = coslut[aa]*i[i2] - sinlut[aa]*r[i2]; aa = aa+a;
			r3 = coslut[aa]*r[i3] + sinlut[aa]*i[i3];
			q3 = coslut[aa]*i[i3] - sinlut[aa]*r[i3];
			f = r[i0]-r2; g = q1-q3; r[i1] = f+g; r[i3] = f-g;
			f = r[i0]+r2; g = r1+r3; r[i0] = f+g; r[i2] = f-g;
			f = i[i0]-q2; g = r1-r3; i[i1] = f-g; i[i3] = f+g;
			f = i[i0]+q2; g = q1+q3; i[i0] = f+g; i[i2] = f-g;
			i0--;
		} while (i0 >= 0);
	}
#else
	float *fptr, *fend;
	long pw, pw3;
	_asm
	{
		mov eax, n
		sub eax, 4
		mov ecx, r
		mov edx, i
		movaps xmm6, neg12
		movaps xmm7, neg01
begit4a:
		movaps xmm0, [ecx+eax*4] ;xmm0: r3 r2 r1 r0
		movaps xmm1, [edx+eax*4] ;xmm1: i3 i2 i1 i0
		movaps xmm2, xmm0        ;xmm2: r3 r2 r1 r0
		movaps xmm3, xmm1        ;xmm3: i3 i2 i1 i0
		shufps xmm2, xmm2, 0x4e  ;xmm2: r1 r0 r3 r2
		shufps xmm3, xmm3, 0x4e  ;xmm3: i1 i0 i3 i2
		xorps xmm2, xmm7         ;xmm2: r1 r0 -r3 -r2
		xorps xmm3, xmm7         ;xmm3: i1 i0 -i3 -i2
		addps xmm2, xmm0         ;xmm2: r1+r3 r0+r2 r1-r3 r0-r2
		addps xmm3, xmm1         ;xmm3: i1+i3 i0+i2 i1-i3 i0-i2
		shufps xmm2, xmm2, 0x72  ;xmm2: r1-r3 r1+r3 r0-r2 r0+r2     2_1 2_3 2_0 2_2
		shufps xmm3, xmm3, 0x72  ;xmm3: i1-i3 i1+i3 i0-i2 i0+i2     3_1 3_3 3_0 3_2
		movaps xmm4, xmm2        ;xmm4: r1-r3 r1+r3 r0-r2 r0+r2     2_1 2_3 2_0 2_2
		unpckhps xmm4, xmm3      ;xmm4: i1-i3 r1-r3 i1+i3 r1+r3     2_1 3_1 2_3 3_3
		movaps xmm5, xmm4        ;xmm5: i1-i3 r1-r3 i1+i3 r1+r3     2_1 3_1 2_3 3_3
		shufps xmm4, xmm4, 0xcc  ;xmm4: i1-i3 r1+r3 i1-i3 r1+r3     3_1 2_3 3_1 2_3
		shufps xmm5, xmm5, 0x99  ;xmm5: r1-r3 i1+i3 r1-r3 i1+i3     2_1 3_3 2_1 3_3
		movlhps xmm2, xmm2       ;xmm2: r0-r2 r0+r2 r0-r2 r0+r2     2_0 2_2 2_0 2_2
		movlhps xmm3, xmm3       ;xmm3: i0-i2 i0+i2 i0-i2 i0+i2     3_0 3_2 3_0 3_2
		xorps xmm4, xmm7         ;xmm4: i1-i3 r1+r3 -(i1-i3) -(r1+r3)    3_1 2_3-3_1-2_3
		xorps xmm5, xmm6         ;xmm5: r1-r3 -(i1+i3) -(r1-r3) i1+i3    2_1-3_3-2_1 3_3
		subps xmm2, xmm4         ;xmm2: (r0-r2)-(i1-i3) (r0+r2)-(r1+r3) (r0-r2)+(i1-i3) (r0+r2)+(r1+r3)
		addps xmm3, xmm5         ;xmm3: (i0-i2)+(r1-r3) (i0+i2)-(i1+i3) (i0-i2)-(r1-r3) (i0+i2)+(i1+i3)
		movaps [ecx+eax*4], xmm2
		movaps [edx+eax*4], xmm3
		sub eax, 4
		jnc short begit4a
	}
	fptr = radix4tab;
	for(pw=4;pw<n;pw<<=2)
	{
		pw3 = pw+pw+pw; fend = &fptr[pw3+pw3]; pw3 = ~pw3;
		_asm
		{
			push ebx
			push esi
			push edi

			mov ecx, n
			mov ebx, pw3
			mov edx, fptr
			mov eax, pw
			sub ecx, 4
			not ebx
begit4:
			and ecx, pw3
			mov esi, r
			mov edi, i
			lea esi, [esi+ecx*4] ;r[i0]:[esi]  r[i1]:[esi+eax*4]  r[i2]:[esi+eax*8]  r[i3]:[esi+ebx*4]
			lea edi, [edi+ecx*4] ;i[i0]:[edi]  i[i1]:[edi+eax*4]  i[i2]:[edi+eax*8]  i[i3]:[edi+ebx*4]

			movaps xmm1, [edx]
			movaps xmm3, [edx+16]
			movaps xmm0, [esi+eax*4] ;r[i1]
			mulps xmm0, xmm1
			movaps xmm2, [edi+eax*4] ;i[i1]
			mulps xmm2, xmm3
			mulps xmm1, [edi+eax*4]  ;i[i1]
			mulps xmm3, [esi+eax*4]  ;r[i1]
			addps xmm0, xmm2         ;      =r1[3] r1[2] r1[1] r1[0]
			subps xmm1, xmm3         ;      =q1[3] q1[2] q1[1] q1[0]

			movaps xmm3, [edx+32]
			movaps xmm2, [esi+eax*8] ;r[i2]
			mulps xmm2, xmm3
			movaps xmm5, [edx+48]
			movaps xmm4, [edi+eax*8] ;i[i2]
			mulps xmm4, xmm5
			mulps xmm3, [edi+eax*8]  ;i[i2]
			mulps xmm5, [esi+eax*8]  ;r[i2]
			addps xmm2, xmm4         ;      =r2[3] r2[2] r2[1] r2[0]
			subps xmm3, xmm5         ;      =q2[3] q2[2] q2[1] q2[0]

			movaps xmm5, [edx+64]
			movaps xmm4, [esi+ebx*4] ;r[i3]
			mulps xmm4, xmm5
			movaps xmm7, [edx+80]
			movaps xmm6, [edi+ebx*4] ;i[i3]
			mulps xmm6, xmm7
			mulps xmm5, [edi+ebx*4]  ;i[i3]
			mulps xmm7, [esi+ebx*4]  ;r[i3]
			addps xmm4, xmm6         ;      =r3[3] r3[2] r3[1] r3[0]
			subps xmm5, xmm7         ;      =q3[3] q3[2] q3[1] q3[0]

			movaps xmm6, [esi]       ;r[i0]
			subps xmm6, xmm2
			movaps xmm7, xmm6
			addps xmm6, xmm1         ;out of registers, so do an extra add
			subps xmm7, xmm1
			add edx, 24*4
			subps xmm6, xmm5
			addps xmm7, xmm5
			cmp edx, fend
			movaps [esi+eax*4], xmm6 ;r[i1]
			movaps [esi+ebx*4], xmm7 ;r[i3]

			addps xmm2, [esi]        ;r[i0]  r2 destroyed (xmm2 free after this block)
			movaps xmm7, xmm0
			addps xmm7, xmm4
			movaps xmm6, xmm2
			addps xmm2, xmm7
			cmovae edx, fptr
			subps xmm6, xmm7
			movaps [esi], xmm2       ;r[i0]
			movaps [esi+eax*8], xmm6 ;r[i2]

			movaps xmm6, [edi]       ;i[i0]
			subps xmm6, xmm3
			movaps xmm7, xmm0
			subps xmm7, xmm4
			movaps xmm2, xmm6
			subps xmm6, xmm7
			sub ecx, 4
			addps xmm2, xmm7
			movaps [edi+eax*4], xmm6 ;i[i1]
			movaps [edi+ebx*4], xmm2 ;i[i3]

			addps xmm3, [edi]        ;i[i0]
			addps xmm1, xmm5
			movaps xmm5, xmm3
			addps xmm3, xmm1
			subps xmm5, xmm1
			movaps [edi], xmm3       ;i[i0]
			movaps [edi+eax*8], xmm5 ;i[i2]

			jnc short begit4

			pop edi
			pop esi
			pop ebx
		}
		fptr = fend;
	}
#endif
}

static void fftradix4odd (float *r, float *i, long nn)
{
#if (FFTRADIX4ASM == 0)
	float *osinlut, *ocoslut;
#endif
	float f, g;
	long z, j;

	nn >>= 1;

#if (FFTRADIX4ASM == 0)
	osinlut = sinlut; sinlut = sinlut2;
	ocoslut = coslut; coslut = coslut2;
#endif
	fftradix4(r,i,nn);
	fftradix4(&r[nn],&i[nn],nn);
#if (FFTRADIX4ASM == 0)
	sinlut = osinlut; coslut = ocoslut;
#endif

#if 0
	for(z=nn-1;z>=0;z--)
	{
		j = z+nn;
		f = r[j]*coslut2[z] + i[j]*sinlut2[z];
		g = i[j]*coslut2[z] - r[j]*sinlut2[z];
		r[j] = r[z]-f; r[z] += f;
		i[j] = i[z]-g; i[z] += g;
	}
#else
	_asm
	{
		push esi
		push edi
		push ebx
		mov ecx, nn
		mov esi, r
		mov edi, i
		mov eax, coslut2
		mov ebx, sinlut2
		lea edx, [ecx+ecx-4]
		sub ecx, 4
beg4o:movaps xmm0, [eax+ecx*4] ;xmm0: c3 c2 c1 c0
		movaps xmm1, [ebx+ecx*4] ;xmm1: s3 s2 s1 s0
		movaps xmm2, xmm0
		movaps xmm3, xmm1
		mulps xmm0, [esi+edx*4]  ;xmm0: c3*r3 c2*r2 c1*r1 c0*r0
		mulps xmm1, [edi+edx*4]  ;xmm1: s3*i3 s2*i2 s1*i1 s0*i0
		mulps xmm2, [edi+edx*4]  ;xmm2: c3*i3 c2*i2 c1*i1 c0*i0
		mulps xmm3, [esi+edx*4]  ;xmm3: s3*r3 s2*r2 s1*r1 s0*r0
		addps xmm0, xmm1         ;xmm0: f3 f2 f1 f0
		subps xmm2, xmm3         ;xmm2: g3 g2 g1 g0
		movaps xmm4, [esi+ecx*4] ;xmm4: r3 r2 r1 r0
		movaps xmm6, [edi+ecx*4] ;xmm6: i3 i2 i1 i0
		movaps xmm5, xmm4        ;xmm5: r3 r2 r1 r0
		movaps xmm7, xmm6        ;xmm7: i3 i2 i1 i0
		subps xmm4, xmm0         ;xmm4: r3-f3 r2-f2 r1-f1 r0-f0
		addps xmm5, xmm0         ;xmm5: i3+f3 i2+f2 i1+f1 i0+f0
		subps xmm6, xmm2         ;xmm6: r3-g3 r2-g2 r1-g1 r0-g0
		addps xmm7, xmm2         ;xmm7: i3+g3 i2+g2 i1+g1 i0+g0
		movaps [esi+edx*4], xmm4
		movaps [esi+ecx*4], xmm5
		movaps [edi+edx*4], xmm6
		sub edx, 4
		movaps [edi+ecx*4], xmm7
		sub ecx, 4
		jnc short beg4o
		pop ebx
		pop edi
		pop esi
	}
#endif
}

	//ai = radians/(2? (or ratio of angle step around full circle)
	//for(i=0;i<num;i++) { sintabptr[i] = sin((double)i*ai*PI*2.0);
static __inline void fcalcsinasm (float *sintabptr, long num, double ai)
{
		//Excellent precision: Based on code from WAVY.C, but stores to float
		//double s, si, f; long i;
		//ai *= PI; f = sin(ai)*2.0; s = 0; si = cos(ai)*f; f *= f;
		//for(i=0;i<num;i++) { sintabptr[i] = s; s += si; si -= s*f; }
	_asm
	{
		fldpi                   ;?
		fmul qword ptr ai       ;ai*?
		fsincos                 ;cos(ai*? sin(ai*?
		fxch st(1)              ;sin(ai*? cos(ai*?
		fadd st, st             ;sin(ai*?*2 cos(ai*?
		fmul st(1), st          ;sin(ai*?*2 cos(ai*?*sin(ai*?*2
		fmul st, st             ;(sin(ai*?*2)?cos(ai*?*sin(ai*?*2
		fldz                    ;0 (sin(ai*?*2)?cos(ai*?*sin(ai*?*2
		mov eax, num            ;s f si
		mov edx, sintabptr
		shl eax, 2
		add edx, eax
		neg eax
 beg: fst dword ptr [eax+edx] ;s f si
		fadd st, st(2)          ;s+si f si
		fld st(1)               ;f s+si f si
		fmul st, st(1)          ;f*(s+si) s+si f si
		fsubp st(3), st         ;s+si f si-f*(s+si)
		add eax, 4
		jnz short beg
		fucompp                 ;st(0)
		fstp st                 ;
	}
}

void initffttables (long m)
{
	long i, j, k, n, no1, no2, no3, a, aa, ainc, pw;
	float *fptr;

	if ((cbl == m) || (m <= 1)) return;
	n = (1<<m); cbl = m; no1 = (n>>1); no2 = (n>>2); no3 = (n>>3);
	oneovern = 1.f / (float)n; twoovern = oneovern+oneovern;

	if (radix4tab) { alignfree(radix4tab); radix4tab = 0; }
	if (bitrev2) { alignfree(bitrev2); bitrev2 = 0; }
	if (bitrev) { alignfree(bitrev); bitrev = 0; }
	if (sinlut) { alignfree(sinlut); sinlut = 0; }
	i = n+no2; if (i < 16) i = 16;
	sinlut = (float *)alignmalloc(i*sizeof(float),16); if (!sinlut) { cbl = 0; return; }
	bitrev = (long *)alignmalloc(no1*sizeof(long),16); if (!bitrev) { cbl = 0; return; }
#if (FFTTEST == 1)
	if (lb0 == lb1)
	{
		printf("   sinlut: %d bytes\n",i*sizeof(float));
		printf("   bitrev: %d bytes\n",no1*sizeof(long));
	}
#endif

	//for(i=0;i<=no2;i++) sinlut[i] = sin((double)i*PI*2.0/(float)n);
	fcalcsinasm(sinlut,no2+1,oneovern);
	for(i=no2+1;i<no1;i++) sinlut[i] = sinlut[no1-i];
	for(i=no1;i<n;i++) sinlut[i] = -sinlut[i-no1];
	if (no2 >= 4)
	{
		coslut = &sinlut[no2];
		for(i=0;i<no2;i++) sinlut[i+n] = sinlut[i];
	}
	else //This ensures coslut is 16-byte aligned (for very small block sizes) :/
	{
		coslut = &sinlut[n];
		for(i=n-1;i>=0;i--) coslut[i] = sinlut[(i+no2)&(n-1)];
	}

	bitrev[0] = 0;
	if (m&1)
	{
		if (n >= 32)
		{
			j = 0;
			for(ainc=(no2>>2),pw=1;ainc>1;ainc>>=2,pw<<=2) j += pw;
			radix4tab = (float *)alignmalloc(j*24*sizeof(float),16); if (!radix4tab) { cbl = 0; return; }
#if (FFTTEST == 1)
			if (lb0 == lb1) printf("radix4tab: %d bytes\n",j*24*sizeof(float));
#endif
			ainc = no2; fptr = radix4tab;
			for(ainc=(no2>>2),pw=1;ainc>1;ainc>>=2,pw<<=2)
				for(i=(ainc<<2),aa=no2-i;aa>=0;aa-=i,fptr+=24)
					for(a=aa,k=0;k<4;k++)
					{
						fptr[k+ 0] = coslut[a]; fptr[k+ 4] = sinlut[a]; j = a+a;
						fptr[k+ 8] = coslut[j]; fptr[k+12] = sinlut[j]; j += a;
						fptr[k+16] = coslut[j]; fptr[k+20] = sinlut[j]; a += ainc;
					}
		}

		fftfunc = fftradix4;
		for(j=1,k=no3;k>0;j<<=2,k>>=2)
			for(i=0;i<j*(4-1);i++) bitrev[i+j] = bitrev[i]+k;

		bitrevreal = bitrev;
	}
	else
	{
		if (sinlut2) { alignfree(sinlut2); sinlut2 = 0; }
		i = no1+no3; if (i < 16) i = 16;
		sinlut2 = (float *)alignmalloc(i*sizeof(float),16); if (!sinlut2) { cbl = 0; return; }
#if (FFTTEST == 1)
		if (lb0 == lb1) printf("  sinlut2: %d bytes\n",i*sizeof(float));
#endif
		if (no3 >= 4)
		{
			coslut2 = &sinlut2[no3];
			for(i=no1+no3-1;i>=0;i--) sinlut2[i] = sinlut[i+i];
		}
		else //This ensures coslut2 is 16-byte aligned (for very small block sizes) :/
		{
			coslut2 = &sinlut2[no1];
			for(i=no1-1;i>=0;i--) { sinlut2[i] = sinlut[i+i]; coslut2[i] = coslut[i+i]; }
		}

		if (n >= 64)
		{
			j = 0;
			for(ainc=(no3>>2),pw=1;ainc>1;ainc>>=2,pw<<=2) j += pw;
			radix4tab = (float *)alignmalloc(j*24*sizeof(float),16); if (!radix4tab) { cbl = 0; return; }
#if (FFTTEST == 1)
			if (lb0 == lb1) printf("radix4tab: %d bytes\n",j*24*sizeof(float));
#endif
			ainc = no3; fptr = radix4tab;
			for(ainc=(no3>>2),pw=1;ainc>1;ainc>>=2,pw<<=2)
				for(i=(ainc<<2),aa=no3-i;aa>=0;aa-=i,fptr+=24)
					for(a=aa,k=0;k<4;k++)
					{
						fptr[k+ 0] = coslut2[a]; fptr[k+ 4] = sinlut2[a]; j = a+a;
						fptr[k+ 8] = coslut2[j]; fptr[k+12] = sinlut2[j]; j += a;
						fptr[k+16] = coslut2[j]; fptr[k+20] = sinlut2[j]; a += ainc;
					}
		}

		fftfunc = fftradix4odd;
		for(j=1,k=no2;k>0;j<<=1,k>>=1)
			for(i=0;i<j*(2-1);i++) bitrev[i+j] = bitrev[i]+k;

		bitrev2 = (long *)alignmalloc(no1*sizeof(long),16); if (!bitrev2) { cbl = 0; return; }
#if (FFTTEST == 1)
		if (lb0 == lb1) printf("  bitrev2: %d bytes\n",no1*sizeof(long));
#endif
		for(i=no2-1;i>=0;i--)
		{
			j = bitrev[((i&0xaaaaaaaa)>>1) + ((i&0x55555555)<<1)];
			bitrev2[i    ] = j;
			bitrev2[i+no2] = j+1;
		}
		bitrevreal = bitrev2;
	}
}

#define WASTE_32K_TO_SPEEDUP_COMPFFT 0
#if (WASTE_32K_TO_SPEEDUP_COMPFFT == 1)
static float oldr[8192], oldi[8192];
#endif

	 //In-place complex-to-complex FFT
	 //   dir=1:forward, dir=-1:inverse
	 //   m: log2(n)
	 //   r/i: real/imag
void compfft (float *r, float *i, long m, long dir)
{
	float f, g, h;
	long j, n, z;

	if (m < 2)
	{
		if (m <= 0) return;
		f = r[1]; g = i[1];
		if (dir > 0)
		{
			r[1] = r[0]-f; r[0] += f;
			i[1] = i[0]-g; i[0] += g;
		}
		else
		{
			r[1] = (r[0]-f)* .5; r[0] = (r[0]+f)* .5;
			i[1] = (i[0]-g)*-.5; i[0] = (i[0]+g)*-.5;
		}
		return;
	}

	initffttables(m+1);
	n = (1<<m);

#if (WASTE_32K_TO_SPEEDUP_COMPFFT == 1)
	if (n <= (sizeof(oldr)/sizeof(oldr[0])))
	{
		memcpy(oldr,r,n<<2);
		memcpy(oldi,i,n<<2);
		if (dir > 0)
		{
			for(z=n-1;z>=0;z--) { j = bitrevreal[z]; r[z] = oldr[j]; i[z] = oldi[j]; }
		}
		else
		{
			g = twoovern; h = -g;
			for(z=n-1;z>=0;z--) { j = bitrevreal[z]; r[z] = oldr[j]*g; i[z] = oldi[j]*h; }
		}
	}
	else
#endif
	{
		if (dir > 0)
		{
			for(z=n-1;z>=0;z--)
			{
				j = bitrev[z];
				if (z < j)
				{
					f = r[z]; r[z] = r[j]; r[j] = f;
					f = i[z]; i[z] = i[j]; i[j] = f;
				}
			}
		}
		else
		{
			g = twoovern; h = -g;
			for(z=0;z<n;z++)
			{
				j = bitrev[z];
				if (z < j)
				{
					f = r[z]*g; r[z] = r[j]*g; r[j] = f;
					f = i[z]*h; i[z] = i[j]*h; i[j] = f;
				}
				else if (z == j) { r[z] *= g; i[z] *= h; }
			}
		}

			//Combined bit reversal requires extra buffer
			//Unfortunately, extra buffer is available in realfft but not compfft :/
		if (fftfunc == fftradix4odd)
		{
			long nn = (n>>1);
			for(z=nn-1;z>=0;z--)
			{
				j = ((z&0xaaaaaaaa)>>1) + ((z&0x55555555)<<1);
				if (j < z)
				{
					f = r[j]; r[j] = r[z]; r[z] = f;
					f = i[j]; i[j] = i[z]; i[z] = f;
					f = r[j+nn]; r[j+nn] = r[z+nn]; r[z+nn] = f;
					f = i[j+nn]; i[j+nn] = i[z+nn]; i[z+nn] = f;
				}
			}
		}
	}

	fftfunc(r,i,n);
}

	//Real value optimized FFT:
	//   m: log2(n)
	//   dir=1:forward, dir=-1:inverse
#pragma warning(disable:4731)
static __declspec(align(16)) float onehalf[4] = {.5f,.5f,.5f,.5f};
void realfft (float *snd, float *r, float *i, long m, long dir)
{
	float f, g, re, ro, ie, io, *sndi;
	long j, n, z, zz, bo2;

	if (((m < 4) && (dir > 0)) || ((m < 3) && (dir < 0)))
	{
		if (dir > 0)
		{
			i[0] = 0;
			switch(m)
			{
				case 0: r[0] = snd[0]; break;
				case 1: r[0] = snd[0]+snd[1]; r[1] = snd[0]-snd[1]; i[1] = 0; break;
				case 2: f = snd[0]+snd[2]; g = snd[1]+snd[3];
					r[0] = f+g; r[1] = snd[0]-snd[2]; r[2] = f-g;
					i[1] = snd[3]-snd[1]; i[2] = 0; break;
				case 3:
					re = snd[0]+snd[4]; ro = snd[2]+snd[6]; f = re+ro;
					ie = snd[1]+snd[5]; io = snd[3]+snd[7]; g = ie+io;
					r[0] = f+g; r[2] = re-ro;
					r[4] = f-g; i[2] = io-ie;
					ie = snd[7]-snd[3]; io = snd[5]-snd[1];
					re = snd[0]-snd[4]; ro = snd[2]-snd[6];
					f = (ie-io)*SQRT2O2; g = (ie+io)*SQRT2O2;
					r[1] = re+f; r[3] = re-f;
					i[1] = g-ro; i[3] = g+ro;
					i[4] = 0; break;
			}
		}
		else
		{
			switch(m)
			{
				case 0: snd[0] = r[0]; break;
				case 1: snd[0] = (r[0]+r[1])*.5f; snd[1] = (r[0]-r[1])*.5f; break;
				case 2: f = (r[0]+r[2])*.25f; re = r[1]*.5f;
						  g = (r[0]-r[2])*.25f; ro = i[1]*.5f;
					snd[0] = f+re; snd[1] = g-ro; snd[2] = f-re; snd[3] = g+ro; break;
			}
		}
		return;
	}

	initffttables(m);

	n = (1<<m); bo2 = (n>>1);
	if (dir > 0)
	{
			//Input: snd(0..n-1), Output: r(0..n\2), i(0..n\2)
#if 0
		for(z=bo2-1;z>=0;z--) { j = (bitrevreal[z]<<1); r[z] = snd[j]; i[z] = snd[j+1]; }
#else
		//for(z=0;z<bo2;z++) printf("%d ",bitrevreal[z]);

		_asm
		{
			push ebx
			push esi
			push edi
			mov ecx, bo2
			mov edx, snd
			mov esi, r
			mov edi, i
			push ebp
			mov ebp, bitrevreal
			lea ebx, [ecx-4]
			shr ecx, 2

begr0:   mov eax, [ebp+ebx*4]
			movss xmm0, [edx+eax*8]
			movss xmm4, [edx+eax*8+4]
			add eax, ecx
			movss xmm1, [edx+eax*8]
			movss xmm5, [edx+eax*8+4]
			add eax, ecx
			movss xmm2, [edx+eax*8]
			movss xmm6, [edx+eax*8+4]
			add eax, ecx
			movss xmm3, [edx+eax*8]
			movss xmm7, [edx+eax*8+4]

			unpcklps xmm0, xmm1
			unpcklps xmm4, xmm5
			unpcklps xmm2, xmm3
			unpcklps xmm6, xmm7
			movlhps xmm0, xmm2
			movlhps xmm4, xmm6
			movaps [esi+ebx*4], xmm0
			movaps [edi+ebx*4], xmm4

			sub ebx, 4
			jns short begr0
			pop ebp
			pop edi
			pop esi
			pop ebx
		}
#endif

		fftfunc(r,i,bo2);

		r[bo2] = 0; i[bo2] = 0;
#if 0
		z = 0; zz = bo2;
		do
		{
			ro = i[z]+i[zz]; io = r[zz]-r[z];
			re = r[z]+r[zz]; ie = i[z]-i[zz];
			r[z] = (re + ro*coslut[z] + io*sinlut[z])*.5f;
			i[z] = (ie + io*coslut[z] - ro*sinlut[z])*.5f;
			r[zz] = re - r[z];
			i[zz] = i[z] - ie;
			z++; zz--;
		} while (z < zz);
		i[z] = -i[z];
#else
		_asm
		{
			push ebx
			push esi
			push edi

			mov eax, coslut
			mov ebx, sinlut
			mov edx, bo2  ;         0 1 2 3 (4) 5 6 7 8
			xor ecx, ecx  ;0 1 2 3  4 5 6 7 (8) 9 10 11 12  13 14 15 16
			sub edx, 3
			mov esi, r
			mov edi, i
			movaps xmm7, onehalf
begrf:
			movups xmm4, [esi+edx*4]
			movups xmm5, [edi+edx*4]
			movaps xmm0, [esi+ecx*4]   ; rz3  rz2  rz1  rz0
			shufps xmm4, xmm4, 0x1b    ;rzz3 rzz2 rzz1 rzz0
			movaps xmm2, [edi+ecx*4]   ; iz3  iz2  iz1  iz0
			shufps xmm5, xmm5, 0x1b    ;izz3 izz2 izz1 izz0
			movaps xmm1, xmm0
			addps xmm0, xmm4           ;xmm0: re
			movaps xmm3, xmm2
			subps xmm4, xmm1           ;xmm4: io
			addps xmm2, xmm5           ;xmm2: ro
			subps xmm3, xmm5           ;xmm3: ie
			movaps xmm1, xmm4          ;xmm1: io
			movaps xmm5, xmm2          ;xmm5: ro

			mulps xmm4, [ebx+ecx*4]    ;xmm4: io*sin
			mulps xmm2, [eax+ecx*4]    ;xmm2: ro*cos
			mulps xmm1, [eax+ecx*4]    ;xmm1: io*cos
			mulps xmm5, [ebx+ecx*4]    ;xmm5: ro*sin
			addps xmm2, xmm0           ;xmm2: re + ro*cos
			addps xmm1, xmm3           ;xmm1: ie + io*cos
			addps xmm2, xmm4           ;xmm2: re + ro*cos + io*sin
			subps xmm1, xmm5           ;xmm1: ie + io*cos - ro*sin

			mulps xmm2, xmm7           ;xmm2: r[z]
			mulps xmm1, xmm7           ;xmm1: i[z]
			movaps [esi+ecx*4], xmm2
			movaps [edi+ecx*4], xmm1
			subps xmm0, xmm2           ;xmm0: r[zz] (backwards)
			subps xmm1, xmm3           ;xmm1: i[zz] (backwards)
			shufps xmm0, xmm0, 0x1b    ;xmm0: r[zz]
			shufps xmm1, xmm1, 0x1b    ;xmm1: i[zz]
			movups [esi+edx*4], xmm0
			movups [edi+edx*4], xmm1
			add ecx, 4
			sub edx, 4
			cmp ecx, edx
			jl short begrf

			pop edi
			pop esi
			pop ebx
		}
		z = (bo2>>1);
		i[z] = -i[z];
#endif
		r[0] -= i[bo2]; i[bo2] = 0;
		r[bo2] -= i[0]; i[0] = 0;
	}
	else
	{
			//Input: r(0..n\2), i(0..n\2), Output: snd(0..n-1)
		sndi = &snd[bo2];
#if 0
		z = 1; zz = bo2-1;
		do
		{
			ro = i[z]+i[zz]; io = r[zz]-r[z];
			re = r[z]+r[zz]; ie = i[zz]-i[z];
			f = ro*coslut[zz] + io*sinlut[zz];
			g = io*coslut[zz] - ro*sinlut[zz];
			snd[ z] = re-f; sndi[ z] = g+ie;
			snd[zz] = re+f; sndi[zz] = g-ie;
			z++; zz--;
		} while (z <= zz);
		sndi[0] = r[0]-r[bo2]; snd[0] = r[0]+r[bo2];
#else
		_asm
		{
			push ebx
			push esi
			push edi

			mov eax, bo2
			mov ebx, sinlut
			mov ecx, snd
			mov edx, sndi
			mov esi, r
			mov edi, i

			push ebp        ;      (0)1 2 3 44 5 6 7(8)
			mov ebp, coslut ;(0)1 2 3 4  5 6 7 88 9 10 11  12 13 14 15(16)
			sub edi, esi               ;edi: i-r
			sub ebp, esi               ;eax: coslut-r
			sub ebx, esi               ;ebx: sinlut-r
			sub ecx, esi               ;ecx: snd-r
			sub edx, esi               ;edx: sndi-r
			lea eax, [eax*4+esi-16]    ;ebp: (zz<<2)+r
			add esi, 4                 ;esi: (z<<2)+r
begri:
			movups xmm0, [edi+esi]     ;xmm0: iz (backwards)
			movaps xmm1, [eax]         ;xmm1: rzz
			movups xmm2, [esi]         ;xmm2: rz (backwards)
			movaps xmm3, [edi+eax]     ;xmm3: izz
			shufps xmm0, xmm0, 0x1b    ;xmm0: iz
			shufps xmm2, xmm2, 0x1b    ;xmm2: rz
			movaps xmm4, xmm0          ;xmm4: iz
			movaps xmm5, xmm1          ;xmm5: rzz
			addps xmm0, xmm3           ;xmm0: iz+izz (ro)
			subps xmm1, xmm2           ;xmm1: rzz-rz (io)
			addps xmm2, xmm5           ;xmm2: rz+rzz (re)
			subps xmm3, xmm4           ;xmm3: izz-iz (ie)
			movaps xmm4, xmm0          ;xmm4: iz+izz (ro)
			movaps xmm5, xmm1          ;xmm5: rzz-rz (io)

			mulps xmm0, [ebp+eax]      ;xmm0: ro*cos
			mulps xmm5, [ebx+eax]      ;xmm5: io*sin
			mulps xmm4, [ebx+eax]      ;xmm4: ro*sin
			mulps xmm1, [ebp+eax]      ;xmm1: io*cos
			addps xmm0, xmm5           ;xmm0: ro*cos + io*sin  f
			subps xmm1, xmm4           ;xmm1: io*cos - ro*sin  g

			movaps xmm4, xmm1          ;xmm4: g
			movaps xmm5, xmm2          ;xmm5: re
			subps xmm1, xmm3           ;xmm1: g-ie (sndi[zz])
			addps xmm2, xmm0           ;xmm2: re+f (snd [zz])
			subps xmm5, xmm0           ;xmm5: re-f (snd [ z])
			addps xmm4, xmm3           ;xmm4: g+ie (sndi[ z])

			movaps [edx+eax], xmm1
			movaps [ecx+eax], xmm2
			shufps xmm5, xmm5, 0x1b
			shufps xmm4, xmm4, 0x1b
			movups [ecx+esi], xmm5
			movups [edx+esi], xmm4
			add esi, 16
			sub eax, 16
			cmp esi, eax
			jl short begri

			pop ebp
			pop edi
			pop esi
			pop ebx
		}
		sndi[0] = r[0]-r[bo2]; snd[0] = r[0]+r[bo2];
#endif

#if 0
		for(z=bo2-1;z>=0;z--) { j = bitrevreal[z]; r[z] = snd[j]; i[z] = sndi[j]; }
#else
		_asm
		{
			push ebx
			push esi
			push edi
			mov eax, bo2
			sub eax, 4
			mov ebx, bitrevreal
			mov edx, snd
			mov esi, r
			mov edi, i
			push ebp
			mov ebp, sndi

begr1:   mov ecx, [ebx+eax*4]
			movss xmm0, [edx+ecx*4]
			movss xmm4, [ebp+ecx*4]
			mov ecx, [ebx+eax*4+4]
			movss xmm1, [edx+ecx*4]
			movss xmm5, [ebp+ecx*4]
			mov ecx, [ebx+eax*4+8]
			movss xmm2, [edx+ecx*4]
			movss xmm6, [ebp+ecx*4]
			mov ecx, [ebx+eax*4+12]
			movss xmm3, [edx+ecx*4]
			movss xmm7, [ebp+ecx*4]

			unpcklps xmm0, xmm1
			unpcklps xmm4, xmm5
			unpcklps xmm2, xmm3
			unpcklps xmm6, xmm7
			movlhps xmm0, xmm2
			movlhps xmm4, xmm6
			movaps [esi+eax*4], xmm0
			movaps [edi+eax*4], xmm4

			sub eax, 4
			jns short begr1
			pop ebp
			pop edi
			pop esi
			pop ebx
		}
#endif

		fftfunc(r,i,bo2);

#if 0
		for(z=0;z<bo2;z++) { snd[z+z] = r[z]*oneovern; snd[z+z+1] = i[z]*oneovern; }
#else
		_asm
		{
			push edi
			mov eax, n
			mov ecx, r
			add eax, eax
			mov edi, i
			mov edx, snd
			add ecx, eax
			add edi, eax
			lea edx, [edx+eax*2]
			neg eax
			movss xmm7, oneovern
			shufps xmm7, xmm7, 0
begreal0:movaps xmm0, [ecx+eax] ;xmm0: r3 r2 r1 r0
			mulps xmm0, xmm7       ;xmm0: r3f r2f r1f r0f
			movaps xmm1, [edi+eax] ;xmm1: i3 i2 i1 i0
			mulps xmm1, xmm7       ;xmm1: i3f i2f i1f i0f
			movaps xmm2, xmm0      ;xmm2: r3f r2f r1f r0f
			unpcklps xmm0, xmm1    ;xmm0: i1f r1f i0f r0f
			movaps [edx+eax*2], xmm0
			unpckhps xmm2, xmm1    ;xmm2: i3f r3f i2f r2f
			movaps [edx+eax*2+16], xmm2
			add eax, 16
			jnc short begreal0
			pop edi
		}
#endif
	}
}

static long curdctblksiz = 0;
static float *dctlut = 0, *dctflut, *dctilut, *tempr, *tempi;

static void initdcttables (long n)
{
	float c, s;
	long i, no2;

	if (n == curdctblksiz) return;
	curdctblksiz = n;

	n = (1<<n); no2 = (n>>1);
	if (dctlut) { alignfree(dctlut); dctlut = 0; }
	i = (((n>>1)+3)&~3)*3 + (((n>>1)+1+3)&~3)*2;
	dctlut = (float *)alignmalloc(i*sizeof(float),16);
	if (!dctlut) { curdctblksiz = 0; return; }
#if (FFTTEST == 1)
	if (lb0 == lb1) printf("   dctlut: %d bytes\n",i*sizeof(float));
#endif
	dctflut = &dctlut[((n>>1)+3)&~3];
	dctilut = &dctflut[((n>>1)+3)&~3];
	tempr = &dctilut[((n>>1)+3)&~3];
	tempi = &tempr[((n>>1)+1+3)&~3];

	fcalcsinasm(dctflut,n,0.25/(float)n); //This is an ugly memory hack... but it IS safe!
	for(i=0;i<no2;i++) dctlut[i] = dctflut[i+i];     //WARNING: Can't fuse loop!
	for(i=0;i<no2;i++) dctflut[i] = dctflut[i+i+1];  //WARNING: Can't fuse loop!
	for(i=0;i<no2;i++) dctilut[i] = .25f/dctflut[i]; //WARNING: Can't fuse loop!
}

	//   //slowdct is equivalent to fastdct, but MUCH SLOWER!
	//void slowdct (float *x, float *y, long logn, long dir)
	//{  long j, k, n; float f, g, s; n = (1<<logn); f = PI / (float)n;
	//   if (dir > 0)
	//      { for(k=0;k<n;k++) { s = 0; for(j=0;j<n;j++) s += cos(((float)j+.5)*(float)k*f)*x[j]; y[k] = s; } }
	//   else
	//   {  g = 2.f / (float)n; y[0] *= .5;
	//      for(j=0;j<n;j++) { s = 0; for(k=0;k<n;k++) s += cos(((float)j+.5)*(float)k*f)*y[k]; x[j] = s*g; }
	//      y[0] *= 2.0; //Hack so it doesn't modify input
	//}  }
void fastdct (float *snd, float *y, long logn, long dir)
{
	long i, j, m, n, no2;
	float a, b, c, d;

	if (logn <= 2)
	{
		if (dir > 0)
		{
			switch (logn)
			{
				case 0: y[0] = snd[0]; break;
				case 1: y[0] = snd[0]+snd[1]; y[1] = (snd[0]-snd[1])*SQRT2O2; break;
				case 2: a = snd[0]+snd[3]; b = snd[1]+snd[2]; y[0] = a+b; y[2] = (a-b)*SQRT2O2;
						  a = snd[0]-snd[3]; b = snd[1]-snd[2];
					y[1] = a*COSPIO8+b*COS3PIO8; y[3] = a*COS3PIO8-b*COSPIO8; break;
			}
		}
		else
		{
			switch (logn)
			{
				case 0: snd[0] = y[0]; break;
				case 1: a = y[0]*.5; b = y[1]*SQRT2O2; snd[0] = a+b; snd[1] = a-b; break;
				case 2: a = y[0]*.25; b = y[2]*(SQRT2O2*.5);
					c = a+b; d = y[1]*(COSPIO8*.5) + y[3]*(COS3PIO8*.5); snd[0] = c+d; snd[3] = c-d;
					c = a-b; d = y[1]*(COS3PIO8*.5) - y[3]*(COSPIO8*.5); snd[1] = c+d; snd[2] = c-d; break;
			}
		}
		return;
	}

	n = (1<<logn); m = n-1; no2 = (n>>1);
	initdcttables(logn);

	if (dir > 0)
	{
#if 0
		for(i=no2-1;i>=0;i--)
		{
			a = (snd[i]+snd[m-i])*.5f;
			b = (snd[i]-snd[m-i])*dctflut[i];
			y[i] = a+b; y[m-i] = a-b;
		}
#else
		_asm
		{
			push esi
			push edi
			mov esi, snd
			mov edi, y
			mov edx, dctflut
			mov eax, no2
			mov ecx, eax
			sub eax, 4
			movaps xmm7, onehalf
begd0:   movaps xmm2, [esi+ecx*4]
			movaps xmm0, [esi+eax*4]
			shufps xmm2, xmm2, 0x1b
			movaps xmm1, xmm0
			addps xmm0, xmm2
			subps xmm1, xmm2
			mulps xmm0, xmm7
			mulps xmm1, [edx+eax*4]
			movaps xmm2, xmm0
			addps xmm0, xmm1
			subps xmm2, xmm1
			movaps [edi+eax*4], xmm0
			shufps xmm2, xmm2, 0x1b
			movaps [edi+ecx*4], xmm2
			add ecx, 4
			sub eax, 4
			jnc short begd0
			pop edi
			pop esi
		}
#endif
		realfft(y,tempr,tempi,logn,1);
		y[n-1] = tempr[no2]*.5f;
		for(i=no2-1;i;i--)
		{
			j = i+i;
			y[j  ] = dctlut[no2-i]*tempr[i] + dctlut[    i]*tempi[i];
			y[j-1] = dctlut[    i]*tempr[i] - dctlut[no2-i]*tempi[i] + y[j+1];
		}
		y[0] = tempr[0];
	}
	else
	{
		tempr[no2] = y[n-1]*2; tempi[no2] = 0;
		for(i=no2-1;i;i--)
		{
			j = i+i; a = y[j-1]-y[j+1];
			tempr[i] = dctlut[no2-i]*y[j] + dctlut[    i]*a;
			tempi[i] = dctlut[    i]*y[j] - dctlut[no2-i]*a;
		}
		tempr[0] = y[0]; tempi[0] = 0;
		realfft(snd,tempr,tempi,logn,-1);
#if 0
		for(i=no2-1;i>=0;i--)
		{
			a = (snd[i]+snd[m-i])*.5f;
			b = (snd[i]-snd[m-i])*dctilut[i];
			snd[i] = a+b; snd[m-i] = a-b;
		}
#else
		_asm
		{
			push esi
			mov esi, snd
			mov edx, dctilut
			mov eax, no2
			mov ecx, eax
			sub eax, 4
			movaps xmm7, onehalf
begd1:   movaps xmm2, [esi+ecx*4]
			movaps xmm0, [esi+eax*4]
			shufps xmm2, xmm2, 0x1b
			movaps xmm1, xmm0
			addps xmm0, xmm2
			subps xmm1, xmm2
			mulps xmm0, xmm7
			mulps xmm1, [edx+eax*4]
			movaps xmm2, xmm0
			addps xmm0, xmm1
			subps xmm2, xmm1
			movaps [esi+eax*4], xmm0
			shufps xmm2, xmm2, 0x1b
			movaps [esi+ecx*4], xmm2
			add ecx, 4
			sub eax, 4
			jnc short begd1
			pop esi
		}
#endif
	}
}

#if 0
float getentropybits (long *frqbuf, long num)
{
	double d;
	long z, t;

	d = 0.0; t = 0;
	for(z=0;z<num;z++)
	{
		if (!frqbuf[z]) continue;
		d -= log((double)frqbuf[z])*(double)frqbuf[z];
		t += frqbuf[z];
	}
	return((float)(log((double)t)*(double)t + d) / log(2));
}
#endif



#ifdef FFTTEST

	//                 FFT.C         (pure radix 2)
	//    radix     *         +         *         +
	// 2: 2        14        25        14        25 :[
	// 3: 4        36        67        44        71   :?
	// 4: 8        88       175       120       187     :)
	// 5: 4       240       435       304       467   :?
	// 6: 2       736      1123       736      1123 :[
	// 7: 8      1216      2435      1728      2627   :?)
	// 8: 2      3968      6019      3968      6019 :[
	// 9: 4      6912     12547      8960     13571   :?
	//10: 8     13824     27907     19968     30211     :)
	//11: 4     33792     61443     44032     66563   :?
	//12: 2     96256    145411     96256    145411 :[
	//13: 8    143360    290819    208896    315395   :?)
	//14: 2    450560    679939    450560    679939 :[
	//15: 4    737280   1343491    966656   1458179   :?
	//16: 8   1409024   2867203   2064384   3112963     :)
	//17: 4   3342336   6094851   4390912   6619139   :?
	//18: 2   9306112  14024707   9306112  14024707 :[
	//19: 8  13369344  27262969  19660800  29622275   :?)
	//20: 2  41418752  62390275  41418752  62390275 :[
	//21: 4  66060288 120586243  87031808 131072003   :?
	//22: 8 123731968 252706819 182452224 274726915     :)
	//23: 4 289406976 528482307 381681664 574619651   :?

	//10/17/2002 P4-2784.31Mhz:
	//             compfft:               realfft:
	//         before:  after:        before:  after:
	//     64     2384    2852  :[       2724    2412  :]
	//     64     2604    2612           3316    2756  :]
	//    128    12400   10708  :]       3288    3432  :[
	//    128    12872   10232  :]       4260    4284
	//    256    13228   15196  :[      13100   11100  :]
	//    256    13696   14044  :[      15856   13336  :]
	//    512    64032   50272  :]      13956   14848  :[
	//    512    64644   47888  :]      20392   20040
	//   1024    63884   70076  :[      61952   48528  :]
	//   1024    65000   66552  :[      75804   59728  :]
	//   2048   341280  229292  :]]     66500   68700  :[
	//   2048   342048  222236  :]]     90584   89168
	//   4096   311820  334280  :[     338216  223816  :]]
	//   4096   312884  319380  :[     388956  265316  :]]
	//   8192  1660424 1044900  :]]    304920  320404  :[
	//   8192  1666316 1010808  :]]    413280  417152  :[
	//  16384  2064008 2118988  :[    1552548  984468  :]]
	//  16384  2097236 2096188        1762848 1170672  :]]
	//  32768 22777504 6248140  :]]   2230584 2313368  :[
	//  32768 22931144 6169900  :]]   2575584 2572780
	//  65536   8:15461472,13198116          14640100/  10000552  forward is invalid sample!
	// 131072                                36232604/  20121748  forward is invalid sample!
	// 262144                               130957352/ 102874480  forward is invalid sample!
	// 524288   8:96320980,599301720        255234312/ 225051028  forward is invalid sample!
	//4194304   8:6690855980,6701608320    2585632084/2365056768  forward is invalid sample!

static __forceinline __int64 rdtsc64 ()
{
	__int64 q;
	_asm
	{
		rdtsc
		mov dword ptr q[0], eax
		mov dword ptr q[4], edx
	}
	return(q);
}

#define isbadfloat(f) (_fpclass(f)&(1+2+4+512))

void main (int argc, char **argv)
{
	__int64 q0, q1, bestq;
	float *snd, *real, *imag, *oreal, *oimag, *snd2;
	long i, lblocksiz, blocksiz, compbytes, testmode;

	if ((argc < 2) || (argc > 3))
	{
		puts("FFT [Function to test] (log2 of block size)                   by Ken Silverman");
		puts("   Function to test: 0=realfft, 1=compfft, 2=fastdct");
		puts("   Leave 2nd param blank to test speed&functionality of block sizes 6-15");
		exit(0);
	}

	switch(atol(argv[1]))
	{
		case 0: testmode = 1; break;
		case 1: testmode = 0; break;
		case 2: testmode = 2; break;
		default: testmode = 1; break;
	}
	if (argc == 3) lb0 = lb1 = min(max(atol(argv[2]),0),23); else { lb0 = 6; lb1 = 15; }

	for(lblocksiz=lb0;lblocksiz<=lb1;lblocksiz++)
	{
		blocksiz = (1<<lblocksiz);

		snd = (float *)alignmalloc(blocksiz*sizeof(float),16); if (!snd) { puts("malloc 1 failed"); exit(0); }
		switch(testmode)
		{
			case 0: case 2: compbytes = blocksiz*sizeof(float); break;
			case 1: compbytes = ((blocksiz>>1)+1)*sizeof(float); break;
		}
		real = (float *)alignmalloc(compbytes,16); if (!real) { puts("malloc 2 failed"); exit(0); }
		imag = (float *)alignmalloc(compbytes,16); if (!imag) { puts("malloc 3 failed"); exit(0); }
		oreal = (float *)alignmalloc(compbytes,16); if (!oreal) { puts("malloc 4 failed"); exit(0); }
		oimag = (float *)alignmalloc(compbytes,16); if (!oimag) { puts("malloc 5 failed"); exit(0); }
		snd2 = (float *)alignmalloc(blocksiz*sizeof(float),16); if (!snd2) { puts("malloc 6 failed"); exit(0); }
		if (lb0 == lb1)
		{
			printf("      snd: %d bytes\n",blocksiz*sizeof(float));
			printf("     real: %d bytes\n",compbytes);
			printf("     imag: %d bytes\n",compbytes);
			printf("    oreal: %d bytes\n",compbytes);
			printf("    oimag: %d bytes\n",compbytes);
			printf("     snd2: %d bytes\n\n",blocksiz*sizeof(float));
		}

			//Generate random time signal data
		for(i=0;i<blocksiz;i++)
			snd[i] = (float)i;
			//snd[i] = (float)(rand()-16384);
			//{ if (i == 1) snd[i] = 1; else snd[i] = 0; }  //(float)(rand()-16384);

		if (lblocksiz < 6)
		{
				//Display time data
			for(i=0;i<blocksiz;i++) printf("time[%d] = %f\n",i,snd[i]);
			puts("");
		}

		bestq = 0x7fffffffffffffff;
		switch (testmode)
		{
			case 0:
				initffttables(lblocksiz+1); //Hack to make table calculation not count in timings
				memcpy(real,snd,blocksiz*sizeof(float)); //remember to allocate real&imag at full size!
				memset(imag,0,blocksiz*sizeof(float));
				if (lb0 == lb1) { q0 = rdtsc64(); compfft(real,imag,lblocksiz,1); q1 = rdtsc64(); bestq = q1-q0; }
				else
				{
					memcpy(oreal,real,compbytes); memcpy(oimag,imag,compbytes);
					for(i=0;i<32;i++)
					{
						memcpy(real,oreal,compbytes); memcpy(imag,oimag,compbytes);
						q0 = rdtsc64(); compfft(real,imag,lblocksiz,1); q1 = rdtsc64();
						if (q1-q0 < bestq) bestq = q1-q0;
					}
				}
				printf("forcompfft 2^%-2d ",lblocksiz);
				break;
			case 1:
				initffttables(lblocksiz); //Hack to make table calculation not count in timings
				if (lb0 == lb1)
				{
					q0 = rdtsc64(); realfft(snd,real,imag,lblocksiz,1); q1 = rdtsc64(); bestq = q1-q0;
				}
				else
				{
					for(i=0;i<32;i++)
					{
						q0 = rdtsc64(); realfft(snd,real,imag,lblocksiz,1); q1 = rdtsc64();
						if (q1-q0 < bestq) bestq = q1-q0;
					}
				}
				printf("forrealfft 2^%-2d ",lblocksiz);
				break;
			case 2:
				initffttables(lblocksiz); //Hack to make table calculation not count in timings
				initdcttables(lblocksiz); //Hack to make table calculation not count in timings
				if (lb0 == lb1) { q0 = rdtsc64(); fastdct(snd,real,lblocksiz,1); q1 = rdtsc64(); bestq = q1-q0; }
				else
				{
					for(i=0;i<32;i++)
					{
						q0 = rdtsc64(); fastdct(snd,real,lblocksiz,1); q1 = rdtsc64();
						if (q1-q0 < bestq) bestq = q1-q0;
					}
				}
				printf("forfastdct 2^%-2d ",lblocksiz);
				break;
		}
		if (fftfunc == fftradix4) printf("radix4");
									else printf("radixM");
		printf(":%12I64dcc",bestq);
		if (lb0 == lb1) puts(""); else printf(" | ");

		if (lblocksiz < 6)
		{
			switch(testmode)
			{
				case 0: case 1:
						//Display frequency data
					for(i=0;i<=(blocksiz>>1);i++) printf("frq[%d]: %f + %fj\n",i,real[i],imag[i]);
					break;
				case 2:
					for(i=0;i<blocksiz;i++) printf("frq[%d]: %f\n",i,real[i]);
					break;
			}
			puts("");
		}

		bestq = 0x7fffffffffffffff;
		switch(testmode)
		{
			case 0:
				initffttables(lblocksiz+1); //Hack to make table calculation not count in timings
				for(i=1;i<(blocksiz>>1);i++) { real[blocksiz-i] = real[i]; imag[blocksiz-i] = -imag[i]; }
				if (lb0 == lb1) { q0 = rdtsc64(); compfft(real,imag,lblocksiz,-1); q1 = rdtsc64(); bestq = q1-q0; }
				else
				{
					memcpy(oreal,real,compbytes); memcpy(oimag,imag,compbytes);
					for(i=0;i<32;i++)
					{
						memcpy(real,oreal,compbytes); memcpy(imag,oimag,compbytes);
						q0 = rdtsc64(); compfft(real,imag,lblocksiz,-1); q1 = rdtsc64();
						if (q1-q0 < bestq) bestq = q1-q0;
					}
				}
				memcpy(snd2,real,blocksiz*sizeof(float)); //remember to allocate real&imag at full size!
				printf("invcompfft 2^%-2d ",lblocksiz);
				break;
			case 1:
				initffttables(lblocksiz); //Hack to make table calculation not count in timings
				if (lb0 == lb1) { q0 = rdtsc64(); realfft(snd2,real,imag,lblocksiz,-1); q1 = rdtsc64(); bestq = q1-q0; }
				else
				{
					memcpy(oreal,real,compbytes); memcpy(oimag,imag,compbytes);
					for(i=0;i<32;i++)
					{
						memcpy(real,oreal,compbytes); memcpy(imag,oimag,compbytes);
						q0 = rdtsc64(); realfft(snd2,real,imag,lblocksiz,-1); q1 = rdtsc64();
						if (q1-q0 < bestq) bestq = q1-q0;
					}
				}
				printf("invrealfft 2^%-2d ",lblocksiz);
				break;
			case 2:
				initffttables(lblocksiz); //Hack to make table calculation not count in timings
				initdcttables(lblocksiz); //Hack to make table calculation not count in timings
				if (lb0 == lb1) { q0 = rdtsc64(); fastdct(snd2,real,lblocksiz,-1); q1 = rdtsc64(); bestq = q1-q0; }
				else
				{
					for(i=0;i<32;i++)
					{
						q0 = rdtsc64(); fastdct(snd2,real,lblocksiz,-1); q1 = rdtsc64();
						if (q1-q0 < bestq) bestq = q1-q0;
					}
				}
				printf("invfastdct 2^%-2d ",lblocksiz);
				break;
		}
		if (fftfunc == fftradix4) printf("radix4");
									else printf("radixM");
		printf(":%12I64dcc\n",bestq);

		if (lblocksiz < 6)
		{
				//Display time data (should be same as original time data)
			for(i=0;i<blocksiz;i++) printf("time[%d] = %f\n",i,snd2[i]);
		}

		for(i=0;i<blocksiz;i++)
			if (((fabs(snd2[i]-snd[i]) > fabs(snd[i]*.01)) && (fabs(snd[i]) > .01)) || (isbadfloat(snd2[i])))
				{ printf("ERROR at %d/%d %f %f\n",i,blocksiz,snd[i],snd2[i]); break; }
		if ((lb0 == lb1) && (i >= blocksiz)) puts("Data good.");

		alignfree(snd2); alignfree(oimag); alignfree(oreal); alignfree(imag); alignfree(real); alignfree(snd);
		if (dctlut) { alignfree(dctlut); dctlut = 0; }
		if (bitrev2) { alignfree(bitrev2); bitrev2 = 0; }
		if (radix4tab) { alignfree(radix4tab); radix4tab = 0; }
		if (bitrev) { alignfree(bitrev); bitrev = 0; }
		if (sinlut) { alignfree(sinlut); sinlut = 0; }
	}
}

#endif

#if 0
!endif
#endif
