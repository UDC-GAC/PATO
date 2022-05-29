#ifndef _TRIPLEX_ALPHABET_HPP_
#define _TRIPLEX_ALPHABET_HPP_

#include "seqan.hpp"

namespace SEQAN_NAMESPACE_MAIN
{

template <typename T = void>
struct TranslateTableTriplex2Ascii_
{
	static char const VALUE[9];
};
template <typename T>
char const TranslateTableTriplex2Ascii_<T>::VALUE[9] = {'A', 'C', 'G', 'T', 'R', 'Y', 'K', 'M', 'N'};


//____________________________________________________________________________
template <typename T = void>
struct TranslateTableDna2Triplex_
{
	static char const VALUE[4];
};
template <typename T>
char const TranslateTableDna2Triplex_<T>::VALUE[4] =
{
	'0', //'A'
	'1', //'C'
	'2', //'G'
	'3'  //'T'
};

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableDna52Triplex_
{
	static char const VALUE[5];
};
template <typename T>
char const TranslateTableDna52Triplex_<T>::VALUE[5] =
{
	'0', //'A'
	'1', //'C'
	'2', //'G'
	'3', //'T'
	'8'  //'N'
};

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableIupac2Triplex_
{
	static char const VALUE[16];
};
template <typename T>
char const TranslateTableIupac2Triplex_<T>::VALUE[16] =
{
	3, //'U'
	3, //'T'
	0, //'A'
	8, //'W' = TA
	1, //'C'
	5, //'Y' = TC
	7, //'M' = AC
	8, //'H' = not-G
	2, //'G'
	6, //'K' = TG
	4, //'R' = AG
	8, //'D' = not-C
	8, //'S' = CG
	8, //'B' = non-A
	8, //'V' = non-T
	8  //'N' = any
};

//____________________________________________________________________________


template <typename T = void>
struct TranslateTableAscii2Triplex_
{
	static char const VALUE[256];
};
template <typename T>
char const TranslateTableAscii2Triplex_<T>::VALUE[256] =
{
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //0
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //1
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //2
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //3

	8,   0,   8,   1,   8,   8,   8,   2,   8,   8,   8,   6,   8,   7,   8,   8, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	8,   8,   4,   8,   3,   3,   8,   8,   8,   5,   8,   8,   8,   8,   8,   8, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,

	8,   0,   8,   1,   8,   8,   8,   2,   8,   8,   8,   6,   8,   7,   8,   8, //6
//   ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	8,   8,   4,   8,   3,   3,   8,   8,   8,   5,   8,   8,   8,   8,   8,   8, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,

	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //8
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //9
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //10
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //11
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //12
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //13
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //14
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8  //15
};

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableByte2Triplex_
{
	static char const VALUE[256];
};
template <typename T>
	char const TranslateTableByte2Triplex_<T>::VALUE[256] = {
	0,   1,   8,   3,   4,   5,   6,   7,   8,   8,   8,   8,   8,   8,   8,   8, //0
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //1
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //2
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //3
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //4
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //5
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //6
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //7
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //8
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //9
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //10
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //11
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //12
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //13
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //14
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8  //15
};

//____________________________________________________________________________


//____________________________________________________________________________

template <typename T = void>
struct TranslateTableDnaRY2Ascii_
{
	static char const VALUE[3];
};
template <typename T>
	char const TranslateTableDnaRY2Ascii_<T>::VALUE[3] = {'R', 'Y', 'N'};

//____________________________________________________________________________
template <typename T = void>
struct TranslateTableDna2DnaRY_
{
	static char const VALUE[4];
};
template <typename T>
	char const TranslateTableDna2DnaRY_<T>::VALUE[4] =
{
	'0', //'A'
	'1', //'C'
	'0', //'G'
	'1'  //'T'
};

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableDna52DnaRY_
{
	static char const VALUE[5];
};
template <typename T>
	char const TranslateTableDna52DnaRY_<T>::VALUE[5] =
{
	'0', //'A'
	'1', //'C'
	'0', //'G'
	'1', //'T'
	'2'  //'N'
};

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableIupac2DnaRY_
{
	static char const VALUE[16];
};
template <typename T>
	char const TranslateTableIupac2DnaRY_<T>::VALUE[16] =
{
	1, //'U'
	1, //'T'
	0, //'A'
	2, //'W' = TA
	1, //'C'
	1, //'Y' = TC
	2, //'M' = AC
	2, //'H' = not-G
	0, //'G'
	2, //'K' = TG
	0, //'R' = AG
	2, //'D' = not-C
	2, //'S' = CG
	2, //'B' = non-A
	2, //'V' = non-T
	2  //'N' = any
};

//____________________________________________________________________________


template <typename T = void>
struct TranslateTableAscii2DnaRY_
{
	static char const VALUE[256];
};
template <typename T>
	char const TranslateTableAscii2DnaRY_<T>::VALUE[256] =
{
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //0
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //1
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //2
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //3

	2,   0,   2,   1,   2,   2,   2,   0,   2,   2,   2,   2,   2,   2,   2,   2, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	2,   2,   0,   2,   1,   1,   2,   2,   2,   1,   2,   2,   2,   2,   2,   2, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,

	2,   0,   2,   1,   2,   2,   2,   0,   2,   2,   2,   2,   2,   2,   2,   2, //6
//   ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	2,   2,   0,   2,   1,   1,   2,   2,   2,   1,   2,   2,   2,   2,   2,   2, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,

	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //8
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //9
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //10
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //11
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //12
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //13
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //14
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2  //15
};

//____________________________________________________________________________


template <typename T = void>
struct TranslateTableByte2DnaRY_
{
	static char const VALUE[256];
};
template <typename T>
	char const TranslateTableByte2DnaRY_<T>::VALUE[256] = {
	0,   1,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //0
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //1
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //2
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //3
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //4
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //5
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //6
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //7
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //8
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //9
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //10
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //11
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //12
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //13
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //14
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   4  //15
};

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableDnaKM2Ascii_
{
	static char const VALUE[3];
};
template <typename T>
	char const TranslateTableDnaKM2Ascii_<T>::VALUE[3] = {'K', 'M', 'N'};

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableDna2DnaKM_
{
	static char const VALUE[4];
};
template <typename T>
	char const TranslateTableDna2DnaKM_<T>::VALUE[4] =
{
	'1', //'A'
	'1', //'C'
	'0', //'G'
	'0'  //'T'
};

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableDna52DnaKM_
{
	static char const VALUE[5];
};
template <typename T>
	char const TranslateTableDna52DnaKM_<T>::VALUE[5] =
{
	'1', //'A'
	'1', //'C'
	'0', //'G'
	'0', //'T'
	'2'  //'N'
};

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableIupac2DnaKM_
{
	static char const VALUE[16];
};
template <typename T>
	char const TranslateTableIupac2DnaKM_<T>::VALUE[16] =
{
	0, //'U'
	0, //'T'
	1, //'A'
	2, //'W' = TA
	1, //'C'
	2, //'Y' = TC
	1, //'M' = AC
	2, //'H' = not-G
	0, //'G'
	0, //'K' = TG
	2, //'R' = AG
	2, //'D' = not-C
	2, //'S' = CG
	2, //'B' = non-A
	2, //'V' = non-T
	2  //'N' = any
};

//____________________________________________________________________________


template <typename T = void>
struct TranslateTableAscii2DnaKM_
{
	static char const VALUE[256];
};
template <typename T>
	char const TranslateTableAscii2DnaKM_<T>::VALUE[256] =
{
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //0
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //1
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //2
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //3

	2,   1,   2,   1,   2,   2,   2,   0,   2,   2,   2,   0,   2,   1,   2,   2, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	2,   2,   2,   2,   0,   0,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,

	2,   1,   2,   1,   2,   2,   2,   0,   2,   2,   2,   0,   2,   1,   2,   2, //6
//   ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	2,   2,   2,   2,   0,   0,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,

	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //8
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //9
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //10
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //11
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //12
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //13
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //14
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2  //15
};

//____________________________________________________________________________


template <typename T = void>
struct TranslateTableByte2DnaKM_
{
	static char const VALUE[256];
};
template <typename T>
	char const TranslateTableByte2DnaKM_<T>::VALUE[256] = {
	0,   1,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //0
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //1
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //2
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //3
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //4
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //5
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //6
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //7
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //8
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //9
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //10
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //11
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //12
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //13
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2, //14
	2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   4  //15
};


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Triplex:
..cat:Alphabets
..summary:Alphabet for DNA (shortened version of Iupac format) with purines (R), pyrimidines (Y) and GT (K), AC (M)  and  'N' character.
..general:Class.SimpleType
..signature:DnaRY
..remarks:
...text:The @Metafunction.ValueSize@ of $DnaRY$ is 9.
The nucleotides are enumerated this way: $'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3, 'R' = 4, 'Y' = 5, 'K' = 6, 'M' = 7$.
The 'N' character ("unkown nucleotide") is encoded by 8.
...text:Objects of type $Triplex$ can be converted from various other types.
An object that has a value not in ${'A','C','G','T','R','Y','K','M'}$ is converted to $'N'$.
...text:$Triplex$ is typedef for $SimpleType<char,Triplex_>$, while $Triplex_$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
*/
struct Triplex_ {};
typedef SimpleType<unsigned char, Triplex_> Triplex;

template <> struct ValueSize< Triplex > { 
	typedef __uint8 Type;
	static const Type VALUE = 9;
};
template <> struct BitsPerValue< Triplex > { 
	typedef __uint8 Type;
    static const Type VALUE = 4;
};

//____________________________________________________________________________

/**
.Spec.DnaRY:
..cat:Alphabets
..summary:Alphabet for DNA grouped in purines (R) and pyrimidines (Y) including 'N' character.
..general:Class.SimpleType
..signature:DnaRY
..remarks:
...text:The @Metafunction.ValueSize@ of $DnaRY$ is 3.
The nucleotides are enumerated this way: $'R' = 0, 'Y' = 1$.
The 'N' character ("unkown nucleotide") is encoded by 2.
...text:Objects of type $DnaRY$ can be converted from various other types.
An object that has a value not in ${'R', 'Y'}$ is converted to $'N'$.
...text:$DnaRY$ is typedef for $SimpleType<char,DnaRY_>$, while $DnaRY_$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
*/
struct DnaRY_ {};
typedef SimpleType<unsigned char, DnaRY_> DnaRY;

template <> struct ValueSize< DnaRY > { 
	typedef __uint8 Type;
	static const Type VALUE = 3;
};
template <> struct BitsPerValue< DnaRY > { 
	typedef __uint8 Type;
    static const Type VALUE = 2;
};

//____________________________________________________________________________


/**
.Spec.DnaKM:
..cat:Alphabets
..summary:Alphabet for DNA grouped in G and T (K) and A and C (M) including 'N' character.
..general:Class.SimpleType
..signature:DnaKM
..remarks:
...text:The @Metafunction.ValueSize@ of $DnaKM$ is 3.
The nucleotides are enumerated this way: $'K' = 0, 'M' = 1$.
The 'N' character ("unkown nucleotide") is encoded by 2.
...text:Objects of type $DnaKM$ can be converted from various other types.
An object that has a value not in ${'K', 'M'}$ is converted to $'N'$.
...text:$DnaKM$ is typedef for $SimpleType<char,DnaKM_>$, while $DnaKM_$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
*/
struct DnaKM_ {};
typedef SimpleType<unsigned char, DnaKM_> DnaKM;

template <> struct ValueSize< DnaKM > { 
	typedef __uint8 Type;
	static const Type VALUE = 3;
};
template <> struct BitsPerValue< DnaKM > { 
	typedef __uint8 Type;
    static const Type VALUE = 2;
};

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign(Ascii & c_target,
				   Triplex const & source)
{
	c_target = TranslateTableTriplex2Ascii_<>::VALUE[source.value];
}
//____________________________________________________________________________
inline void assign(Ascii & c_target,
				   DnaRY const & source)
{
	c_target = TranslateTableDnaRY2Ascii_<>::VALUE[source.value];
}
//____________________________________________________________________________

inline void assign(Ascii & c_target,
				   DnaKM const & source)
{
	c_target = TranslateTableDnaKM2Ascii_<>::VALUE[source.value];
}

//////////////////////////////////////////////////////////////////////////////
//Triplex (3 letters)

template <>
struct CompareType<Triplex, __uint8> { typedef Triplex Type; };
inline void assign(Triplex & target, __uint8 c_source)
{
	target.value = TranslateTableByte2Triplex_<>::VALUE[c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Triplex, Ascii> { typedef Triplex Type; };
inline void assign(Triplex & target, Ascii c_source)
{
	target.value = TranslateTableAscii2Triplex_<>::VALUE[(unsigned char)c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Triplex, Unicode> { typedef Triplex Type; };
inline void assign(Triplex & target, Unicode c_source)
{
	target.value = TranslateTableAscii2Triplex_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Triplex, Dna> { typedef Triplex Type; };
inline void assign(Triplex & target, Dna const & c_source)
{
	target.value = TranslateTableDna2Triplex_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Triplex, Dna5> { typedef Triplex Type; };
inline void assign(Triplex & target, Dna5 const & c_source)
{
	target.value = TranslateTableDna52Triplex_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Triplex, Iupac> { typedef Triplex Type; };
inline void assign(Triplex & target, Iupac const & source)
{
	target.value = TranslateTableIupac2Triplex_<>::VALUE[source.value];
}

//////////////////////////////////////////////////////////////////////////////
//DnaRY (3 letters)

template <>
struct CompareType<DnaRY, __uint8> { typedef DnaRY Type; };
inline void assign(DnaRY & target, __uint8 c_source)
{
	target.value = TranslateTableByte2DnaRY_<>::VALUE[c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<DnaRY, Ascii> { typedef DnaRY Type; };
inline void assign(DnaRY & target, Ascii c_source)
{
	target.value = TranslateTableAscii2DnaRY_<>::VALUE[(unsigned char)c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<DnaRY, Unicode> { typedef DnaRY Type; };
inline void assign(DnaRY & target, Unicode c_source)
{
	target.value = TranslateTableAscii2DnaRY_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<DnaRY, Dna> { typedef DnaRY Type; };
inline void assign(DnaRY & target, Dna const & c_source)
{
	target.value = TranslateTableDna2DnaRY_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<DnaRY, Dna5> { typedef DnaRY Type; };
inline void assign(DnaRY & target, Dna5 const & c_source)
{
	target.value = TranslateTableDna52DnaRY_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<DnaRY, Iupac> { typedef DnaRY Type; };
inline void assign(DnaRY & target, Iupac const & source)
{
	target.value = TranslateTableIupac2DnaRY_<>::VALUE[source.value];
}

//////////////////////////////////////////////////////////////////////////////
//DnaKM (3 letters)

template <>
struct CompareType<DnaKM, __uint8> { typedef DnaKM Type; };
inline void assign(DnaKM & target, __uint8 c_source)
{
	target.value = TranslateTableByte2DnaKM_<>::VALUE[c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<DnaKM, Ascii> { typedef DnaKM Type; };
inline void assign(DnaKM & target, Ascii c_source)
{
	target.value = TranslateTableAscii2DnaKM_<>::VALUE[(unsigned char)c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<DnaKM, Unicode> { typedef DnaKM Type; };
inline void assign(DnaKM & target, Unicode c_source)
{
	target.value = TranslateTableAscii2DnaKM_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<DnaKM, Dna> { typedef DnaKM Type; };
inline void assign(DnaKM & target, Dna const & c_source)
{
	target.value = TranslateTableDna2DnaKM_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<DnaKM, Dna5> { typedef DnaKM Type; };
inline void assign(DnaKM & target, Dna5 const & c_source)
{
	target.value = TranslateTableDna52DnaKM_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<DnaKM, Iupac> { typedef DnaKM Type; };
inline void assign(DnaKM & target, Iupac const & source)
{
	target.value = TranslateTableIupac2DnaKM_<>::VALUE[source.value];
}

//////////////////////////////////////////////////////////////////////////////
//Complement Functor

template <typename T = void>
struct TranslateTableTriplex2TriplexComplement_
{
	static char const VALUE[9];
};
template <typename T>
	char const TranslateTableTriplex2TriplexComplement_<T>::VALUE[9] = {'T', 'G', 'C', 'A', 'Y', 'R', 'M', 'K'};

//____________________________________________________________________________
template <typename T = void>
struct TranslateTableDnaRY2DnaRYComplement_
{
	static char const VALUE[3];
};
template <typename T>
	char const TranslateTableDnaRY2DnaRYComplement_<T>::VALUE[3] = {'Y', 'R', 'N'};

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableDnaKM2DnaKMComplement_
{
	static char const VALUE[3];
};
template <typename T>
	char const TranslateTableDnaKM2DnaKMComplement_<T>::VALUE[3] = {'M', 'K', 'N'};


//////////////////////////////////////////////////////////////////////////////
//test for match (Ns never match)
inline bool isMatch(Triplex val1, Triplex val2)
{
	if (val1.value == 8 || val2.value == 8) // 'N'
		return false;
	if (val1.value == val2.value)
		return true;
	else 
		return false;
}
	
template <typename TAlphabet>
inline bool isMatch(TAlphabet val1, TAlphabet val2)
{
	if (val1 == val2)
		return true;
	else 
		return false;
}

//////////////////////////////////////////////////////////////////////////////
//Repeat mask

template <>
inline bool _repeatMaskValue(Triplex const &val) 
{
	return val == unknownValue<Triplex>(); // 'N'
	
}
	
template <>
inline bool _repeatMaskValue(DnaRY const &val)
{
	return val == unknownValue<DnaRY>(); // 'N'
}

template <>
inline bool _repeatMaskValue(DnaKM const &val)
{
	return val == unknownValue<DnaKM>(); // 'N'
}
	
//////////////////////////////////////////////////////////////////////////////
//typedefs

typedef String<Triplex, Alloc<void> > TriplexString;
typedef String<DnaRY, Alloc<void> > DnaRYString;
typedef String<DnaKM, Alloc<void> > DnaKMString;

typedef ModView< FunctorComplement<Triplex> >	ModComplementTriplex;
typedef ModView< FunctorComplement<DnaRY> >	ModComplementDnaRY;
typedef ModView< FunctorComplement<DnaKM> >	ModComplementDnaKM;

typedef ModifiedString<TriplexString, ModView< FunctorComplement<Triplex> > >		TriplexStringComplement;
typedef ModifiedString<DnaRYString, ModView< FunctorComplement<DnaRY> > >		DnaRYStringComplement;
typedef ModifiedString<DnaKMString, ModView< FunctorComplement<DnaKM> > >		DnaKMStringComplement;

typedef ModifiedString<
			ModifiedString<	TriplexString, ModView< FunctorComplement<Triplex> > >,
			ModReverse
		>	TriplexStringReverseComplement;

typedef ModifiedString<
			ModifiedString<	DnaRYString, ModView< FunctorComplement<DnaRY> > >,
			ModReverse
		>	DnaRYStringReverseComplement;

typedef ModifiedString<
			ModifiedString<	DnaKMString, ModView< FunctorComplement<DnaKM> > >,
			ModReverse
		>	DnaKMStringReverseComplement;


template <>
struct FunctorComplement<Triplex> : public ::std::unary_function<Triplex,Triplex>
{
    inline Triplex operator()(Triplex x) const {
		return TranslateTableTriplex2TriplexComplement_<>::VALUE[x.value];
	}
};

template <>
struct FunctorComplement<DnaRY> : public ::std::unary_function<DnaRY,DnaRY>
{
    inline DnaRY operator()(DnaRY x) const {
		return TranslateTableDnaRY2DnaRYComplement_<>::VALUE[x.value];
	}
};


template <>
struct FunctorComplement<DnaKM> : public ::std::unary_function<DnaKM,DnaKM>
{
    inline DnaKM operator()(DnaKM x) const {
		return TranslateTableDnaKM2DnaKMComplement_<>::VALUE[x.value];
	}
};

}

#endif
