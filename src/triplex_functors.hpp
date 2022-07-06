#ifndef TRIPLEX_FUNCTORS_HPP
#define TRIPLEX_FUNCTORS_HPP

#include "triplex_alphabet.hpp"

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
// Triplex motif functors
//////////////////////////////////////////////////////////////////////////////
/**
 * mask all purines (G,A) in a sequences with 'R' and pyrimidins (C,T/U) as 'Y'
 */
struct FunctorRYFilter : public ::std::unary_function<Triplex,Triplex>
{
	inline Triplex operator()(Triplex x) const {
		if ((x == 'G') || (x == 'g') || (x == 'A') || (x == 'a') || (x == 'R'))
			return 'R';
		else if ((x == 'C') || (x == 'c') || (x == 'T') || (x == 't') || (x == 'U') || (x == 'u') || (x == 'Y'))
			return 'Y';
		else
			return 'N';
	}
};

/**
 * mask all pyrimidins (C,T/U) as 'Y'
 */
struct FunctorGAYFilter : public ::std::unary_function<Triplex,Triplex>
{
	inline Triplex operator()(Triplex x) const {
		if ((x == 'G') || (x == 'g'))
			return 'G';
		else if ((x == 'A') || (x == 'a'))
			return 'A';
		else if (x == 'R')
			return 'R';
		else if ((x == 'C') || (x == 'c') || (x == 'T') || (x == 't') || (x == 'U') || (x == 'u') || (x == 'Y'))
			return 'Y';
		else
			return 'N';
	}
};

/**
 * mask all purines (G,A) in a sequences with 'R'
 */
struct FunctorCTRFilter : public ::std::unary_function<Triplex,Triplex>
{
	inline Triplex operator()(Triplex x) const {
		if ((x == 'C') || (x == 'c'))
			return 'C';
		else if ((x == 'T') || (x == 't') || (x == 'U') || (x == 'u'))
			return 'T';
		else if (x == 'Y')
			return 'Y';
		else if ((x == 'G') || (x == 'g') || (x == 'A') || (x == 'a') || (x == 'R'))
			return 'R';
		else
			return 'N';
	}
};

/**
 * mask all (C,A) in a sequences with 'M'
 */
struct FunctorGTMFilter : public ::std::unary_function<Triplex,Triplex>
{
	inline Triplex operator()(Triplex x) const {
		if ((x == 'G') || (x == 'g'))
			return 'G';
		else if ((x == 'T') || (x == 't') || (x == 'U') || (x == 'u'))
			return 'T';
		else if (x == 'K')
			return 'K';
		else if ((x == 'C') || (x == 'c') || (x == 'A') || (x == 'a') || (x == 'M'))
			return 'M';
		else
			return 'N';
	}
};


/**
 * mask all (T/U,G) in a sequences with 'K' and all (A,C) with 'M'. Covers both DNA and RNA.
 */
struct FunctorKMFilter : public ::std::unary_function<Triplex,Triplex>
{
	inline Triplex operator()(Triplex x) const {
		if ((x == 'C') || (x == 'c') || (x == 'A') || (x == 'a') || (x == 'M'))
			return 'M';
		else if ((x == 'G') || (x == 'g') || (x == 'T') || (x == 't') || (x == 'U') || (x == 'u') || (x == 'K'))
			return 'K';
		else
			return 'N';
	}
};

/**
 * Translate the TC motif into a corresponding triplex target site and mask all remaining
 * character with 'N'
 */
struct FunctorTCMotif : public ::std::unary_function<Triplex,Triplex>
{
	inline Triplex operator()(Triplex x) const {
		if ((x == 'C') || (x == 'c') )
			return 'G';
		else if ((x == 'T') || (x == 't') || (x == 'U') || (x == 'u'))
			return 'A';
		else
			return 'N';
	}
};

/**
 * Translate the TC motif into a corresponding triplex target site and mask all remaining
 * character with 'N'
 */
struct FunctorGAMotif : public ::std::unary_function<Triplex,Triplex>
{
	inline Triplex operator()(Triplex x) const {
		if ((x == 'G') || (x == 'g'))
			return 'G';
		else if ((x == 'A') || (x == 'a') )
			return 'A';
		else
			return 'N';
	}
};

/**
 * Translate the GT motif into a corresponding triplex target site and mask all remaining
 * character with 'N'
 */
struct FunctorGTMotif : public ::std::unary_function<Triplex,Triplex>
{
	inline Triplex operator()(Triplex x) const {
		if ((x == 'G') || (x == 'g') )
			return 'G';
		else if ((x == 'T') || (x == 't') || (x == 'U') || (x == 'u'))
			return 'A';
		else
			return 'N';
	}
};

/**
 * Translate the TC motif into a corresponding triplex target site and mask all remaining
 * character with 'N'
 */
struct FunctorTCMotifPretty : public ::std::unary_function<Triplex,char>
{
	inline char operator()(Triplex x) const {
		if ((x == 'C') || (x == 'c') )
			return 'C';
		else if ((x == 'T') || (x == 't') || (x == 'U') || (x == 'u'))
			return 'T';
		else if ((x == 'G') || (x == 'g') )
			return 'g';
		else if ((x == 'A') || (x == 'a') )
			return 'a';
		else
			return 'n';
	}
};

/**
 * Translate the TC motif into a corresponding triplex target site and mask all remaining
 * character with 'N'
 */
struct FunctorGAMotifPretty : public ::std::unary_function<Triplex,char>
{
	inline char operator()(Triplex x) const {
		if ((x == 'G') || (x == 'g'))
			return 'G';
		else if ((x == 'A') || (x == 'a') )
			return 'A';
		else if ((x == 'C') || (x == 'c') )
			return 'c';
		else if ((x == 'T') || (x == 't') )
			return 't';
		else
			return 'n';
	}
};

/**
 * Translate the GT motif into a corresponding triplex target site and mask all remaining
 * character with 'N'
 */
struct FunctorGTMotifPretty : public ::std::unary_function<Triplex,char>
{
	inline char operator()(Triplex x) const {
		if ((x == 'G') || (x == 'g') )
			return 'G';
		else if ((x == 'T') || (x == 't') || (x == 'U') || (x == 'u'))
			return 'T';
		else if ((x == 'C') || (x == 'c') )
			return 'c';
		else if ((x == 'A') || (x == 'a') )
			return 'a';
		else
			return 'n';
	}
};

/**
 * Translate the TC motif into a corresponding triplex target site and mask all remaining
 * character with 'N'
 */
struct FunctorTCMotifOutput : public ::std::unary_function<Triplex,char>
{
	inline char operator()(Triplex x) const {
		if ((x == 'C') || (x == 'c') )
			return 'C';
		else if ((x == 'T') || (x == 't') || (x == 'U') || (x == 'u'))
			return 'T';
		else if ((x == 'G') || (x == 'g') )
			return 'G';
		else if ((x == 'A') || (x == 'a') )
			return 'A';
		else
			return 'N';
	}
};

/**
 * Translate the TC motif into a corresponding triplex target site and mask all remaining
 * character with 'N'
 */
struct FunctorGAMotifOutput : public ::std::unary_function<Triplex,char>
{
	inline char operator()(Triplex x) const {
		if ((x == 'G') || (x == 'g'))
			return 'G';
		else if ((x == 'A') || (x == 'a') )
			return 'A';
		else if ((x == 'C') || (x == 'c') )
			return 'C';
		else if ((x == 'T') || (x == 't') )
			return 'T';
		else
			return 'N';
	}
};

/**
 * Translate the GT motif into a corresponding triplex target site and mask all remaining
 * character with 'N'
 */
struct FunctorGTMotifOutput : public ::std::unary_function<Triplex,char>
{
	inline char operator()(Triplex x) const {
		if ((x == 'G') || (x == 'g') )
			return 'G';
		else if ((x == 'T') || (x == 't') || (x == 'U') || (x == 'u'))
			return 'T';
		else if ((x == 'C') || (x == 'c') )
			return 'C';
		else if ((x == 'A') || (x == 'a') )
			return 'A';
		else
			return 'N';
	}
};



/**
 * Masks non-purine characters of a putative triplex target site 
 * with 'Y' (Different character than masking for TFOs!)
 */
struct FunctorTTSMotif : public ::std::unary_function<Triplex,Triplex>
{
	inline Triplex operator()(Triplex x) const {
		if ((x == 'G') || (x == 'g'))
			return 'G';
		else if ((x == 'A') || (x == 'a') )
			return 'A';
		else
			return 'Y';
	}
};


/**
 * Masks non-purine characters of a putative triplex target site 
 * with 'Y' (Different character than masking for TFOs!)
 */
struct FunctorTTSMotifPretty : public ::std::unary_function<Triplex,char>
{
	inline char operator()(Triplex x) const {
		if ((x == 'G') || (x == 'g'))
			return 'G';
		else if ((x == 'A') || (x == 'a') )
			return 'A';
		else if ((x == 'C') || (x == 'c') )
			return 'c';
		else if ((x == 'T') || (x == 't') )
			return 't';
		else
			return 'n';
	}
};

/**
 * Masks non-purine characters of a putative triplex target site 
 * with 'Y' (Different character than masking for TFOs!)
 */
struct FunctorTTSMotifOutput : public ::std::unary_function<Triplex,char>
{
	inline char operator()(Triplex x) const {
		if ((x == 'G') || (x == 'g'))
			return 'G';
		else if ((x == 'A') || (x == 'a') )
			return 'A';
		else if ((x == 'C') || (x == 'c') )
			return 'C';
		else if ((x == 'T') || (x == 't') )
			return 'T';
		else
			return 'N';
	}
};


/**
 * Masks non-pyrimidine characters of a putative triplex target site 
 * with 'Y' (Different character than masking for TFOs!)
 */
struct FunctorTTSMotifCompl : public ::std::unary_function<Triplex,Triplex>
{
	inline Triplex operator()(Triplex x) const {
		if ((x == 'C') || (x == 'c'))
			return 'G';
		else if ((x == 'T') || (x == 't') )
			return 'A';
		else
			return 'Y';
	}
};

/**
 * Masks non-pyrimidine characters of a putative triplex target site 
 * with 'Y' (Different character than masking for TFOs!)
 */
struct FunctorTTSMotifComplPretty : public ::std::unary_function<Triplex,char>
{
	inline char operator()(Triplex x) const {
		if ((x == 'C') || (x == 'c'))
			return 'G';
		else if ((x == 'T') || (x == 't') )
			return 'A';
		else if ((x == 'G') || (x == 'g') )
			return 'c';
		else if ((x == 'A') || (x == 'a') )
			return 't';
		else
			return 'n';
	}
};

/**
 * Masks non-pyrimidine characters of a putative triplex target site 
 * with 'Y' (Different character than masking for TFOs!)
 */
struct FunctorTTSMotifComplOutput : public ::std::unary_function<Triplex,char>
{
	inline char operator()(Triplex x) const {
		if ((x == 'C') || (x == 'c'))
			return 'G';
		else if ((x == 'T') || (x == 't') )
			return 'A';
		else if ((x == 'G') || (x == 'g') )
			return 'C';
		else if ((x == 'A') || (x == 'a') )
			return 'T';
		else
			return 'N';
	}
};

}

#endif
