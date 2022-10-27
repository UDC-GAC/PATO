// ==========================================================================
//                                triplexator
// ==========================================================================
// Copyright (c) 2011,2012, Fabian Buske, UQ
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Fabian Buske or the University of Queensland nor 
//       the names of its contributors may be used to endorse or promote products 
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Fabian Buske <fbuske@uq.edu.au>
// ==========================================================================

#ifndef TRIPLEX_FUNCTORS_HPP
#define TRIPLEX_FUNCTORS_HPP

#include "triplex_alphabet.hpp"

namespace seqan
{

/**
 * mask all purines (G,A) in a sequences with 'R' and pyrimidins (C,T/U) as 'Y'
 */
struct FunctorRYFilter : public std::function<Triplex(Triplex)>
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
