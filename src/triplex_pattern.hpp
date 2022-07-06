#ifndef TRIPLEX_PATTERN_HPP
#define TRIPLEX_PATTERN_HPP

#include <seqan/basic.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>

#include "triplex_functors.hpp"

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
// ModStringTriplex
//////////////////////////////////////////////////////////////////////////////


/**
.Spec.ModStringTriplex:
..cat:ModStringTriplexs
..summary:Container class for a triplex pattern. Holds the TFO segment and TTS segment respectively
..general:Class.ModStringTriplex
..signature:ModStringTriplex<THost, TString>
..param.THost:Type of the host sequence.
..param.TString:The TTS pattern string
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..param.parallel:Boolean indicating if the TFO pattern is parallel or antiparallel to the purines in the TTS
*/

///.Metafunction.Host.param.T.type:Class.ModStringTriplex


template <typename THost_, typename TString>
class ModStringTriplex
{
public:
	typedef typename Value<TString>::Type			TAlphabet;
	typedef typename Host<ModStringTriplex>::Type 	THost;
	typedef typename Infix<THost>::Type				TSegment;
	typedef unsigned								TId;
	typedef Pair<TId, typename Position<ModStringTriplex<THost_, TString> >::Type >	TDuplicate;
	typedef ::std::vector<TDuplicate>				TDuplicates;	

	typedef ModifiedString<TSegment, ModView< FunctorTCMotif > >		TtcMotif;
	typedef ModifiedString<TSegment, ModView< FunctorGTMotif > >		TgtMotif;
	typedef ModifiedString< ModifiedString<TSegment, ModView< FunctorGTMotif > >,ModReverse>  TgtMotifRev;
	typedef ModifiedString< ModifiedString<TSegment, ModView< FunctorGAMotif > >,ModReverse>  	TgaMotif;
	typedef ModifiedString<TSegment, ModView< FunctorTTSMotif > >		TttsMotif;
	typedef ModifiedString< ModifiedString<TSegment, ModView< FunctorTTSMotifCompl > >, ModReverse>	TttsMotifRevComp;
	
	TString mask_string;					// a string which masks the host segment
	bool parallel;							// whether the motif string is parallel or antiparallel oriented to the host
	TSegment segment;						// the host segment masked by the string
	TId seqNo;								// the sequence id
	int copies;								// number of copies
	bool isTFO;								// indicates whether this is a TFO (=true) or a TTS (=false)
	char motif;								// the motif of the triplex
	double score;
	TDuplicates duplicates;
	
	void _updateMaskString(){
		
		if (!isTFO){
			if (motif == '+')
				mask_string = TttsMotif(segment);
			else
				mask_string = TttsMotifRevComp(segment);
		} else {
			if (parallel){
				if (motif == 'M')
					mask_string = TgtMotif(segment);
				else if (motif == 'Y')
					mask_string = TtcMotif(segment);
			} else {
				if (motif == 'M')
					mask_string = TgtMotifRev(segment);
				else if (motif == 'R')
					mask_string = TgaMotif(segment);
			}
		}
	}
	
//____________________________________________________________________________

public:

	ModStringTriplex():
		mask_string(""),
		segment()
	{
	}
	
	
	ModStringTriplex(typename Parameter_<THost>::Type _host, 
					 bool _parallel_orientation, 
					 unsigned _seqNo, 
					 bool _isTFO, 
					 char _motif
					 ):
		parallel(_parallel_orientation),
		segment(_host),
		seqNo(_seqNo),
		isTFO(_isTFO),
		motif(_motif)
	{
		_updateMaskString();
		copies = -1;
	}
	
	ModStringTriplex(TSegment _segment,  
					 bool _parallel_orientation, 
					 unsigned _seqNo, 
					 bool _isTFO, 
					 char _motif
					 ):
		parallel(_parallel_orientation),
		segment(_segment),
		seqNo(_seqNo),
		isTFO(_isTFO),
		motif(_motif)
	{
		_updateMaskString();
		copies = -1;
	}
	
	
	ModStringTriplex(typename Parameter_<THost>::Type _host, 
					 typename Position<THost>::Type _begin_index, 
					 typename Position<THost>::Type _end_index,  
					 bool _parallel_orientation, 
					 unsigned _seqNo, 
					 bool _isTFO, 
					 char _motif
					 ):
		parallel(_parallel_orientation),
		segment(_host, _begin_index, _end_index),
		seqNo(_seqNo),
		isTFO(_isTFO),
		motif(_motif)
	{
		_updateMaskString();
		copies = -1;
	}

	ModStringTriplex(typename Parameter_<THost>::Type _host, 
					 typename Position<THost>::Type _begin_index, 
					 typename Position<THost>::Type _end_index,  
					 bool _parallel_orientation, 
					 unsigned _seqNo, 
					 bool _isTFO, 
					 char _motif, 
					 unsigned _copies):
		parallel(_parallel_orientation),
		segment(_host, _begin_index, _end_index),
		seqNo(_seqNo),
		isTFO(_isTFO),
		motif(_motif),
		copies(_copies)
	{
		_updateMaskString();
	}
	
	ModStringTriplex(typename Parameter_<THost>::Type _host, 
					 typename Iterator<THost, Standard>::Type _begin, 
					 typename Iterator<THost, Standard>::Type _end,  
					 bool _parallel_orientation, 
					 unsigned _seqNo, 
					 bool _isTFO, 
					 char _motif, 
					 int _copies):
		parallel(_parallel_orientation),
		segment(_host, _begin, _end),
		seqNo(_seqNo),
		isTFO(_isTFO),
		motif(_motif),
		copies(_copies)
	{
		_updateMaskString();
	}

	~ ModStringTriplex()
	{
	}

	template <typename TSource>
	inline ModStringTriplex &
	operator = (TSource const & source)
	{
		assign(*this, source);
		return *this;
	}

//	inline ModStringTriplex & operator = (ModStringTriplex const & source)
//	{
//
////		assign(*this, source);
//		this = source;
//		return *this;
//	}

	template <typename TPos>
	inline typename Reference<ModStringTriplex>::Type
	operator [] (TPos pos)
	{
		return value(*this, pos);
	}

	template <typename TPos>
	inline typename Reference<ModStringTriplex const>::Type
	operator [] (TPos pos) const
	{
		return value(*this, pos);
	}
};

///.Function.value.param.container.type:Class.String

template <typename THost, typename TString, typename TPos>
inline typename Reference< ModStringTriplex<THost, TString> >::Type
value(ModStringTriplex<THost, TString> & me, TPos pos)
{
	return *(begin(me.mask_string, Standard()) + pos);
}


template <typename THost, typename TString, typename TPos>
inline typename Reference< ModStringTriplex<THost, TString> const >::Type
value(ModStringTriplex<THost, TString> const & me, TPos pos)
{
	return *(begin(me.mask_string, Standard()) + pos);
}

//////////////////////////////////////////////////////////////////////////////

///Function.host.param.object.type:Class.ModStringTriplex

template <typename THost_, typename TString>
inline typename Parameter_<THost_>::Type
host(ModStringTriplex<THost_, TString> & me)
{
	return host(me.segment);
}

template <typename THost_, typename TString>
inline typename Parameter_<THost_>::Type
host(ModStringTriplex<THost_, TString> const & me)
{
	return host(me.segment);
}


//////////////////////////////////////////////////////////////////////////////
// stream operators
//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename THost, typename TString>
inline TStream &
operator << (TStream & target,
		ModStringTriplex<THost, TString > const & source)
{
	write(target, ttsString(source));
	return target;
}


//____________________________________________________________________________

///.Function.begin.param.object.type:Class.ModStringTriplex

template <typename THost_, typename TString>
inline typename Iterator<ModStringTriplex<THost_,TString>, Standard>::Type
begin(ModStringTriplex<THost_,TString>& me, Standard)
{
	return begin(me.mask_string);
}
template <typename THost_, typename TString>
inline typename Iterator<ModStringTriplex<THost_,TString> const, Standard>::Type
begin(ModStringTriplex<THost_,TString> const & me, Standard)
{
	return begin(me.mask_string);
}
	
//____________________________________________________________________________

///.Function.merge.param.object.type:Class.ModStringTriplex
// @TODO: merge with respect to reverse complement	
template <typename THost_, typename TString>
bool merge(ModStringTriplex<THost_,TString> & m1, ModStringTriplex<THost_,TString> & m2){
	if (m1.seqNo != m2.seqNo || m1.motif != m2.motif || m1.parallel != m2.parallel)
		return false;
	else {
		setBeginPosition(m1, std::min(beginPosition(m1),beginPosition(m2)));
		setEndPosition(m1, std::max(endPosition(m1),endPosition(m2)));
		m1._updateMaskString();
		return true;
	}
}


//____________________________________________________________________________

///.Function.beginPosition.param.object.type:Class.ModStringTriplex

template <typename THost_, typename TString>
inline typename Position<ModStringTriplex<THost_,TString> >::Type
beginPosition(ModStringTriplex<THost_,TString> & me)
{
	return beginPosition(me.segment);
}
template <typename THost_, typename TString>
inline typename Position<ModStringTriplex<THost_,TString> const>::Type
beginPosition(ModStringTriplex<THost_,TString> const & me)
{
	return beginPosition(me.segment);
}

//____________________________________________________________________________

///.Function.isParallel.param.object.type:Class.ModStringTriplex

template <typename THost_, typename TString>
inline bool
isParallel(ModStringTriplex<THost_,TString> & me)
{
	return me.parallel;
}	
	
template <typename THost_, typename TString>
inline bool
isParallel(ModStringTriplex<THost_,TString> const & me)
{
	return me.parallel;
}	

//____________________________________________________________________________


/**
 .Function.addDuplicate:
	..summary:Adds a new object to the duplicate container
	..cat:Dependent Objects
	..signature:addDuplicate(object, id, position)
	..param.object:The object that will be associated with a new duplicate.
	...type:Class.ModStringTriplex
	..param.id:the sequence identifier of the duplicate.
	..param.position:the position the duplicate starts in the sequence identified by id.
	...text:Note id must be of the same type as the id used in object.
	..see:Function.getDuplicates
	..see:Function.getDuplicateAt
	*/
	
template <typename THost_, typename TString, typename TId, typename TPos>
inline void
addDuplicate(ModStringTriplex<THost_,TString> & me, TId const & seqnr, TPos const & pos)
{
	typename ModStringTriplex<THost_,TString>::TDuplicate d(seqnr, pos);
	appendValue(me.duplicates, d);
}	

//____________________________________________________________________________

/**
 .Function.getDuplicates:
	..summary:Returns the container of duplicates indicated by sequence id and start position
	..cat:Dependent Objects
	..signature:getDuplicates(object)
	..param.object:The object the duplicate list is return from.
	...type:Class.ModStringTriplex
	..returns: a container of duplicate pairs <id, position>
	..see:Function.addDuplicate
	..see:Function.getDuplicateAt
	*/

template <typename THost_, typename TString>
inline typename ModStringTriplex<THost_,TString>::TDuplicates
getDuplicates(ModStringTriplex<THost_,TString> & me)
{
	return me.duplicates;
}
	
template <typename THost_, typename TString>
inline typename ModStringTriplex<THost_,TString>::TDuplicates
getDuplicates(ModStringTriplex<THost_,TString> const & me)
{
	return me.duplicates;
}
	
//____________________________________________________________________________

/**
 .Function.getDuplicateAt:
	..summary:Returns the duplicate at the requested index
	..cat:Dependent Objects
	..signature:getDuplicateAt(object, index)
	..param.object:The object the duplicate is return from.
	..param.index:the index in the container that specifies this duplicate
	...type:Class.ModStringTriplex
	..returns: a duplicate pair<id, position>
	..see:Function.addDuplicate
	..see:Function.getDuplicates
	*/

template <typename THost_, typename TString, typename TIndex>
inline typename ModStringTriplex<THost_,TString>::TDuplicate
getDuplicateAt(ModStringTriplex<THost_,TString> & me, TIndex const & pos)
{
	return value(me.duplicates, pos);
}

template <typename THost_, typename TString, typename TIndex>
inline typename ModStringTriplex<THost_,TString>::TDuplicate
getDuplicateAt(ModStringTriplex<THost_,TString> const & me, TIndex const & pos)
{
	return value(me.duplicates, pos);
}
	
//____________________________________________________________________________

/**
 .Function.isTFO:
	..summary:Indicates whether the object is of type tfo
	..cat:Dependent Objects
	..signature:isTFO(object)
	..param.object:The object
	...type:Class.ModStringTriplex
	..returns:true if object is of type TFO, false otherwise
	*/

template <typename THost_, typename TString>
inline bool
isTFO(ModStringTriplex<THost_,TString> & me)
{
	return me.isTFO;
}	

template <typename THost_, typename TString>
inline bool
isTFO(ModStringTriplex<THost_,TString> const & me)
{
	return me.isTFO;
}	
//____________________________________________________________________________

/**
 .Function.getMotif:
	..summary:returns the motif type of this object
	..cat:Dependent Objects
	..signature:getMotif(object)
	..param.object:The object
	...type:Class.ModStringTriplex
	..returns:motif type (either R,Y,M for TFOs or +,- for TTSs)
	*/

template <typename THost_, typename TString>
inline char
getMotif(ModStringTriplex<THost_,TString> & me)
{
	return me.motif;
}
template <typename THost_, typename TString>
inline char
getMotif(ModStringTriplex<THost_,TString> const & me)
{
	return me.motif;
}
	
//____________________________________________________________________________

///.Function.getSequenceNo.param.object.type:Class.ModStringTriplex
	
template <typename THost_, typename TString>
inline typename Id<ModStringTriplex<THost_,TString> >::Type
getSequenceNo(ModStringTriplex<THost_,TString> & me)
{
	return me.seqNo;
}
	
	
template <typename THost_, typename TString>
inline typename Id<ModStringTriplex<THost_,TString> >::Type
getSequenceNo(ModStringTriplex<THost_,TString> const & me)
{
	return me.seqNo;
}

//____________________________________________________________________________

/**
 .Function.setSequenceNo:
	..summary:Assigns a new sequence no to the object.
	..cat:Dependent Objects
	..signature:setSequenceNo(object, new_seqno)
	..param.object:An object.
	...type:Spec.ModStringTriplex
	...type:Spec.ModStringTriplex
	..param.new_seqno:new sequence no 
	*/

template <typename THost_, typename TString>
inline void
setSequenceNo(ModStringTriplex<THost_, TString> & me, typename Id<ModStringTriplex<THost_, TString> >::Type new_seqno)
{
	me.seqNo = new_seqno;
}
	
//____________________________________________________________________________

/**
.Function.setBegin:
..summary:Sets begin of object in host.
..cat:Dependent Objects
..signature:setBegin(object, new_begin)
..param.object:An object.
...type:Spec.ModStringTriplex
...type:Spec.ModStringTriplex
..param.new_begin:iterator to the new first item in $host(object)$ that belongs of $object$.
...type:Metafunction.Iterator
..see:Function.begin
..see:Function.beginPosition
*/
template <typename THost_, typename TString, typename TIterator>
inline void
setBegin(ModStringTriplex<THost_, TString> & me, TIterator new_begin)
{
	me.data_begin_position = new_begin - begin(host(me));//, Standard());
}


//____________________________________________________________________________

/**
.Function.setBeginPosition:
..summary:Sets begin position of object in host.
..cat:Dependent Objects
..signature:setBeginPosition(object, new_begin)
..param.object:An object.
...type:Spec.ModStringTriplex
...type:Spec.ModStringTriplex
..param.new_begin:position of the new first item in $host(object)$ that belongs of $object$.
...type:Metafunction.Position
..see:Function.begin
..see:Function.beginPosition
..see:Function.setBegin
*/

template <typename THost_, typename TString_, typename TPosition_>
inline void
setBeginPosition(ModStringTriplex<THost_, TString_> & me, TPosition_ new_begin)
{
	setBeginPosition(me.segment,new_begin);
}

//____________________________________________________________________________

///.Function.end.param.object.type:Class.ModStringTriplex

template <typename THost_, typename TString_>
inline typename Iterator<ModStringTriplex<THost_, TString_>, Standard>::Type
end(ModStringTriplex<THost_, TString_> & me,
	Standard)
{
	return end(me.mask_string);
}
template <typename THost_, typename TString_>
inline typename Iterator<ModStringTriplex<THost_, TString_> const, Standard>::Type
end(ModStringTriplex<THost_, TString_> const & me,
	Standard)
{
	return end(me.mask_string);
}

//____________________________________________________________________________

///.Function.endPosition.param.object.type:Class.ModStringTriplex

template <typename THost_, typename TString>
inline typename Position<ModStringTriplex<THost_, TString> >::Type
endPosition(ModStringTriplex<THost_, TString> & me)
{
	return endPosition(me.segment);
}
template <typename THost_, typename TString>
inline typename Position<ModStringTriplex<THost_, TString> >::Type
endPosition(ModStringTriplex<THost_, TString> const & me)
{
	return endPosition(me.segment);
}

//____________________________________________________________________________

///.Function.score.param.object.type:Class.ModStringTriplex

template <typename THost_, typename TString>
inline double
score(ModStringTriplex<THost_, TString> & me)
{
	return me.score;
}
template <typename THost_, typename TString>
inline double
score(ModStringTriplex<THost_, TString> const & me)
{
	return me.score;
}

//____________________________________________________________________________

/**
 .Function.setScore:
	..summary:Assigns a new score to the object
	..cat:Dependent Objects
	..signature:setScore(object, score)
	..param.object:The object
	..param.score:The new score assigned to the object	 
	...type:Class.ModStringTriplex
	*/

template <typename THost_, typename TString_, typename TScore>
inline void
setScore(ModStringTriplex<THost_, TString_> & me, TScore score)
{
	me.score = score;
}

//____________________________________________________________________________

///.Function.getSegment.param.object.type:Class.ModStringTriplex

template <typename THost, typename TString>
inline typename Infix<THost>::Type
getSegment(ModStringTriplex<THost, TString> & me)
{
	return me.segment;
}
template <typename THost, typename TString>
inline typename Infix<THost>::Type const
getSegment(ModStringTriplex<THost, TString> const & me)
{
	return me.segment;
}

//____________________________________________________________________________

/**
 .Function.guanineRate:
	..summary:Returns the proportion of guanines in the target sequence available for triplex formation
	..cat:Dependent Objects
	..signature:errorString(object)
	..param.object:The object
	..returns:guanine rate normalized over feature length 
	...type:Class.ModStringTriplex
	*/

template <typename THost, typename TString>
inline double
guanineRate(ModStringTriplex<THost, TString> & me)
{
	double guanines = 0.;
	for (unsigned i=0; i<length(me); ++i){
		if (value(me.mask_string,i)=='G' || value(me.mask_string,i)=='g') ++guanines;
	}
	return guanines/length(me);
}

template <typename THost, typename TString>
inline double 
guanineRate(ModStringTriplex<THost, TString> const & me)
{
	double guanines = 0.;
	for (unsigned i=0; i<length(me); ++i){
		if (value(me.mask_string,i)=='G' || value(me.mask_string,i)=='g') ++guanines;
	}
	return guanines/length(me);
}	


//____________________________________________________________________________

/**
 .Function.errorString:
	..summary:Returns a string indicating which position violates triplex rules
	..cat:Dependent Objects
	..signature:errorString(object)
	..param.object:The object
	..returns:String of mismatches (errors) 
	...Note:similar to a cigar string all errors are concatenated,
	o - stands for a mismatch in the oligo (TFO)
	d - stands for a mismatch in the duplex (TTS)
	...type:Class.ModStringTriplex
	*/

template <typename THost, typename TString>
inline CharString
errorString(ModStringTriplex<THost, TString> & me)
{
	std::ostringstream errors;	
	
	CharString ps = prettyString(me);
	if (me.isTFO){
		for (unsigned i=0; i<length(ps); ++i){
			if (! isupper(value(ps,i))) errors << 'o' << i;
		}
	} else {
		if (me.motif == '+'){
			for (unsigned i=0; i<length(ps); ++i){
				if (! isupper(value(ps,i))) errors << 'd' << i;
			}
		} else{
			for (unsigned i=length(ps); i>0; --i){
				if (! isupper(value(ps,i-1))) errors << 'd' << (length(ps)-i);
			}			
		}	
	}
	return errors.str();
}
template <typename THost, typename TString>
inline CharString const
errorString(ModStringTriplex<THost, TString> const & me)
{
	std::ostringstream errors;	
	CharString ps = prettyString(me);
	if (me.isTFO){
		for (unsigned i=0; i<length(ps); ++i){
			if (! isupper(value(ps,i))) errors << 'o' << i;
		}
	} else {
		if (me.motif == '+'){
			for (unsigned i=0; i<length(ps); ++i){
				if (! isupper(value(ps,i))) errors << 'd' << i;
			}
		} else{
			for (unsigned i=length(ps); i>0; --i){
				if (! isupper(value(ps,i-1))) errors << 'd' << (length(ps)-i);
			}			
		}	
	}
	return errors.str();
}	

//____________________________________________________________________________

/**
 .Function.ttsString:
	..summary:Returns the target string (TTS)
	..cat:Dependent Objects
	..signature:ttsString(object)
	..param.object:The object
	..returns:target string
	...type:Class.ModStringTriplex
	*/

template <typename THost, typename TString>
inline TString
ttsString(ModStringTriplex<THost, TString> & me)
{
	return me.mask_string;
}
template <typename THost, typename TString>
inline TString const
ttsString(ModStringTriplex<THost, TString> const & me)
{
	return me.mask_string;
}

//____________________________________________________________________________

/**
 .Function.tfoString:
	..summary:Returns the TFO
	..cat:Dependent Objects
	..signature:tfoString(object)
	..param.object:The object
	..returns:tfo string
	...Note: string is in 5'-3' orientation
	...type:Class.ModStringTriplex
	*/
template <typename THost, typename TString>
inline THost
tfoString(ModStringTriplex<THost, TString> & me)
{
	THost tempSeq2(me.segment);
	return tempSeq2;
}
template <typename THost, typename TString>
inline THost const
tfoString(ModStringTriplex<THost, TString> const & me)
{
	THost tempSeq2(me.segment);
	return tempSeq2;
}

//____________________________________________________________________________

/**
 .Function.outputString:
	..summary:Prepares a string that mimics the triplex feature space for output
	..cat:Dependent Objects
	..signature:outputString(object)
	..param.object:The object
	..returns:string for output
	...Note: The string is masked according to the triplex motif
	...type:Class.ModStringTriplex
	..see:prettyString
	*/

template <typename THost, typename TString>
inline CharString
outputString(ModStringTriplex<THost, TString> & me)
{
	typedef ModifiedString<THost, ModView< FunctorTTSMotifOutput > >  TttsMotifOutput;
	typedef ModifiedString< ModifiedString<THost, ModView< FunctorTTSMotifComplOutput > >, ModReverse>	TttsMotifRevCompOutput;
	typedef ModifiedString<THost, ModView< FunctorTCMotifOutput > >  TtcMotifOutput;
	typedef ModifiedString<THost, ModView< FunctorGTMotifOutput > >  TgtMotifOutput;
	typedef ModifiedString<THost, ModView< FunctorGAMotifOutput > >  TgaMotifOutput;
	
	if (!me.isTFO){
		if (me.motif == '+'){
			return TttsMotifOutput(me.segment);
		} else {
			return TttsMotifRevCompOutput(me.segment);
		}
	} else {
		if (me.motif == 'M')
			return TgtMotifOutput(me.segment);
		else if (me.motif == 'Y')
			return TtcMotifOutput(me.segment);
		else if (me.motif == 'R')
			return TgaMotifOutput(me.segment);
		else
			return "";
	}
}
template <typename THost, typename TString>
inline CharString const
outputString(ModStringTriplex<THost, TString> const & me)
{
	typedef ModifiedString<THost, ModView< FunctorTTSMotifOutput > >  TttsMotifOutput;
	typedef ModifiedString< ModifiedString<THost, ModView< FunctorTTSMotifComplOutput > >, ModReverse>	TttsMotifRevCompOutput;
	typedef ModifiedString<THost, ModView< FunctorTCMotifOutput > >  TtcMotifOutput;
	typedef ModifiedString<THost, ModView< FunctorGTMotifOutput > >  TgtMotifOutput;
	typedef ModifiedString<THost, ModView< FunctorGAMotifOutput > >  TgaMotifOutput;
	
	if (!me.isTFO){
		if (me.motif == '+'){
			return TttsMotifOutput(me.segment);
		} else {
			return TttsMotifRevCompOutput(me.segment);
		}
	} else {
		if (me.motif == 'M')
			return TgtMotifOutput(me.segment);
		else if (me.motif == 'Y')
			return TtcMotifOutput(me.segment);
		else if (me.motif == 'R')
			return TgaMotifOutput(me.segment);
		else
			return "";
	}
}
	
//____________________________________________________________________________

/**
 .Function.prettyString:
	..summary:Prepares a string that mimics the triplex feature space for output
	..cat:Dependent Objects
	..signature:prettyString(object)
	..param.object:The object
	..returns:string for output
	...Note: Same as outputString but mismatches are converted into small letters to make
	them more obvious
	...type:Class.ModStringTriplex
	..see:outputString
	*/


template <typename THost, typename TString>
inline CharString
prettyString(ModStringTriplex<THost, TString> & me)
{
	typedef typename Infix<THost>::Type	TSegment;
	typedef ModifiedString<TSegment, ModView< FunctorTTSMotifPretty > >  TttsMotifOutput;
	typedef ModifiedString< ModifiedString<TSegment, ModView< FunctorTTSMotifComplPretty > >, ModReverse>	TttsMotifRevCompOutput;
	typedef ModifiedString<TSegment, ModView< FunctorTCMotifPretty > >  TtcMotifOutput;
	typedef ModifiedString<TSegment, ModView< FunctorGTMotifPretty > >  TgtMotifOutput;
	typedef ModifiedString<TSegment, ModView< FunctorGAMotifPretty > >  TgaMotifOutput;
	
	if (!me.isTFO){
		if (me.motif == '+'){
			return TttsMotifOutput(me.segment);
		} else {
			return TttsMotifRevCompOutput(me.segment);
		}
	} else {
		if (me.motif == 'M')
			return TgtMotifOutput(me.segment);
		else if (me.motif == 'Y')
			return TtcMotifOutput(me.segment);
		else if (me.motif == 'R')
			return TgaMotifOutput(me.segment);
		else
			return "";
	}
}

template <typename THost, typename TString>
inline CharString const
prettyString(ModStringTriplex<THost, TString> const & me)
{
	typedef typename Infix<THost>::Type	TSegment;
	typedef ModifiedString<TSegment, ModView< FunctorTTSMotifPretty > >  TttsMotifOutput;
	typedef ModifiedString< ModifiedString<TSegment, ModView< FunctorTTSMotifComplPretty > >, ModReverse>	TttsMotifRevCompOutput;
	typedef ModifiedString<TSegment, ModView< FunctorTCMotifPretty > >  TtcMotifOutput;
	typedef ModifiedString<TSegment, ModView< FunctorGTMotifPretty > >  TgtMotifOutput;
	typedef ModifiedString<TSegment, ModView< FunctorGAMotifPretty > >  TgaMotifOutput;
	
	if (!me.isTFO){
		if (me.motif == '+'){
			return TttsMotifOutput(me.segment);
		} else {
			return TttsMotifRevCompOutput(me.segment);
		}
	} else {
		if (me.motif == 'M')
			return TgtMotifOutput(me.segment);
		else if (me.motif == 'Y')
			return TtcMotifOutput(me.segment);
		else if (me.motif == 'R')
			return TgaMotifOutput(me.segment);
		else
			return "";
	}
}
	
/**
.Function.setEnd:
..summary:Sets end of object in host.
..cat:Dependent Objects
..signature:setEnd(object, new_end)
..param.object:An object.
...type:Spec.ModStringTriplex
...type:Spec.ModStringTriplex
..param.new_end:Iterator behind the last item in $host(object)$ belongs of $object$.
...type:Metafunction.Iterator
..see:Function.end
..see:Function.endPosition
..see:Function.setBegin
*/

template <typename THost_, typename TString, typename TIterator>
inline void
setEnd(ModStringTriplex<THost_, TString> & me, TIterator new_end)
{
	me.data_end_position = new_end - begin(host(me));//, Standard());
}

//____________________________________________________________________________


/**
.Function.setEndPosition:
..summary:Sets begin position of object in host.
..cat:Dependent Objects
..signature:setEndPosition(object, new_end)
..param.object:An object.
...type:Spec.ModStringTriplex
...type:Spec.ModStringTriplex
..param.new_end:position behind the last item in $host(object)$ that belongs of $object$.
...type:Metafunction.Position
..see:Function.end
..see:Function.endPosition
..see:Function.setBeginPosition
..see:Function.setEnd
*/

template <typename THost_, typename TString_, typename TPosition_>
inline void
setEndPosition(ModStringTriplex<THost_, TString_> & me, TPosition_ new_end)
{
	setEndPosition(me.segment,new_end);
//	me.data_end_position = new_end;
}

//____________________________________________________________________________

///.Function._setLength.param.object.type:Class.ModStringTriplex

template <typename THost_, typename TString>
inline void
_setLength(
	ModStringTriplex<THost_, TString> & me,
	typename Size<THost_>::Type new_length)
{
	me.data_end_position = me.data_begin_position + new_length;
}

	
//____________________________________________________________________________

/**
.Function.setHost:
..summary:Sets the host of an object.
..cat:Dependent Objects
..signature:setHost(object, host)
..param.object:The object that will get a new host.
...type:Class.ModStringTriplex
..param.host:The new host.
..remarks:After this operation, $object$ depends on $host$.
...text:Note that setting the host can invalidate $object$.
For example, if one changes the host of a @Class.ModStringTriplex@ object, it is possible
that begin- and end-position of the ModStringTriplex does not fit into the new host sequence.
..see:Function.host
*/
template <typename THost_, typename TString>
inline void
setHost(ModStringTriplex<THost_, TString> & me, typename Parameter_<THost_>::Type _host)
{
	me.data_host = _toPointer(_host);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.length.param.object.type:Class.ModStringTriplex

template <typename THost, typename TString>
inline typename Size<ModStringTriplex<THost, TString> const>::Type
length(ModStringTriplex<THost, TString> const & me)
{
	return endPosition(me) - beginPosition(me);
}

template <typename THost, typename TString>
inline TString
substr(ModStringTriplex<THost, TString> & me, int beginPos, int endPos)
{
	typedef typename Iterator<TString, Rooted>::Type	TIter;

	::std::string tmp_sub;
	tmp_sub.reserve(endPos-beginPos);

	TIter it = begin(me.mask_string, Rooted());
	TIter itEnd = begin(me.mask_string, Rooted());
	goFurther(it, beginPos);
	goFurther(itEnd, endPos);
	for(;it != itEnd;++it){
		tmp_sub+=*(it);
	}
	TString tmp_tstring(tmp_sub);
	return tmp_tstring;
}

//////////////////////////////////////////////////////////////////////////////
// METAFUNCTIONS
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Class.String

//////////////////////////////////////////////////////////////////////////////

template <typename THost_, typename TString>
struct Value<ModStringTriplex<THost_, TString> >
{
	typedef typename Value<TString>::Type Type;
};
template <typename THost_, typename TString>
struct Value<ModStringTriplex<THost_, TString> const >:
	public Value<ModStringTriplex<THost_, TString> >
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Host.param.T.type:Class.String

template <typename THost_, typename TString>
struct Host<ModStringTriplex<THost_, TString> >
{
	typedef THost_ Type;
};
template <typename THost_, typename TString>
struct Host<ModStringTriplex<THost_, TString> const>:
	public Host<ModStringTriplex<THost_, TString> >
{
};


//////////////////////////////////////////////////////////////////////////////

///.Metafunction.IsSequence.param.T.type:Class.String

template <typename THost_, typename TString>
struct IsSequence<ModStringTriplex<THost_, TString> > {
	typedef True Type;
	enum { VALUE = true };
};

	


//////////////////////////////////////////////////////////////////////////////

///.Metafunction.orientation.param.T.type:Class.String

template <typename THost_, typename TString>
inline char
orientationString(ModStringTriplex<THost_, TString> & me)
{
	if (me.parallel)
		return 'P';
	else
		return 'A';
}


template <typename THost_, typename TString>
inline char
orientationString(ModStringTriplex<THost_, TString> const & me)
{
	if (me.parallel)
		return 'P';
	else
		return 'A';
}
	
	
//____________________________________________________________________________

///.Function.begin.param.object.type:Class.ModStringTriplex

template <typename THost_, typename TString>
inline int
duplicates(ModStringTriplex<THost_,TString> & me)
{
	return me.copies;
}
template <typename THost_, typename TString>
inline int
duplicates(ModStringTriplex<THost_,TString> const & me)
{
	return me.copies;
}

template <typename THost_, typename TString, typename TPos>
inline void
duplicates(ModStringTriplex<THost_,TString> & me, TPos & copies)
{
	me.copies = copies;
}


//////////////////////////////////////////////////////////////////////////////

///.Metafunction.GetValue.param.T.type:Class.ModStringTriplex

template <typename THost_, typename TString>
struct GetValue<ModStringTriplex<THost_, TString> >
{
	typedef typename GetValue<TString>::Type Type;
};

template <typename THost_, typename TString>
struct GetValue<ModStringTriplex<THost_, TString> const >
{
	typedef typename GetValue<TString const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Class.ModStringTriplex

template <typename THost_, typename TString>
struct Iterator<ModStringTriplex<THost_, TString>, Rooted>
{
	typedef ModStringTriplex<THost_, TString> TSequence;
	typedef typename Iterator<TString, Standard>::Type TIterator;
	typedef Iter<TSequence, AdaptorIterator<TIterator> > Type;
};
template <typename THost_, typename TString>
struct Iterator<ModStringTriplex<THost_, TString> const, Rooted>
{
	typedef ModStringTriplex<THost_, TString> const TSequence;
	typedef typename Iterator<TString const, Standard>::Type TIterator;
	typedef Iter<TSequence, AdaptorIterator<TIterator> > Type;
};

template <typename THost_, typename TString>
struct Iterator<ModStringTriplex<THost_, TString>, Standard>:
	Iterator<TString, Standard>
{
};
template <typename THost_, typename TString>
struct Iterator<ModStringTriplex<THost_, TString> const, Standard>:
	Iterator<TString, Standard>
{
};



//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Size.param.T.type:Class.ModStringTriplex

template <typename THost_, typename TString>
struct Size<ModStringTriplex<THost_, TString> >
{
	typedef typename Size<THost_>::Type Type;
};

template <typename THost_, typename TString>
struct Size<ModStringTriplex<THost_, TString> const >
{
	typedef typename Size<THost_>::Type Type;
};

template <typename THost_, typename TString>
struct Position<ModStringTriplex<THost_, TString> >
{
	typedef typename Position<THost_>::Type Type;
};

template <typename THost_, typename TString>
struct Position<ModStringTriplex<THost_, TString> const >
{
	typedef typename Position<THost_>::Type Type;
};
	
template <typename THost_, typename TString>
struct Id<ModStringTriplex<THost_, TString> >
{
	typedef typename ModStringTriplex<THost_, TString>::TId Type;
};
template <typename THost_, typename TString>
struct Id<ModStringTriplex<THost_, TString> const>
{
	typedef typename ModStringTriplex<THost_, TString>::TId Type;
};

}

#endif
