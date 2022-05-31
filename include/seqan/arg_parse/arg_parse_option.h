// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_OPTION_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_OPTION_H_

#include <string>
#include <vector>
#include <seqan/arg_parse/arg_parse_argument.h>

namespace seqan {

// ----------------------------------------------------------------------------
// Class ArgParseOption
// ----------------------------------------------------------------------------

/**
.Class.ArgParseOption
..base:Class.ArgParseArgument
..cat:Miscellaneous
..summary:Stores information for a specific command line option.
..signature:ArgParseOption
..remarks:A @Class.ArgParseOption@ object can be added to a @Class.ArgumentParser@ via @Function.addOption@.
..include:seqan/arg_parse.h
*/

/**
.Memfunc.ArgParseOption#ArgParseOption
..class:Class.ArgParseOption
..summary:Constructor
..signature:ArgParseOption (shortName, longName, helpText [, argument])
..param.shortName:A std::string containing the short-name option identifier (e.g. $"h"$ for the $-h/--help$ option).
Although not suggested the short-name can contain more than 1 character.
...remarks:Note that the leading "-" is not passed.
..param.longName:A std::string containing the long-name option identifier (e.g. $"help"$ for the $-h/--help$ option).
...remarks:Note that the leading "--" is not passed.
..param.helpText:A std::string containing the help text associated with this option.
..param.argument:A @Class.ArgParseArgument@ for the option (e.g., an integer argument).
...type:Class.ArgParseArgument
*/

///.Function.isListArgument.param.argument.type:Class.ArgParseOption
///.Function.isStringArgument.param.argument.type:Class.ArgParseOption
///.Function.isIntegerArgument.param.argument.type:Class.ArgParseOption
///.Function.isDoubleArgument.param.argument.type:Class.ArgParseOption
///.Function.isInputFileArgument.param.argument.type:Class.ArgParseOption
///.Function.isOutputFileArgument.param.argument.type:Class.ArgParseOption
///.Function.setMinValue.param.argument.type:Class.ArgParseOption
///.Function.setMaxValue.param.argument.type:Class.ArgParseOption
///.Function.setValidValues.param.argument.type:Class.ArgParseOption
///.Function.getArgumentValue.param.argument.type:Class.ArgParseOption
///.Function.getArgumentValues.param.argument.type:Class.ArgParseOption
///.Function.hasValue.param.argument.type:Class.ArgParseOption
///.Function.isSet.param.argument.type:Class.ArgParseOption
///.Function.numberOfAllowedValues.param.argument.type:Class.ArgParseOption

class ArgParseOption
    : public ArgParseArgument
{
public:
    // ----------------------------------------------------------------------------
    // Members to specify the names of the ArgParseOption
    // ----------------------------------------------------------------------------
    std::string         shortName;     // short option name
    std::string         longName;      // long option name

    // ----------------------------------------------------------------------------
    // Members representing type, content and restrictions of the ArgParseOption
    // ----------------------------------------------------------------------------
    bool                _isFlag;      // true if this a bool option, that has no
                                      // argument we will internally represent it as a
                                      // string option set to either "true" or "false"
    bool                _isRequired;  // true if this ArgParseOption must be set
    bool                _isHidden;    // true if this ArgParseOption should not be
                                      // shown on the command line

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------
    ArgParseOption(std::string const & _shortName,
                   std::string const & _longName,
                   std::string const & _help,
                   ArgumentType argumentType,
                   bool isListArgument = false,
                   std::string const & argumentLabel = "",
                   unsigned numberOfValues = 1) :
        ArgParseArgument(argumentType, isListArgument, argumentLabel, numberOfValues),
        shortName(_shortName),
        longName(_longName),
        _isFlag(false),
        _isRequired(false),
        _isHidden(false)
    {
        _helpText = _help;    
    }

    template <typename TValue>
    ArgParseOption(std::string const & _shortName,
                   std::string const & _longName,
                   std::string const & _help,
                   ArgumentType argumentType,
                   bool isListArgument,
                   std::string const & argumentLabel,
                   unsigned numberOfValues,
                   TValue const & _default) :
        ArgParseArgument(argumentType, isListArgument, argumentLabel, numberOfValues),
        shortName(_shortName),
        longName(_longName),
        _isFlag(false),
        _isRequired(false),
        _isHidden(false)
    {
        _helpText = _help;  

        std::stringstream strm;
        strm << _default;
        appendValue(defaultValue, strm.str());
    }

    ArgParseOption(std::string const & _shortName,
                   std::string const & _longName,
                   std::string const & _help) :
        ArgParseArgument(ArgParseArgument::STRING),
        shortName(_shortName),
        longName(_longName),
        _isFlag(true),
        _isRequired(false),
        _isHidden(false)
    {
        defaultValue.push_back("false");
        setValidValues(*this, "true false");
        _helpText = _help;
    }

};

// ----------------------------------------------------------------------------
// Function isStringArgument()
// ----------------------------------------------------------------------------

///.Function.isStringArgument.param.argument.type:Class.ArgParseOption

inline bool isStringArgument(ArgParseOption const & me)
{
    return isStringArgument(me) && !me._isFlag;
}

// ----------------------------------------------------------------------------
// Function isBooleanOption()
// ----------------------------------------------------------------------------

/**
.Function.isBooleanOption
..summary:Returns whether option is a switch.
..cat:Miscellaneous
..signature:isBooleanOption(option)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..returns:$true$ if the option is a switch.
..include:seqan/arg_parse.h
*/

inline bool isBooleanOption(ArgParseOption const & me)
{
    return me._isFlag;
}

// ----------------------------------------------------------------------------
// Function isVisible()
// ----------------------------------------------------------------------------

/**
.Function.isVisible
..summary:Returns whether option is visible on the help screen. Default is true.
..cat:Miscellaneous
..signature:isHiddenOption(option)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..returns:$true$ if the option is shown on the help screen.
..include:seqan/arg_parse.h
*/

inline bool isVisible(ArgParseOption const & me)
{
    return me._isHidden;
}

// ----------------------------------------------------------------------------
// Function hideOption()
// ----------------------------------------------------------------------------

/**
.Function.hideOption
..summary:Hides the ArgParseOption from the help screen.
..cat:Miscellaneous
..signature:hideOption(option [, hide])
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..param.hide:The new visibility of the option. Default is false.
...type:nolink:bool
..include:seqan/arg_parse.h
*/

inline void hideOption(ArgParseOption & me, bool hide = true)
{
    me._isHidden = hide;
}

// ----------------------------------------------------------------------------
// Function isRequired()
// ----------------------------------------------------------------------------

/**
.Function.isRequired
..summary:Returns whether the option is mandatory.
..cat:Miscellaneous
..signature:isRequired(option)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..returns:$true$ if the option is mandatory.
..include:seqan/arg_parse.h
*/

inline bool isRequired(ArgParseOption const & me)
{
    return me._isRequired;
}

// ----------------------------------------------------------------------------
// Function setRequired()
// ----------------------------------------------------------------------------

/**
.Function.setRequired
..summary:Sets whether or not the option is mandatory.
..cat:Miscellaneous
..signature:setRequired(option, required)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..param.required:The new required value of the option.
...type:nolink:bool
..include:seqan/arg_parse.h
*/

inline void setRequired(ArgParseOption & me, bool required)
{
    me._isRequired = required;
}

// ----------------------------------------------------------------------------
// Function getArgumentLabel()
// ----------------------------------------------------------------------------

///.Function.getArgumentLabel.param.argument.type:Class.ArgParseOption

inline std::string const getArgumentLabel(ArgParseOption const & me)
{
    if(isBooleanOption(me))
        return "";
    else
        return getArgumentLabel(static_cast<ArgParseArgument>(me));
}

// ----------------------------------------------------------------------------
// Helper Function _writeOptName()
// ----------------------------------------------------------------------------

template <typename TStream>
inline void _writeOptName(TStream & target, ArgParseOption const & me)
{
    //IOREV _notio_ irrelevant for iorev
    _streamWrite(target, empty(me.shortName) ? "" : "-");
    _streamWrite(target, me.shortName);
    _streamWrite(target, (empty(me.shortName) || empty(me.longName)) ? "" : ", ");
    if (!empty(me.longName))
    {
        _streamWrite(target, "--");
        _streamWrite(target, me.longName);
    }
}

// ----------------------------------------------------------------------------
// Function write()                                           ArgParseOption
// ----------------------------------------------------------------------------

/**
.Function.write
..summary:Writes the basic information about the @Class.ArgParseOption@ to the provided stream.
..cat:Miscellaneous
..signature:write(stream,option)
..param.stream:The target stream.
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..include:seqan/arg_parse.h
*/

template <typename TStream>
inline void write(TStream & target, ArgParseOption const & me)
{
    //IOREV _nodoc_ this specialization is not documented
    _streamPut(target, '\t');
    _writeOptName(target, me);
    _streamPut(target, '\t');
    _streamPut(target, '\t');
    _streamWrite(target, me._helpText);
}

// ----------------------------------------------------------------------------
// operator<<()                                               ArgParseOption
// ----------------------------------------------------------------------------

template <typename TStream>
inline TStream & operator<<(TStream & target, ArgParseOption const & source)
{
    //IOREV _nodoc_ this specialization is not documented
    write(target, source);
    return target;
}

} // namespace seqan

#endif // SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_OPTION_H_
