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

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_DOC_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_DOC_H_

#include <seqan/arg_parse/argument_parser.h>

#include <seqan/misc/tool_doc.h>

namespace seqan {

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function getAppName()
// --------------------------------------------------------------------------

/**
.Function.getAppName
..summary:Get tool name of @Class.ArgumentParser@ object.
..cat:Miscalleneous
..signature:getAppName(parser)
..param.doc:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..returns:Tool name of argument parser object.
...type:std::string
..include:seqan/arg_parse.h
*/

inline CharString const & getAppName(ArgumentParser const & parser)
{
    return getName(parser._toolDoc);
}

// ----------------------------------------------------------------------------
// Helper Function _parseAppName()
// ----------------------------------------------------------------------------

inline void _parseAppName(ArgumentParser & parser, std::string const & candidate)
{
    //IOREV _notio_ irrelevant for io-revision
    int i = length(candidate) - 1;

    for (; i >= 0; --i)
        if (candidate[i] == '\\' || candidate[i] == '/')
            break;

    setName(parser._toolDoc, candidate.substr(i + 1));
}

// ----------------------------------------------------------------------------
// Helper Function _addLine()
// ----------------------------------------------------------------------------

/**
.Function.addLine:
..summary:Adds a line of text to the help output of the @Class.ArgumentParser@ in the block of
@Class.ArgParseOption@s.
..cat:Miscellaneous
..signature:addLine(parser, text)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.text:A line of text that will be added to the help output.
...type:Shortcut.CharString
..include:seqan/arg_parse.h
*/

template <typename TString>
inline void addLine(ArgumentParser & me, TString const & line)
{
    addOption(me, ArgParseOption("", "", line));
}

// ----------------------------------------------------------------------------
// Function addSection()
// ----------------------------------------------------------------------------

/**
.Function.addSection:
..summary:Begins a new section of @Class.ArgParseOption@ the help output of
the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addSection(parser, text)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.text:A section header that will be added to the help output.
...type:Shortcut.CharString
..include:seqan/arg_parse.h
..example.code:
ArgumentParser parser;

[...] // init parser

addSection(parser, "In-/Output-Options");
addOption("i", ... );
addOption("o", ... );

addSection(parser, "Other Options");
addOption("x", ... );
*/

template <typename TString>
inline void addSection(ArgumentParser & me, TString const & line)
{
    addLine(me, "");
    addLine(me, line);
}

// ----------------------------------------------------------------------------
// Function addUsageLine()
// ----------------------------------------------------------------------------

/**
.Function.addUsageLine:
..summary:Adds a line of text to the usage output of the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addUsageLine(parser, text)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.text:A text line that will be added to the usage output.
..include:seqan/arg_parse.h
*/

inline void addUsageLine(ArgumentParser & me, std::string const & line)
{
    me._usageText.push_back(line);
}

// ----------------------------------------------------------------------------
// Helper Function _addUsage()
// ----------------------------------------------------------------------------

void _addUsage(ToolDoc & toolDoc, ArgumentParser const & me)
{
    for (unsigned i = 0; i < length(me._usageText); ++i)
    {
        std::string text = "\\fB";
        append(text, getAppName(me));
        append(text, "\\fP ");
        append(text, me._usageText[i]);
        addText(toolDoc, text, false);
    }
}

// ----------------------------------------------------------------------------
// Function addDescription()
// ----------------------------------------------------------------------------

/**
.Function.addDescription
..summary:Appends a description paragraph to the @Class.ArgumentParser@ documentation.
..cat:Miscellaneous
..signature:addDescription(parser, text)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.text:The description paragraph.
..returns:$void$
..include:seqan/arg_parse.h
*/

inline void addDescription(ArgumentParser & me, std::string const & description)
{
    addText(me._description, description);
}

// ----------------------------------------------------------------------------
// Function setAppName()
// ----------------------------------------------------------------------------

/**
.Function.setAppName
..summary:Sets application name of @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:setAppName(parser, appName)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.appName:The name of the application.
..returns:$void$
..include:seqan/arg_parse.h
*/

inline void setAppName(ArgumentParser & me, std::string const & name)
{
    setName(me._toolDoc, name);
}

// ----------------------------------------------------------------------------
// Function setShortDescription()
// ----------------------------------------------------------------------------

/**
.Function.setShortDescription
..summary:Sets short description test of @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:setShortDescription(parser, text)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.text:The short description text.
..returns:$void$
..include:seqan/arg_parse.h
*/

inline void setShortDescription(ArgumentParser & me, std::string const & description)
{
    setShortDescription(me._toolDoc, description);
}

// ----------------------------------------------------------------------------
// Function setVersion()
// ----------------------------------------------------------------------------

/**
.Function.setVersion
..summary:Sets version string of @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:setVersion(parser, versionString)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.versionString:The version string to set.
..returns:$void$
..include:seqan/arg_parse.h
*/

inline void setVersion(ArgumentParser & me, std::string const & versionString)
{
    setVersion(me._toolDoc, versionString);
    if(!hasOption(me, "version")) addOption(me, ArgParseOption("", "version", "Display version information"));
}

// --------------------------------------------------------------------------
// Function getVersion()                                              ToolDoc
// --------------------------------------------------------------------------

/**
.Function.getVersion
..cat:Miscalleneous
..summary:Get version string from @Class.ArgumentParser@ object.
..signature:getVersion(parser)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..returns:Date string.
...type:Shortcut.CharString
..include:seqan/misc/tool_doc.h
*/

inline CharString const & getVersion(ArgumentParser const & me)
{
    return getVersion(me._toolDoc);
}

// ----------------------------------------------------------------------------
// Function setDate()
// ----------------------------------------------------------------------------

/**
.Function.setDate
..summary:Sets date string of @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:setDate(parser, date)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.date:The date string.
..returns:$void$
..include:seqan/arg_parse.h
*/

inline void setDate(ArgumentParser & me, std::string const & date)
{
    setDate(me._toolDoc, date);
}

// ----------------------------------------------------------------------------
// Function addTextSection()
// ----------------------------------------------------------------------------

/**
.Function.addTextSection
..summary:Adds a text section to the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addTextSection(parser, title)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.title:The section title.
..returns:$void$
..remarks:This will result in an additional section heading to be printed.
..include:seqan/arg_parse.h
*/

inline void addTextSection(ArgumentParser & me, std::string const & title)
{
    addSection(me._toolDoc, title);
}

// ----------------------------------------------------------------------------
// Function addTextSubSection()
// ----------------------------------------------------------------------------

/**
.Function.addTextSubSection
..summary:Adds a text subsection to the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addTextSubSection(parser, title)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.title:The subsection title.
..returns:$void$
..remarks:This will result in an additional subsection heading to be printed.
..include:seqan/arg_parse.h
*/

inline void addTextSubSection(ArgumentParser & me, std::string const & title)
{
    addSubSection(me._toolDoc, title);
}

// ----------------------------------------------------------------------------
// Function addText()
// ----------------------------------------------------------------------------

/**
.Function.addText
..summary:Appends a text paragraph to the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addText(parser, text)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.text:The content of the text.
..returns:$void$
..include:seqan/arg_parse.h
*/

inline void addText(ArgumentParser & me, std::string const & text)
{
    addText(me._toolDoc, text);
}

// ----------------------------------------------------------------------------
// Function addListItem()
// ----------------------------------------------------------------------------

/**
.Function.addListItem
..summary:Appends a list item to the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addListItem(parser, item, description)
..description:
This method adds a list item to the parser's output.
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.item:The item text.
..param.description:The description text.
..returns:$void$
..include:seqan/arg_parse.h
*/

inline void addListItem(ArgumentParser & me, std::string const & item, std::string const & description)
{
    addListItem(me._toolDoc, item, description);
}

// ----------------------------------------------------------------------------
// Function printShortHelp()
// ----------------------------------------------------------------------------

/**
.Function.printShortHelp
..summary:Prints a short help message for the parser to a stream
..cat:Miscellaneous
..signature:printShortHelp(parser[, stream])
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.stream:Target stream (e.g. $std::cerr$).
..include:seqan/arg_parse.h
*/

inline void printShortHelp(ArgumentParser const & me, std::ostream & stream)
{
    // TODO: maybe we can get this a bit prettier
    ToolDoc shortDoc(me._toolDoc);
    clearEntries(shortDoc);

    _addUsage(shortDoc, me);

    std::stringstream shortHelp;
    shortHelp << "Try '" << getAppName(me) << " --help' for more information.\n";
    addText(shortDoc, shortHelp.str());

    print(stream, shortDoc, "txt");
}

// ----------------------------------------------------------------------------
// Function printVersion()
// ----------------------------------------------------------------------------

/**
.Function.printVersion
..summary:Prints the version information of the parser to a stream.
..cat:Miscellaneous
..signature:printVersion(parser[, stream])
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.stream:Target std::ostream (e.g. $std::cerr$).
...default: $std::cerr$
..include:seqan/arg_parse.h
*/

inline void printVersion(ArgumentParser const & me, std::ostream & stream)
{
    stream << getAppName(me) << " version " << getVersion(me) << std::endl;
}

inline void printVersion(ArgumentParser const & me)
{
    printVersion(me, std::cerr);
}

// ----------------------------------------------------------------------------
// Function printHelp()
// ----------------------------------------------------------------------------

/**
.Function.printHelp
..summary:Prints the complete help message for the parser to a stream.
..cat:Miscellaneous
..signature:printHelp(parser[, stream][, format])
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.stream:Target std::ostream (e.g. $std::cerr$).
...default: $std::cerr$
..param.format:Format to print, one of "html", "man", "txt".
..include:seqan/arg_parse.h
*/

inline void printHelp(ArgumentParser const & me, std::ostream & stream, CharString const & format)
{
    ToolDoc toolDoc(me._toolDoc);
    clearEntries(toolDoc);  // We will append me._toolDoc later.

    // Build synopsis section.
    addSection(toolDoc, "Synopsis");
    _addUsage(toolDoc, me);

    // Add description to tool documentation.
    append(toolDoc, me._description);

    // Add options to description section.
    for (unsigned i = 0; i < length(me.optionMap); ++i)
    {
        ArgParseOption const & opt = me.optionMap[i];
        if (empty(opt.shortName) && empty(opt.longName))  // this is not an option but a text line
        {
            if (empty(opt._helpText))  // TODO(holtgrew): Should go away in future.
                continue;  // Skip empty lines.

            // Is command line parser section, maps to ToolDoc subsection.
            std::string title = opt._helpText;
            append(title, ":");
            addSubSection(toolDoc, title);
        }
        else if (!isVisible(opt))
        {
            // Build list item term.
            std::string term;
            if (!empty(opt.shortName))
            {
                term = "\\fB-";
                append(term, opt.shortName);
                append(term, "\\fP");
            }
            if (!empty(opt.shortName) && !empty(opt.longName))
                append(term, ", ");
            if (!empty(opt.longName))
            {
                append(term, "\\fB--");
                append(term, opt.longName);
                append(term, "\\fP");
            }
            // Get arguments, autogenerate if necessary.
            std::string arguments = getArgumentLabel(opt);

            // Write arguments to term line -> only exception, boolean flags
            if (!empty(arguments))
            {
                // Tokenize argument names.
                std::istringstream iss(toCString(arguments));
                std::vector<std::string> tokens;
                std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
                          std::back_inserter<std::vector<std::string> >(tokens));
                // Append them, formatted in italic.
                for (unsigned i = 0; i < length(tokens); ++i)
                {
                    append(term, " \\fI");
                    append(term, tokens[i]);
                    append(term, "\\fP");
                }
            }

            // Add list item.
            addListItem(toolDoc, term, opt._helpText);
        }
    }

    append(toolDoc, me._toolDoc);
    print(stream, toolDoc, format);
}

inline void printHelp(ArgumentParser const & me, std::ostream & stream)
{
    printHelp(me, stream, "txt");
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_DOC_H_
