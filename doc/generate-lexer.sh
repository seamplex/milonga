#!/bin/sh
# 
# remember to tell pygmentize it has a new lexer
#
# $ ./generate-lexer.sh > lexer/milongalexer/lexer.py 
# $ cd lexer
# $ python setup.py install --user
#

. ./keywords.sh 

cat << EOF 
import re

from pygments.lexer import RegexLexer, include, bygroups, using, \\
     this, combined, ExtendedRegexLexer
from pygments.token import Error, Punctuation, Literal, Token, \\
     Text, Comment, Operator, Keyword, Name, String, Number, Generic

__all__ = ['MilongaLexer']

class MilongaLexer(RegexLexer):
    """
    Lexer for milonga input files.

    *Not in pygments official distribution
    """
    name = 'milonga'
    aliases = ['milonga']
    filenames = ['*.mil']
    mimetypes = ['text/plain']

    tokens = {
        'root': [
EOF

echo "            # Keywords"
echo -n "             (r'\\\\b("
echo -n $UPPER | sed s/\ /\|/g
echo ")\s*\\\\b',"
echo "             Keyword),"
echo             
echo "            # Reserved variables"
echo -n "             (r'\\\\b("
echo -n $VARS | sed s/\ /\|/g
echo ")\s*\\\\b',"
echo "             Name.Builtin),"
echo
echo "            # Built-in functions, functionals and vectorfunctions"
echo -n "             (r'\\\\b("
echo -n $FUNCS | sed s/\ /\|/g
echo ")\s*\\\\b',"
echo "             Name.Function),"

cat << EOF

            # Numbers
            (r'([0-9]*\.[0-9]*|[0-9]+)(e[+-]?[0-9]+)?', Number),
            (r' .[0-9]+', Number),

            # Operators
            (r'[=!+-/\*\^]', Operator),

            # Punctuation
            (r'[&\|\(\)\[\]\{\}<>@\\\\$,:;%\"\']', Punctuation),

            # Comment (but not if escaped)
            (r'^#.*?\n', Comment),
            (r'\#', Text),
            (r'\\\#', Text),
            (r'\\\\\#', Text),
            (r'.#.*?\n', Comment),

            (r'[a-z][a-z0-9_]*', Text),
            (r'[A-Z][A-Z0-9_]*', Text),
            (r' ', Text),
            (r' .*', Text),

        ],

    }

EOF