UPPER1=`grep strcasecmp ../src/parser.c ../src/wasora/parser.c ../src/wasora/mesh/parser.c | sed -r 's/[^"]*("[^"]*")?/ \1/g;s/" +"/\n/g;s/ *"//g'| mawk '$1 in p{next}{p[$1];print}' | mawk '/^[A-Z]/'`
UPPER2=`grep keywords   ../src/parser.c ../src/wasora/parser.c ../src/wasora/mesh/parser.c | sed -r 's/[^"]*("[^"]*")?/ \1/g;s/" +"/\n/g;s/ *"//g'| mawk '$1 in p{next}{p[$1];print}' | mawk '/^[A-Z]/'`
LOWER=`grep strcasecmp  ../src/parser.c ../src/wasora/parser.c ../src/wasora/mesh/parser.c | sed -r 's/[^"]*("[^"]*")?/ \1/g;s/" +"/\n/g;s/ *"//g'| mawk '$1 in p{next}{p[$1];print}' | mawk '/^[a-z]/' | sort`
VARS=`grep variable     ../src/init.c ../src/wasora/init.c ../src/wasora/mesh/init.c       | grep -v "computing" | sed -r 's/[^"]*("[^"]*")?/ \1/g;s/" +"/\n/g;s/ *"//g'| mawk '$1 in p{next}{p[$1];print}' | sort`
FUNCS=`cat ../src/wasora/builtin.h                                                  | sed -r 's/[^"]*("[^"]*")?/ \1/g;s/" +"/\n/g;s/ *"//g'| mawk '$1 in p{next}{p[$1];print}' | sort`
UPPER=`echo $UPPER1 $UPPER2 | sort`
