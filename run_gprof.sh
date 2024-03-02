URL=https://raw.githubusercontent.com/jrfonseca/gprof2dot/master/gprof2dot.py

FILE=gprof2dot.py
if [ -f $FILE ]; then
   echo ""
else
   echo "File $FILE does not exist. Downloading..."
   wget $URL
   chmod a+x $FILE
fi

gprof $1 | ./gprof2dot.py -ws -n0.01 -e0.01 | dot -Tpng -o "$1"_graph.png
