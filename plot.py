#! /usr/bin/env python
import sys
import codecs
import re
import collections
from pylab import *
print "Hello World"

def func(a, filt):
	if filt == 1:#median filter
		b = list(a);
		b.sort();
		ret = b[(len(b)-1)//2];
	elif filt == 2:#moving average (simple)
		L = len(a);
		ret = 0;
		for i in a:
			ret = ret+i;
		ret = ret/L;
	return ret;

if len(sys.argv) != 2:
	print "Wrong number of arguments";
	print "arg1: inputfile"
	raise SystemExit(1);

finput = open(sys.argv[1]);


outbuf = [];
i = 0;
for line in finput:
	outbuf.append(float(line));
	
	
Y = array(outbuf);
X = arange(len(outbuf));

plot(X,Y);

show();

























