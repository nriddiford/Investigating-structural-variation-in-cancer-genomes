#!/usr/bin/env python

import sys
import os
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-f",
    "--cnv_file",
    dest="cnv_in",
    help="CNV calls")

parser.add_option("-o",
    "--out_file",
    dest="out_file",
    help="Output")

(options, args) = parser.parse_args()

if not options.cnv_in:
    parser.error('No input file provided')

f = open(options.cnv_in,'r')

if options.out_file:
    out = open(options.out_file, 'w')
    print "Writing filtered calls to '%s'" % options.out_file

else:
    base_name = (os.path.splitext(options.cnv_in)[0])
    out_base = base_name.split('_')[0]
    window = base_name.split('_')[1]
    outfile = "_".join([out_base, window, "filtered.txt"])
    out = open(outfile,'w')
    print "Writing filtered calls to '%s'" % outfile


filtered_count = 0
cnv_count = 0

i = 1
for l in f:
    parts = l.rstrip().split()

    if i == 1 and parts[0] == 'cnv': continue

    cnv, chromosome, start, end, size, log2, pvalue = parts

    size = int(size)
    log2 = abs(float(log2))

    cnv_out = '\t'.join(parts) + '\n'

    if size > 500 and size < 1000 and log2 >= 2:
        cnv_count += 1
        out.write(cnv_out)

    if size >= 1000 and log2 >= 0.8:
        cnv_count += 1
        out.write(cnv_out)


    elif size > 10000 and abs(float(log2)) >= 0.53:
        print "Big cnv: %s" % parts
        cnv_count += 1
        out.write(cnv_out)

    else:
        filtered_count += 1

print "%s filtered, %s kept from %s" % (filtered_count, cnv_count, out_base)

f.close()
