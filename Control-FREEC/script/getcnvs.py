#!/usr/bin/python

import sys, os
from optparse import OptionParser


def get_sig(CNVs):
    sig_calls = []
    with open(CNVs, 'r') as cnv_in:
        all_count = 0
        sig_count = 0
        for l in cnv_in:
            parts = l.rstrip().split('\t')
            if parts[0] == 'chr':
                sig_calls.append(parts)
                continue
            all_count += 1
            sig1, sig2 = float(parts[5]), float(parts[6])
            if sig1 <= 0.01 or sig2 <= 0.01:
                sig_calls.append(parts)
                sig_count += 1

        filtered_calls = all_count - sig_count
        print("Filtered %s calls" % filtered_calls)
        print("%s significant calls" % sig_count)

    return(sig_calls)


def write_cnvs(sig_calls, cnv_out):
    with open(cnv_out, 'w') as cnv_out_file:
        for l in sig_calls:
            cnv_out_file.write('\t'.join(l) + "\n")


def write_gff(sig_calls, gff_out):
    up, down = "#25CAA2", "#CA254D"
    with open(gff_out, 'w') as cnv_out_gff:
        cnv_out_gff.write("##gff-version 3\n")
        cnv_out_gff.write("#track name=\"$id Freec\" gffTags=on\n")

        for l in sig_calls:
            chrom, start,	end, cn, status, wilco, kolmo = l
            if chrom == 'chr':
                continue

            start = int(start) + 1
            end = int(end) + 1

            if status == 'gain':
                info_field = "Name=" + status + ";" + "copy_no=" + cn + ";" + "colour=" + up + ";" + "wilco_p=" + wilco + ";" + "kolmo_p=" + kolmo + ";"
                vals = [ chrom, 'control-freec', 'DUP', start, end, '.', '+', ".", info_field ]
                gff_fields = '\t'.join(str(v) for v in vals)
                cnv_out_gff.write(gff_fields + "\n")

            elif status == 'loss':
                info_field = "Name=" + status + ";" + "copy_no=" + cn + ";" + "colour=" + down + ";" + "wilco_p=" + wilco + ";" + "kolmo_p=" + kolmo + ";"
                vals = [ chrom, 'control-freec', 'DEL', start, end, '.', '+', ".", info_field ]
                gff_fields = '\t'.join(str(v) for v in vals)
                cnv_out_gff.write(gff_fields + "\n")

            else:
                print("Unknown copy number status: %s" % cn)


def get_args():
    parser = OptionParser()

    parser.add_option("-i", \
                    "--CNV_file", \
                    dest="CNV_file",
                    action="store",
                    help="CNV file with p-values ", \
                        metavar="FILE")

    parser.add_option("-o", \
                      "--out_file", \
                      dest="out_file",
                      action="store",
                      help="Outfile for significant CNVs ", \
                      metavar="FILE")

    parser.add_option("-g", \
                      "--gff_file", \
                      dest="gff_file",
                      action="store",
                      help="Outfile for significant CNVs gff ", \
                      metavar="FILE")

    options, args = parser.parse_args()

    if options.CNV_file is None:
      parser.print_help()
      print

    return(options, args)


def main():
    options, args = get_args()

    if options.CNV_file is not None:
        cnv_in = options.CNV_file
        file_name = (os.path.basename(cnv_in))
        sample = file_name.split('.')[0]
        gff_out = sample + "_" + "control-freec.gff3"
        cnv_out = sample + "_" + "sig_cnvs.txt"

        sig_calls = get_sig(cnv_in)
        write_cnvs(sig_calls, cnv_out)
        write_gff(sig_calls, gff_out)


if __name__ == "__main__":
    sys.exit(main())
