#!/usr/bin/env python
import sys, re
import pysam
import argparse

def complement(base):

    return {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'N': 'N',
        'D': 'D',
        '.': '.',
    }[base]

def reverse_complement(seq):
    
    return ''.join([complement(base) for base in reversed(seq)])

def main_to_mr(args):
    
    samfile = pysam.AlignmentFile(args.i, "rb")

    qname2reads = {}
    for x in samfile:
        if x.is_unmapped:
            continue
        if x.is_qcfail:
            continue
        if x.is_duplicate:
            continue
        if x.is_secondary:
            continue
        
        if x.query_name not in qname2reads:
            qname2reads[x.query_name] = x
            continue
        y = qname2reads[x.query_name]
        if x.is_read1 and y.is_read2:
            r1 = x
            r2 = y
        elif x.is_read2 and y.is_read1:
            r1 = y
            r2 = x
        else:
            sys.stderr.write('multiple mapping detected for %s, skip' % x.query_name)
            continue

        qname2reads.pop(x.query_name, None)
        if r1.tid != r2.tid:
            continue

        md1 = r1.get_tag('MD') if r1.has_tag('MD') else None
        md2 = r2.get_tag('MD') if r2.has_tag('MD') else None
        xm1 = r1.get_tag('XM') if r1.has_tag('XM') else None
        xm2 = r2.get_tag('XM') if r2.has_tag('XM') else None
        bs1 = r1.get_tag('XG') if r1.has_tag('XG') else None
        bs2 = r2.get_tag('XG') if r2.has_tag('XG') else None
        nm1 = r1.get_tag('NM') if r1.has_tag('NM') else None
        nm2 = r2.get_tag('NM') if r2.has_tag('NM') else None

        # get sequence and qual
        s1 = ''
        n1 = ''
        q1 = ''
        qpos = 0
        for ct, cl in r1.cigar:
            if ct == 0:     # match
                s1 += r1.seq[qpos:qpos+cl]
                q1 += r1.qual[qpos:qpos+cl]
            elif ct == 1:   # insertion
                qpos += cl
            elif ct == 2:   # deletion
                s1 += 'N'*cl
                q1 += 'B'*cl
            elif ct == 4:
                qpos += cl

        s2 = ''
        n2 = ''
        q2 = ''
        qpos = 0
        for ct, cl in r2.cigar:
            if ct == 0:     # match
                s2 += r2.seq[qpos:qpos+cl]
                q2 += r2.qual[qpos:qpos+cl]
            elif ct == 1:   # insertion
                qpos += cl
            elif ct == 2:   # deletion
                s2 += 'N'*cl
                q2 += 'B'*cl
            elif ct == 4:
                qpos += cl

        # get mismatch string, for bismark
        if (md1 is not None and md2 is not None and xm1 is not None and
            xm2 is not None and bs1 is not None and bs2 is not None):
            qpos = 0
            for m in re.finditer(r'(\d+)([\^ATCGN]+)', md1):
                skip = int(m.group(1))
                qpos += skip
                n1 += '.'*skip
                c = m.group(2)
                if c[0] == '^':
                    n1 += 'D'+'.'*(len(c)-2)    # only record first base as deletion
                elif bs1 == 'CT' and c == 'C' and s1[qpos] == 'T':
                    n1 += '.'
                    qpos += len(c)
                elif bs1 == 'GA' and c == 'G' and s1[qpos] == 'A':
                    n1 += '.'
                    qpos += len(c)
                else:
                    n1 += c
                    qpos += len(c)

            m = re.match(r'.*?(\d+)$', md1)
            n1 += '.'*int(m.group(1))

            qpos = 0
            for m in re.finditer(r'(\d+)([\^ATCGN]+)', md2):
                skip = int(m.group(1))
                qpos += skip
                n2 += '.'*skip
                c = m.group(2)
                if c[0] == '^':
                    n2 += 'D'+'.'*(len(c)-2)    # only record first base as deletion
                elif bs2 == 'CT' and c == 'C' and s2[qpos] == 'T':
                    n2 += '.'
                    qpos += len(c)
                elif bs2 == 'GA' and c == 'G' and s2[qpos] == 'A':
                    n2 += '.'
                    qpos += len(c)
                else:
                    n2 += c
                    qpos += len(c)

            m = re.match(r'.*?(\d+)$', md2)
            n2 += '.'*int(m.group(1))

        # merge paired reads
        if r1.is_reverse:
            strand = '-'
            rbeg = r2.reference_start
            rend = r1.reference_end
            rlen = rend - rbeg
            if rlen > args.maxrlen or rlen < args.k: # skip improper pair and SV jump
                continue
            s = reverse_complement(s1)[:rlen]
            q = ''.join(reversed(q1))[:rlen]
            n = reverse_complement(n1)[:rlen]
            gap = r1.reference_start - r2.reference_end
            if gap > 0:
                s += 'N'*gap
                s += reverse_complement(s2)
                q += 'B'*gap
                q += ''.join(reversed(q2))
                if len(n) > 0: n += 'N'*gap
                n += reverse_complement(n2)
            else:
                s += reverse_complement(s2)[-gap:]
                n += reverse_complement(n2)[-gap:]
                q += ''.join(reversed(q2))[-gap:]
        else:
            strand = '+'
            rbeg = r1.reference_start
            rend = r2.reference_end
            rlen = rend - rbeg
            if rlen > args.maxrlen or rlen < args.k: # skip improper pair and SV jump
                continue
            s = s1[:rlen]
            q = q1[:rlen]
            n = n1[:rlen]
            gap = r2.reference_start - r1.reference_end
            if gap > 0:
                s += 'N'*gap
                s += s2
                q += 'B'*gap
                q += q2
                if len(n) > 0: n += 'N'*gap
                n += n2
            else:
                s += s2[-gap:]
                q += q2[-gap:]
                n += n2[-gap:]

        if len(n) > 0:
            nm = len(n)-n.count('.')-n.count('N')
        elif nm1 is not None and nm2 is not None:
            nm = nm1 + nm2
            
        args.o.write('%s\t%d\t%d\tFRAG:%s\t%d\t%s\t%s\t%s\n' % (samfile.getrname(r1.tid), rbeg, rend, r1.query_name, nm, strand, s, q))
        if args.v > 0:
            print
            print bs1, bs2
            print n
            print s
            print len(s)
            print len(n)
            print rend - rbeg
            print r1
            print r2

        assert(rend - rbeg == len(s))
        assert(len(s) == len(n) or len(n) == 0)

def main_to_methylKit(args):

    out = open(args.o,"w") if args.o is not None else sys.stdout
    out.write("chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n")
    for line in args.i:
        fields = line.strip().split('\t')
        strand = "F" if fields[5]=="C" else "R"
        out.write("%s.%s\t%s\t%s\t%s\t%d\t%1.2f\t%1.2f\n" % 
                  (fields[0], fields[2], fields[0], fields[2], strand,
                   int(fields[4]), float(fields[3])*100, (1-float(fields[3]))*100))
    
    return

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Python scripts for Biscuits')
    subparsers = parser.add_subparsers()
    
    parser_to_mr = subparsers.add_parser('to_mr', help='convert bam to mr file for methpipe (automatically determine bam input type and parse indels and clippings)')
    parser_to_mr.add_argument('-i', required=True, help='input bam')
    parser_to_mr.add_argument('-o', type=argparse.FileType('w'), help='output', default=sys.stdout)
    parser_to_mr.add_argument('-v', type=int, default=0, help='verbose level')
    parser_to_mr.add_argument('-l', '--maxrlen', type=int, default=1000, help='maximum template length [1000]')
    parser_to_mr.add_argument('-k', type=int, default=40, help='minimum template length [40]')
    parser_to_mr.add_argument('-m', type=int, default=10000000, help='number of total unpaired reads cached in memory before flushing')
    parser_to_mr.set_defaults(func=main_to_mr)

    
    parser_to_methylKit = subparsers.add_parser('to_methylKit', help='convert biscuit vcf2bed output to methylKit format')
    parser_to_methylKit.add_argument('-i', type=argparse.FileType('r'), default='-', help='input table')
    parser_to_methylKit.add_argument('-o', help='output', default=None)
    parser_to_methylKit.set_defaults(func=main_to_methylKit)
    
    args = parser.parse_args()
    try:
        args.func(args)
    except IOError as e:
        sys.exit()
