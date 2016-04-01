#!/usr/bin/env python
"""
The MIT License

Copyright (c) 2015,2016 by Van Andel Research Institute

Wanding Zhou (wanding.zhou@vai.org)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""
import os, sys, re

from subprocess import Popen, PIPE, check_call

def opengz(fn,m='r'):
    
    if fn.endswith('.gz'):
        import gzip
        fh = gzip.open(fn,m)
    else:
        fh = open(fn,m)

    return fh

def wrap(line):

    lw = 80
    # fields = line.strip().split('\t')
    # yield '\t'.join(fields[:4])
    # yield ' '*3+'\t'.join(fields[4:6])
    # line = '\t'.join(fields[6:])

    fields = line.strip().split('\t')[::-1]
    first = True
    outfields = []
    while fields:
        outfields.append(fields.pop())
        if len(''.join(outfields)) + 2*(len(outfields)-1) > 80 or len(fields) == 0:
            yield '  '.join(outfields)
            outfields = ['  ']
    
if len(sys.argv) == 1:
    indir = 'biscuit.wiki'
    outdir = 'outwiki'
else:
    indir = sys.argv[1]
    outdir = sys.argv[2]

if os.path.isdir(indir):
    ifns = os.listdir(indir)
    fromdir = True
    if not os.path.exists(outdir):
        os.mkdir(outdir)
else:
    ifns = [indir]
    fromdir = False

def calloutput(line, capture=True):

    line = line.strip('$').strip()
    ar = re.split(r'[\'"]', line.strip())

    B = []
    for i, e in enumerate(ar):
        e = e.strip()
            
        if i % 2 == 0:
            if len(e)>0:
                B.extend(e.split())
        else:
            B.append(e)

    if B[0] == 'biscuit':
        B[0] = 'bin/biscuit'

    p = Popen(B, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    out = out.strip()
    err = err.strip()
    if p.returncode != 0 and capture:
        print err
        raise Exception('Run error!')

    return out, err, p.returncode

def calloutput_out(line, capture=True):

    out,err,code = calloutput(line, capture=capture)
    return out.split('\n')

def calloutput_err(line, capture=True):

    out,err,code = calloutput(line, capture=capture)
    return err.split('\n')

prog='biscuit-develop'

# get version
version = ''
for line in calloutput_err(prog, capture=False):
    m = re.match(r'Version: ([\d\.]*)', line)
    if m:
        version = m.group(1)

raw_input("confirm version: %s" % version)

def linesub(line, sub):

    fields = []
    for f in line.split():
        if f in sub:
            fields.append(sub[f])
        elif f[0] == '>' and f[1:] in sub:
            fields.append('>'+sub[f[1:]])
        else:
            fields.append(f)

    return ' '.join(fields)

for ifn in ifns:
    if not ifn.endswith(".md"):
        continue

    if fromdir:
        ifh = open(os.path.join(indir, ifn))
        ofn = os.path.join(outdir, ifn)
    else:
        ifh = open(ifn)
        ofn = outdir
        
    ofh = open(ofn, "w")
    
    result = ''
    tofill = False
    infill = False
    sub = {}
    for line in ifh:

        # line = line.replace('@VERSION', version)
        line = line.rstrip()
        line2 = linesub(line, sub)
        if line.startswith('$$$ '): # non-capture run without print
            calloutput(line2)
            ofh.write(line+'\n')

        elif line.startswith('$$ '): # non-capture run
            print
            print '======'+line2+'====='
            calloutput(line2)
            ofh.write(line+'\n')

        elif line.startswith('$$+ '): # non-capture run with pipe inside, use shell
            print
            print '======'+line2+'====='
            check_call(line2[4:], shell=True)
            ofh.write(line+'\n')

        elif line.startswith('@@ '):       # udpate substitution
            pairs = line[3:].split(':')
            sub[pairs[0].strip()] = pairs[1].strip()
            # sub1 = eval(line[3:])
            # for k,v in sub1.iteritems():
            #     sub[k] = v
            ofh.write(line+'\n')

        elif line.startswith('##compare'): # compare non-capture runs
            m = re.match('##compare (\S*) vs (\S*)', line)
            if not m:
                raise Exception('Invalid compare')
            sys.stdout.write('compare %s vs %s ... ' % (m.group(1), m.group(2)))
            if not os.path.exists(m.group(1)):
                raise Exception("file %s nonexistent" % m.group(1))
            if not os.path.exists(m.group(2)):
                cc = raw_input("file %s nonexistent, first-time run? going to create (y or n):" % m.group(2))
                if cc == 'y': 
                    check_call(['mkdir','-p',os.path.dirname(m.group(2))])
                    check_call(['cp',m.group(1),m.group(2)])
                else:
                    sys.exit(1)
                
            res, err, code = calloutput('zdiff %s %s' % (m.group(1), m.group(2)), capture=False)
            if code not in [0,1]:
                print err
                raise Exception('Run error!')
            if ''.join(res):
                raw_input('different')
            else:
                print('same')

            ofh.write(line+'\n')

        elif line.startswith('$ '):    # capture run

            ofh.write(line+'\n')
            print
            print '======'+line2+'======'
            result, err, code = calloutput(line2)
            print '\n'.join([rr for rr in result.split('\n') if not rr.startswith('#')])
            tofill = True
            infill = False
            oldfill = ''

        elif line.startswith('$+ '):    # pure file display

            ofh.write(line+'\n')
            result = ''
            for i,l in enumerate(opengz(line2[3:])):
                if i>5:
                    break
                result += l
            tofill = True
            infill = False
            oldfill = ''

        elif line.startswith('```text') and tofill: # load old output in a capture run
            infill = True
            tofill = False
            ofh.write(line+'\n')

        elif line.startswith('```') and infill: # compare capture run

            newfill = ''
            for rr in result.split('\n'):
                if (len(rr.strip()) == 0
                    or rr.startswith('#')):
                    continue

                for r in wrap(rr):
                    newfill += r+'\n'

            if newfill != oldfill:
                print '\n+++ OLD +++'
                print oldfill
                print '\n+++ NEW +++'
                print newfill
                raw_input("Difference ...")
            else:
                print "\nSame!\n"

            ofh.write(newfill)
            ofh.write('```\n')

            infill = False

        elif infill:
            oldfill += line

        else:

            ofh.write(line+'\n')

    ifh.close()
    ofh.close()
