# first karyotype file
import glob

# I run this notebook from the working directory hence I don't need to
# specify that path; instead of path/to/CircosPlot_FgFo/data/genomes/
# I can just say data/genomes/
indir	= '/data/bigexpansion/jcerca/013_ScalesiaGenome/03_circosplots/03_maskedSubgenomes_A_vs_B/00_data/'
outdir  = '/data/bigexpansion/jcerca/013_ScalesiaGenome/03_circosplots/03_maskedSubgenomes_A_vs_B/02_circos/01_karyotype/'
def karyotype_file_from_fasta(fasta_fname, out_fname, color = 'grey'):
    '''
    This function takes a fastafile (fasta_fname) and "translates" this
    to a karyotype file (out_fname). Grey is the default color for chromosomes.
    '''
    outfile = open(out_fname, 'w')
    l      = 0 # counter for the sequence length
    cindex = 1 # counter for the chromosomes, note that this only works
               # because the fastafile is ordered!
    cid    = ''
    for line in open(fasta_fname).readlines():
        if line[0] == '>':
            '''
            This is the fastaheader. We need to store the name in this header
            (we assign it to the variable cid), while we count the length of
            the subsequent lines (the sequence of this chromosome), until we
            encounter the next header
            '''
            if l > 0:
                '''
                If l > 0 it means it is not the first sequence.
                We first have to process the previous sequence before we can move on
                with this one
                '''
                #process previous sequence
                outfile.write('chr - '+cid+' '+str(cindex)+' 0 '+' '+str(l)+' '+color+'\n')
                #make everything ready for this one
                l = 0
                cindex += 1

            '''If you encounter a header, store it in cid.
            '''
            cid = line[1:].strip().split()[0]
            print(cid, cindex)

        else:
            '''
            If it is not header it is part of the sequence, add the length of this line.
            str.strip() takes of any newlines, spaces etc. so we don't count those.
            '''
            l += len(line.strip())

    '''
    When you've reached the end of the file, you still need to process the last header.
    We do this here.
    '''
    if l > 0:
        outfile.write('chr - '+cid+' '+str(cindex)+' 0 '+' '+str(l)+' '+color+'\n')

    outfile.close()

colors = ['green', 'orange']
filenames = glob.glob(indir+'*.fasta')

if len(colors) != len(filenames): #check if your list of colors is long enough
    print("I need more colors, I have ", len(colors), " color(s) for ", len(filenames), " files")
else:
    for i, fasta_fname in enumerate(glob.glob(indir+'*.fasta')):
        print(fasta_fname.split('/')[-1].split('.fasta')[0],colors[i])
        out_fname = outdir + fasta_fname.split('/')[-1].split('.fasta')[0]+'.txt'
        karyotype_file_from_fasta(fasta_fname, out_fname, color = colors[i])
