import glob

def coords_to_links(coords_fname, links_fname, chr2color = None, min_length = 0):
    linksfile = open(links_fname, 'w')
    lines     = open(coords_fname).readlines()[5:] #skip header
    for line in lines:
        columns = line.strip().split()

        if len(columns) != 17:
            print('Unexpected fileformat, can not continue!')
            break

        chrA  = columns[15]
        chrB  = columns[16]
        fromA = columns[0]
        toA   = columns[1]
        fromB = columns[3]
        toB   = columns[4]

        lenA  = int(columns[6])
        lenB  = int(columns[7])
        if lenA >= min_length and lenB >= min_length:
            out = chrA+' '+fromA+' '+toA+' '+chrB+' '+fromB+' '+toB

            if chr2color != None:
                if chrA in chr2color:
                    linksfile.write(out+' color='+chr2color[chrA]+'\n')
                elif chrB in chr2color:
                    linksfile.write(out+' color='+chr2color[chrB]+'\n')
            else:
                linksfile.write(out+'\n')
    linksfile.close()


def karyotype_file_from_fasta(fasta_fname, out_fname, chr2color = None):

    outfile = open(out_fname, 'w')
    l      = 0 #counter for the sequence length
    cindex = 1 #counter for the chromosomes
    cid    = ''
    color  = 'grey'
    for line in open(fasta_fname).readlines():
        if line[0] == '>':
            if l > 0:
                outfile.write('chr - '+cid+' '+str(cindex)+' 0 '+' '+str(l)+' '+color+'\n')
                l = 0
                cindex += 1
            cid = line[1:].strip().split()[0]
            if chr2color != None and cid in chr2color:
                color = chr2color[cid]
            print(cid, cindex)

        else:
            l += len(line.strip())

    #last one:
    if l > 0:
        outfile.write('chr - '+cid+' '+str(cindex)+' 0 '+' '+str(l)+' '+color+'\n')

    outfile.close()

# write karyotype files from fasta
chr2color = {}
chr2color['ScDrF4C_6'] = 'lblue'
chr2color['ScDrF4C_1633'] = 'blue'
chr2color['ScDrF4C_1634'] = 'lred'
chr2color['ScDrF4C_9'] = 'red'
chr2color['ScDrF4C_4'] = 'lyellow'
chr2color['ScDrF4C_25'] = 'yellow'
chr2color['ScDrF4C_21'] = 'lpurple'
chr2color['ScDrF4C_19'] = 'purple'
chr2color['ScDrF4C_10'] = 'lgreen'
chr2color['ScDrF4C_17'] = 'green'
chr2color['ScDrF4C_14'] = 'lgrey'
chr2color['ScDrF4C_24'] = 'grey'
chr2color['ScDrF4C_30'] = 'lorange'
chr2color['ScDrF4C_5'] = 'orange'
chr2color['ScDrF4C_18'] = 'black'
chr2color['ScDrF4C_15'] = 'lgrey'
chr2color['ScDrF4C_16'] = 'grey'

for i, fasta_fname in enumerate(glob.glob('/data/bigexpansion/jcerca/013_ScalesiaGenome/03_circosplots/03_maskedSubgenomes_A_vs_B/00_data/*fasta')):

    out_fname = '/data/bigexpansion/jcerca/013_ScalesiaGenome/03_circosplots/03_maskedSubgenomes_A_vs_B/02_circos/04_coloring/' + fasta_fname.split('/')[-1].split('.fasta')[0]+'.color-coded.txt'
    karyotype_file_from_fasta(fasta_fname, out_fname, chr2color = chr2color)

coords_fname = '/data/bigexpansion/jcerca/013_ScalesiaGenome/03_circosplots/03_maskedSubgenomes_A_vs_B/01_alignment/Scalesia_maskedsubgenomes.promer_def.coords'
links_fname  = '/data/bigexpansion/jcerca/013_ScalesiaGenome/03_circosplots/03_maskedSubgenomes_A_vs_B/02_circos/04_coloring/filteredLinks_ColorCoded_ScalA_vs_ScalB.colour.txt'
coords_to_links(coords_fname, links_fname, chr2color = chr2color, min_length = 1500)
