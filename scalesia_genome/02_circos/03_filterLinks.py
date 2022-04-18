def coords_to_links(coords_fname, links_fname, min_length = 0):
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
            linksfile.write(chrA+' '+fromA+' '+toA+' '+chrB+' '+fromB+' '+toB+'\n')
    linksfile.close()

coords_fname = '/data/bigexpansion/jcerca/013_ScalesiaGenome/03_circosplots/03_maskedSubgenomes_A_vs_B/01_alignment/Scalesia_maskedsubgenomes.promer_def.coords'
min_length   = 1550
links_fname  = '/data/bigexpansion/jcerca/013_ScalesiaGenome/03_circosplots/03_maskedSubgenomes_A_vs_B/02_circos/03_Filteredlinks/filteredLinks_ScalA_vs_ScalB.min_length'+str(min_length)+'bp.txt'
coords_to_links(coords_fname, links_fname, min_length = min_length)
