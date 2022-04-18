def coords_to_links(coords_fname, links_fname):
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

        linksfile.write(chrA+' '+fromA+' '+toA+' '+chrB+' '+fromB+' '+toB+'\n')
    linksfile.close()
    print(coords_fname,": Done")

coords_fname ='/data/bigexpansion/jcerca/013_ScalesiaGenome/03_circosplots/03_maskedSubgenomes_A_vs_B/01_alignment/Scalesia_maskedsubgenomes.promer_def.coords'
links_fname  = '/data/bigexpansion/jcerca/013_ScalesiaGenome/03_circosplots/03_maskedSubgenomes_A_vs_B/02_circos/02_links/masked_ScalA_vs_B.links'
coords_to_links(coords_fname, links_fname)
