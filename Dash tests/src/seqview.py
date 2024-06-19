def basicView(seq):
    from src.utils import find_orf,translate
    import matplotlib.pyplot as plt
    
    #break seq into 60 nt lines
    ntlines = [seq[n:n+60] for n in range(0,len(seq),60)]
    ntlines.reverse()
    #find the orf coordinates
    orf_coords = find_orf(seq)
    #translte the orf
    aa_seq = translate(seq[orf_coords[0]:orf_coords[1]])
    #pad the aa sequence with zeros to fill the UTRs, and add double zeros
    #between aas to align with nt seq
    aa_seq_spaced = ' ' * orf_coords[0] + ''.join([aa+'  ' for aa in aa_seq]) +' ' * (len(seq)-orf_coords[1])
    aalines = [aa_seq_spaced[n:n+60] for n in range(0,len(aa_seq_spaced),60)]
    aalines.reverse()
    #plot both lines
    fig,ax = plt.subplots(figsize=(8,len(ntlines)*0.6))
    for lineno, ntline in enumerate(ntlines):
        ax.text(0,lineno,ntlines[lineno],fontsize=9,fontname='Liberation Mono')
        if lineno == 0:
            ax.text(0.75,lineno,str(len(seq)),fontname='Liberation Mono',fontweight='bold')
        else:
            ax.text(0.75,lineno,str((len(ntlines))*60-lineno*60),fontname='Liberation Mono',fontweight='bold')
        ax.text(0,lineno-0.3,aalines[lineno],fontsize=9,fontname='Liberation Mono',color='green',fontweight='bold')
    ax.set_ylim((0,len(ntlines)))
    ax.axis('off')
    return fig,ax

def basicStats(seq):
    from CodOpY.misc import codon_choices
    import numpy as np
    from utils import makeCodedict

    stats = {}

    #calculate GC content
    stats['GC Content'] = sum([1 for nt in seq if nt in ['G','C']])/len(seq)

    #calculate codon optimality