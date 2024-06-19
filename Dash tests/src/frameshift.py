def mapFrameshifts(seq, plot=True):
    
    import matplotlib.pyplot as plt
    from src.utils import find_orf
    
    orf_coords = find_orf(seq)
    orf_seq = seq[orf_coords[0]:orf_coords[1]]
    codonSeq = [orf_seq[n:n+3] for n in range(0,len(seq),3)]
    slipperySites = []
    slipperySeqs = []
    for index, codon in enumerate(codonSeq):
        if codonSeq[index] == 'TTT':
            if codonSeq[index+1][0] in ['T','C']:
                slipperySites.append((index+1)*3-2+orf_coords[0])
                slipperySeqs.append('TTT' + codonSeq[index+1][0])
    if plot:       
        y = [1.2 for x in slipperySites]
        fig,ax = plt.subplots(figsize=(8,2))
        ax.scatter(slipperySites,y,c='red',marker='v')
        texty=1.3
        if len(slipperySites)==0:
            ax.text(len(seq)/2,1.5,'No slippery sites detected.',color='green',
                     fontweight='bold',horizontalalignment='center')
        else:
            for idx,site in enumerate(slipperySites):
                if idx > 0:
                    if ((site-slipperySites[idx-1])<0.1*len(seq)) and (texty<2.2):
                        texty+=0.3
                    else:
                        texty=1.3
                ax.text(site,texty,slipperySeqs[idx],horizontalalignment='center',fontsize=8)
                ax.text(site,texty+0.15,str(slipperySites[idx]),horizontalalignment='center',fontsize=8,fontweight='bold')
        ax.plot((0, len(seq)),(1,1),linewidth=2,color='grey')
        ax.plot(find_orf(seq),(1,1),linewidth=5,color='black')
        ax.text(0,0.75,'0',horizontalalignment='center',fontsize=7)
        ax.text(len(seq),0.75,str(len(seq)),horizontalalignment='center',fontsize=7)
        ax.set_ylim((0.5,2.5))
        ax.axis('off')
        return fig,ax 
    else:
        return slipperySites