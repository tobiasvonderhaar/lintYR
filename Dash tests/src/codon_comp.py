def plotCodonCounts(seq, genCode=1, mode='abs'):
    #has 'mode' been set to an allowed value?
    if mode not in ['abs','rel']:
        raise ValueError('mode must be either abs or rel')
        
    from src.utils import makeReverseCodedict, countCodons, find_orf
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    
    #generate and tabulate codon counts or frequencies
    codonDict = makeReverseCodedict(genCode=genCode)
    orf_coords = find_orf(seq)
    orf=seq[orf_coords[0]:orf_coords[1]]
    counts = countCodons(orf)
    countDict = {}
    if mode == 'rel':
        freqDict = {}
    freqDict={}
    sortedCodonDict = dict(sorted(codonDict.items()))
    for aa in sortedCodonDict.keys():
        countDict[aa] = []
        for codon in codonDict[aa]:
            countDict[aa].append(counts[codon])
        if mode == 'rel':
            freqDict[aa] = [v / sum(countDict[aa]) for v in countDict[aa]]
    #prepare lists with labels and values in plotting order
    if mode == 'rel':
        plotDict = freqDict
    else:
        plotDict = countDict
    codons, values, aas,aaxpos = [],[],[],[]
    lastxposref = 0
    for aa in sortedCodonDict.keys():
        aas.append(aa)
        aaxpos.append(((lastxposref + (lastxposref + len(codonDict[aa]))) / 2) - 1)
        lastxposref += len(codonDict[aa])
        codons += codonDict[aa]
        values += plotDict[aa]
    #plot figure
    fig,ax = plt.subplots(figsize = (8,3))
    #draw bars
    ax.bar(codons, values,color='black',width=0.5)
    top = ax.get_ylim()[1]
    ax.set_xlim(-1,64)
    ax.set_ylim((0,top*1.2))
    #prepare patches to identify aa areas
    rbase = -0.5
    for idx, aa in enumerate(aas):
        #only draw pactches for every second box
        if idx%2 == 1:
            rect = plt.Rectangle((rbase,0),\
                                 len(plotDict[aa]),top*1.2,\
                                     facecolor='grey',alpha=0.5,zorder=0)
            ax.add_patch(rect)
        rbase += len(plotDict[aa])
    #format other aspects of plot
    ax.tick_params(axis='x', labelrotation = 90, labelsize=6)
    if mode == 'rel':
        ax.set_ylabel('Frequency')
    else:
        ax.set_ylabel('Count')
    ax.set_xlabel('Codon')
    ax.set_title('Codon usage')
    #add amino acid labels
    aaDict = dict(zip(aas,aaxpos))
    for aa, pos in aaDict.items():
        ax.text(pos, top*1.1, aa, fontsize=8)
    return fig,ax