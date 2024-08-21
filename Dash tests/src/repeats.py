def reportNtRepeats(seq,complexity_limit=20):
    
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from .utils import find_orf
    
    #construct a list of all 6-mer nucleotides
    nts = ['A','C','G','T']
    k6mers = [n1+n2+n3+n4+n5+n6 for n1 in nts for n2 in nts for n3 in nts for n4 in nts for n5 in nts for n6 in nts]
    #break seq into initial (6) mers
    mersize = 6
    seqmers = [seq[n:n+mersize] for n in range(len(seq)-mersize)]
    #set count of the initial mers
    countmers = {mersize:dict()}
    for mer in k6mers:
        if seqmers.count(mer) > 1:
            countmers[6][mer] = seqmers.count(mer)
    #extend repeated 6-mers and continue
    extend_further = True
    while extend_further:
        extend_further = False
        mersize += 1
        countmers[mersize] = dict()
        seqmers = [seq[n:n+mersize] for n in range(len(seq)-mersize)]
        for mer in countmers[mersize-1]:
            for nt in nts:
                plusmer = mer + nt
                if seqmers.count(plusmer) > 1:
                    countmers[mersize][plusmer] = seqmers.count(plusmer)
                    extend_further = True
    #map reportable repeats
    xs,lengths,motifs,complexities = [],[],[],[]
    complexity_limit = 20 # product of k and counts
    for k, kdict in countmers.items():
        for merseq,count in kdict.items():
            if min((0.07*k)**2,1)*k*count > complexity_limit:
                subseq = seq
                xs.append([])
                lengths.append(k)
                motifs.append(merseq)
                complexities.append(min((0.07*k)**2,1)*k*count)
                while merseq in subseq:
                    xs[-1].append(subseq.index(merseq))
                    subseq = subseq[xs[-1][-1]+1:]
    #assess whether any of the motifs is present within a logner motif,
    #if so, discard
    clean_xs,clean_lengths,clean_motifs,clean_complexities = [],[],[],[]
    for idx, motif in enumerate(motifs):
        is_submotif = False
        for n in range(idx+1,len(motifs)):
            if motif in motifs[n]:
                is_submotif = True
        if not is_submotif:
            clean_xs.append(xs[idx])
            clean_lengths.append(lengths[idx])
            clean_motifs.append(motifs[idx])
            clean_complexities.append(complexities[idx])
    #if more than three entries disaply ony the three most complex
    if len(clean_xs) > 3:
        final_xs,final_lengths,final_motifs,final_complexities = [],[],[],[]
        ordered_complexities = clean_complexities.copy()
        ordered_complexities.sort()
        selected_complexities = ordered_complexities[-3:]
        for idx,c in enumerate(clean_complexities):
            if c in selected_complexities:    
                final_xs.append(clean_xs[idx])
                final_lengths.append(clean_lengths[idx])
                final_motifs.append(clean_motifs[idx])
                final_complexities.append(clean_complexities[idx])
    else:
        final_xs = clean_xs
        final_lengths = clean_lengths
        final_motifs = clean_motifs
        final_complexities = clean_complexities
        #plot results
        usecmap=cm.coolwarm
    fig,ax = plt.subplots(figsize=(8,2.5))
    #display the basic sequence diagram
    ax.plot((0, len(seq)),(1,1),linewidth=2,color='grey')
    ax.plot(find_orf(seq),(1,1),linewidth=5,color='black')
    ax.text(0,0.95,'0',horizontalalignment='center',fontsize=7)
    ax.text(len(seq),0.95,str(len(seq)),horizontalalignment='center',fontsize=7)
    ax.set_ylim((0.7,1.35))
    if len(final_xs) == 0:
        ax.text(len(seq)/2,1.15,'No notable repeat sequences detected.',color='green',
                         fontweight='bold',horizontalalignment='center')
    for index, xvals in enumerate(final_xs):
        xvals.sort()
        yvals = [1.1+index*0.1]*len(xvals)
        plt.scatter(xvals,yvals,marker='|',color = usecmap(final_complexities[index]/40))
        plt.text(xvals[0],yvals[0]-0.05,final_motifs[index],fontsize=6)

    ax.text(len(seq)/2,0.75,'Repeat complexity',horizontalalignment='center',fontsize=9)
    ax.text(0,0.85,'Low',fontsize=9,fontweight='bold',verticalalignment='center')
    ax.text(len(seq),0.85,'High',fontsize=9,fontweight='bold',
            horizontalalignment='right',verticalalignment='center')
    for p in range(int(len(seq)*0.1),int(len(seq)*0.9),int(len(seq)/10)):
        plt.scatter(p,0.85,color=usecmap(p/len(seq)))
    ax.axis('off')
    return fig,ax