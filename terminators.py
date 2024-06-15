def map_typeII_transcription_terminators(seq):
    import pandas as pd
    motifs = pd.read_csv('Termination sites.csv')
    motif_seq,motif_start,motif_stop,motif_strength,motif_context = [],[],[],[],[]
    for row in motifs.iterrows():
        motif = row[1]['Sequence']
        eff = row[1]['Efficiency']
        searchstart=0
        while motif in seq[searchstart:]:
            motif_seq.append(motif)
            motif_pos = seq[searchstart:].index(motif)+searchstart
            motif_start.append(motif_pos) 
            motif_stop.append(motif_pos+len(motif))
            motif_strength.append(eff)
            motif_context.append(seq[motif_pos+len(motif):motif_pos+len(motif)+8])
            searchstart = motif_pos + 1
    return pd.DataFrame({'Sequence':motif_seq,'From':motif_start,'To':motif_stop,'Strength':motif_strength,
                         'Context':motif_context})

def plotTypeIITranscriptionTerminators(seq):
    
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import pandas as pd
    from terminators import map_typeII_transcription_terminators
    from utils import find_orf, context_score
    
    terminators = map_typeII_transcription_terminators(seq).sort_values(by='From')
        
    usecmap=cm.coolwarm
    fig,ax = plt.subplots(figsize=(8,2.5))
    if terminators.shape[0] == 0:
        ax.text((find_orf(seq)[0]+find_orf(seq)[1])/2,1.1,'No Type II terminators found',
               horizontalalignment='center',fontweight='bold',color='green')
    else:
        unitsize=len(seq)/10
        labelexistsat=[-1*unitsize]
        for row in terminators.iterrows():
            newlabelpoint = row[1]['From']
            corescore = row[1]['Strength']/66.7
            contextscore = context_score(row[1]['Context'])
            totalscore = (corescore+contextscore)/2
            color = usecmap(totalscore)
            if newlabelpoint - labelexistsat[-1] > unitsize:
                lineheight=0.25
            else:
                lineheight-=0.07
            ax.plot((newlabelpoint,newlabelpoint),(1,1+lineheight),color=color,linewidth=1)
            ax.text(newlabelpoint,1+lineheight+0.025,row[1]['From'],fontsize=5,fontweight='bold')
            ax.text(newlabelpoint,1+lineheight,row[1]['Sequence'],fontsize=5,fontweight='bold')
            ax.text(newlabelpoint,1+lineheight-0.025,row[1]['Context'],fontsize=5)
            labelexistsat.append(newlabelpoint)
    
    ax.plot((0, len(seq)),(1,1),linewidth=2,color='grey')
    ax.plot(find_orf(seq),(1,1),linewidth=5,color='black')
    ax.text(0,0.95,'0',horizontalalignment='center',fontsize=7)
    ax.text(len(seq),0.95,str(len(seq)),horizontalalignment='center',fontsize=7)
    ax.text(len(seq)/2,0.75,'Termination Risk',horizontalalignment='center',fontsize=9)
    ax.text(0,0.85,'Low',fontsize=9,fontweight='bold',verticalalignment='center')
    ax.text(len(seq),0.85,'High',fontsize=9,fontweight='bold',
            horizontalalignment='right',verticalalignment='center')
    for p in range(int(len(seq)*0.1),int(len(seq)*0.9),int(len(seq)/10)):
        plt.scatter(p,0.85,color=usecmap(p/len(seq)))
    ax.set_ylim((0.7,1.35))
    ax.axis('off')
    return fig,ax
    
#=======================================================================

def plotStemLoops(seq,limit=-7.5,window=25):
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from utils import find_orf,context_score
    from seqfold import dg
    #determine local secondary structure
    sary = []
    for wdw in [seq[n:n+window] for n in range(0,len(seq)-window)]:
        sary.append(dg(wdw, temp = 37.0))
    #plot the location of local secondary structure lower than the set limit
    rel_unit=len(seq)/20
    usecmap=cm.coolwarm
    fig,ax = plt.subplots(figsize=(8,2.5))
    ax.plot((0, len(seq)),(1,1),linewidth=2,color='grey')
    ax.plot(find_orf(seq),(1,1),linewidth=5,color='black')
    ax.text(0,0.95,'0',horizontalalignment='center',fontsize=7)
    ax.text(len(seq),0.95,str(len(seq)),horizontalalignment='center',fontsize=7)
    if len(sary) == 0:
        ax.text(len(seq)/2,1.1,'No local secondary structure below ' + str(limit) + ' detected.',
               horizontalalignment='center',fontweight='bold',fontsize=8,color='green')
    else:
        for x,val in enumerate(sary):
            if val < limit:
                #is there a lower point within 20 nt?
                if val == min(sary[x-10:x+10]):
                    ax.text(x,1.05,'||',horizontalalignment='center',fontsize=7)
                    ax.text(x,1.2,r'$\Delta$G = '+str(val),horizontalalignment='center',
                            fontsize=7,color=usecmap((val-limit)/(-11-limit)))
                    ax.text(x,1.15,'Context = ',horizontalalignment='center',fontsize=7)
                    ax.text(x,1.1,seq[x+25:x+33],horizontalalignment='center',
                            fontsize=7,color=usecmap(context_score(seq[x+25:x+33])))
    ax.text(len(seq)/2,0.75,'Termination Risk',horizontalalignment='center',fontsize=9)
    ax.text(0,0.85,'Low',fontsize=9,fontweight='bold',verticalalignment='center')
    ax.text(len(seq),0.85,'High',fontsize=9,fontweight='bold',
            horizontalalignment='right',verticalalignment='center')
    for p in range(int(len(seq)*0.1),int(len(seq)*0.9),int(len(seq)/10)):
        plt.scatter(p,0.85,color=usecmap(p/len(seq)))
    ax.set_ylim((0.7,1.35))
    ax.axis('off')
    return fig,ax