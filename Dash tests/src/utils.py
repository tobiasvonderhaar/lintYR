def context_score(seq_context):
    w = len(seq_context)
    max_score = sum([1-(i*1/w) for i,nuc in enumerate(seq_context)])
    return sum([1-(i*1/w) for i,nuc in enumerate(seq_context) if nuc == 'T'])/max_score

#================================================================

def find_orf(seq):
    RNA=False
    if 'U' in seq:
        RNA=True
        seq = seq.replace('U','T')
    start = seq.index('ATG')
    codon_seq = [seq[start:][n:n+3] for n in range(0,len(seq[start:]),3)]
    stops=[]
    try:
        stops.append(codon_seq.index('TAA')*3+3+start)
    except:
        pass
    try:
        stops.append(codon_seq.index('TGA')*3+3+start)
    except:
        pass
    try:
        stops.append(codon_seq.index('TAG')*3+3+start)
    except:
        pass
    if len(stops)>=1:
        stop=min(stops)
    else:
        stop=len(seq)
    return (start,stop)

#================================================================

def countCodons(seq):
    #convert sequence DNA (in case RNA) and uppercase
    seq = seq.upper()
    seq = seq.replace('U','T')
    #divide sequence into full length codons
    codonSeq = [seq[n:n+3] for n in range(0,len(seq),3) if len(seq[n:n+3]) == 3]
    #make a reference list of codons
    codons = []
    for nt1 in ['A','C','G','T']:
        for nt2 in ['A','C','G','T']:
            for nt3 in ['A','C','G','T']:
                codons.append(nt1+nt2+nt3)
    counts = {}
    for codon in codons:
        counts[codon] = codonSeq.count(codon)
    return counts

#==============================================================================

def makeCodedict(genCode=1):
    '''
    Returns a code dictionary specifying the amino acids encoded by each codon.

    Parameters
    ==========
    gen_code : int
        The genetic code to be used. Values follow the definitions used by
        NCBI at
        www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes


    Returns
    dict
        A dictionary specifying the 64 codons and their encoded amino acids in
        one letter code (with *=Stop).
    '''

    codedict = {'AAA':'K','AAC':'N','AAG':'K','AAT':'N','ACA':'T','ACC':'T',
                'ACG':'T','ACT':'T','AGA':'R','AGC':'S','AGG':'R','AGT':'S',
                'ATA':'I','ATC':'I','ATG':'M','ATT':'I','CAA':'Q','CAC':'H',
                'CAG':'Q','CAT':'H','CCA':'P','CCC':'P','CCG':'P','CCT':'P',
                'CGA':'R','CGC':'R','CGG':'R','CGT':'R','CTA':'L','CTC':'L',
                'CTG':'L','CTT':'L','GAA':'E','GAC':'D','GAG':'E','GAT':'D',
                'GCA':'A','GCC':'A','GCG':'A','GCT':'A','GGA':'G','GGC':'G',
                'GGG':'G','GGT':'G','GTA':'V','GTC':'V','GTG':'V','GTT':'V',
                'TAA':'*','TAC':'Y','TAG':'*','TAT':'Y','TCA':'S','TCC':'S',
                'TCG':'S','TCT':'S','TGA':'*','TGC':'C','TGG':'W','TGT':'C',
                'TTA':'L','TTC':'F','TTG':'L','TTT':'F'}
    if genCode == 2: #Vertebrate mitochondrial
        codedict['AGA'] = '*'
        codedict['AGG'] = '*'
        codedict['ATA'] = 'M'
        codedict['TGA'] = 'W'
    elif genCode == 3: #Yeast mitochondrial
        codedict['ATA'] = 'M'
        codedict['TGA'] = 'W'
        codedict['CTA'] = 'T'
        codedict['CTC'] = 'T'
        codedict['CTG'] = 'T'
        codedict['CTT'] = 'T'
    elif genCode == 4: #Protozoan, Coelenterate MitoMito, Mycopl./Spiropl. nuclear
        codedict['TGA'] = 'W'
    elif genCode == 5: #Invertebrate Mitochondrial
        codedict['AGA'] = 'S'
        codedict['AGG'] = 'S'
        codedict['ATA'] = 'M'
        codedict['TGA'] = 'W'
    elif genCode == 6: #Ciliate Nuclear
        codedict['TAA'] = 'Q'
        codedict['TAG'] = 'Q'
    elif genCode == 9: #Echinoderm and Flatworm Mitochondrial
        codedict['AAA'] = 'N'
        codedict['AGA'] = 'S'
        codedict['AGG'] = 'S'
        codedict['TGA'] = 'W'
    elif genCode == 10: #Euplotid Nuclear
        codedict['TGA'] = 'C'
    elif genCode == 12: #Yeast CTG Clade nuclear
        codedict['CTG'] = 'S'
    elif genCode == 13: #Ascidian Mitochondrial
        codedict['AGA'] = 'G'
        codedict['AGG'] = 'G'
        codedict['ATA'] = 'M'
        codedict['TGA'] = 'W'
    elif genCode == 14: #Alternative Flatworm Mitochondrial
        codedict['AAA'] = 'N'
        codedict['AGA'] = 'S'
        codedict['AGG'] = 'S'
        codedict['TAA'] = 'Y'
        codedict['TGA'] = 'W'
    elif genCode == 16: #Chlorophycean Mitochondrial
        codedict['TAG'] = 'L'
    elif genCode == 21: #Trematode Mitochondrial
        codedict['AAA'] = 'N'
        codedict['AGA'] = 'S'
        codedict['AGG'] = 'S'
        codedict['ATA'] = 'M'
        codedict['TGA'] = 'W'
    elif genCode == 22: #Scenedesmus obliquus Mitochondrial
        codedict['TAG'] = 'L'
        codedict['TCA'] = '*'
    elif genCode == 24: #Rhabdopleuridae Mitochondrial
        codedict['AGA'] = 'S'
        codedict['AGG'] = 'K'
        codedict['TGA'] = 'W'
    elif genCode == 25: #Gracilibacteria 
        codedict['ATA'] = 'G'
    elif genCode == 26: #Pachysolen tannophilus nuclear 
        codedict['CUG'] = 'A'
    elif genCode == 29: #Mesodinium nuclear
        codedict['TAA'] = 'Y'
        codedict['TAG'] = 'Y'
    elif genCode == 30: #Peritrich nuclear
        codedict['TAA'] = 'E'
        codedict['TAG'] = 'E'
    elif genCode == 33: #Cephalodiscidae Mitochondrial
        codedict['TAA'] = 'Y'
        codedict['TAG'] = 'W'
        codedict['AGA'] = 'S'
        codedict['AGG'] = 'K'
    return codedict

#==============================================================================

def translate(dnaSeq,genCode=1):
    '''
    Returns an amino acid sequence translated for an input DNA 
    sequence.

    Parameters
    ==========
    seq : str
        The DNA sequence to be translated.
    genCode : int
        The genetic code to be used. The standard genetic code is 1, variant 
        codes follow the numbering system used at 
        www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes.
    
    Returns
    str
        The amino acid sequence translated from seq.
    '''
    from src.utils import makeCodedict
    codedict = makeCodedict(genCode=genCode)
    dnaSeq = dnaSeq.replace('U','T').upper()
    codonSeq = [dnaSeq[n:n+3] for n in range(0,len(dnaSeq),3)]
    try:
        proteinSeq = ''.join([codedict[c] for c in codonSeq if len(c) == 3])
    except:
        raise ValueError('The DNA sequence cannot be translated.')

    return proteinSeq

#==============================================================================

def makeReverseCodedict(genCode=1):
    from src.utils import makeCodedict
    forwardDict = makeCodedict(genCode=genCode)
    reverseDict = {}
    for k,v in forwardDict.items():
        if v in reverseDict.keys():
            reverseDict[v].append(k)
        else:
            reverseDict[v] = [k]
    return reverseDict