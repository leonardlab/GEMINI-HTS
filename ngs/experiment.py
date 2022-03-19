
from .classifier import NGSClassifier
from .descriptor import NGSDescriptor

from segments.experiment import Experiment
from segments.util import reverse_complement


class NGSExperiment(Experiment):

    def aliases(self):
        return [ "ngs" ]

    def fragmentDescriptor(self):
        return NGSDescriptor().fragmentDescriptor()

    def classifier(self, descriptor):
        return NGSClassifier(descriptor)

    @property
    def target(self):
        return NGSDescriptor().target

    def barcode(self, r1s, r2s):
        # return the barcode; only used if barcode pruning. return
        # what it would be on a match; it'll only be a prune skip if
        # there's a matching one already in the cache.
        return r1s.characters[:5] + "_" + reverse_complement(r2s.characters[:5])

    def fragmentCacheKey(self, frag):
        # using this key looks like it changes a very few size=1
        # deletions and a couple of muts. seems to be similar minimal
        # impact as r1r2 cache key.
        return None
        #return frag.seq.characters[5:-5]

    def r1r2CacheKey(self, r1s, r2s):
        # looks like a speedup of ~2x using this cache key (saving
        # aligning every fragment), but it will change the results
        # very slightly (single-muts only). ok for dev runs, but
        # should probably not use it for final numbers.  
        #return None
        return r1s.characters[5:] + "_" + r2s.characters[5:]

    def addExtraReads(self, reads):
        ngsd = NGSDescriptor()
        target = ngsd.fp + ngsd.target + ngsd.rp
        genome = _parsedHumanGenomeFragment()
        tidx = genome.find(target)
        assert(-1 != tidx)
        left, right = genome[:tidx], genome[tidx + len(target):]
        assert(left + target + right == genome)
        assert(target not in left)
        assert(target not in right)
        if 0:
            reads.add("cxcr4_left", left, 16)
            reads.add("cxcr4_right", right, 16)
            reads.add("cxcr4_left_rc", reverse_complement(left), 16)
            reads.add("cxcr4_right_rc", reverse_complement(right), 16)


def _parsedHumanGenomeFragment():
    glines = _human_CXCR4_genome.split('\n')
    g = ''
    for line in glines:
        bit = ''.join(line.strip().split(' ')[1:]).upper()
        if bit:
            g += bit
    assert(8747 == len(g))
    return g

_human_CXCR4_genome = '''
        1 caattctgaa tcctgccttt tgcacttaat gtttcataag tatttcccca tgtcactaaa
       61 aattcttcca aataacattc acgatgtcca tatggaattt cagatgtgga tgaaccaaaa
      121 tcttgtcaac tattccacta acagtggtta tttagggatg ttcagacatt tcactattta
      181 aaaaaaaatg tttccacaaa tacctttgtg gcataagttt ttatgagtgg agttactgtt
      241 ctgaagttcc tgctgaatag aaaatgcttt ccagtgaggc tgtcccaagc cacattccca
      301 tcagtgacaa gcgagagaca gctggtcttt tcaaatccgg agaccaaata ttatctttga
      361 aaaaaaatgg atttttgcct aatttggtag tcaccaaata gcatctcatt gttcttttaa
      421 ttatctgctt ccttttagta gagatcccta aaaagatctg aaaggagtct tcagataaag
      481 gaaggagctt tcttttgtct gtctacaatc aacaaatatt tattatgcaa accattttgc
      541 tccgagtttt ctcctctttc cctttttgga cagatttggg agatctcacc tttcaggttt
      601 tagacatcgt gcagggagga gttttgaggt agggtgcagc ttacggtcca ggataaaaca
      661 tactgattct gccactacca ggctttgtga aaagcaagtc atgaaaacgc tctgaaattc
      721 tagaccttca gtagatagga tctaccgtgt ctataaaaat atgaagatcc ttaagtttta
      781 ttaaagattc gaaaaaagta aaagtgtttt tacggtttta ttttcatttt tatttcttac
      841 cgttatcgtt tattataaag gatattataa aggatacaga tgaagagata cgtaatgcaa
      901 ggcctgtgag aaggggcgtg gagcttccga aacctcttcc agccaccacc ctccaagaac
      961 ctggagtttc tttttttttt tttaattcta caaatgtaat attagaattg attttatctg
     1021 gccattagtg tgtgtcctaa ctcgttcgtt tctgagagtc ccatctcccg gcccgggata
     1081 tcatctttcc tgtgtcagtg aaagtgcaga gtagatgaga acctttaacc accaacatta
     1141 gggaggggtc ccagacaaag ggggtaagtc atgctctgta gagaaaaggt tccctgcctc
     1201 cgaactacct ctggaacact ccagtaaatg tttcctcttt tgatatagaa aagagggatc
     1261 gtgtgtagag tgcagtctgg gcaatccctc tcctcgggac catttcgggg taggggcctc
     1321 tggggtccgt gtcgcgacgc acgcgcctcg gtcccagcta tctccgcagc gggccacccc
     1381 gcctgcggac gcagtttctc ggccccgccc cacactcgct cccccgcccc acccagtctc
     1441 cgcgccggag ggaagtggcg cgagggggaa agcactgtct gcgcgcccac tgcaaacctc
     1501 agccagtctg agatcgcttt aaacgtctga cccccacccc cactccgccc cgcccagttc
     1561 ttcaacctaa tttctgattc gtgccaaagc ttgtcctctg ctcaaaatcg tggaagacgc
     1621 cgagtatggg gaccgaagac ctgggttcaa gcccggcttg gaatccctgc ccatccctgg
     1681 catttcatct ctccgggctt atttgctggt ttctccgaat gcgggccttg tctggttcac
     1741 gctggatccc caacgcctag aacagtgcgt ggcacgcagt tcgtccttct ataaatatcg
     1801 gactaaatgc atctctgtga tggtaatacc cacacggtgt tgtgagaatg aatgagtgat
     1861 tctgtgcaag ttcctagtga tctgttacaa aaagtactgg tcgctaaatt actcttataa
     1921 taaagcatac ttttaggata ataaagcact attcgcgaat tggttaccgc tattatgaaa
     1981 ttactgagca atacatatct acatctgatc agtctccaga attatgccaa atcctacctt
     2041 cttctgaaag tatctcctaa ttatctgcac ctgaccctag tgatgctgtg aatgtgcaag
     2101 tatagctaca tcctccgaag gaaaggatct ttactccttt tacctcctga atgggctgcg
     2161 tctgctgaaa gcgcggggga atgggcggtt ggaagcttgg ccctacttcc agcattgccg
     2221 cctactggtt gggttactcc agcaagtcac tccccttccc tgggcctcag tgtctctact
     2281 gtagcattcc caggtctgga attccatcca ctttagcaag gatggacgcg ccacagagag
     2341 acgcgttcct agcccgcgct tcccacctgt cttcaggcgc atcccgcttc cctcaaactt
     2401 aggaaatgcc tctgggaggt cctgtccggc tccggactca ctaccgacca cccgcaaaca
     2461 gcagggtccc ctgggcttcc caagccgcgc acctctccgc cccgcccctg cgccctcctt
     2521 cctcgcgtct gcccctctcc cccaccccgc cttctccctc cccgccccag cggcgcatgc
     2581 gccgcgctcg gagcgtgttt ttataaaagt ccggccgcgg ccagaaactt cagtttgttg
     2641 gctgcggcag caggtagcaa agtgacgccg agggcctgag tgctccagta gccaccgcat
     2701 ctggagaacc agcggttacc atggagggga tcagtgtaag tccagtttca acctgctttg
     2761 tcataaatgt acaaacgttt gaacttagag cgcagcccct ctccgagcgg gcagaagcgg
     2821 ccaggacatt ggaggtaccc gtactccaaa aaagggtcac cgaaaggagt tttcttgacc
     2881 atgcctatat agtgcgggtg ggtggggggg gagcaggatt ggaatctttt tctctgtgag
     2941 tcgaggagaa acgactggaa agagcgttcc agtggctgca tgtgtctccc ccttgagtcc
     3001 cgccgcgcgc ggcggcttgc acgctgtttg caaacgtaag aacattctgt gcacaagtgc
     3061 agagaaggcg tgcgcgctgc ctcgggactc agaccaccgg tctcttcctt ggggaagcgg
     3121 ggatgtcttg gagcgagtta cattgtctga atttagaggc ggagggcggc gtgcctgggc
     3181 tgagttccca ggaggagatt gcgcccgctt taacttcggg gttaagcgcc tggtgactgt
     3241 tcttgacact gggtgcgtgt ttgttaaact ctgtgcggcc gacggagctg tgccagtctc
     3301 ccagcacagt aggcagaggg cgggagaggc gggtggaccc accgcgccga tcctctgagg
     3361 ggatcgagtg gtggcagcag ctaggagttg atccgcccgc gcgctttggg tttgaggggg
     3421 aaaccttccc gccgtccgaa gcgcgcctct tccccacggc cgcgagtggg tcctgcagtt
     3481 cgagagtttg gggtcgtgca gaggtcagcg gagtggtttg acctcccctt tgacaccgcg
     3541 cagctgccag ccctgagatt tgcgctccgg ggataggagc gggtacgggg tgaggggcgg
     3601 gggcggttaa gaccgcacct gggctgccag gtcgccgccg cgaagactgg caggtgcaag
     3661 tggggaaacc gtttggctct ctccgagtcc agttgtgatg tttaaccgtc ggtggtttcc
     3721 agaaaccttt tgaaaccctc ttgctaggga gtttttggtt tcctgcagcg gcgcgcaatt
     3781 caaagacgct cgcggcggag ccgcccagtc gctccccagc accctgtggg acagagcctg
     3841 gcgtgtcgcc cagcggagcc cctgcagcgc tgcttgcggg cggttggcgt gggtgtagtg
     3901 ggcagccgcg gcggcccggg gctggacgac ccggcccccc gcgtgcccac cgcctggagg
     3961 cttccagctg cccacctccg gccgggttaa ctggatcagt ggcggggtaa tgggaaacca
     4021 cccgggagag tgaggaaatg aaacttgggg cgaggaccac gggtgcagac cccgttacct
     4081 tctccaccca ggaaaatgcc ccgctcccta acgtcccaaa cgcgccaagt gataaacacg
     4141 aggatggcaa gagacccaca caccggagga gcgcccgctt gggggaggag gtgccgtttg
     4201 ttcattttct gacactcccg cccaatatac cccaagcacc gaagggcctt cgttttaaga
     4261 ccgcattctc tttacccact acaagttgct tgaagcccag aatggtttgt atttaggcag
     4321 gcgtgggaaa attaagtttt tgcgccttag gagaatgagt ctttgcaacg cccccgccct
     4381 ccccccgtga tcctcccttc tcccctcttc cctccctggg cgaaaaactt cttacaaaaa
     4441 gttaatcact gcccctccta gcagcaccca ccccaccccc cacgccgcct gggagtggcc
     4501 tctttgtgtg tatttttttt ttcctcctaa ggaaggtttt ttttcttccc tctagtgggc
     4561 ggggcagagg agttagccaa gatgtgactt tgaaaccctc agcgtctcag tgcccttttg
     4621 ttctaaacaa agaattttgt aattggttct accaaagaag gatataatga agtcactatg
     4681 ggaaaagatg gggaggagag ttgtaggatt ctacattaat tctcttgtgc ccttagccca
     4741 ctacttcaga atttcctgaa gaaagcaagc ctgaattggt tttttaaatt gctttaaaaa
     4801 atttttttaa ctgggttaat gcttgctgaa ttggaagtga atgtccattc ctttgcctct
     4861 tttgcagata tacacttcag ataactacac cgaggaaatg ggctcagggg actatgactc
     4921 catgaaggaa ccctgtttcc gtgaagaaaa tgctaatttc aataaaatct tcctgcccac
     4981 catctactcc atcatcttct taactggcat tgtgggcaat ggattggtca tcctggtcat
     5041 gggttaccag aagaaactga gaagcatgac ggacaagtac aggctgcacc tgtcagtggc
     5101 cgacctcctc tttgtcatca cgcttccctt ctgggcagtt gatgccgtgg caaactggta
     5161 ctttgggaac ttcctatgca aggcagtcca tgtcatctac acagtcaacc tctacagcag
     5221 tgtcctcatc ctggccttca tcagtctgga ccgctacctg gccatcgtcc acgccaccaa
     5281 cagtcagagg ccaaggaagc tgttggctga aaaggtggtc tatgttggcg tctggatccc
     5341 tgccctcctg ctgactattc ccgacttcat ctttgccaac gtcagtgagg cagatgacag
     5401 atatatctgt gaccgcttct accccaatga cttgtgggtg gttgtgttcc agtttcagca
     5461 catcatggtt ggccttatcc tgcctggtat tgtcatcctg tcctgctatt gcattatcat
     5521 ctccaagctg tcacactcca agggccacca gaagcgcaag gccctcaaga ccacagtcat
     5581 cctcatcctg gctttcttcg cctgttggct gccttactac attgggatca gcatcgactc
     5641 cttcatcctc ctggaaatca tcaagcaagg gtgtgagttt gagaacactg tgcacaagtg
     5701 gatttccatc accgaggccc tagctttctt ccactgttgt ctgaacccca tcctctatgc
     5761 tttccttgga gccaaattta aaacctctgc ccagcacgca ctcacctctg tgagcagagg
     5821 gtccagcctc aagatcctct ccaaaggaaa gcgaggtgga cattcatctg tttccactga
     5881 gtctgagtct tcaagttttc actccagcta acacagatgt aaaagacttt tttttatacg
     5941 ataaataact tttttttaag ttacacattt ttcagatata aaagactgac caatattgta
     6001 cagtttttat tgcttgttgg atttttgtct tgtgtttctt tagtttttgt gaagtttaat
     6061 tgacttattt atataaattt tttttgtttc atattgatgt gtgtctaggc aggacctgtg
     6121 gccaagttct tagttgctgt atgtctcgtg gtaggactgt agaaaaggga actgaacatt
     6181 ccagagcgtg tagtgaatca cgtaaagcta gaaatgatcc ccagctgttt atgcatagat
     6241 aatctctcca ttcccgtgga acgtttttcc tgttcttaag acgtgatttt gctgtagaag
     6301 atggcactta taaccaaagc ccaaagtggt atagaaatgc tggtttttca gttttcagga
     6361 gtgggttgat ttcagcacct acagtgtaca gtcttgtatt aagttgttaa taaaagtaca
     6421 tgttaaactt acttagtgtt atgttctgat ttctgttgac attcttttgg ctagtagaag
     6481 acaaaagtaa tacatttatg gtatgcaaag cactatccta ggtatttcat tgtaatattt
     6541 tacttacccc ttatcacaac tctgatagat tctgcttctg ttactaatta cattttatag
     6601 aagaggaaac ggaggcacag aaagcctaag taacttggtt aaaggcatgt agtaagtatc
     6661 aaatcctgta ttttaaacca ggtaacatga cttaacgaat ctgaagcctt caccacttta
     6721 aattcaaatg gaagtttaga aatggccagc cagcacctat ttgtatgaaa ggtcatcttt
     6781 cagaggataa gcatgtataa agaagaaaag gtatgcagtc gtgtttggat tttactccac
     6841 catccacttg tgaaacccag gtctgtgcaa tgccagacgg tgtgtgcttt cctcatccag
     6901 tatcctcagt gtagataacc atcactccct tttcacagac aagagaactg agattcagag
     6961 actttccata cattgcactt tcaagggggc aaagccaaga actaattctg tttattgttc
     7021 cagctcttgc tcttaactct tacctactat tgcccttcag aacacctggg cataagtcaa
     7081 ctgaactgct aataaagaaa gccaaaagtg aatgttttct tcataaaatt aaccatgacc
     7141 aaaatactcc tcttgtaata tcttctatgc aaatctcaac acttttattc ttaaactatc
     7201 gcaacaccta gcacctcctc aaggactcag ccaagcagct acaagttaat actgatattt
     7261 gttagagtca gaaggaaggt ccactgaagc aagctccctg ttgctcacat tttgcacaag
     7321 attttggaga cttatgtaac cacccgttgc tattaacacg accattgtgc aagccccagg
     7381 ctcttgagta aatttcagct ttggtttcta tttaaagata atttctaaac tctagccata
     7441 cctacctcac attggaacac aaacagggta cactccaggc atgcactcag ataataagta
     7501 ggatataatt acgacaatat ttggtctact tttagtaatt gtttctggca cagaaaatcc
     7561 attttggagg aaaaattgca atgccttatc tttctgaggc aaatcacatt tgttcaaggc
     7621 aaattataga tcctgtgaag ggaaataact taattactta aaatagaatc caatttggct
     7681 gtacattttt gctgccgtct atggatctgg ggtaattcaa agtggtattc atattctact
     7741 tgaggacaca attagatttc agataggaaa ttatcttgag gtttcttggt tttccctgag
     7801 aagcctaatt ggatcaccct tcatttaagc atagttttac atgcactctc tcaaaggctt
     7861 agtcttaaag ccacaaccat tgagacagac ttcacttgaa ccctctctat aaatatttat
     7921 tctccgggag acaatagaag aaatccttgg aaggcatgct ttttctttct catcttggct
     7981 tgaaacctcc ttaccccaga ttcctctcct ttaccgtgga gtcacaacaa aaggaactga
     8041 gccaaaacaa aattcccagt gtcaccagtc ttaatggata tttcattctc ccttggaaca
     8101 aagatggaat agcttttttt ccaaaagaaa aacaagcctt ggctctctcc ctgccccaaa
     8161 agggtgcccc ccacccccat cattctctgt cccaaccctg ccatgttaga gcgtctccaa
     8221 agccttccct gtgtcgtggt ttgtctgaca atgtggggaa acccagtctg ctggccagcc
     8281 cttgcatgaa gtagctgatt gttccctctc ctcatccctt atgaatgggg cccttgaagt
     8341 tcagtcatgt agattcagtt gtataatgaa agctaaaata tttaaattgt atgcatgctg
     8401 ccaataacag catacatctg acatctaact tattaataac attaagcctg caactagggg
     8461 ggaaagtgga tgttttttct tgcaaagcct ttgttttcct aaaatgacac ttgaaaattt
     8521 atctccccct actgcaggct tcccagcccc cttttataat tatgcttaaa ttaaaataat
     8581 gattctggga tactcttttg gggagatacc ctacaggctt tattttaata attgaactaa
     8641 gtgtttgtga ctttctccta gatattgtca aatattaaat aaaggctcca taaacaattg
     8701 agctgtctta ttcccagata atacccattt aggaggggca aggatcc
'''
