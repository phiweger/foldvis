### TODO

- ENA API wrapper to get the AF2 proteins by ID
- [ ] `fold.to_dataframe()` for easy integration with `sklearn` and `altair`
- [ ] consurf https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4987940/
- [ ] positive and negative control
    - https://elifesciences.org/articles/22709
    - transferrin receptor (science paper iron piracy)
    - pyoverdin/ ~ receptor paper (printed out)
- [ ] plot annotation for complexes should have an option to either use a common color scale (plddt) or separate ones (position), ie one for each chain
- [ ] COSMIS score of functional constraint https://www.nature.com/articles/s41467-022-30936-x#Sec12 -- delta simple conservation?
- [ ] explore dnamic stuff:
    - https://github.com/prody/ProDy
    - https://biologicalmodeling.org/coronavirus/home
    - membrane modelling https://github.com/pstansfeld/MemProtMD/
- [ ] in complexes, only the first seq is annotated -- have some nice way to support complexes or multichain (toggle off, delta color, mark interface, ...)
- [ ] TMAlign better? also, what is pulcra? https://twitter.com/thesteinegger/status/1517172868127150080?s=20&t=SCQJV1iKb7zuMytzKG3BQg
- [ ] replace pdb parsing with https://github.com/haddocking/pdb-tools
- [ ] when interactdome is not available, show active site prediction from pfam (there is a hmmer option to ptrdict those, too, I think)
- [ ] streamlit? https://github.com/napoles-uach/streamlit_3dmol
- [ ] deposit `3Dmol.js` code locally instead of using CDN
- [ ] install using `pip`
- [ ] allow display of binding sites, like [here](https://merenlab.org/2020/07/22/interacdome/)
- [ ] add [state sequence](https://github.com/steineggerlab/foldseek/issues/15) representing the structure
- [ ] return a contact map, and common protein stats 
- [ ] these are not of type "Structure" -- https://github.com/phiweger/foldvis/blob/main/foldvis/utils.py#L122
- [ ] check that `foldseek` available
- [ ] allow model selection
- [ ] download pdb structures
- [ ] get aligned coords like in https://www.rcsb.org/alignment, see https://www.rcsb.org/docs/tools/pairwise-structure-alignment
- [ ] scoring residue conserveation https://pubmed.ncbi.nlm.nih.gov/12112692/
- [ ] mimicry example https://www.pnas.org/doi/full/10.1073/pnas.2021785118 -- 7jtl and 4hr9; why no algignment w/ foldseek but Tmalign in https://www.rcsb.org/alignment? Would we find this alignment with the states only? Make sure mode 2 and min 0.1 really work on some test cases

surface fingerprints

https://twitter.com/ilvailvail/status/1519808480332369921?s=20&t=7OCi-h3hFU5wnYer_yhtcw

https://github.com/phiweger/foldvis/issues/1


```
foldseek createdb 4hr9.pdb queryDB
foldseek lndb queryDB_h queryDB_ss_h
foldseek convert2fasta queryDB_ss queryDB_ss.fasta

foldseek createdb 7jtl.pdb queryDB
foldseek lndb queryDB_h queryDB_ss_h
foldseek convert2fasta queryDB_ss queryDB_ss.fasta

DPDDADDPDDPDCCPPQFFWDWDWQADPQWVVRTDTATDTPDQAGQDPVRDGDNQKGKDFDKDWDWIWGDVVGPDPGHTDIDIDIGGDTIDIDRDD

DAEAEEEEAAQDKAQDDDPADPLWDKWKWFAFAQDPPGDTHTDPDDDDDDPDDDDDPSPQQVDVVPGMGRPDADPRGWMKMWIGNDPVSPDTDIHTYHYDYD
```


- [ ] get test examples

```
import requests, base64
r = requests.get('https://mmtf.rcsb.org/v1.0/full/5lgo')
view = py3Dmol.view()
view.addModel(base64.b64encode(r.content).decode(),'mmtf')
view.addUnitCell()
view.zoomTo()
```

- [ ] extract the CA carbon atoms of some chain and renumber

```
# https://colab.research.google.com/github/pb3lab/ibm3202/blob/master/tutorials/lab02_molviz.ipynb#scrollTo=jFANtPwvF2GK

#Do you want to know how many chains are contained in the PDB
#You just downloaded? grep can help you! 
#(remember to check your directory)
!grep 'COMPND'.*'CHAIN' 6ANE.pdb

#How many residues has each chain? Now we can use awk!
!awk '$1=="ATOM" && $3=="CA" && $5=="A" {print $0}' 6ANE.pdb | wc -l
```

- [ ] interesting methods to implement here too? (electrostatics, solvability) https://github.com/Electrostatics
- [ ] electrostatic coloring:

> Implement coloring an isosurface by a different volume data (e.g. isosurface drawn from electron density, colored by electrostatic potential). -- release 1.6.2, https://github.com/3dmol/3Dmol.js/releases

- http://bits.csb.pitt.edu/surfacemaker/
- https://www.biostars.org/p/299342/
- https://www.youtube.com/watch?v=WC-r53vBLvM

- [ ] disulfide bond prediction
- [ ] allow a binary mark for just the mapped domain, even though interacdome does not predict binding
- [ ] Analysis (distance, distance from surface, ...) https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
- [ ] select chain:
    - https://stackoverflow.com/questions/11685716/how-to-extract-chains-from-a-pdb-file
    - https://stackoverflow.com/questions/25677167/how-to-extract-all-chains-from-a-pdb-file
    - https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/dealing-with-coordinates

