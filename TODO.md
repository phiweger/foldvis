### TODO

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
