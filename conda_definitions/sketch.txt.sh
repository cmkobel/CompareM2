
0c811240a515a9f4fea202f66f386e5c_.yaml  #Python 3.10.9   name: prokka

7747644331801b571bf2c3eb25e862b1_.yaml  #Python 3.11.0   name: r

#py 3.7
cat 948ee5f83efea16322d34bf146e28f62_.yaml  #Python 3.7.12   name: busco
cat a064302a62bbba97189271f7653c8ced_.yaml  #Python 3.7.12   name: mashtree
cat d13e2a550a90b04eb5c1744285ce8d27_.yaml  #Python 3.7.12   name: abricate
cat d38d6bd1aadd2a8eab503be05d1c1636_.yaml  #Python 3.7.12   name: mlst


# All currently use python 3.7.12
name: py_3_7_tools
channels:
  - defaults
  - conda-forge
  - bioconda
dependencies:
  - busco=5.4.4
  - mashtree
  - abricate
  - mlst





# py 3.8
cat 5510677d8a9035ea306e157b25016c12_.yaml  #Python 3.8.15   name: gtdbtk
cat f1498541056ec3a47322283507e25cdb_.yaml  #Python 3.8.15   name: checkm2_conda


name: py_3_8_tools
channels:
  - conda-forge
  - bioconda
dependencies:
  - gtdbtk
dependencies:
  - checkm2
  - gtdbtk



