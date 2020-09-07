[![Build Status](https://travis-ci.org/sandialabs/RetroSynth.svg?branch=master)](https://travis-ci.org/sandialabs/RetSynth)

# RetSynth

The overaching goal of RetSynth is to streamline the arduous and complex step of selecting enzyme/reactions pairs to produce a target compound for bioengineering microbial organisms. Additionally, the Gene compatability component of this software identifies optimal gene sequences for added gene/enzymes based on the host organisms codon bias.

## Documentation

See documentation at http://sandialabs.github.io/RetroSynth/

## Build

The difficult part of ensuring that BioRetroSynth can run is installing the non-python dependencies which include:
	
    GNU/GLPK 	 Download from the website http://ftp.gnu.org/gnu/glpk/
	    
    GraphViz     Download from the website http://graphviz.org/ or using MacPorts

```bash
git clone https://github.com/sandialabs/RetSynth.git
pip install -r requirements.txt
python setup.py install
```

### Dependencies
-------------
RetSynth is tested to work under Python 2.7 and 3

* glpk==0.4.6
* pulp==1.6.8
* cobra==0.14.1
* bs4
* pygraphviz==1.3.1
* beautifulsoup4
* python-libsbml-experimental==5.10.0
* pubchempy==1.0.4
* openpyxl==3.0.4
* tqdm==4.47.0
* scipy==1.5.1
* openbabel==2.4.0
* soupsieve==2.0.1
* filelock==3.0.12
* distlib==0.3.1
* jsonschema==3.2.0
* matplotlib==3.2.2
* graphviz==0.14
* lxml
* bio

## License

BSD - 3-Clause Copyright 2017 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
