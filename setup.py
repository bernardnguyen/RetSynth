from setuptools import setup, find_packages

setup(name="RetSynth_gc",
      version="2.1.0",
      description="A retrosynthetic tool that can identify enzyme/reaction pairs that when added \
                   to a desired organism would produce a target chemical compound with the added feature that \
                    we identify optimal gene sequences",
      author=["Leanne Whitmore", "Bernard Nguyen", "Corey M. Hudson"],
      author_email=["coreymhudson@gmail.com",
                    "leanne382@gmail.com",
                    "bernguy@sandia.gov"
                    ],
      platforms=["linux",
                 "osx", 
                 "windows", 
                 "cygwin"
                 ],
      license="BSD 3 clause",
      url="https://github.com/sandialabs/RetSynth",
      classifiers=['Development Status :: 3 - Alpha',
                   'Intended Audience :: Developers, AI Researchers, Bioinformaticists, bioengineers',
		   'Topic :: integer linear programming :: flux balance analysis',
   		   'License :: BSD 3 clause',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 2.6',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3',
		   ],
      keywords='bioengineering, integer linear programming',
      test_suite="tests",
      packages=find_packages(),
      install_requires=[
          'argparse',
          'openbabel==2.4.0',
          'soupsieve==2.0.1',
          'filelock==3.0.12',
          'distlib==0.3.1',
          'jsonschema==3.2.0',
          'matplotlib==3.2.2',
          'graphviz==0.14',
          'glpk==0.4.6',
          'pulp==1.6.8',
          'cobra==0.14.1',
          'bs4',
          'beautifulsoup4',
          'python-libsbml-experimental==5.10.0',
          'pubchempy==1.0.4',
          'openpyxl==3.0.4',
          'tqdm==4.47.0',
          'scipy==1.5.1'
      ],
      include_package_data=True,
      zip_safe=False,
      scripts=[
          'rsgc/rs_gc.py',
      ]
)
