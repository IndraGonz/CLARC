from setuptools import setup, find_packages

setup(
    name='CLARC',
    version='1.1.8',
    packages=find_packages(),
    scripts=['scripts/acccog_blastn.sh'],
    entry_points={
        'console_scripts': [
            'clarc = CLARC.__init__:main'
        ]
    },
    install_requires=[],
    author='Indra Gonzalez Ojeda',
    author_email='igonzalezojeda@g.harvard.edu',
    description='Package to cluster COG definitions based on a subpopulation of a general pangenome analysis',
    url='https://github.com/IndraGonz/CLARC',
)
