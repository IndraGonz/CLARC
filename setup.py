from setuptools import setup, find_packages

setup(
    name='CLARC',
    version='1.0.1',  # Adding the linkage calculation
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'clarc = CLARC.__init__:main'
        ]
    },
    install_requires=[
        'numpy',
        'scipy',
        'pandas'
    ],
    author='Indra Gonzalez Ojeda',
    author_email='igonzalezojeda@g.harvard.edu',
    description='Package to cluster COG definitions based on a subpopulation of a general pangenome analysis (in Roary so far)',
    url='tbd',
)
