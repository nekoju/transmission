from setuptools import setup, find_packages

setup(
        name='transmission',
        version='0.0.1',
        py_modules=['transmission'],
        description=('Tools for inferring symbiont transmission mode from'
                     'metagenomic data'),
        url='http://github.com/mpjuers/transmission',
        author='Mark Juers',
        author_email='mpjuers@indiana.edu',
        license='GPL>=3',
        packages=find_packages(),
        include_pakage_data=True,
        install_requires=[
            'click',
            'msprime',
            'numpy',
            'rpy2',
            'tqdm'
            ],
        entry_points='''
        [console_scripts]
        transmission-priorgen=transmission.transmission:simulate_prior_stats
        ''',
        zip_safe=False
    )
