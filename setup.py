from setuptools import setup

setup(
        name='transmission',
        version='0.1',
        description=('tools for inferring symbiont transmission mode from host'
                     'metagenomic data'),
        url='http://github.com/mpjuers/transmission',
        author='Mark Juers',
        author_email='mpjuers@indiana.edu',
        license='GPL>=3',
        packages=['transmission'],
        install_requires=[
            'msprime',
            'numpy'
            ],
        zip_safe=False
        )
