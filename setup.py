from setuptools import setup

setup(
    name                 = 'QRePS',
    version              = '1.1.0',
    description          = '''A utility for statistical comparison of protein identification results in proteomics and calculation of quantitative metrics for total cellular response''',
    long_description     = (''.join(open('README.md').readlines())),
    long_description_content_type="text/markdown",
    author               = 'Elizaveta Kazakova',
    author_email         = 'kazakova.em@phystech.edu',
    install_requires     = ['numpy', 'pandas', 'statsmodels', 'scipy', 'matplotlib', 'seaborn', 'scikit-learn', 'requests'],
    classifiers          = ['Intended Audience :: Science/Research',
                            'Programming Language :: Python :: 3',
                            'Topic :: Scientific/Engineering :: Bio-Informatics',
                            'Topic :: Scientific/Engineering :: Chemistry',
                            'Topic :: Scientific/Engineering :: Physics'],
    license              = 'License :: OSI Approved :: Apache Software License',
    packages             = ['QRePS'],
    entry_points         = {'console_scripts': ['qreps=QRePS.stattest_metrics:main']},
    )
