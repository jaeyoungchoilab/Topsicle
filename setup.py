from setuptools import setup, find_packages

setup(
    name='Topsicle',
    version='1.0.0',
    packages=find_packages(),
    install_requires=[
        biopython=='1.75',
        matplotlib=='3.0.0',
        matplotlib-inline=='0.1.6',
        numpy=='1.19.5',
        pandas=='2.2.0',
        ruptures=='1.1.9',
        seaborn=='0.11.2',
    ],
    description="Telomere boundary identification in long read sequencing - Topsicle",
    author="Dr. Jae Young Choi",
    author_email="jaeyoung.choi@ku.edu",
    url="https://github.com/jychoilab",  # GitHub URL
    python_requires='>=3.6',

)
