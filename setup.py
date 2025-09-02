from setuptools import setup, find_packages

setup(
    name='Topsicle',
    version='1.0.0',
    packages=find_packages(),
    install_requires=[
        "biopython>=1.75",
        "cython>=0.29.21",
        "matplotlib>=3.3.4",
        "matplotlib-inline>=0.1.6",
        "numpy>=1.22.4",
        "pandas>=2.2.0",
        "ruptures==1.1.9",
        "seaborn>=0.11.2"
    ],
    entry_points={
        'console_scripts': [
            'topsicle = Topsicle.main:main'
        ]
    },
    description="Telomere boundary identification in long read sequencing - Topsicle",
    author="Dr. Jae Young Choi and Linh Nguyen",
    author_email="jaeyoung.choi@ku.edu and nguyen.linh.1010@ku.edu",
    url="https://github.com/jaeyoungchoilab/Topsicle",  
    python_requires='>=3.6',

)
