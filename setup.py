from setuptools import setup, find_packages

setup(
    name='glycogenius',
    version='0.1.0',
    author='Hector Franco Loponte',
    author_email='hectorfloponte@gmail.com',
    description='GlycoGenius is a python script that aims to be an all-in-one solution for data analysis of glycomics data.',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: GNU GPL v3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)