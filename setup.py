from setuptools import setup, find_packages

long_description_from_file = ""
with open("README.md", "r", encoding="utf-8") as f:
    for lines in f:
        if lines[0] != "!":
            long_description_from_file+= lines
    f.close()

setup(
    name='glycogenius',
    version='0.3.14',
    author='Hector Franco Loponte',
    author_email='hectorfloponte@gmail.com',
    description='GlycoGenius is an all-in-one solution for data analysis of glycomics data.',
    long_description=long_description_from_file,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    install_requires=["pandas", "scipy", "pyteomics",
                     "dill", "numpy", "lxml",
                     "openpyxl", "setuptools"],
    entry_points={
        'console_scripts': [
            'glycogenius = glycogenius:glycogenius',
        ]
    }
)