from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name='glycogenius',
    version='0.2.2',
    author='Hector Franco Loponte',
    author_email='hectorfloponte@gmail.com',
    description='GlycoGenius is an all-in-one solution for data analysis of glycomics data.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="currently unavailable",
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    install_requires=["pandas", "scipy", "pyteomics",
                     "dill", "numpy"],
    entry_points={
        'console_scripts': [
            'glycogenius = glycogenius:glycogenius',
        ]
    }
)