from setuptools import setup, find_packages

setup(
    name='glycogenius',
    version='0.1.0',
    author='Hector Franco Loponte',
    author_email='hectorfloponte@gmail.com',
    description='GlycoGenius is an all-in-one solution for data analysis of glycomics data.',
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