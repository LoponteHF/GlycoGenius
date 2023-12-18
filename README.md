![GlycoGenius logo](./logo.png "GlycoGenius")

# GlycoGenius - Glycomics Data Analysis Tool

GlycoGenius is a python script that aims to be an all-in-one solution for data analysis of glycomics data.

Glycobiologists analyzing glycans' mass spectrometer data usually rely on several different tools to perform different tasks on different parts of their workflow and, in between the use of all those different tools, there's usually a lot of manual work that has to be done at least for data curating, which could also be considered cherry-picking in the research milieu.

With that in mind, this tool aims to put all the usual workflow for glycomics in a single place.

In order to do that, this tool is able to do several different tasks:
- Create glycans libraries based on user input, which can be monosaccharides numbers or specific glycans;
- Automatically identify noise level in samples;
- Process the spectra data and creates refined extracted ion chromatograms (EICs) for each glycan analyzed;
- Peak-picks multiple peaks in a single EIC, which allows the identification and possibly quantification of plausible isomers;
- Provides scorings of isotopic distribution peaks and chromatogram peak curve fitting based on relation and correlation;
- Identify PPM differences between theoretical mass and identified mz;
- Calculates signal-to-noise ratio;
- Identifies MS2 glycans' fragments and assign them to its respective precursor.

## Installation
~~~
	pip install glycogenius
~~~
## Usage

1. Export your MS data to an MzXML or MzML file;
2. Type 'glycogenius' in the terminal;
3. Follow instructions;
   - You can analyze directly on the CLI;
   - You can export a parameters file for advanced executions;
4. If you exported the parameters file, pipeline it to glycogenius after setting it up.
   - ie. in terminal type:
~~~
        cat .\glycogenius_parameters.ini | glycogenius
~~~
## Credits

Pyteomics:

> Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013) “Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software Prototyping in Proteomics”, Journal of The American Society for Mass Spectrometry, 24(2), 301–304. DOI: 10.1007/s13361-012-0516-6

> Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018) “Pyteomics 4.0: five years of development of a Python proteomics framework”, Journal of Proteome Research. DOI: 10.1021/acs.jproteome.8b00717

Dill for Python:

> M.M. McKerns, L. Strand, T. Sullivan, A. Fang, M.A.G. Aivazis, "Building a framework for predictive science", Proceedings of the 10th Python in Science Conference, 2011; http://arxiv.org/pdf/1202.1056

> wMichael McKerns and Michael Aivazis, "pathos: a framework for heterogeneous computing", 2010- ;	https://uqfoundation.github.io/project/pathos

Numpy:

> Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2. (Publisher link).

SciPy:

> Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.

Pandas:

> The pandas development team, Pandas, Zenoddo, Feb 2020, DOI:10.5281/zenodo.3509134

## License

This project is licensed under [GNU GPLv3 or later](https://spdx.org/licenses/GPL-3.0-or-later.html)