from setuptools import setup, find_packages

MAJOR = 1
MINOR = 0
PATCH = 0
VERSION = "{}.{}.{}".format(MAJOR, MINOR, PATCH)

with open("PeriProt/version.py", "w") as version_file:
    version_file.write("__version__ = '{}'\n".format(VERSION))


setup(
    name='PeriProt',
    version=VERSION,
    url='https://github.com/Moutoussamy/PeriProt',
    license='GPL3',
    author='Emmanuel Edouard Moutoussamy',
    author_email='e.e.moutoussamy@gmail.com',
    description='PeriProt: Tool box to analyse peripheral membrane protein structure and MD',
    platforms=["Linux", "Solaris", "Mac OS-X", "darwin", "Unix", "win32"],
    install_requires=['argparse',
                      'matplotlib',
                      'numpy',
                      'pandas',
                      'MDAnalysis ==  0.20.1'],

    entry_points={'console_scripts':['PeriProt=PeriProt.PeriProt:main']},

    packages=find_packages(),
)