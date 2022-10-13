from setuptools import setup, find_packages

from pysname import __version__

setup(
    name='pysname',
    version=__version__,

    url='https://github.com/mjakuba/pysname',
    author='Michael Jakuba',
    author_email='mjakuba@whoi.edu',

    packages=find_packages(),
)
