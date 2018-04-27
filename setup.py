from setuptools import setup#, find_packages

setup(
    name="hifan--mseutils",
    version="0.0.0",
    packages=["mseutils"],
    entry_points={
      'console_scripts': [
          'mseutils = mseutils.__main__:main',
          'csv2mgf = mseutils.csv2mgf:main'
      ]
    },
    )