from setuptools import setup, find_packages

setup(name='pyreactors',
      version='0.1',
      description='Python package for modeling reactor antineutrino spectra.',
      packages=find_packages(),
      # packages=['pyspectrum'],
      author_email='metalsera94@gmail.com',
      zip_safe=False,
      install_requires=[
          'numpy', 'pandas', 'tqdm'
      ],
      )
