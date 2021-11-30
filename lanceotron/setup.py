from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
   name='lanceotron',
   version='1.0.1',
   python_requires='>3.6', 
   description='Command-line interface to the lanceotron deep learning peak caller',
   long_description=long_description,
   long_description_content_type="text/markdown",
   url="http://chrisbcole.me/lanceotron/",
   project_urls={
      "Bug Tracker": "https://github.com/Chris1221/lanceotron/issues",
      "Source": "https://github.com/Chris1221/lanceotron",
   },
   author='Chris Cole, Lance Hentges',
   author_email='ccole@well.ox.ac.uk',
   packages=['lanceotron', 'lanceotron.static'],  
   include_package_data=True,
   keywords="deep learning, peak calling, keras",
   install_requires = [
      "pyBigWig",
      "sklearn",
      "numpy",
      "pandas",
      "tensorflow>2",
      "scipy",
      "natsort",
      "tqdm"
   ],
   entry_points={
        'console_scripts': [
            'lanceotron = lanceotron.cli:cli',
        ],
    }
)
