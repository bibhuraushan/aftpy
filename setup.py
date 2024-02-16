from setuptools import setup
# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

setup(
    name='aftpy',
    version='1.0.1',
    packages=['aftpy'],
    include_package_data=True,
    package_data={"style": ["bkj_style.mplstyle", "list_of_files.csv"]},
    url='https://github.com/bibhuraushan/aftpy/',
    license='MIT',
    author='Bibhuti Kumar Jha',
    author_email='bibhuraushan1@gmail.com',
    description='Python package to read and download AFTmap data.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords=['Sun', 'dataloader', 'AFT', 'SFT'],  # Keywords that define your package best
    install_requires=[  # I get to this in a second
        'matplotlib',
        'numpy',
        'astropy',
        'h5py',
        'sunpy',
        'pandas'
    ],
    classifiers=[
        'License :: OSI Approved :: MIT License',  # Again, pick a license,
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
)
