from setuptools import setup, find_packages

setup(
    name='ModiFinder',
    version='1.05',
    packages=find_packages(),
    install_requires=[
        'cairosvg',
        'rdkit',
        'requests',
        'xlsxwriter',
        'prettytable',
        'dash',
        'numpy',
        'tqdm',
        'IPython',
        'furl',
        'dash-bootstrap-components',
        'pillow',
        'openpyxl',
        'openpyxl-image-loader',
        'msbuddy',
        'lightgbm',
        'seaborn'
    ],
    author='Reza Shahneh',
    author_email='mzare008@ucr.edu',
    description='ModiFinder package',
)