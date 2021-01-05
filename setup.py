from setuptools import setup

setup(
    name='binf',
    version='1.0',
    description='Bioinformatic toolkit',
    classifiers=[
        'Programming Language :: Python :: 3.7',
    ],
    url='https://gitlab.mff.cuni.cz/dorcakot/bioinformatic-toolbox',
    install_requires=[
        'biopython',
        'numpy',
        'matplotlib'
    ],
    include_package_data=True,
    zip_safe=False,
    packages=[
        'binf',
    ],
    entry_points={
        'console_scripts': [
            'binf=binf.main:main',
        ],
    },
)
