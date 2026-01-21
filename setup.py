from setuptools import setup, find_packages
import os

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="SCiMS",  
    version="1.0.0",  
    author="Hanh Tran",
    description="SCiMS: Sex Calling in Metagenomic Sequencing",
    long_description=long_description,
    long_description_content_type="text/markdown",  
    url="https://github.com/hanhntran/SCiMS", 
    packages=find_packages(),  # Finds all packages (scims, tests, etc.)
    package_data={
        'scims': ['data/training_data_hmp_1000x_normalizedXY.txt'],
        # Optionally, if you want to include tests:
        'tests': ['*'],
    },
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
    ],
    python_requires=">=3.9",
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "setuptools",
    ],
    entry_points={
        "console_scripts": [
            "scims = scims.__main__:main",
        ],
    },
)
