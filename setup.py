"""
A setuptools based setup module.

See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

import io
from os.path import dirname
from os.path import join

from setuptools import find_packages
from setuptools import setup


def read(*names, **kwargs):
    """Read file."""
    with io.open(
        join(dirname(__file__), *names), encoding=kwargs.get("encoding", "utf8")
    ) as fh:
        return fh.read()


setup(
    name="blobtoolkit-pipeline",
    version="3.0.0",
    description="blobtoolkit pipeline",
    long_description="blobtoolkit pipeline",
    long_description_content_type="text/markdown",
    url="https://github.com/blobtoolkit/pipeline",
    author="blobtoolkit",
    author_email="blobtoolkit@genomehubs.org",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3 :: Only",
    ],
    keywords="bioinformatics",
    packages=find_packages(where="."),
    python_requires=">=3.6, <4",
    install_requires=[
        "blobtoolkit==3.0.0",
        "defusedxml==0.7.1",
        "docopt==0.6.2",
        "GitPython==3.1.26",
        "PyYAML==6.0",
        "requests==2.27.1",
        "snakemake==6.15.3",
        "tolkein==0.4.0",
        "ujson==5.1.0",
    ],
    extras_require={
        "dev": ["pycodestyle>=2.6.0", "pydocstyle>=5.0.2", "pylint>=2.5.3"],
        "test": [],
    },
    entry_points={
        "console_scripts": ["btk = btk:main",],
        "btk.subcmd": ["pipeline = pipeline:main",],
        "pipeline.subcmd": [
            "data = pipeline.lib.window_stats:main",
            "run = pipeline.lib.unchunk_blast:main",
            "chunk-fasta = pipeline.lib.chunk_fasta:main",
            "count-busco-genes = pipeline.lib.count_busco_genes:main",
            "extract-busco-genes = pipeline.lib.extract_busco_genes:main",
            "generate-config = pipeline.lib.generate_config:main",
            "generate-static-images = pipeline.lib.generate_static_images:main",
            "transfer-completed = pipeline.lib.transfer_completed:main",
            "unchunk-blast = pipeline.lib.unchunk_blast:main",
            "window-stats = pipeline.lib.window_stats:main",
        ],
    },
    project_urls={
        "Bug Reports": "https://github.com/blobtoolkit/pipeline/issues",
        "Source": "https://github.com/blobtoolkit/pipeline",
    },
    include_package_data=True,
    zip_safe=True,
)
