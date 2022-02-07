"""
A setuptools based setup module.

See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

import io
from os import makedirs
from os import remove
from os.path import abspath
from os.path import dirname
from os.path import isfile
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
    name="blobtoolkit-pipeline",  # Required
    version="3.0.0",
    description="blobtoolkit pipeline",  # Optional
    long_description="blobtoolkit pipeline",  # Optional
    long_description_content_type="text/markdown",
    # long_description="%s\n%s"
    # % (
    #     re.compile("^.. start-badges.*^.. end-badges", re.M | re.S).sub(
    #         "", read("README.rst")
    #     ),
    #     re.sub(":[a-z]+:`~?(.*?)`", r"``\1``", read("CHANGELOG.rst")),
    # ),
    # long_description_content_type="text/x-rst",  # Optional (see note above)
    url="https://github.com/blobtoolkit/pipeline",  # Optional
    # This should be your name or the name of the organization which owns the
    # project.
    author="blobtoolkit",  # Optional
    # This should be a valid email address corresponding to the author listed
    # above.
    author_email="blobtoolkit@genomehubs.org",  # Optional
    # Classifiers help users find your project by categorizing it.
    #
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 5 - Production/Stable",
        # Indicate who your project is intended for
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        # Pick your license as you wish
        "License :: OSI Approved :: MIT License",
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate you support Python 3. These classifiers are *not*
        # checked by 'pip install'. See instead 'python_requires' below.
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3 :: Only",
    ],
    # This field adds keywords for your project which will appear on the
    # project page. What does your project relate to?
    #
    # Note that this is a list of additional keywords, separated
    # by commas, to be used to assist searching for the distribution in a
    # larger catalog.
    keywords="bioinformatics",  # Optional
    # When your source code is in a subdirectory under the project root, e.g.
    # `src/`, it is necessary to specify the `package_dir` argument.
    # package_dir={"": "."},  # Optional
    # You can just specify package directories manually here if your project is
    # simple. Or you can use find_packages().
    #
    # Alternatively, if you just want to distribute a single Python file, use
    # the `py_modules` argument instead as follows, which will expect a file
    # called `my_module.py` to exist:
    #
    #   py_modules=["my_module"],
    #
    packages=find_packages(where="."),  # Required
    # Specify which Python versions you support. In contrast to the
    # 'Programming Language' classifiers above, 'pip install' will check this
    # and refuse to install the project if the version does not match. See
    # https://packaging.python.org/guides/distributing-packages-using-setuptools/#python-requires
    python_requires=">=3.6, <4",
    # This field lists other packages that your project depends on to run.
    # Any package you put here will be installed by pip when your project is
    # installed, so they must be valid existing projects.
    #
    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
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
    ],  # Optional
    # List additional groups of dependencies here (e.g. development
    # dependencies). Users will be able to install these using the "extras"
    # syntax, for example:
    #
    #   $ pip install sampleproject[dev]
    #
    # Similar to `install_requires` above, these must be valid existing
    # projects.
    extras_require={  # Optional
        "dev": ["pycodestyle>=2.6.0", "pydocstyle>=5.0.2", "pylint>=2.5.3"],
        "test": [],
    },
    entry_points={
        "console_scripts": ["window-stats = pipeline.lib.window_stats:main",],
        # "pipeline.subcmd": [
        #     "add = blobtools.lib.add:cli",
        #     "create = blobtools.lib.add:cli",
        #     "filter = blobtools.lib.filter:cli",
        #     "host = blobtools.lib.host:cli",
        #     "remove = blobtools.lib.remove:cli",
        #     "replace = blobtools.lib.add:cli",
        #     "validate = blobtools.lib.validate:cli",
        #     "view = blobtools.lib.view:cli",
        # ],
    },
    project_urls={
        "Bug Reports": "https://github.com/blobtoolkit/pipeline/issues",
        "Source": "https://github.com/blobtoolkit/pipeline",
    },
    include_package_data=True,
    zip_safe=True,
)
