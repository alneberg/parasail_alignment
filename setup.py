#!/usr/bin/env python
"""Python setup.py for parasail_alignment package"""

import io
import os
from setuptools import find_packages, setup


def read(*paths, **kwargs):
    """Read the contents of a text file safely.
    >>> read("project_name", "VERSION")
    '0.1.0'
    >>> read("README.md")
    ...
    """

    content = ""
    with io.open(
        os.path.join(os.path.dirname(__file__), *paths),
        encoding=kwargs.get("encoding", "utf8"),
    ) as open_file:
        content = open_file.read().strip()
    return content


def read_requirements(path):
    return [
        line.strip()
        for line in read(path).split("\n")
        if not line.startswith(('"', "#", "-", "git+"))
    ]


setup(
    name="parasail_alignment",
    version=read("parasail_alignment", "VERSION"),
    description="Mostly exploratory code to try things out with parasail",
    url="https://github.com/alneberg/parasail_alignment/",
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    author="Johannes Alneberg",
    packages=find_packages(exclude=["tests", ".github"]),
    install_requires=read_requirements("requirements.txt"),
    entry_points={
        "console_scripts": ["parasail_alignment = parasail_alignment.__main__:main"]
    },
    extras_require={"test": read_requirements("requirements-test.txt")},
)