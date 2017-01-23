#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
from glob import glob

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = [
    # TODO: put package requirements here
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='uncle_PSL',
    version='0.1.0',
    description="A BLAT to SAM converter.",
    long_description=readme,
    author="Botond Sipos",
    author_email='Botond.Sipos@nanoporetech.com',
    url='',
    packages=[
        'uncle_PSL',
    ],
    package_dir={'uncle_PSL':
                 'uncle_PSL'},
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords='uncle_PSL',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
    tests_require=test_requirements,
    scripts=[x for x in glob('scripts/*.py') if x != 'scripts/__init__.py']
)
