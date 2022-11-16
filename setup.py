"""Setup ARCTIC3D."""
from setuptools import find_packages, setup

with open("requirements.txt") as f:
    required = f.read().splitlines()

setup(
    name="arctic3d",
    license="Apache License 2.0",
    version="0.0.0",
    author="",
    description="",
    author_email="",
    include_package_data=True,
    packages=find_packages("src"),
    package_dir={"": "src"},
    classifiers=[],
    python_requires=">=3.10, <4",
    install_requires=required,
    entry_points={
        "console_scripts": [
            "arctic3d=arctic3d.cli:maincli",
            "arctic3d_resclust=arctic3d.cli_resclust:maincli",
            "arctic3d_localise=arctic3d.cli_localise:maincli",
        ],
    },
)
