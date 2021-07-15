import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "materials_visualization",
    version = "0.0.1",
    author = "Kevin Whitham",
    author_email = "kevin.whitham@gmail.com",
    description = "Routines to visualize atomic models in jupyter",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "",
    install_requires=[""],
    classifiers = [
            "Programming Language :: Python :: 3",
            "License :: OSI Approvied :: GNU General PUblic License v3 or later (GPLv3+)",
            "Operating System :: OS Independent",
    ],

    packages = setuptools.find_packages(),
    python_requires = ">=3",
)

