import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="peglit",
    version="1.0.1",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://www.github.com/sshen8/peglit",
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    install_requires=["pandas", "tqdm", "matplotlib", "Levenshtein", "scipy", "numpy", "sklearn"],
    entry_points={
        "console_scripts": [
            "peglit=peglit.peglit:main",
            "peglit.score=peglit.score:main",
            "peglit.inspect=peglit.inspect:main",
        ],
    },
)
