# RDF READER

## User Guide

### Installation

``` bash
pip install rdfreader
```

### Basic Usage

``` python
from rdfreader import RDFParser

rdf_file_name = "reactions.rdf"

with open(rdf_file_name, "r") as rdf_file:

    # create a RDFParser object, this is a generator that yields Reaction objects
    rdfreader = RDFParser(
        rdf_file,
        except_on_invalid_molecule=False,  # will return None instead of raising an exception if a molecule is invalid
        except_on_invalid_reaction=False,  # will return None instead of raising an exception if a reaction is invalid 
    )

    for rxn in rdfreader:
        if rxn is None:
            continue # the parser failed to read the reaction, go to the next one
  
        # rxn is a Reaction object, it is several attributes, including:
        print(rxn.smiles) # reaction SMILES string
        print(rxn.properties) # a dictionary of properties extracted from the RXN record
        
        reactants = rxn.reactants # a list of Molecule objects
        products = rxn.products
        solvents = rxn.solvents 
        catalysts = rxn.catalysts 
 
        # Molecule objects have several attributes, including:
        print(reactants[0].smiles)
        print(reactants[0].properties) # a dictionary of properties extracted from the MOL record (often empty)
        reactants[0].rd_mol # an RDKit molecule object
```

## Developer Guide

The project is managed and packaged using [poetry](https://python-poetry.org/docs/#installation).

### Installation

``` bash
git clone https://github.com/deepmatterltd/rdfreader
poetry install  # create a virtual environment and install the project dependencies
pre-commit install  # install pre-commit hooks, these mostly manage codestyle
```

### Contributions

Contributions are welcome via the [fork and pull request model](https://docs.github.com/en/get-started/quickstart/contributing-to-projects).

Before you commit changes, ensure these pass the hooks installed by pre-commit. This should be run automatically on each commit if you have run `pre-commit install`, but can be run manually from the terminal with `pre-commit run`.

### Releases

Releases are managed by GitHub releases/workflow. The version number in the pyproject file should ideally be kept up to date to the current release but is ignored by the release workflow.

To release a new version:
- Push the changes to GitHub.
- Use the github website to create a release. Tag the commit to be released with a version number, e.g. v1.2.3. The tag should be in v*.*.* format and be compatible with version numbering in poetry/pypi (the version number will be taken from the tag with the leading "v" stripped).
- When the release is published, a github workflow will run, build a wheel and publish it to PyPI.

### Example Data

You can find example data in the `test/resources directory`. `spresi-100.rdf` contains 100 example records from SPRESI.