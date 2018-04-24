# Draft metadata proposal for VGP

This is a first draft of a metadata proposal meant to track important information about the raw data and files.
It is not currently written as a schema, but as an example document with comments and questions scattered throughout.
See [FAANG](https://github.com/FAANG/faang-metadata) and [HumanCellAtlas](https://github.com/HumanCellAtlas/metadata-schema) for more formalised schemas along similar lines that we could possibly aim for but would require time and dedicated resource.

Proposal would be to have/manage these YAML files in a github repo (e.g. `vgp-metadata`) with a `species/<species>/metadata.yaml` file structure mirroring the AWS bucket.
Advantages would be:

 * If we organise a well defined schema and validator the data could be validated via continuous integration with travis or the like.
 * Gives a place where internal and external users can point out mistakes and suggest changes or additions.
 * Gives an automatic history of changes.

Adding the files to AWS could probably be set up to happen via a git-hook?
But, we could just leave the AWS bucket for the data and the github repo for the metadata.
Do we want these files to be the source of populating the website/DB-backend?
Strategy would probably be OK for the ~260 ordinal species.
Beyond that we will hopefully have obtained some dedicated resource to develop a proper metadata solution.
