# nf-variant-call

Genotyping of Metamist sequencing groups: CRAM -> single-sample reblocked gVCF with GATK HaplotypeCaller,
followed by gVCF QC. Both are published to and registered in each SG's own dataset.

Documentation:
- [docs/DESIGN.md](docs/DESIGN.md): design brief, flow, outputs, resources, status and next steps
- [GLOSSARY.md](GLOSSARY.md): domain terms, shared with nf-align and nf-combine
- [docs/RISKS.md](docs/RISKS.md): what may bite later
- [docs/adr/](docs/adr/): decisions with real trade-offs
