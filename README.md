# mev-sctk-doubletfinder

This repository contains a WDL-format Cromwell-compatible workflow for executing the DoubletFinder algorithem (https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30073-0) for single-cell RNA-seq data as provided through the Single-Cell Toolkit (https://github.com/compbiomed/singleCellTK).

Outputs include a file containing the barcodes of likely doublets/multiplets *and* a count matrix subset where those cell barcodes have been removed.

To use, simply fill in the the `inputs.json` with the path to the single-cell counts file and submit to a Cromwell job runner.

Alternatively (if you do not want to use Cromwell), you can pull the docker image (https://github.com/web-mev/mev-sctk-doubletfinder/pkgs/container/mev-sctk-doubletfinder), start the container, and run: 

```
Rscript /opt/software/doubletfinder_qc.R \
    -f <path to raw counts tab-delimited file> \
    -o <prefix (string) for the subsetted count matrix filename> \
    -d <prefix (string) for the multiplet barcodes filename>
```