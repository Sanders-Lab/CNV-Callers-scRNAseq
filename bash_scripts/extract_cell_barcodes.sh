#!/bin/bash

samtools view BAM_TNBC1.bam | grep -oEh "CB:Z:[[:alnum:]]{16}-[[:digit:]]{1}" > tnbc1_cell_barcodes.txt
