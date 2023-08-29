## Manifest 

Directory where data are stored: `{CONSTRAINT_TOOLS_DATA}/khurana`

`pacbio-deleted-enhancers.khurana-scores.hg19.csv` was renamed from `trioTypeDF.csv`,
and represents 21 novel enhancers completely deleted in a homozygous fashion in the three individuals sequenced by Chaisson et al.

`disease-enhancers.khurana-scores.hg19.csv` was renamed from `diseaseTypeDF.csv`,
and represents the 90 "disease enhancers" that Duo et al
extracted from Zhang et al 2018 ( https://academic.oup.com/nar/article/46/D1/D78/4559115 )

Also see: https://mail.google.com/mail/u/0/#inbox/QgrcJHsNjBschvpNVgQXNsdnxpDwSvJrLHb

## Workflow

1. [DONE] Convert enhancer data to bed format (for later transformation to hg38): 
`download-process-data/khurana/convert-to-bed.ipynb`
2. [DONE] Use the web interface to translate hg19 coordinates of enhancers to grch38 coordinates (liftover): 
https://genome.ucsc.edu/cgi-bin/hgLiftOver
3. [DONE] Merge hg38 coordinates with Khurana scores: `download-process-data/khurana/create-hg38-enhancers-and-khurana-scores.ipynb`
4. [DONE] Intersect enhancers with Chen windows and Chen scores: `experiments/germline-model/chen-et-al-2022/intersect-khurana-enhancers-with-chen-windows.sh`
5. [DONE] Reproduce Fig 4a of Khurana paper and extend to Chen scores: `experiments/germline-model/chen-et-al-2022/khurana.1.ipynb`


