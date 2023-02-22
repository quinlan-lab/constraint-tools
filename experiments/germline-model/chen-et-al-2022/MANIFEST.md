1. Chen + McHale scores on Chen windows, created using merge_chen_zscores_with_mchale_zscores.ipynb
```
${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.bed
```
2. The amount of overlap of Chen windows with things such as enhancers and exons: 
```
${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.overlapAmounts.bed
```
3. Classify Chen windows as overlapping enhancers, exons: 
```
${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale-enhancer-exon.bed
```
4. Additionally, label each window with a replication-timing score: 
```
${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale-enhancer-exon-replicationTiming.bed
```
5. Label Chen windows with their z-score quantiles: 
```
${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale-quantiles.bed
```
6. Add Clinvar and Mendelian variant counts to Chen windows: 
```
${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale-enhancer-exon-clinvar.bed
${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale-enhancer-exon-mendelian-variants.bed
```
7. Intersect noncoding SVs with Chen windows, created using intersect-noncoding-svs-with-windows.sh : 
```
{CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/{source}-noncoding-svs-chen-mchale.bed
{CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/{source}-noncoding-svs-chen-mchale-enhancer-exon.bed
```
