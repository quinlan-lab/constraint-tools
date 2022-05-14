-- https://www.biostars.org/p/93011/#402520
-- http://uswest.ensembl.org/info/docs/api/core/core_schema.html

SELECT 
  seq_region.name AS chromosome,
  exon.seq_region_start AS exon_start,
  exon.seq_region_end AS exon_end,
  xref.display_label AS gene_name,
  exon_transcript.rank AS exon_rank,
  exon.seq_region_strand AS strand,
  exon.is_constitutive AS exon_is_constitutive,
  gene.biotype as gene_biotype
FROM gene 
-- https://www.w3schools.com/sql/sql_join_inner.asp
INNER JOIN exon_transcript ON exon_transcript.transcript_id = gene.canonical_transcript_id
INNER JOIN exon ON exon.exon_id = exon_transcript.exon_id
INNER JOIN seq_region ON seq_region.seq_region_id = gene.seq_region_id
-- https://www.w3schools.com/sql/sql_join_left.asp
LEFT JOIN xref ON xref.xref_id = gene.display_xref_id
-- WHERE seq_region.name='20' AND exon.seq_region_start > 63400000 AND exon.seq_region_end < 63500000
ORDER BY seq_region.name, exon.seq_region_start;

