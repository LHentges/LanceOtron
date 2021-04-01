To calculate the average coverage for the peak caller exclusive regions, the following DeepTools commands were carried out:

computeMatrix reference-point -S /public/lhentges/benchmarking_datasets/CTCF_spleen_ENCFF656CCY.bw /public/lhentges/benchmarking_datasets/control_spleen_ENCFF305CPV.bw /public/lhentges/benchmarking_datasets/DNase_spleen_ENCFF518TQH.bw -R ../peak_calls/2-LanceOtron-with-input_ChIP-seq_CTCF_spleen_3-MACS2-with-input_exclusive-peaks.bed ../peak_calls/3-MACS2-with-input_ChIP-seq_CTCF_spleen_2-LanceOtron-with-input_exclusive-peaks.narrowPeak --referencePoint center -a 1000 -b 1000 -out CTCF-spleen_LoT-and-MACS2_matrix.tab.gz

plotProfile -m CTCF-spleen_LoT-and-MACS2_matrix.tab.gz -out ../results/CTCF-spleen_LoT-and-MACS2_average_coverage.png --samplesLabel "CTCF" "control" "Open Chromatin" --regionsLabel "LanceOtron only" "MACS2 only" --plotType=heatmap
