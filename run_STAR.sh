#!/bash/bin
for fileName in $(cat sampleList); do
STAR --genomeDir /home/chienlab/Workspace/refs/STAR_hg38 \
      --genomeLoad LoadAndKeep \                                                                                     
      --readFilesIn ${fileName}_1.fq.gz ${fileName}_2.fq.gz \                                                                      
      --twopassMode Basic \                                                                                                      
      --outReadsUnmapped None \                                                                                                  
      --chimSegmentMin 12 \                                                                                                    
      --chimJunctionOverhangMin 12 \                                                                                           
      --alignSJDBoverhangMin 10 \                                                                                              
      --alignMatesGapMax 100000 \                                                                                             
      --alignIntronMax 100000 \                                                                                                
      --chimSegmentReadGapMax 3 \                                                                                    
      --alignSJstitchMismatchNmax 5 -1 5 5 \
      --runThreadN 6 \                                                                                                           
      --outSAMstrandField intronMotif \
      --outSAMunmapped Within \
      --outSAMtype BAM Unsorted \
      --outSAMattrRGline ID:GRPundef \
      --chimMultimapScoreRange 10 \
      --chimMultimapNmax 10 \
      --chimNonchimScoreDropMin 10 \
      --peOverlapNbasesMin 12 \
      --peOverlapMMp 0.1 \
      --chimOutJunctionFormat 1; # required as of STAR v2.6.1
done
STAR --genomeLoad Remove

